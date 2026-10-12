#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "CandidateDBReader.h"
#include "Reporter.h"
#include "TaxonomyWrapper.h"
#include "DBReader.h"
#include "GroupGenerator.h"
#include "Util.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

constexpr size_t MAX_REPORT_ROWS = 50;
constexpr size_t OUTPUT_CHUNK = 65536;


constexpr double ABUNDANCE_EXPONENT = 1.0;
// Floor on the iteration count so a single lucky step cannot report
// convergence before the abundances have moved at all.
constexpr int EM_MIN_ITER = 3;

constexpr double SQUAREM_STEP_MIN = -20.0;
constexpr double SQUAREM_STEP_MAX = 1.0;
// An aggressive step is only kept if it does not lower the likelihood; otherwise
// it is halved back towards -1 this many times before giving up on it.
constexpr int SQUAREM_MAX_BACKTRACK = 4;

void normalize(std::vector<double> &x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += x[i];
    }
    if (sum > 0.0) {
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] /= sum;
        }
    }
}

void addScoreToClade(
    TaxonomyWrapper *taxonomy,
    TaxID taxId,
    double score,
    std::unordered_map<TaxID, double> &cladeScoreSums)
{
    if (taxId == 0 || !taxonomy->nodeExists(taxId)) {
        return;
    }
    const TaxonNode *taxon = taxonomy->taxonNode(taxId);
    cladeScoreSums[taxon->taxId] += score;
    while (taxon->parentTaxId != taxon->taxId && taxonomy->nodeExists(taxon->parentTaxId)) {
        taxon = taxonomy->taxonNode(taxon->parentTaxId);
        cladeScoreSums[taxon->taxId] += score;
    }
}

void addAverageScoresToMap(
    const std::unordered_map<TaxID, double> &cladeScoreSums,
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    std::unordered_map<TaxID, double> &avgScores)
{
    avgScores.reserve(avgScores.size() + cladeScoreSums.size());
    for (const auto &taxScore : cladeScoreSums) {
        auto countIt = cladeCounts.find(taxScore.first);
        if (countIt == cladeCounts.end() || countIt->second.cladeCount == 0) {
            continue;
        }
        avgScores[taxScore.first] = taxScore.second / static_cast<double>(countIt->second.cladeCount);
    }
}

// ---------------------------------------------------------------------------
// Candidate DB -> EM observations
// ---------------------------------------------------------------------------

// One (query, candidate species) pair: the unit the E-step assigns responsibility to.
// A query is one read, or, with --read-groupings 1, one read group.
struct Observation {
    uint32_t spIdx;   // index into CandidateSet::speciesList
    float weight;     // score ^ --score-power
};

struct CandidateSet {
    // Observations grouped by EM query: query q owns [queryOffsets[q], queryOffsets[q+1]).
    // Without read grouping, query q is candidate DB entry q.
    std::vector<Observation> obs;
    std::vector<uint32_t> queryOffsets;
    std::vector<TaxID> speciesList;                 // compact index -> species taxId
    std::unordered_map<TaxID, uint32_t> sp2idx;
    // How many reads each query stands for in the M-step. Empty means 1 for every
    // query: the per-read EM, or --group-weight 0.
    std::vector<uint32_t> queryMass;
    size_t readsWithCandidates = 0;    // candidate DB entries with >= 1 candidate
    size_t queriesWithCandidates = 0;  // EM queries with >= 1 observation
    double totalMass = 0.0;            // R: summed mass of the queries with observations

    size_t queryCount() const {
        return queryOffsets.empty() ? 0 : queryOffsets.size() - 1;
    }

    double massOf(size_t q) const {
        return queryMass.empty() ? 1.0 : static_cast<double>(queryMass[q]);
    }
};

bool loadObservations(
    CandidateDBReader &reader,
    size_t entryCount,
    const LocalParameters &par,
    CandidateSet &cs)
{
    const double scorePower = par.scorePower;
    cs.queryOffsets.reserve(entryCount + 1);
    cs.queryOffsets.push_back(0);

    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (reader.getByIndex(i, entry, 0)) {
            for (const SpeciesCandidate &candidate : entry.candidates) {
                auto inserted = cs.sp2idx.emplace(candidate.speciesId,
                                                  static_cast<uint32_t>(cs.speciesList.size()));
                if (inserted.second) {
                    cs.speciesList.push_back(candidate.speciesId);
                }
                const double weight = std::pow(static_cast<double>(candidate.idScore), scorePower);
                Observation observation;
                observation.spIdx = inserted.first->second;
                observation.weight = static_cast<float>(weight);
                cs.obs.push_back(observation);
            }
            if (!entry.candidates.empty()) {
                ++cs.readsWithCandidates;
            }
        }
        // queryOffsets indexes obs with 32 bits; refuse rather than silently wrap.
        if (cs.obs.size() > static_cast<size_t>(UINT32_MAX)) {
            std::cerr << "Error: more than " << UINT32_MAX << " (read, candidate) pairs in the candidate DB."
                      << std::endl;
            return false;
        }
        cs.queryOffsets.push_back(static_cast<uint32_t>(cs.obs.size()));
    }
    cs.queriesWithCandidates = cs.readsWithCandidates;
    cs.totalMass = static_cast<double>(cs.readsWithCandidates);
    return !cs.speciesList.empty();
}

// ---------------------------------------------------------------------------
// Read groups (--read-groupings 1)
// ---------------------------------------------------------------------------
//
// The grouping step writes two text files into its output directory:
//   groupMap : one row per read, "<read id>\t<group id>". Read ids are 1-based
//              positions in the read file(s); group id 0 = the read is in no group.
//   groups   : one row per group, "<group id>\t<member read id>\t..." (any order).
// create-candidates keys its entries by 0-based read position, so read id r is
// candidate DB key r - 1. Both steps must have seen the same reads in the same
// order; the row count and the groups/groupMap cross-check below catch the
// obvious mismatches, nothing can catch a reordered input.

constexpr uint32_t GROUP_UNSET = UINT32_MAX;

struct ReadGrouping {
    std::vector<uint32_t> groupOfKey;   // candidate DB key -> group id (0 = ungrouped)
    size_t groupedReads = 0;
    size_t groupCount = 0;              // distinct non-zero group ids
};

// Parses the unsigned integers on `line` (tab or space separated, trailing
// separators allowed). Returns false on any other token.
bool parseUnsignedRow(const std::string &line, std::vector<uint64_t> &values) {
    values.clear();
    size_t i = 0;
    const size_t n = line.size();
    while (i < n) {
        while (i < n && (line[i] == '\t' || line[i] == ' ' || line[i] == '\r')) {
            ++i;
        }
        if (i >= n) {
            break;
        }
        if (!std::isdigit(static_cast<unsigned char>(line[i]))) {
            return false;
        }
        uint64_t value = 0;
        while (i < n && std::isdigit(static_cast<unsigned char>(line[i]))) {
            value = value * 10 + static_cast<uint64_t>(line[i] - '0');
            ++i;
        }
        values.push_back(value);
    }
    return true;
}

void openLargeTextFile(std::ifstream &in, std::vector<char> &buffer, const std::string &path) {
    buffer.resize(16 * 1024 * 1024);
    in.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    in.open(path);
}

bool loadGroupMap(const std::string &path, size_t entryCount, ReadGrouping &grouping) {
    std::ifstream in;
    std::vector<char> ioBuffer;
    openLargeTextFile(in, ioBuffer, path);
    if (!in.is_open()) {
        std::cerr << "Error: could not open the read-to-group map " << path << "." << std::endl;
        return false;
    }
    std::cout << "Read-to-group map     : " << path << std::endl;

    grouping.groupOfKey.assign(entryCount, GROUP_UNSET);
    grouping.groupedReads = 0;
    size_t rows = 0;
    std::string line;
    std::vector<uint64_t> values;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!parseUnsignedRow(line, values) || values.size() != 2) {
            std::cerr << "Error: malformed row " << (rows + 1) << " in " << path
                      << " (expected \"<read id>\\t<group id>\"): " << line << std::endl;
            return false;
        }
        ++rows;
        const uint64_t readId = values[0];
        const uint64_t groupId = values[1];
        if (readId == 0 || readId > entryCount) {
            std::cerr << "Error: " << path << " names read " << readId << ", but the candidate DB has "
                      << entryCount << " entries (read ids are 1-based)." << std::endl;
            std::cerr << "       The grouping and create-candidates must be run on the same reads in the same order."
                      << std::endl;
            return false;
        }
        if (groupId >= GROUP_UNSET) {
            std::cerr << "Error: group id " << groupId << " in " << path << " is out of range." << std::endl;
            return false;
        }
        uint32_t &slot = grouping.groupOfKey[readId - 1];
        if (slot != GROUP_UNSET) {
            std::cerr << "Error: read " << readId << " appears twice in " << path << "." << std::endl;
            return false;
        }
        slot = static_cast<uint32_t>(groupId);
        if (groupId != 0) {
            ++grouping.groupedReads;
        }
    }
    if (rows != entryCount) {
        std::cerr << "Error: " << path << " has " << rows << " reads but the candidate DB has " << entryCount
                  << " entries." << std::endl;
        std::cerr << "       Both must come from the same read file(s) in the same order." << std::endl;
        return false;
    }
    // rows == entryCount with distinct in-range ids means every key was set.

    std::vector<uint32_t> ids;
    ids.reserve(grouping.groupedReads);
    for (const uint32_t g : grouping.groupOfKey) {
        if (g != 0) {
            ids.push_back(g);
        }
    }
    std::sort(ids.begin(), ids.end());
    grouping.groupCount = static_cast<size_t>(std::unique(ids.begin(), ids.end()) - ids.begin());
    return true;
}

// The groups file repeats what groupMap says; reading it only guards against the
// two files coming from different runs.
bool checkGroupsFile(const std::string &path, size_t entryCount, const ReadGrouping &grouping) {
    std::ifstream in;
    std::vector<char> ioBuffer;
    openLargeTextFile(in, ioBuffer, path);
    if (!in.is_open()) {
        std::cout << "Warning: " << path << " is not readable; the group map is used unchecked." << std::endl;
        return true;
    }
    std::cout << "Group member list     : " << path << std::endl;

    size_t groups = 0;
    size_t members = 0;
    std::string line;
    std::vector<uint64_t> values;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (!parseUnsignedRow(line, values) || values.size() < 2) {
            std::cerr << "Error: malformed row " << (groups + 1) << " in " << path
                      << " (expected \"<group id>\\t<read id>...\"): " << line << std::endl;
            return false;
        }
        ++groups;
        const uint64_t groupId = values[0];
        for (size_t i = 1; i < values.size(); ++i) {
            const uint64_t readId = values[i];
            if (readId == 0 || readId > entryCount) {
                std::cerr << "Error: group " << groupId << " in " << path << " lists read " << readId
                          << ", but the candidate DB has " << entryCount << " entries." << std::endl;
                return false;
            }
            if (grouping.groupOfKey[readId - 1] != groupId) {
                std::cerr << "Error: read " << readId << " is in group " << groupId << " according to " << path
                          << " but in group " << grouping.groupOfKey[readId - 1]
                          << " according to the group map. The two files are from different runs." << std::endl;
                return false;
            }
            ++members;
        }
    }
    if (groups != grouping.groupCount || members != grouping.groupedReads) {
        std::cerr << "Error: " << path << " lists " << groups << " groups over " << members
                  << " reads, but the group map has " << grouping.groupCount << " groups over "
                  << grouping.groupedReads << " reads. The two files are from different runs." << std::endl;
        return false;
    }
    return true;
}

// Side information about the EM queries when reads are grouped.
struct GroupedQueries {
    uint32_t groupCount = 0;                   // queries [0, groupCount) are read groups, ascending group id
    std::vector<uint32_t> groupIds;            // group query -> group id
    std::vector<uint32_t> readQuery;           // candidate DB index -> EM query
    std::vector<uint32_t> queryReads;          // EM query -> reads it stands for (1 for a single read)
    std::vector<uint32_t> groupReadsWithCand;  // group query -> members with >= 1 candidate
    std::vector<float> obsMeanScore;           // mean idScore of obs[0 .. queryOffsets[groupCount])
    size_t maxGroupReads = 0;
};

// Pools the candidates of every read group into one EM query. The score of a
// group for a species is the sum of its members' idScores for that species divided
// by the group's read count (members without candidates count in the divisor), and
// the observation weight is that score ^ --score-power, exactly as for a single
// read. Reads in no group become single-read queries after the group queries.
bool loadGroupedObservations(
    CandidateDBReader &reader,
    size_t entryCount,
    const LocalParameters &par,
    const ReadGrouping &grouping,
    CandidateSet &cs,
    GroupedQueries &gq)
{
    const double scorePower = par.scorePower;

    // Group ids -> compact query indices, ascending by group id.
    {
        std::vector<uint32_t> ids;
        ids.reserve(grouping.groupedReads);
        for (const uint32_t g : grouping.groupOfKey) {
            if (g != 0) {
                ids.push_back(g);
            }
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        gq.groupIds.swap(ids);
    }
    gq.groupCount = static_cast<uint32_t>(gq.groupIds.size());
    gq.readQuery.assign(entryCount, 0);
    gq.queryReads.assign(gq.groupCount, 0);
    gq.groupReadsWithCand.assign(gq.groupCount, 0);

    // Pass 1: every read's raw (species, idScore) pairs in candidate DB order, and
    // the EM query each read belongs to.
    std::vector<Observation> raw;     // weight holds the raw idScore here
    std::vector<uint32_t> rawOffsets;
    rawOffsets.reserve(entryCount + 1);
    rawOffsets.push_back(0);

    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        const bool readable = reader.getByIndex(i, entry, 0);
        uint32_t groupId = 0;
        if (readable) {
            if (entry.queryId >= grouping.groupOfKey.size()) {
                std::cerr << "Error: candidate DB key " << entry.queryId << " has no row in the group map ("
                          << grouping.groupOfKey.size() << " reads)." << std::endl;
                return false;
            }
            groupId = grouping.groupOfKey[entry.queryId];
        }

        uint32_t q;
        if (groupId != 0) {
            q = static_cast<uint32_t>(std::lower_bound(gq.groupIds.begin(), gq.groupIds.end(), groupId)
                                      - gq.groupIds.begin());
            ++gq.queryReads[q];
        } else {
            // Ungrouped (or unreadable) entries are single-read queries.
            q = static_cast<uint32_t>(gq.queryReads.size());
            gq.queryReads.push_back(1);
        }
        gq.readQuery[i] = q;

        if (readable) {
            for (const SpeciesCandidate &candidate : entry.candidates) {
                auto inserted = cs.sp2idx.emplace(candidate.speciesId,
                                                  static_cast<uint32_t>(cs.speciesList.size()));
                if (inserted.second) {
                    cs.speciesList.push_back(candidate.speciesId);
                }
                Observation observation;
                observation.spIdx = inserted.first->second;
                observation.weight = candidate.idScore;
                raw.push_back(observation);
            }
            if (!entry.candidates.empty()) {
                ++cs.readsWithCandidates;
                if (q < gq.groupCount) {
                    ++gq.groupReadsWithCand[q];
                }
            }
        }
        if (raw.size() > static_cast<size_t>(UINT32_MAX)) {
            std::cerr << "Error: more than " << UINT32_MAX << " (read, candidate) pairs in the candidate DB."
                      << std::endl;
            return false;
        }
        rawOffsets.push_back(static_cast<uint32_t>(raw.size()));
    }
    if (cs.speciesList.empty()) {
        return false;
    }

    const size_t queryCount = gq.queryReads.size();
    for (uint32_t g = 0; g < gq.groupCount; ++g) {
        gq.maxGroupReads = std::max<size_t>(gq.maxGroupReads, gq.queryReads[g]);
    }

    // Pass 2: bucket the raw pairs by query (counting sort keeps the read order
    // inside a query, so the result does not depend on the thread count).
    std::vector<uint32_t> bucketOffsets(queryCount + 1, 0);
    for (size_t i = 0; i < entryCount; ++i) {
        bucketOffsets[gq.readQuery[i] + 1] += rawOffsets[i + 1] - rawOffsets[i];
    }
    for (size_t q = 0; q < queryCount; ++q) {
        bucketOffsets[q + 1] += bucketOffsets[q];
    }
    std::vector<Observation> pooled(raw.size());
    {
        std::vector<uint32_t> cursor(bucketOffsets.begin(), bucketOffsets.end() - 1);
        for (size_t i = 0; i < entryCount; ++i) {
            uint32_t &dst = cursor[gq.readQuery[i]];
            for (uint32_t j = rawOffsets[i]; j < rawOffsets[i + 1]; ++j) {
                pooled[dst++] = raw[j];
            }
        }
    }
    raw.clear();
    raw.shrink_to_fit();
    rawOffsets.clear();
    rawOffsets.shrink_to_fit();

    // Pass 3: inside each query, merge the pairs of one species into its mean
    // idScore over the query's reads. Sorting by species index first makes the
    // merge a single sweep and the observation order canonical.
    std::vector<uint32_t> mergedLen(queryCount, 0);
#ifdef OPENMP
#pragma omp parallel for schedule(dynamic, 4096) shared(pooled, bucketOffsets, mergedLen, gq, queryCount)
#endif
    for (size_t q = 0; q < queryCount; ++q) {
        const uint32_t begin = bucketOffsets[q];
        const uint32_t end = bucketOffsets[q + 1];
        if (begin == end) {
            continue;
        }
        std::sort(pooled.begin() + begin, pooled.begin() + end,
                  [](const Observation &a, const Observation &b) { return a.spIdx < b.spIdx; });
        const double reads = static_cast<double>(gq.queryReads[q]);
        uint32_t w = begin;
        for (uint32_t j = begin; j < end;) {
            const uint32_t s = pooled[j].spIdx;
            double sum = 0.0;
            while (j < end && pooled[j].spIdx == s) {
                sum += static_cast<double>(pooled[j].weight);
                ++j;
            }
            pooled[w].spIdx = s;
            pooled[w].weight = static_cast<float>(sum / reads);   // mean idScore, in [0, 1]
            ++w;
        }
        mergedLen[q] = w - begin;
    }

    // Pass 4: compact, and turn the mean scores into EM weights. The means of the
    // group queries are kept for the group table.
    cs.queryOffsets.assign(queryCount + 1, 0);
    gq.obsMeanScore.clear();
    size_t dst = 0;
    for (size_t q = 0; q < queryCount; ++q) {
        cs.queryOffsets[q] = static_cast<uint32_t>(dst);
        const uint32_t begin = bucketOffsets[q];
        const uint32_t len = mergedLen[q];
        for (uint32_t k = 0; k < len; ++k) {
            const float mean = pooled[begin + k].weight;
            const uint32_t spIdx = pooled[begin + k].spIdx;
            pooled[dst].spIdx = spIdx;
            pooled[dst].weight = static_cast<float>(std::pow(static_cast<double>(mean), scorePower));
            if (q < gq.groupCount) {
                gq.obsMeanScore.push_back(mean);
            }
            ++dst;
        }
    }
    cs.queryOffsets[queryCount] = static_cast<uint32_t>(dst);
    pooled.resize(dst);
    pooled.shrink_to_fit();
    cs.obs.swap(pooled);

    if (par.groupWeight != 0) {
        cs.queryMass = gq.queryReads;
    } else {
        cs.queryMass.clear();
    }
    cs.queriesWithCandidates = 0;
    cs.totalMass = 0.0;
    for (size_t q = 0; q < queryCount; ++q) {
        if (cs.queryOffsets[q + 1] > cs.queryOffsets[q]) {
            ++cs.queriesWithCandidates;
            cs.totalMass += cs.massOf(q);
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Per-species length normalization (representative genome size)
// ---------------------------------------------------------------------------

// Reads <DB dir>/species2genomeSize.tsv, which build writes one row per species:
//   <internal species id>\t<original taxID>\t<representative genome length (nt)>
// The representative is a single assembly (the one named for the species, else the
// most contiguous one near the median length), so the size is that of one genome.
// Unlike the species' summed distinct k-mers it does not grow with the number of
// assemblies the DB holds for the species.
bool loadGenomeSizes(
    const std::string &dbDir,
    std::unordered_map<TaxID, uint64_t> &sp2genomeSize)
{
    const std::string fileName = dbDir + "/species2genomeSize.tsv";
    std::ifstream in(fileName);
    if (!in.is_open()) {
        std::cout << "Warning: could not open " << fileName << "." << std::endl;
        return false;
    }
    std::cout << "Genome sizes          : " << fileName << std::endl;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        const size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) {
            continue;
        }
        const size_t tab2 = line.find('\t', tab1 + 1);
        if (tab2 == std::string::npos) {
            continue;
        }
        // Exception-free parsing (the project is built with -fno-exceptions).
        const std::string spIdStr = line.substr(0, tab1);
        const std::string sizeStr = line.substr(tab2 + 1);
        char *endPtr = nullptr;
        const unsigned long long spId = std::strtoull(spIdStr.c_str(), &endPtr, 10);
        if (endPtr == spIdStr.c_str()) {
            continue;
        }
        const unsigned long long genomeSize = std::strtoull(sizeStr.c_str(), &endPtr, 10);
        if (endPtr == sizeStr.c_str()) {
            continue;
        }
        sp2genomeSize[static_cast<TaxID>(spId)] = static_cast<uint64_t>(genomeSize);
    }
    return true;
}

// lengthFactor is the 1/L of P(read | species) = score^p / L: a read is one of the
// ~L start positions of the species' genome, so a larger genome makes any given
// read less likely and does not attract reads on size alone. L is the
// representative genome size; needs <DB dir>/species2genomeSize.tsv.
std::vector<double> buildLengthFactors(
    const std::string &dbDir,
    const CandidateSet &cs)
{
    std::vector<double> lengthFactor(cs.speciesList.size(), 1.0);

    std::unordered_map<TaxID, uint64_t> sp2genomeSize;
    if (!loadGenomeSizes(dbDir, sp2genomeSize)) {
        std::cout << "Warning: no genome sizes are available; "
                  << "using a length factor of 1.0 for every species." << std::endl;
        return lengthFactor;
    }

    std::vector<size_t> missingIdx;
    std::vector<double> knownFactors;
    knownFactors.reserve(cs.speciesList.size());
    for (size_t i = 0; i < cs.speciesList.size(); ++i) {
        const auto it = sp2genomeSize.find(cs.speciesList[i]);
        if (it == sp2genomeSize.end() || it->second == 0) {
            missingIdx.push_back(i);
            continue;
        }
        lengthFactor[i] = 1.0 / static_cast<double>(it->second);
        knownFactors.push_back(lengthFactor[i]);
    }

    // Unknown genome size: a 0.0 factor would silently delete the species, and 1.0
    // would weigh it ~genome-size times above every sized species, letting it take
    // every read it competes for. The median factor of the sized species is neutral.
    // With no sized species at all, every factor stays 1.0 (no normalization).
    double fillFactor = 1.0;
    if (!knownFactors.empty()) {
        const size_t mid = knownFactors.size() / 2;
        std::nth_element(knownFactors.begin(), knownFactors.begin() + mid, knownFactors.end());
        fillFactor = knownFactors[mid];
    }
    for (const size_t i : missingIdx) {
        lengthFactor[i] = fillFactor;
    }

    std::cout << "Length factor         : 1/genomeSize"
              << " (" << missingIdx.size() << "/" << cs.speciesList.size() << " species without a size";
    if (!missingIdx.empty()) {
        std::cout << ", given the median factor " << fillFactor;
    }
    std::cout << ")" << std::endl;
    return lengthFactor;
}


// ---------------------------------------------------------------------------
// EM
// ---------------------------------------------------------------------------

constexpr double LL_FLOOR = 1e-300;

// Mass-weighted mean log-likelihood: a query standing for n reads counts n times.
double meanLogLikelihood(
    const CandidateSet &cs,
    const std::vector<double> &theta,
    const std::vector<double> &lengthFactor,
    double alpha)
{
    const size_t queryCnt = cs.queryCount();
    double total = 0.0;
#ifdef OPENMP
#pragma omp parallel for schedule(static) reduction(+:total) \
        shared(cs, theta, lengthFactor, alpha, queryCnt)
#endif
    for (size_t q = 0; q < queryCnt; ++q) {
        const uint32_t begin = cs.queryOffsets[q];
        const uint32_t end = cs.queryOffsets[q + 1];
        if (begin == end) {
            continue;
        }
        double denom = 0.0;
        for (uint32_t j = begin; j < end; ++j) {
            const uint32_t s = cs.obs[j].spIdx;
            denom += static_cast<double>(cs.obs[j].weight)
                     * std::pow(theta[s], alpha)
                     * lengthFactor[s];
        }
        total += cs.massOf(q) * std::log(std::max(denom, LL_FLOOR));
    }
    return (cs.totalMass > 0.0)
        ? total / cs.totalMass
        : 0.0;
}

// One EM update. Returns the mean log-likelihood of the INPUT theta, which the
// per-query denominators give away for free.
double emStep(
    const CandidateSet &cs,
    const std::vector<double> &theta,
    const std::vector<double> &lengthFactor,
    const std::vector<double> &priorTerm,
    double alpha,
    std::vector<double> &out)
{
    out = priorTerm;
    const size_t queryCnt = cs.queryCount();
    double llTotal = 0.0;

#ifdef OPENMP
#pragma omp parallel shared(cs, theta, lengthFactor, out, alpha, queryCnt, llTotal)
#endif
    {
        std::vector<double> local(out.size(), 0.0);
        std::vector<double> weights;
        double llLocal = 0.0;
#ifdef OPENMP
#pragma omp for schedule(static)
#endif
        for (size_t q = 0; q < queryCnt; ++q) {
            const uint32_t begin = cs.queryOffsets[q];
            const uint32_t end = cs.queryOffsets[q + 1];
            if (begin == end) {
                continue;
            }
            weights.clear();
            double denom = 0.0;
            for (uint32_t j = begin; j < end; ++j) {
                const uint32_t s = cs.obs[j].spIdx;
                const double weight = static_cast<double>(cs.obs[j].weight)
                                      * std::pow(theta[s], alpha)
                                      * lengthFactor[s];
                weights.push_back(weight);
                denom += weight;
            }
            const double mass = cs.massOf(q);
            llLocal += mass * std::log(std::max(denom, LL_FLOOR));
            if (denom <= 0.0) {
                continue;
            }
            for (uint32_t j = begin; j < end; ++j) {
                local[cs.obs[j].spIdx] += mass * (weights[j - begin] / denom);
            }
        }
#ifdef OPENMP
#pragma omp critical
#endif
        {
            for (size_t s = 0; s < local.size(); ++s) {
                out[s] += local[s];
            }
            llTotal += llLocal;
        }
    }

    normalize(out);
    return (cs.totalMass > 0.0)
        ? llTotal / cs.totalMass
        : 0.0;
}

std::vector<double> runEm(
    const CandidateSet &cs,
    const std::vector<double> &lengthFactor,
    const LocalParameters &par)
{
    const size_t speciesCnt = cs.speciesList.size();

    // Uniform Dirichlet pseudocount, constant across iterations and species.
    // Expressed as a fraction of the total read mass rather than as an absolute
    // count, so it does not need retuning per sample:
    //   priorTerm[s] = phi * R / S,  total prior mass = phi * R
    // With phi = 0 the update is plain EM.
    const double uniformTerm =
        (speciesCnt > 0)
            ? std::max(0.0, par.priorMass) * cs.totalMass
                  / static_cast<double>(speciesCnt)
            : 0.0;
    const std::vector<double> priorTerm(speciesCnt, uniformTerm);

    std::vector<double> x0(speciesCnt, 1.0 / static_cast<double>(speciesCnt));
    std::vector<double> x1(speciesCnt, 0.0);
    std::vector<double> x2(speciesCnt, 0.0);
    std::vector<double> xNew(speciesCnt, 0.0);

    const time_t start = time(nullptr);
    int iterationsRun = 0;
    bool converged = false;
    int backtracks = 0;     // aggressive steps halved back towards -1
    int fallbacks = 0;      // iterations that ended up at the plain double-EM point
    std::vector<double> r(speciesCnt);
    std::vector<double> v(speciesCnt);
    for (int iter = 0; iter < par.maxIter; ++iter) {
        iterationsRun = iter + 1;

        // SQUAREM: two EM steps, then extrapolate along their difference.
        // Both calls hand back the log-likelihood of their input for free.
        emStep(cs, x0, lengthFactor, priorTerm, ABUNDANCE_EXPONENT, x1);
        const double llx1 = emStep(cs, x1, lengthFactor, priorTerm, ABUNDANCE_EXPONENT, x2);

        double normR = 0.0, normV = 0.0;
        for (size_t i = 0; i < speciesCnt; ++i) {
            r[i] = x1[i] - x0[i];
            v[i] = x2[i] - x1[i] - r[i];
            normR += r[i] * r[i];
            normV += v[i] * v[i];
        }
        normR = std::sqrt(normR);
        normV = std::sqrt(normV);
        double accel = (normV == 0.0) ? -1.0 : -(normR / normV);
        accel = std::max(SQUAREM_STEP_MIN, std::min(SQUAREM_STEP_MAX, accel));

        bool accepted = false;
        for (int attempt = 0; attempt < SQUAREM_MAX_BACKTRACK; ++attempt) {
            bool leftSimplex = false;
            for (size_t i = 0; i < speciesCnt; ++i) {
                xNew[i] = x0[i] - 2.0 * accel * r[i] + accel * accel * v[i];
                if (xNew[i] < 0.0) {
                    leftSimplex = true;
                    break;
                }
            }
            if (!leftSimplex) {
                normalize(xNew);
                if (meanLogLikelihood(cs, xNew, lengthFactor, ABUNDANCE_EXPONENT) >= llx1) {
                    accepted = true;
                    break;
                }
            }
            if (accel >= -1.0) {
                break;              // already at the plain double-EM step
            }
            accel = std::min(-1.0, (accel - 1.0) / 2.0);
            ++backtracks;
        }
        if (!accepted) {
            xNew = x2;
            accel = -1.0;
            ++fallbacks;
        }

        double change = 0.0;
        for (size_t i = 0; i < speciesCnt; ++i) {
            change = std::max(change, std::fabs(xNew[i] - x0[i]));
        }
        x0.swap(xNew);

        std::cout << "  iteration " << iterationsRun << ": step=" << accel
                  << " change=" << change << std::endl;
        if (change < par.tol && iterationsRun >= EM_MIN_ITER) {
            std::cout << "Converged after       : " << iterationsRun << " iterations" << std::endl;
            converged = true;
            break;
        }
    }
    std::cout << "EM steps              : " << (2 * iterationsRun) << " (2 per iteration)"
              << ", SQUAREM backtracks " << backtracks
              << ", double-EM fallbacks " << fallbacks << "/" << iterationsRun << std::endl;
    std::cout << "EM time               : " << double(time(nullptr) - start) << " s" << std::endl;
    if (!converged) {
        std::cout << "Warning: hit --max-iter (" << par.maxIter << ") without converging to --tol "
                  << par.tol << "." << std::endl;
    }
    return x0;
}

void printAbundanceDistribution(const std::vector<double> &theta) {
    if (theta.empty()) {
        return;
    }
    std::vector<double> values = theta;
    std::sort(values.begin(), values.end());
    auto quantile = [&values](double q) {
        return values[static_cast<size_t>(q * static_cast<double>(values.size() - 1))];
    };
    const std::ios_base::fmtflags flags = std::cout.flags();
    const std::streamsize precision = std::cout.precision();
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Abundance distribution: min=" << values.front()
              << " p50=" << quantile(0.50)
              << " p90=" << quantile(0.90)
              << " p99=" << quantile(0.99)
              << " max=" << values.back() << std::endl;
    std::cout.flags(flags);
    std::cout.precision(precision);
}

// ---------------------------------------------------------------------------
// Reclassification + output
// ---------------------------------------------------------------------------

// Turns a query's unnormalized per-species scores into a taxon: species are taken
// best-first until their normalized posteriors sum to cumulativeThreshold, and the
// LCA of those species is the answer. `selected` receives them best-first.
// Returns false when nothing can be assigned.
bool resolvePosterior(
    std::vector<std::pair<TaxID, double>> &scores,
    double denom,
    double cumulativeThreshold,
    TaxonomyWrapper *taxonomy,
    std::vector<TaxID> &selected,
    TaxID &classification,
    double &posteriorSum)
{
    classification = 0;
    posteriorSum = 0.0;
    selected.clear();
    if (denom <= 0.0 || scores.empty()) {
        return false;
    }

    for (auto &entryScore : scores) {
        entryScore.second /= denom;
    }
    std::sort(scores.begin(), scores.end(),
              [](const std::pair<TaxID, double> &a, const std::pair<TaxID, double> &b) {
                  if (a.second != b.second) {
                      return a.second > b.second;
                  }
                  return a.first < b.first; // deterministic tie-break
              });

    double sum = 0.0;
    for (size_t i = 0; i < scores.size() && sum < cumulativeThreshold; ++i) {
        sum += scores[i].second;
        selected.push_back(scores[i].first);
    }
    if (selected.empty()) {
        return false;
    }

    const TaxonNode *lca = taxonomy->LCA(selected);
    if (lca == nullptr) {
        selected.clear();
        return false;
    }
    classification = lca->taxId;
    posteriorSum = sum;
    return true;
}

void reclassifyRead(
    const CandidateDBEntry &entry,
    const std::unordered_map<TaxID, double> &theta,
    const std::unordered_map<TaxID, double> &lengthFactor,
    double minTheta,
    double cumulativeThreshold,
    double scorePower,
    TaxonomyWrapper *taxonomy,
    Query &query,
    std::vector<std::pair<TaxID, double>> &scratchScores,
    std::vector<TaxID> &scratchSpecies)
{
    query.name = entry.queryName;
    query.queryLength = static_cast<int>(entry.queryLength);
    query.queryLength2 = 0;
    query.classification = 0;
    query.idScore = 0.0f;
    query.eValue = -1.0f;
    query.taxCnt.clear();

    scratchScores.clear();
    double denom = 0.0;
    for (const SpeciesCandidate &candidate : entry.candidates) {
        auto thetaIt = theta.find(candidate.speciesId);
        if (thetaIt == theta.end() || thetaIt->second < minTheta) {
            continue;
        }
        auto lfIt = lengthFactor.find(candidate.speciesId);
        const double lf = (lfIt == lengthFactor.end()) ? 1.0 : lfIt->second;
        const double weight = std::pow(static_cast<double>(candidate.idScore), scorePower);
        const double score = thetaIt->second * weight * lf;
        if (score <= 0.0) {
            continue;
        }
        denom += score;
        scratchScores.emplace_back(candidate.speciesId, score);
    }

    TaxID classification = 0;
    double sum = 0.0;
    if (!resolvePosterior(scratchScores, denom, cumulativeThreshold, taxonomy,
                          scratchSpecies, classification, sum)) {
        return;
    }
    query.classification = classification;
    query.idScore = static_cast<float>(sum);

    for (const SpeciesCandidate &candidate : entry.candidates) {
        if (std::find(scratchSpecies.begin(), scratchSpecies.end(), candidate.speciesId)
                == scratchSpecies.end()) {
            continue;
        }
        if (candidate.speciesId == scratchSpecies[0]) {
            query.eValue = static_cast<float>(std::exp(static_cast<double>(candidate.logE)));
        }
        for (const auto &taxCount : candidate.taxCnt) {
            query.taxCnt[taxCount.first] += static_cast<int>(taxCount.second);
        }
    }
}

// Per-query result of the grouped reclassification.
struct QueryAssignments {
    std::vector<TaxID> classification;   // 0 = unclassified
    std::vector<float> score;            // summed posterior of the selected species
    std::vector<uint32_t> selectedCount; // selected species of query q live at
    std::vector<TaxID> selected;         //   selected[queryOffsets[q] .. + selectedCount[q]), best first
};

// Classifies every EM query (group or single read) from its pooled observations,
// the grouped counterpart of reclassifyRead. Every member read then inherits the
// result, which is the point of grouping: a read whose own candidates are
// ambiguous, or absent, is labelled by the reads it overlaps with.
void assignQueries(
    const CandidateSet &cs,
    const std::vector<double> &theta,
    const std::vector<double> &lengthFactor,
    double minTheta,
    double cumulativeThreshold,
    TaxonomyWrapper *taxonomy,
    QueryAssignments &qa)
{
    const size_t queryCnt = cs.queryCount();
    qa.classification.assign(queryCnt, 0);
    qa.score.assign(queryCnt, 0.0f);
    qa.selectedCount.assign(queryCnt, 0);
    qa.selected.assign(cs.obs.size(), 0);

#ifdef OPENMP
#pragma omp parallel shared(cs, theta, lengthFactor, minTheta, cumulativeThreshold, taxonomy, qa, queryCnt)
#endif
    {
        std::vector<std::pair<TaxID, double>> scores;
        std::vector<TaxID> selected;
#ifdef OPENMP
#pragma omp for schedule(dynamic, 1024)
#endif
        for (size_t q = 0; q < queryCnt; ++q) {
            const uint32_t begin = cs.queryOffsets[q];
            const uint32_t end = cs.queryOffsets[q + 1];
            if (begin == end) {
                continue;
            }
            scores.clear();
            double denom = 0.0;
            for (uint32_t j = begin; j < end; ++j) {
                const uint32_t s = cs.obs[j].spIdx;
                if (theta[s] < minTheta) {
                    continue;
                }
                const double score = theta[s] * static_cast<double>(cs.obs[j].weight) * lengthFactor[s];
                if (score <= 0.0) {
                    continue;
                }
                denom += score;
                scores.emplace_back(cs.speciesList[s], score);
            }
            TaxID classification = 0;
            double sum = 0.0;
            if (!resolvePosterior(scores, denom, cumulativeThreshold, taxonomy, selected, classification, sum)) {
                continue;
            }
            qa.classification[q] = classification;
            qa.score[q] = static_cast<float>(sum);
            qa.selectedCount[q] = static_cast<uint32_t>(selected.size());
            std::copy(selected.begin(), selected.end(), qa.selected.begin() + begin);
        }
    }
}

// Fills one read's output row from its query's assignment. The taxon and score
// are the query's; the E-value and match counts stay the read's own, restricted
// to the selected species, as in the per-read path. A read without a candidate
// for the top species keeps E-value -1.
void fillGroupedRead(
    const CandidateDBEntry &entry,
    uint32_t q,
    const CandidateSet &cs,
    const QueryAssignments &qa,
    Query &query)
{
    query.name = entry.queryName;
    query.queryLength = static_cast<int>(entry.queryLength);
    query.queryLength2 = 0;
    query.classification = 0;
    query.idScore = 0.0f;
    query.eValue = -1.0f;
    query.taxCnt.clear();

    const TaxID classification = qa.classification[q];
    if (classification == 0) {
        return;
    }
    query.classification = classification;
    query.idScore = qa.score[q];

    const TaxID *selectedBegin = qa.selected.data() + cs.queryOffsets[q];
    const TaxID *selectedEnd = selectedBegin + qa.selectedCount[q];
    for (const SpeciesCandidate &candidate : entry.candidates) {
        if (std::find(selectedBegin, selectedEnd, candidate.speciesId) == selectedEnd) {
            continue;
        }
        if (candidate.speciesId == *selectedBegin) {
            query.eValue = static_cast<float>(std::exp(static_cast<double>(candidate.logE)));
        }
        for (const auto &taxCount : candidate.taxCnt) {
            query.taxCnt[taxCount.first] += static_cast<int>(taxCount.second);
        }
    }
}

// One row per read group: what the EM saw as the group's query (its pooled
// candidate scores) and what it decided.
void writeGroupTable(
    const std::string &fileName,
    const CandidateSet &cs,
    const GroupedQueries &gq,
    const QueryAssignments &qa,
    TaxonomyWrapper *taxonomy)
{
    std::ofstream out(fileName);
    if (!out.is_open()) {
        std::cerr << "Error: could not open " << fileName << " for writing." << std::endl;
        return;
    }
    out << "#group_id\tn_reads\tn_reads_with_candidates\ttaxID\tscore\trank\tspeciesId:mean_idScore\n";

    std::vector<std::pair<float, TaxID>> ranked;
    for (uint32_t g = 0; g < gq.groupCount; ++g) {
        out << gq.groupIds[g] << '\t' << gq.queryReads[g] << '\t' << gq.groupReadsWithCand[g] << '\t';
        const TaxID classification = qa.classification[g];
        if (classification != 0) {
            const TaxonNode *node = taxonomy->taxonNode(classification);
            out << taxonomy->getOriginalTaxID(classification) << '\t' << qa.score[g] << '\t'
                << taxonomy->getString(node->rankIdx);
        } else {
            out << "0\t0\t-";
        }
        out << '\t';

        ranked.clear();
        for (uint32_t j = cs.queryOffsets[g]; j < cs.queryOffsets[g + 1]; ++j) {
            ranked.emplace_back(gq.obsMeanScore[j], cs.speciesList[cs.obs[j].spIdx]);
        }
        std::sort(ranked.begin(), ranked.end(),
                  [](const std::pair<float, TaxID> &a, const std::pair<float, TaxID> &b) {
                      if (a.first != b.first) {
                          return a.first > b.first;
                      }
                      return a.second < b.second;
                  });
        // Formatted by hand so the fixed precision does not leak into the stream
        // state used by the other columns.
        char token[64];
        for (const auto &species : ranked) {
            std::snprintf(token, sizeof(token), "%d:%.4f ",
                          static_cast<int>(taxonomy->getOriginalTaxID(species.second)),
                          static_cast<double>(species.first));
            out << token;
        }
        out << '\n';
    }
    out.close();
}

} // namespace

// em
//
//   metabuli em <candidate DB> <DB dir> <output directory> <job ID> [options]
//
// Outputs (<out> = <output directory>/<job ID>, matching classify-candidates):
//   <out>_classifications.tsv  per-read reclassification, same columns as classify-candidates
//   <out>_report.tsv           taxonomy report over those reclassified reads
//   <out>_EM_report.tsv        taxonomy report over the raw EM abundances
//   <out>_abundance.tsv        per-species EM abundance and length factor
//   <out>_groups.tsv           (--read-groupings 1) per-group pooled scores and assignment
//
// The likelihood weight of assigning query q to species s is
//   score(q,s)^p  *  theta(s)  *  lengthFactor(s)
// where a query is a read, score its idScore, and lengthFactor(s) = 1 / the size of
// the species' representative genome (species2genomeSize.tsv), so that
// score^p * lengthFactor plays P(read | s). With --read-groupings 1 a query
// is a read group and score(q,s) is the mean idScore of its reads for s (reads
// without a candidate for s contribute 0); reads in no group stay single queries.
// The M-step adds a uniform pseudocount:
//   theta~(s) = sum_q mass(q) * posterior(q,s) + phi * R / S
// with mass(q) = 1 for a read and the read count for a group (--group-weight 1),
// R the summed mass of the queries with candidates.
//
// The groups live in <output directory>/grouping (groupMap, groups): em builds them there
// by grouping the reads of --query-file against --common-kmer-db, or, without those two
// flags, reuses what is already there (an earlier run, or a `metabuli grouping` result
// copied in).
int emAlgorithm(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    // The built-in grouping shares the `grouping` command's defaults; flags still win.
    setGroupGenerationDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string candidateDb = par.filenames[0];
    const std::string dbDir = par.filenames[1];
    const std::string outDir = par.filenames[2];
    const std::string jobId = par.filenames[3];

    // --- Validate inputs ---
    if (!FileUtil::fileExists(candidateDb.c_str())) {
        std::cerr << "Error: candidate DB " << candidateDb << " is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::fileExists((candidateDb + ".index").c_str())) {
        std::cerr << "Error: candidate DB index " << candidateDb << ".index is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        std::cerr << "Error: database directory " << dbDir << " is not found." << std::endl;
        return 1;
    }

    // --- Validate the EM knobs before doing any work ---
    if (par.lcaPosterior <= 0.0 || par.lcaPosterior > 1.0) {
        std::cerr << "Error: --lca-posterior must be in (0.0, 1.0]." << std::endl;
        return 1;
    }
    if (par.minReads < 0) {
        std::cerr << "Error: --min-reads must be >= 0." << std::endl;
        return 1;
    }

    // --- Read grouping inputs ---
    // The groups live in <output directory>/grouping:
    //   --query-file <R1[,R2]> --common-kmer-db <db>   group the reads here, (re)writing that directory
    //   neither                                        reuse the groupMap/groups already there
    const bool useGroups = par.readGroupings != 0;
    const std::string groupingDir = outDir + "/grouping";
    const std::string groupMapFileName = groupingDir + "/groupMap";
    const std::string groupsFileName = groupingDir + "/groups";
    std::vector<std::string> groupingReadFiles;
    bool groupHere = false;
    if (useGroups) {
        if (!par.queryFile.empty() || !par.commonKmerDb.empty()) {
            if (par.queryFile.empty() || par.commonKmerDb.empty()) {
                std::cerr << "Error: grouping the reads needs both --query-file <R1[,R2]> and --common-kmer-db <dir>."
                          << std::endl;
                return 1;
            }
            groupingReadFiles = Util::split(par.queryFile, ",");
            if (groupingReadFiles.empty() || groupingReadFiles.size() > 2) {
                std::cerr << "Error: --query-file takes one read file, or two comma-separated files for paired-end."
                          << std::endl;
                return 1;
            }
            for (const std::string &readFile : groupingReadFiles) {
                if (!FileUtil::fileExists(readFile.c_str())) {
                    std::cerr << "Error: read file " << readFile << " (from --query-file) is not found." << std::endl;
                    return 1;
                }
            }
            // Two files = paired-end, one file = single-end.
            par.seqMode = (groupingReadFiles.size() == 2) ? 2 : 1;
            std::string error;
            if (!prepareGroupGeneration(par, par.commonKmerDb, error, true)) {
                std::cerr << "Error: " << error << std::endl;
                return 1;
            }
            groupHere = true;
        } else if (!FileUtil::fileExists(groupMapFileName.c_str())) {
            std::cerr << "Error: --read-groupings 1 needs the read groups in " << groupingDir << ". Either" << std::endl;
            std::cerr << "         pass --query-file <R1[,R2]> --common-kmer-db <dir> to group the reads now, or" << std::endl;
            std::cerr << "         put the groupMap and groups of a `metabuli grouping` run there." << std::endl;
            return 1;
        }
    }

    if (!FileUtil::directoryExists(outDir.c_str())) {
        FileUtil::makeDir(outDir.c_str());
    }

    const int threadCount = par.threads <= 0 ? 1 : par.threads;
#ifdef OPENMP
    omp_set_num_threads(threadCount);
#endif

    if (groupHere) {
        std::cout << "Grouping reads        : " << par.queryFile << std::endl;
        std::cout << "Common k-mer DB       : " << par.commonKmerDb << std::endl;
        std::cout << "Grouping output       : " << groupingDir << std::endl;
        std::cout << "Grouping k-mers       : syncmer " << par.syncmer << ", s-mer length " << par.smerLen
                  << " (from the DB's kmer_params when present)" << std::endl;
        std::cout << "==================== read grouping ====================" << std::endl;
        const time_t groupingStart = time(nullptr);
        if (!runReadGrouping(par, groupingReadFiles, par.commonKmerDb, groupingDir)) {
            std::cerr << "Error: read grouping did not produce " << groupMapFileName << "." << std::endl;
            return 1;
        }
        std::cout << "==================== read grouping done (" << double(time(nullptr) - groupingStart)
                  << " s) ====================" << std::endl;
    } else if (useGroups) {
        std::cout << "Read groups reused    : " << groupingDir
                  << " (pass --query-file and --common-kmer-db to group again)" << std::endl;
    }

    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, "");
    if (taxonomy == nullptr) {
        std::cerr << "Error: could not load the taxonomy from " << dbDir << "." << std::endl;
        return 1;
    }

    CandidateDBReader reader(candidateDb, threadCount);
    if (!reader.open(DBReader<unsigned int>::SORT_BY_ID)) {
        std::cerr << "Error: failed to open candidate DB " << candidateDb << std::endl;
        delete taxonomy;
        return 1;
    }
    const size_t entryCount = reader.size();

    std::cout << "Candidate DB          : " << candidateDb << std::endl;
    std::cout << "Candidate DB entries  : " << entryCount << std::endl;
    std::cout << "Score power (p)       : " << par.scorePower << std::endl;
    std::cout << "Uniform prior (phi)   : " << par.priorMass << std::endl;
    std::cout << "Abundance exponent    : " << ABUNDANCE_EXPONENT << std::endl;
    std::cout << "SQUAREM               : on (step clamped to ["
              << SQUAREM_STEP_MIN << ", " << SQUAREM_STEP_MAX << "])" << std::endl;
    std::cout << "Read grouping         : " << (useGroups ? "on" : "off") << std::endl;
    if (useGroups) {
        std::cout << "Group weight          : "
                  << (par.groupWeight != 0 ? "read count of the group" : "1 per group") << std::endl;
    }
    std::cout << "--------------------" << std::endl;

    // --- Pass 1: candidate DB (+ read groups) -> EM observations ---
    CandidateSet cs;
    GroupedQueries gq;
    if (useGroups) {
        ReadGrouping grouping;
        if (!loadGroupMap(groupMapFileName, entryCount, grouping)
            || !checkGroupsFile(groupsFileName, entryCount, grouping)) {
            reader.close();
            delete taxonomy;
            return 1;
        }
        std::cout << "Read groups           : " << grouping.groupCount << std::endl;
        std::cout << "Reads in a group      : " << grouping.groupedReads << "/" << entryCount << std::endl;
        if (!loadGroupedObservations(reader, entryCount, par, grouping, cs, gq)) {
            std::cout << "No species candidates found in " << candidateDb << ". Nothing to do." << std::endl;
            reader.close();
            delete taxonomy;
            return 1;
        }
        std::cout << "Group size            : mean "
                  << (gq.groupCount > 0
                          ? static_cast<double>(grouping.groupedReads) / static_cast<double>(gq.groupCount)
                          : 0.0)
                  << ", max " << gq.maxGroupReads << std::endl;
        std::cout << "Reads with candidates : " << cs.readsWithCandidates << "/" << entryCount << std::endl;
        std::cout << "EM queries            : " << cs.queryCount() << " (" << gq.groupCount << " groups + "
                  << (cs.queryCount() - gq.groupCount) << " single reads), "
                  << cs.queriesWithCandidates << " with candidates" << std::endl;
        // Always a whole number: every query mass is a read count (or 1).
        std::cout << "Read mass (R)         : " << static_cast<unsigned long long>(cs.totalMass) << std::endl;
    } else {
        if (!loadObservations(reader, entryCount, par, cs)) {
            std::cout << "No species candidates found in " << candidateDb << ". Nothing to do." << std::endl;
            reader.close();
            delete taxonomy;
            return 1;
        }
        std::cout << "Reads with candidates : " << cs.readsWithCandidates << "/" << entryCount << std::endl;
    }
    std::cout << "Distinct species      : " << cs.speciesList.size() << std::endl;
    std::cout << "Observations          : " << cs.obs.size() << std::endl;

    // --- Length factor ---
    const std::vector<double> lengthFactorByIdx = buildLengthFactors(dbDir, cs);
    // --- EM ---
    std::cout << "--------------------" << std::endl;
    std::cout << "Running EM ..." << std::endl;
    const std::vector<double> theta = runEm(cs, lengthFactorByIdx, par);
    printAbundanceDistribution(theta);

    // --- Reclassification threshold: --min-reads expected reads ---
    double minTheta = 0.0;
    if (par.minReads > 0 && cs.totalMass > 0.0) {
        minTheta = static_cast<double>(par.minReads) / cs.totalMass;
        std::cout << "Min. expected reads   : " << par.minReads
                  << " => theta >= " << minTheta << std::endl;
    }

    // Species-keyed views for the per-read pass.
    std::unordered_map<TaxID, double> thetaBySp;
    std::unordered_map<TaxID, double> lengthFactorBySp;
    thetaBySp.reserve(cs.speciesList.size());
    lengthFactorBySp.reserve(cs.speciesList.size());
    for (size_t i = 0; i < cs.speciesList.size(); ++i) {
        thetaBySp[cs.speciesList[i]] = theta[i];
        lengthFactorBySp[cs.speciesList[i]] = lengthFactorByIdx[i];
    }

    // --- Output wiring ---
    const std::string outPrefix = outDir + "/" + jobId;
    const std::string classificationFileName = outPrefix + "_classifications.tsv";
    const std::string reportFileName = outPrefix + "_report.tsv";
    const std::string emReportFileName = outPrefix + "_EM_report.tsv";
    const std::string abundanceFileName = outPrefix + "_abundance.tsv";
    const std::string groupFileName = outPrefix + "_groups.tsv";

    Reporter reporter(par, taxonomy, reportFileName);
    reporter.setReadClassificationFileName(classificationFileName);
    // Only ReportType::EM is emitted here; the reclassified-read report is the
    // ReportType::Default one above, so both EM slots can share a name.
    reporter.setEmReportFileNames(emReportFileName, emReportFileName);

    // --- Pass 2 (grouped): classify every EM query once ---
    QueryAssignments qa;
    if (useGroups) {
        std::cout << "--------------------" << std::endl;
        std::cout << "Classifying read groups ..." << std::endl;
        assignQueries(cs, theta, lengthFactorByIdx, minTheta, par.lcaPosterior, taxonomy, qa);
        size_t classifiedGroups = 0;
        for (uint32_t g = 0; g < gq.groupCount; ++g) {
            classifiedGroups += (qa.classification[g] != 0);
        }
        std::cout << "Groups classified     : " << classifiedGroups << "/" << gq.groupCount << std::endl;
        writeGroupTable(groupFileName, cs, gq, qa, taxonomy);
    }

    // --- Pass 3: reclassify every read and write the classification file ---
    std::cout << "--------------------" << std::endl;
    std::cout << "Reclassifying reads ..." << std::endl;
    reporter.openReadClassificationFile();

    std::unordered_map<TaxID, unsigned int> reclassifiedTaxCounts;
    std::unordered_map<TaxID, double> cladeScoreSums;
    size_t reclassifiedCnt = 0;
    size_t unreadableCnt = 0;

    std::vector<Query> chunk;
    for (size_t chunkStart = 0; chunkStart < entryCount; chunkStart += OUTPUT_CHUNK) {
        const size_t chunkEnd = std::min(entryCount, chunkStart + OUTPUT_CHUNK);
        chunk.assign(chunkEnd - chunkStart, Query());

#ifdef OPENMP
#pragma omp parallel shared(reader, chunk, chunkStart, chunkEnd, thetaBySp, lengthFactorBySp, \
        minTheta, taxonomy, par, useGroups, cs, gq, qa)
#endif
        {
            int threadIdx = 0;
#ifdef OPENMP
            threadIdx = omp_get_thread_num();
#endif
            CandidateDBEntry entry;
            std::vector<std::pair<TaxID, double>> scratchScores;
            std::vector<TaxID> scratchSpecies;
#ifdef OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
            for (size_t i = chunkStart; i < chunkEnd; ++i) {
                if (!reader.getByIndex(i, entry, threadIdx)) {
                    continue;
                }
                if (useGroups) {
                    fillGroupedRead(entry, gq.readQuery[i], cs, qa, chunk[i - chunkStart]);
                } else {
                    reclassifyRead(entry, thetaBySp, lengthFactorBySp, minTheta,
                                   par.lcaPosterior, par.scorePower, taxonomy,
                                   chunk[i - chunkStart], scratchScores, scratchSpecies);
                }
            }
        }

        // writeReadClassification stops at the first entry with an empty name, so
        // drop unreadable entries instead of letting one truncate the chunk.
        const size_t beforeCompaction = chunk.size();
        chunk.erase(std::remove_if(chunk.begin(), chunk.end(),
                                   [](const Query &query) { return query.name.empty(); }),
                    chunk.end());
        unreadableCnt += beforeCompaction - chunk.size();

        for (const Query &query : chunk) {
            if (query.classification != 0) {
                reclassifiedTaxCounts[query.classification] += 1;
                addScoreToClade(taxonomy, query.classification, query.idScore, cladeScoreSums);
                ++reclassifiedCnt;
            }
        }
        reporter.writeReadClassification(chunk);
    }
    reporter.closeReadClassificationFile();
    std::cout << "Reads reclassified    : " << reclassifiedCnt << "/" << entryCount << std::endl;
    if (unreadableCnt > 0) {
        std::cout << "Warning: " << unreadableCnt
                  << " candidate DB entries could not be read and are absent from the output." << std::endl;
    }

    // --- Report over the reclassified reads (classify-candidates format) ---
    const std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    std::unordered_map<TaxID, TaxonCounts> reclassifiedCladeCounts =
        taxonomy->getCladeCounts(reclassifiedTaxCounts, parentToChildren);
    std::unordered_map<TaxID, double> avgScores;
    addAverageScoresToMap(cladeScoreSums, reclassifiedCladeCounts, avgScores);
    reporter.writeReportFile(static_cast<int>(entryCount), reclassifiedCladeCounts,
                             avgScores, ReportType::Default);

    // --- Report over the raw EM abundances ---
    // theta is a distribution over the read mass that had at least one candidate,
    // so the remainder is reported as unclassified.
    std::unordered_map<TaxID, unsigned int> emTaxCounts;
    unsigned long long explained = 0;
    for (size_t i = 0; i < cs.speciesList.size(); ++i) {
        const unsigned int count = static_cast<unsigned int>(theta[i] * cs.totalMass);
        if (count == 0) {
            continue;
        }
        emTaxCounts[cs.speciesList[i]] = count;
        explained += count;
    }
    if (static_cast<unsigned long long>(entryCount) > explained) {
        emTaxCounts[0] = static_cast<unsigned int>(entryCount - explained);
    }
    std::unordered_map<TaxID, TaxonCounts> emCladeCounts =
        taxonomy->getCladeCounts(emTaxCounts, parentToChildren);
    reporter.writeReportFile(static_cast<int>(entryCount), emCladeCounts, ReportType::EM);

    // --- Per-species abundance table ---
    {
        std::vector<size_t> order(cs.speciesList.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&theta, &cs](size_t a, size_t b) {
            if (theta[a] != theta[b]) {
                return theta[a] > theta[b];
            }
            return cs.speciesList[a] < cs.speciesList[b];
        });

        std::ofstream abundanceFile(abundanceFileName);
        if (!abundanceFile.is_open()) {
            std::cerr << "Error: could not open " << abundanceFileName << " for writing." << std::endl;
        } else {
            abundanceFile << "#speciesId\tabundance\testimated_reads\tlength_factor\n";
            for (const size_t i : order) {
                abundanceFile << taxonomy->getOriginalTaxID(cs.speciesList[i]) << '\t'
                              << theta[i] << '\t'
                              << (theta[i] * cs.totalMass) << '\t'
                              << lengthFactorByIdx[i] << '\n';
            }
            abundanceFile.close();
        }

        std::cout << "--------------------" << std::endl;
        std::cout << "speciesId\tabundance\testimated_reads" << std::endl;
        for (size_t r = 0; r < order.size() && r < MAX_REPORT_ROWS; ++r) {
            const size_t i = order[r];
            std::cout << taxonomy->getOriginalTaxID(cs.speciesList[i]) << '\t'
                      << theta[i] << '\t'
                      << (theta[i] * cs.totalMass) << std::endl;
        }
        if (order.size() > MAX_REPORT_ROWS) {
            std::cout << "... (" << (order.size() - MAX_REPORT_ROWS) << " more species)" << std::endl;
        }
    }

    reader.close();
    delete taxonomy;

    std::cout << "--------------------" << std::endl;
    std::cout << "Reclassification      : " << classificationFileName << std::endl;
    std::cout << "Report                : " << reportFileName << std::endl;
    std::cout << "EM abundance report   : " << emReportFileName << std::endl;
    std::cout << "EM abundance table    : " << abundanceFileName << std::endl;
    if (useGroups) {
        std::cout << "Read group table      : " << groupFileName << std::endl;
    }
    std::cout << "EM completed." << std::endl;
    return 0;
}
