#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "CandidateDBReader.h"
#include "Reporter.h"
#include "TaxonomyWrapper.h"
#include "DBReader.h"

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

// One (read, candidate species) pair: the unit the E-step assigns responsibility to.
struct Observation {
    uint32_t spIdx;   // index into CandidateSet::speciesList
    float weight;     // idScore ^ --score-power
};

struct CandidateSet {
    // Observations grouped by read: read r owns [readOffsets[r], readOffsets[r+1]).
    std::vector<Observation> obs;
    std::vector<uint32_t> readOffsets;
    std::vector<TaxID> speciesList;                 // compact index -> species taxId
    std::unordered_map<TaxID, uint32_t> sp2idx;
    size_t readsWithCandidates = 0;
};

bool loadObservations(
    CandidateDBReader &reader,
    size_t entryCount,
    const LocalParameters &par,
    CandidateSet &cs)
{
    const double scorePower = par.scorePower;
    cs.readOffsets.reserve(entryCount + 1);
    cs.readOffsets.push_back(0);

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
        // readOffsets indexes obs with 32 bits; refuse rather than silently wrap.
        if (cs.obs.size() > static_cast<size_t>(UINT32_MAX)) {
            std::cerr << "Error: more than " << UINT32_MAX << " (read, candidate) pairs in the candidate DB."
                      << std::endl;
            return false;
        }
        cs.readOffsets.push_back(static_cast<uint32_t>(cs.obs.size()));
    }
    return !cs.speciesList.empty();
}

// ---------------------------------------------------------------------------
// Per-species length normalization (unique k-mer count)
// ---------------------------------------------------------------------------

bool loadUniqKmerCounts(
    const std::string &dbDir,
    TaxonomyWrapper *taxonomy,
    std::vector<uint32_t> &sp2uniqKmerCnt)
{
    const std::string cacheFileName = dbDir + "/sp2uniqKmerCnt";
    if (FileUtil::fileExists(cacheFileName.c_str())) {
        std::ifstream cacheFile(cacheFileName);
        if (!cacheFile.is_open()) {
            std::cout << "Warning: could not open " << cacheFileName << "." << std::endl;
            return false;
        }
        std::cout << "Unique k-mer counts   : " << cacheFileName << std::endl;
        TaxID taxId = 0;
        uint32_t count = 0;
        while (cacheFile >> taxId >> count) {
            if (taxId >= 0 && static_cast<size_t>(taxId) < sp2uniqKmerCnt.size()) {
                sp2uniqKmerCnt[taxId] = count;
            }
        }
        return true;
    }

    const std::string infoFileName = dbDir + "/info";
    if (!FileUtil::fileExists(infoFileName.c_str())) {
        std::cout << "Warning: neither " << cacheFileName << " nor " << infoFileName
                  << " is available." << std::endl;
        return false;
    }

    const size_t taxIdCount = static_cast<size_t>(taxonomy->getMaxTaxID()) + 1;
    constexpr TaxID SPECIES_UNRESOLVED = -1;
    std::vector<TaxID> taxId2speciesId(taxIdCount, SPECIES_UNRESOLVED);

    const size_t totalKmers = FileUtil::getFileSize(infoFileName) / sizeof(TaxID);
    std::cout << "Unique k-mer counts   : counting " << totalKmers << " k-mers in "
              << infoFileName << " ..." << std::endl;

    ReadBuffer<TaxID> readBuffer(infoFileName, 16 * 1024 * 1024);
    const time_t countStart = time(nullptr);
    const size_t reportStep = (totalKmers > 0) ? (totalKmers / 10 + 1) : 0;
    size_t nextReport = (reportStep > 0) ? reportStep : SIZE_MAX;
    size_t processed = 0;
    TaxID taxId = 0;
    while ((taxId = readBuffer.getNext()) != 0) {
        if (taxId > 0 && static_cast<size_t>(taxId) < taxIdCount) {
            TaxID speciesTaxId = taxId2speciesId[taxId];
            if (speciesTaxId == SPECIES_UNRESOLVED) {
                speciesTaxId = taxonomy->getTaxIdAtRank(taxId, "species");
                taxId2speciesId[taxId] = speciesTaxId;
            }
            if (speciesTaxId > 0 && static_cast<size_t>(speciesTaxId) < sp2uniqKmerCnt.size()) {
                ++sp2uniqKmerCnt[speciesTaxId];
            }
        }
        if (++processed >= nextReport) {
            std::cout << "  " << (100 * processed / totalKmers) << "% ("
                      << double(time(nullptr) - countStart) << " s)" << std::endl;
            nextReport += reportStep;
        }
    }
    std::cout << "  done in " << double(time(nullptr) - countStart) << " s" << std::endl;

    std::ofstream cacheOut(cacheFileName);
    if (!cacheOut.is_open()) {
        std::cout << "Warning: could not write " << cacheFileName << "; the next run will have to count "
                  << infoFileName << " again." << std::endl;
        return true;
    }
    for (size_t i = 0; i < sp2uniqKmerCnt.size(); ++i) {
        if (sp2uniqKmerCnt[i] > 0) {
            cacheOut << i << " " << sp2uniqKmerCnt[i] << "\n";
        }
    }
    cacheOut.close();
    std::cout << "Cached to             : " << cacheFileName
              << " (later runs skip this scan)" << std::endl;
    return true;
}

// lengthFactor divides the likelihood weight so that a species owning more
// unique k-mers does not attract reads on size alone. Fixed at 1/uniqKmerCnt;
// needs <DB dir>/sp2uniqKmerCnt or <DB dir>/info.
std::vector<double> buildLengthFactors(
    const std::string &dbDir,
    TaxonomyWrapper *taxonomy,
    const CandidateSet &cs)
{
    std::vector<double> lengthFactor(cs.speciesList.size(), 1.0);

    std::vector<uint32_t> sp2uniqKmerCnt(static_cast<size_t>(taxonomy->getMaxTaxID()) + 1, 0);
    if (!loadUniqKmerCounts(dbDir, taxonomy, sp2uniqKmerCnt)) {
        std::cout << "Warning: no unique k-mer counts are available; "
                  << "using a length factor of 1.0 for every species." << std::endl;
        return lengthFactor;
    }

    size_t missing = 0;
    for (size_t i = 0; i < cs.speciesList.size(); ++i) {
        const TaxID sp = cs.speciesList[i];
        const uint32_t count = (sp >= 0 && static_cast<size_t>(sp) < sp2uniqKmerCnt.size())
            ? sp2uniqKmerCnt[sp]
            : 0;
        if (count == 0) {
            // Unknown k-mer count: a 0.0 factor would silently delete the species,
            // so leave it unnormalized instead.
            ++missing;
            lengthFactor[i] = 1.0;
            continue;
        }
        lengthFactor[i] = 1.0 / static_cast<double>(count);
    }
    std::cout << "Length factor         : 1/uniqKmerCnt"
              << " (" << missing << "/" << cs.speciesList.size() << " species without a count)"
              << std::endl;
    return lengthFactor;
}


// ---------------------------------------------------------------------------
// EM
// ---------------------------------------------------------------------------

constexpr double LL_FLOOR = 1e-300;

double meanLogLikelihood(
    const CandidateSet &cs,
    const std::vector<double> &theta,
    const std::vector<double> &lengthFactor,
    double alpha)
{
    const size_t readCnt = cs.readOffsets.size() - 1;
    double total = 0.0;
#ifdef OPENMP
#pragma omp parallel for schedule(static) reduction(+:total) \
        shared(cs, theta, lengthFactor, alpha, readCnt)
#endif
    for (size_t r = 0; r < readCnt; ++r) {
        const uint32_t begin = cs.readOffsets[r];
        const uint32_t end = cs.readOffsets[r + 1];
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
        total += std::log(std::max(denom, LL_FLOOR));
    }
    return (cs.readsWithCandidates > 0)
        ? total / static_cast<double>(cs.readsWithCandidates)
        : 0.0;
}

// One EM update. Returns the mean log-likelihood of the INPUT theta, which the
// per-read denominators give away for free.
double emStep(
    const CandidateSet &cs,
    const std::vector<double> &theta,
    const std::vector<double> &lengthFactor,
    const std::vector<double> &priorTerm,
    double alpha,
    std::vector<double> &out)
{
    out = priorTerm;
    const size_t readCnt = cs.readOffsets.size() - 1;
    double llTotal = 0.0;

#ifdef OPENMP
#pragma omp parallel shared(cs, theta, lengthFactor, out, alpha, readCnt, llTotal)
#endif
    {
        std::vector<double> local(out.size(), 0.0);
        std::vector<double> weights;
        double llLocal = 0.0;
#ifdef OPENMP
#pragma omp for schedule(static)
#endif
        for (size_t r = 0; r < readCnt; ++r) {
            const uint32_t begin = cs.readOffsets[r];
            const uint32_t end = cs.readOffsets[r + 1];
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
            llLocal += std::log(std::max(denom, LL_FLOOR));
            if (denom <= 0.0) {
                continue;
            }
            for (uint32_t j = begin; j < end; ++j) {
                local[cs.obs[j].spIdx] += weights[j - begin] / denom;
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
    return (cs.readsWithCandidates > 0)
        ? llTotal / static_cast<double>(cs.readsWithCandidates)
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
            ? std::max(0.0, par.priorMass) * static_cast<double>(cs.readsWithCandidates)
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
    if (denom <= 0.0 || scratchScores.empty()) {
        return;
    }

    for (auto &entryScore : scratchScores) {
        entryScore.second /= denom;
    }
    std::sort(scratchScores.begin(), scratchScores.end(),
              [](const std::pair<TaxID, double> &a, const std::pair<TaxID, double> &b) {
                  if (a.second != b.second) {
                      return a.second > b.second;
                  }
                  return a.first < b.first; // deterministic tie-break
              });

    scratchSpecies.clear();
    double sum = 0.0;
    for (size_t i = 0; i < scratchScores.size() && sum < cumulativeThreshold; ++i) {
        sum += scratchScores[i].second;
        scratchSpecies.push_back(scratchScores[i].first);
    }
    if (scratchSpecies.empty()) {
        return;
    }

    const TaxonNode *lca = taxonomy->LCA(scratchSpecies);
    if (lca == nullptr) {
        return;
    }
    query.classification = lca->taxId;
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
//
// The likelihood weight of assigning read q to species s is
//   idScore(q,s)^p  *  theta(s)  *  lengthFactor(s)
// and the M-step adds a uniform pseudocount:
//   theta~(s) = sum_q posterior(q,s) + phi * R / S
int emAlgorithm(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
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

    if (!FileUtil::directoryExists(outDir.c_str())) {
        FileUtil::makeDir(outDir.c_str());
    }

    const int threadCount = par.threads <= 0 ? 1 : par.threads;
#ifdef OPENMP
    omp_set_num_threads(threadCount);
#endif

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
    std::cout << "--------------------" << std::endl;

    // --- Pass 1: candidate DB -> EM observations ---
    CandidateSet cs;
    if (!loadObservations(reader, entryCount, par, cs)) {
        std::cout << "No species candidates found in " << candidateDb << ". Nothing to do." << std::endl;
        reader.close();
        delete taxonomy;
        return 1;
    }
    std::cout << "Reads with candidates : " << cs.readsWithCandidates << "/" << entryCount << std::endl;
    std::cout << "Distinct species      : " << cs.speciesList.size() << std::endl;
    std::cout << "Observations          : " << cs.obs.size() << std::endl;

    // --- Length factor ---
    const std::vector<double> lengthFactorByIdx = buildLengthFactors(dbDir, taxonomy, cs);
    // --- EM ---
    std::cout << "--------------------" << std::endl;
    std::cout << "Running EM ..." << std::endl;
    const std::vector<double> theta = runEm(cs, lengthFactorByIdx, par);
    printAbundanceDistribution(theta);

    // --- Reclassification threshold: --min-reads expected reads ---
    double minTheta = 0.0;
    if (par.minReads > 0 && cs.readsWithCandidates > 0) {
        minTheta = static_cast<double>(par.minReads)
            / static_cast<double>(cs.readsWithCandidates);
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

    Reporter reporter(par, taxonomy, reportFileName);
    reporter.setReadClassificationFileName(classificationFileName);
    // Only ReportType::EM is emitted here; the reclassified-read report is the
    // ReportType::Default one above, so both EM slots can share a name.
    reporter.setEmReportFileNames(emReportFileName, emReportFileName);

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
        minTheta, taxonomy, par)
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
                reclassifyRead(entry, thetaBySp, lengthFactorBySp, minTheta,
                               par.lcaPosterior, par.scorePower, taxonomy,
                               chunk[i - chunkStart], scratchScores, scratchSpecies);
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
    // theta is a distribution over the reads that had at least one candidate, so
    // the remainder is reported as unclassified.
    std::unordered_map<TaxID, unsigned int> emTaxCounts;
    unsigned long long explained = 0;
    for (size_t i = 0; i < cs.speciesList.size(); ++i) {
        const unsigned int count = static_cast<unsigned int>(
            theta[i] * static_cast<double>(cs.readsWithCandidates));
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
                              << (theta[i] * static_cast<double>(cs.readsWithCandidates)) << '\t'
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
                      << (theta[i] * static_cast<double>(cs.readsWithCandidates)) << std::endl;
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
    std::cout << "EM completed." << std::endl;
    return 0;
}
