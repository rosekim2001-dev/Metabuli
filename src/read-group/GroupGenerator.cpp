#include "GroupGenerator.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"
#include "Kmer.h"
#include <sys/stat.h>
#include <queue>
#include <iomanip>
#include <functional>

// One stream head in mergeGraph's k-way merge: the relation plus the subGraph_* index it
// came from, so the merge knows which stream to refill after popping.
struct MergeHeapEntry {
    Relation rel;
    size_t stream;
};

// Min-heap ordering for MergeHeapEntry. std::priority_queue pops the "largest" element, so
// this comparator answers "does a come after b". Relation::operator< compares (id1,id2)
// only; equal keys tie-break on stream index to keep the pop order deterministic.
struct MergeHeapGreater {
    bool operator()(const MergeHeapEntry & a, const MergeHeapEntry & b) const {
        if (a.rel < b.rel) { return false; }
        if (b.rel < a.rel) { return true; }
        return a.stream > b.stream;
    }
};

using MergeHeap = std::priority_queue<MergeHeapEntry, std::vector<MergeHeapEntry>, MergeHeapGreater>;

// Human-readable byte size for disk-usage logs.
static std::string humanBytes(uint64_t bytes) {
    const char* unit[] = {"B", "KB", "MB", "GB", "TB"};
    double v = static_cast<double>(bytes);
    int i = 0;
    while (v >= 1024.0 && i < 4) { v /= 1024.0; ++i; }
    char buf[64];
    snprintf(buf, sizeof(buf), "%.2f %s", v, unit[i]);
    return std::string(buf);
}

// Min, median and max of `values`, by partial_sort so the whole vector is not ordered: the
// caller's vectors run to hundreds of thousands of entries and only three order statistics are
// wanted. Takes a copy because it reorders. Returns {0,0,0} on an empty input.
struct OrderStats { uint64_t min; uint64_t median; uint64_t max; };
static OrderStats orderStats(std::vector<uint64_t> values) {
    if (values.empty()) {
        return {0, 0, 0};
    }
    const size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    const uint64_t median = values[mid];
    const auto minmax = std::minmax_element(values.begin(), values.end());
    return {*minmax.first, median, *minmax.second};
}

// Sum sizes of existing files in `paths`, log "<label>: N files, X freed", then remove them.
// Non-existent paths are skipped (not counted).
static void reportAndRemoveFiles(const std::vector<std::string>& paths, const std::string& label) {
    uint64_t total = 0;
    size_t count = 0;
    for (const std::string& p : paths) {
        struct stat st;
        if (stat(p.c_str(), &st) == 0) {
            total += static_cast<uint64_t>(st.st_size);
            ++count;
            std::remove(p.c_str());
        }
    }
    std::cout << "[disk] " << label << ": " << count << " files, "
              << humanBytes(total) << " freed" << std::endl;
}

GroupGenerator::GroupGenerator(LocalParameters & par) : par(par) {
    commonKmerDB = par.filenames[1 + (par.seqMode == 2)];
    outDir       = par.filenames[2 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    kmerFormat = par.kmerFormat;
    printLog = par.printLog;
    commonKmerSpan = par.commonKmerSpan;
    // 0 means "follow --threads", which is what every run before --partitions existed did.
    // Clamped to at least 1: routeOf takes id1 % partitionCnt, and --threads is itself >= 1.
    partitionCnt = (par.partitions > 0) ? par.partitions : par.threads;
    if (partitionCnt < 1) { partitionCnt = 1; }
    edgeMode = par.edgeMode;

    geneticCode = new RegularGeneticCode();
    queryIndexer = new QueryIndexer(par);
    queryIndexer->setKmerLen(12);
    kmerExtractor = new KmerExtractor(par, geneticCode, kmerFormat);
}

GroupGenerator::~GroupGenerator() {
    delete queryIndexer;
    delete kmerExtractor;
    delete geneticCode;
}

void GroupGenerator::startGroupGeneration(const LocalParameters &par) {  
    Buffer<Kmer> queryKmerBuffer;
    Buffer<std::pair<uint32_t, uint32_t>> matchBuffer; // seq id, pos
    vector<Query> queryList;

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;

    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;

        // new code
        if (tries == 1) {
                cout << "Indexing query file ...";
        }
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();

        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "Done" << endl;
            cout << "Total number of sequences: " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
        }

        // Set up kseq
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = nullptr;
        if (par.seqMode == 2) { kseq2 = KSeqFactory(par.filenames[1].c_str()); }

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        for (size_t splitIdx = 0; splitIdx < queryReadSplit.size(); splitIdx++) {
            // Allocate memory for query list
            queryList.clear();
            queryList.resize(queryReadSplit[splitIdx].end - queryReadSplit[splitIdx].start);

            // Allocate memory for query k-mer buffer
            queryKmerBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            queryKmerBuffer.init();
            matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt);
            matchBuffer.init();

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2);

            filterCommonKmers(queryKmerBuffer, matchBuffer, commonKmerDB, processedReadCnt);
            time_t t = time(nullptr);
            writeKmers(queryKmerBuffer, processedReadCnt);
            cout << "Writing query k-mer file: " << double(time(nullptr) - t) << " s" << endl;
            processedReadCnt += queryReadSplit[splitIdx].readCnt;
            cout << "The number of processed sequences: " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
            cout << "-----------------------------------" << endl;
        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        } 
    }   

    makeSubGraph(processedReadCnt);   

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<uint32_t> queryGroupInfo;
    queryGroupInfo.resize(processedReadCnt + 1, 0); 

    if (printLog) {
        mergeGraph_one(processedReadCnt);
    } else {
        mergeGraph();

        // An edge weight counts the k-mers two reads share, so weight / kmersPerRead is the
        // fraction of a read the two overlap by. That makes the core threshold a property of read
        // geometry: one ratio means the same overlap on every dataset, which an absolute weight
        // does not. There is no fallback -- --min-overlap-ratio is validated to be > 0 before the
        // run starts, so a threshold that cannot be formed is a hard error, not a switch to some
        // other rule.
        const double kmersPerRead = (processedReadCnt > 0)
            ? static_cast<double>(totalFilteredKmers) / static_cast<double>(processedReadCnt)
            : 0.0;
        const int coreThr = static_cast<int>(par.minOverlapRatio * kmersPerRead + 0.5);
        if (coreThr < 2) {
            cerr << "Error: core threshold resolves to " << coreThr << " (overlap ratio "
                 << par.minOverlapRatio << " x " << kmersPerRead << " k-mers/read)." << endl;
            cerr << "       It must be at least 2, or every shared k-mer would form a core edge"
                 << " and no weak band would exist." << endl;
            cerr << "       Raise --min-overlap-ratio, or check that the common k-mer filter did"
                 << " not remove nearly every k-mer." << endl;
            exit(EXIT_FAILURE);
        }

        // Weak band lower bound, as a fraction of the core threshold. This was briefly an
        // absolute count of shared k-mers (--min-edge), adopted because 0.3333 reproduces an
        // absolute 5 only where the core is 15. Measurement reversed the argument: an absolute 10
        // is 10/33 = 0.30 of the core on CAMI2 strain-madness but 10/26 = 0.38 on
        // species-exclusion, so it is the absolute form that means different things on different
        // data. And the axis is sharp -- on strain-madness at 16 threads a bound of 11 gives
        // species purity 0.932 while 10 gives 0.763, one k-mer of band width for 0.169 of purity.
        // The clamp and the warning are kept from the absolute version; they were the good part.
        const int requestedWeakThr =
            static_cast<int>(par.weakBandRatio * static_cast<float>(coreThr) + 0.5f);
        int weakThr = requestedWeakThr;
        if (weakThr < 1) { weakThr = 1; }
        if (weakThr >= coreThr) { weakThr = coreThr - 1; }

        cout << "Core threshold: " << coreThr << " (overlap ratio " << par.minOverlapRatio
             << " x " << kmersPerRead << " k-mers/read)" << endl;
        cout << "Weak band: (" << weakThr << ", " << coreThr << "] (ratio "
             << par.weakBandRatio << " x core " << coreThr << "); Phase 3 floor "
             << weakThr << endl;

        // Where that threshold falls in the weights it is cutting. The rule is unchanged for the
        // star edge mode on purpose: a star weight is not the number of k-mers two reads share --
        // an edge exists only where one of the two was the hub -- so ratio x k-mers-per-read means
        // something different there. Correcting the formula would make the correction itself the
        // thing under test and blur the only comparison that decides whether the mode is worth
        // having, which is the two modes' purity-recall curves. --min-overlap-ratio absorbs the
        // shift instead; this line is how far it has to travel.
        {
            uint64_t weightTotal = 0;
            uint64_t belowCore = 0;
            uint64_t inBand = 0;
            for (size_t b = 0; b < edgeWeightHist.size(); ++b) {
                const uint64_t n = edgeWeightHist[b];
                if (n == 0) { continue; }
                weightTotal += n;
                // Bucket b covers [2^b, 2^(b+1)); classify it by its lower bound, which is exact
                // for the powers of two and off by less than one bucket elsewhere.
                const size_t lo = (size_t(1) << b);
                if (lo < static_cast<size_t>(coreThr)) {
                    belowCore += n;
                    if (lo > static_cast<size_t>(weakThr)) { inBand += n; }
                }
            }
            cout << "[edges] weight distribution: " << weightTotal << " edges, "
                 << (weightTotal ? 100.0 * static_cast<double>(belowCore)
                                       / static_cast<double>(weightTotal) : 0.0)
                 << "% below the core threshold, "
                 << (weightTotal ? 100.0 * static_cast<double>(inBand)
                                       / static_cast<double>(weightTotal) : 0.0)
                 << "% in the weak band" << endl;
        }
        if (requestedWeakThr >= coreThr) {
            cout << "[WARN] --weak-band-ratio " << par.weakBandRatio << " puts the bound at "
                 << requestedWeakThr << ", at or above the core threshold " << coreThr
                 << "; clamped to " << weakThr << ". The weak band would otherwise be empty,"
                 << " leaving nothing for Phase 2 or Phase 3 to work with." << endl;
        }

        // Phase 1: form core groups with strong edges
        makeGroups(coreThr, processedReadCnt, groupInfo, queryGroupInfo);

        // Phase 2: merge units joined by several independent weak links
        mergeBySupport(coreThr, weakThr, processedReadCnt, groupInfo, queryGroupInfo);

        // Phase 3: link Phase-1 singletons among themselves with weak edges
        {
            std::vector<bool> isSingleton(processedReadCnt + 1, false);
            for (uint32_t i = 1; i <= processedReadCnt; i++) {
                if (queryGroupInfo[i] == 0) isSingleton[i] = true;
            }
            makeGroupsPhase3(weakThr, processedReadCnt, isSingleton, groupInfo, queryGroupInfo);
        }

        saveGroupsToFile(groupInfo, queryGroupInfo);

        // relations_*.bin are reparsed every refinement iteration but are no longer
        // needed once grouping is saved. Indices: relations_0 .. relations_{2*threads}.
        std::vector<std::string> relationFiles;
        for (int i = 0; i <= partitionCnt * 2; i++) {
            relationFiles.push_back(outDir + "/relations_" + std::to_string(i) + ".bin");
        }
        reportAndRemoveFiles(relationFiles, "relations");
    }
}

void GroupGenerator::filterCommonKmers(Buffer<Kmer> & qKmers,
                                       Buffer<std::pair<uint32_t, uint32_t>> & matchBuffer,
                                       const string & commonKmerDB,
                                       size_t processedReadCnt) {
    string gtdbListDB;
    std::string diffIdxFileName = commonKmerDB +"/diffIdx";
    std::string diffIdxSplitFileName = commonKmerDB + "/split";

    size_t blankCnt = std::find_if(qKmers.buffer,
                                   qKmers.buffer + qKmers.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - qKmers.buffer;

    size_t queryKmerNum = qKmers.startIndexOfReserve - blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    // Filter out meaningless target splits
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    std::vector<QueryKmerSplit> querySplits;
    size_t quotient = queryKmerNum / par.threads;
    size_t remainder = queryKmerNum % par.threads;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < par.threads; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        uint64_t queryAA = qKmers.buffer[startIdx].value;
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= diffIdxSplits.data[j].ADkmer) {
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j - (j != 0)]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
        }
        startIdx = endIdx + 1;
    }
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);

    time_t beforeFilter = time(nullptr);
    std::cout << "Common k-mer searching : " << std::flush;
    #pragma omp parallel default(none), shared(matchBuffer, commonKmerDB, querySplits, qKmers, cout)
    {
        Buffer<std::pair<uint32_t, uint32_t>> localMatches(1024 * 1024 * 2);  // 16 Mb <queryID, pos>
        KmerDbReader * kmerDbReader
            = new KmerDbReader(commonKmerDB + "/diffIdx", 1024 * 1024, 1024 * 1024);
        std::vector<std::pair<uint32_t, uint32_t>> tempMatches;  
        bool hasOverflow = false;
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            kmerDbReader->setReadPosition(querySplits[i].diffIdxSplit);
            uint64_t tKmer = kmerDbReader->next();
            Kmer qKmer(UINT64_MAX, 0);
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the AA matches if queries are identical
                if (qKmer.value == qKmers.buffer[j].value) {
                    if (unlikely(!localMatches.afford(tempMatches.size()))) {
                        if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                            hasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                    memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                           sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
                    for (size_t k = 0; k < tempMatches.size(); k++) {
                        localMatches.buffer[posToWrite + k].first = qKmers.buffer[j].qInfo.sequenceID;
                        localMatches.buffer[posToWrite + k].second = qKmers.buffer[j].qInfo.pos;
                    }
                    continue;
                }
                tempMatches.clear();
                // Get next query, and start to find
                qKmer = qKmers.buffer[j];

                // Skip target k-mers lexiocographically smaller
                while (!kmerDbReader->isCompleted() && qKmer.value > tKmer) {
                    tKmer = kmerDbReader->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!kmerDbReader->isCompleted() && qKmer.value == tKmer) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, (uint32_t) qKmer.qInfo.pos);
                    tKmer = kmerDbReader->next();                                      
                }

                if (unlikely(!localMatches.afford(tempMatches.size()))) {
                    if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                        hasOverflow = true;
                        break;
                    }
                }

                size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                       sizeof(std::pair<uint32_t, uint32_t>) * tempMatches.size());
            } // End of one split

            // Move matches in the local buffer to the shared buffer
            if (!Buffer<std::pair<uint32_t, uint32_t>>::moveSmallToLarge(&localMatches, &matchBuffer)) {
                hasOverflow = true;
            }
        } // End of omp for (Iterating for splits)
    } // End of omp parallel
    std::cout << time(nullptr) - beforeFilter << " s (" << matchBuffer.startIndexOfReserve << " matches found)" << endl;

    time_t here = time(nullptr);
    SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve);    
    std::cout << "Sorting matches        : " << double(time(nullptr) - here) << " s" << std::endl;
 
    // Sort query k-mers by <seqID, pos>
    time_t firstSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer + blankCnt, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQKmerByIdAndPos);
    std::cout << "Query k-mer sorting (1): " << double(time(nullptr) - firstSort) << " s" << std::endl;

    // Filter neighbor k-mers
    here = time(nullptr);
    std::cout << "Filtering common k-mers: " << std::flush;
    size_t storePos = blankCnt;
    size_t lookingPos = blankCnt;
    size_t matchIdx = 0;
    while (lookingPos < qKmers.startIndexOfReserve) {
        if (matchIdx < matchBuffer.startIndexOfReserve) {
            // copy
            if (qKmers.buffer[lookingPos].qInfo.sequenceID < matchBuffer.buffer[matchIdx].first){
                qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
            }
            // next target check
            else if(qKmers.buffer[lookingPos].qInfo.sequenceID > matchBuffer.buffer[matchIdx].first){
                matchIdx++;
            }
            // same seq
            else{
                // copy
                if (int64_t(qKmers.buffer[lookingPos].qInfo.pos) < int(matchBuffer.buffer[matchIdx].second) - commonKmerSpan){
                    qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
                }
                // next target check
                else if(int(matchBuffer.buffer[matchIdx].second) + commonKmerSpan < int64_t(qKmers.buffer[lookingPos].qInfo.pos)){
                    matchIdx++;
                }
                // pass
                else{
                    lookingPos++;
                }
            }            
        }
        else{
            qKmers.buffer[storePos++] = qKmers.buffer[lookingPos++];
        }
    }
    qKmers.startIndexOfReserve = size_t(storePos);
    cout << double(time(nullptr) - here) << " s" << endl;

    // Per-read counts over the k-mers that just survived, before the sort below reorders them.
    // The copy loop preserved the <sequenceID, pos> order, so the last kept entry carries the
    // largest id in this split and one resize covers the whole range. Split-local ids become
    // global by adding processedReadCnt, the same offset writeKmers stamps into the files.
    if (storePos > blankCnt) {
        const size_t maxGlobalId =
            static_cast<size_t>(qKmers.buffer[storePos - 1].qInfo.sequenceID) + processedReadCnt;
        if (kmerCntPerRead.size() <= maxGlobalId) {
            kmerCntPerRead.resize(maxGlobalId + 1, 0);
        }
        for (size_t i = blankCnt; i < storePos; ++i) {
            kmerCntPerRead[static_cast<size_t>(qKmers.buffer[i].qInfo.sequenceID)
                           + processedReadCnt]++;
        }
    }
    
    // sort buffer by kmer
    time_t secondSort = time(nullptr);
    SORT_PARALLEL(qKmers.buffer, qKmers.buffer + qKmers.startIndexOfReserve, Kmer::compareQueryKmer);
    secondSort = time(nullptr) - secondSort;    
    cout << "Query k-mer sorting (2): " << double(secondSort) << " s" << endl;
    cout << "Filtered k-mer number  : " << storePos - blankCnt << endl;
    // Accumulate across splits; this point is outside the OpenMP region above, so no race.
    totalFilteredKmers += storePos - blankCnt;
}

void GroupGenerator::writeKmers(Buffer<Kmer>& queryKmerBuffer,
                                size_t processedReadCnt) {
    size_t blankCnt = std::find_if(queryKmerBuffer.buffer,
                                   queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.qInfo.sequenceID != 0;}
                                  ) - queryKmerBuffer.buffer;                                        
    size_t queryKmerNum = queryKmerBuffer.startIndexOfReserve - blankCnt;

    // Make k-mer boundaries based on actual distribution
    if (!boundariesInitialized) {
        size_t quotient = queryKmerNum / par.threads;
        size_t remainder = queryKmerNum % par.threads;
        size_t idx = blankCnt;
        for (size_t i = 0; i < par.threads - 1; i++) {
            idx = idx + quotient - (i == 0);
            if (remainder > 0) {
                idx++;
                remainder--;
            }
            if (idx >= queryKmerNum + blankCnt) {
                std::cout << "Warning: endIdx exceeded queryKmerNum, adjusting to max." << std::endl;
                idx = queryKmerNum + blankCnt - 1;
            }
            kmerBoundaries.emplace_back(queryKmerBuffer.buffer[idx].value);
        }
        boundariesInitialized = true;
    }

    // Make query k-mer ranges for each split
    std::vector<std::pair<size_t, size_t>> queryKmerRanges = getKmerRanges(queryKmerBuffer, blankCnt);
    #pragma omp parallel default(none), shared(queryKmerBuffer, queryKmerRanges, processedReadCnt)
    {
        size_t threadId = omp_get_thread_num();
        size_t startIdx = queryKmerRanges[threadId].first;
        size_t endIdx = queryKmerRanges[threadId].second;
        WriteBuffer<uint16_t> diffBuffer(this->outDir + "/kmer_delta_" + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        WriteBuffer<uint32_t> infoBuffer(this->outDir + "/kmer_info_"  + to_string(this->numOfSplits) + "_" + to_string(threadId), 1024 * 1024);
        uint64_t lastKmer = 0;
        for (size_t i = startIdx; i < endIdx; i++) {
            queryKmerBuffer.buffer[i].qInfo.sequenceID += processedReadCnt;
            uint32_t id = static_cast<uint32_t>(queryKmerBuffer.buffer[i].qInfo.sequenceID);
            infoBuffer.write(&id);
            IndexCreator::getDiffIdx(lastKmer, queryKmerBuffer.buffer[i].value, diffBuffer);
        }
        infoBuffer.flush();
        diffBuffer.flush();
    }
    this->numOfSplits++;
}

std::vector<std::pair<size_t, size_t>> GroupGenerator::getKmerRanges(const Buffer<Kmer> & kmerBuffer, 
                                                                     size_t offset) {
    // Custom comparator to find the first k-mer that IS GREATER THAN the value (>)
    auto value_less_than_kmer = [](uint64_t val, const Kmer& kmer) {
        return val < kmer.value;
    };

    std::vector<std::pair<size_t, size_t>> ranges;
    size_t startIdx = offset;

    // For every boundary, find the first element > it. This marks the end of the current range.
    for (const uint64_t& boundary : kmerBoundaries) {
        // Find the first element that is > boundary.
        auto it_cut = std::upper_bound(kmerBuffer.buffer + startIdx,
                                       kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                                       boundary, 
                                       value_less_than_kmer);

        // The distance gives the index relative to the beginning of the vector
        size_t cut_idx = std::distance(kmerBuffer.buffer, it_cut);

        ranges.push_back({startIdx, cut_idx});
        startIdx = cut_idx;
    }

    // Add the final range for all elements > the last boundary
    ranges.push_back({startIdx, kmerBuffer.startIndexOfReserve});
    return ranges;
}

void GroupGenerator::scanMHistogram(size_t processedReadCnt,
                                    std::vector<uint64_t> & mHistKmers,
                                    std::vector<uint64_t> & mHistPairs) {
    mHistKmers.assign(M_HIST_BUCKETS, 0);
    mHistPairs.assign(M_HIST_BUCKETS, 0);

    // Same partitioning as the emit pass: thread t owns k-mer value range t, so the per-thread
    // histograms are over disjoint k-mers and their sum is order-independent. Two runs with
    // different --threads therefore produce the same histogram, and the same cap.
    #pragma omp parallel num_threads(par.threads)
    {
        const int threadIdx = omp_get_thread_num();
        std::vector<uint64_t> localKmers(M_HIST_BUCKETS, 0);
        std::vector<uint64_t> localPairs(M_HIST_BUCKETS, 0);
        uint64_t localValues = 0;

        scanKmerRuns(threadIdx, processedReadCnt, localValues,
                     [&](uint64_t, const std::vector<uint32_t> & queryIds) {
            const size_t m = queryIds.size();
            if (m < 2) {
                return;
            }
            const size_t bucket = mHistBucket(m);
            localKmers[bucket]++;
            localPairs[bucket] += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
        });

        #pragma omp critical
        {
            for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
                mHistKmers[b] += localKmers[b];
                mHistPairs[b] += localPairs[b];
            }
        }
    }
}

void GroupGenerator::makeSubGraph(size_t processedReadCnt) {
    cout << "Connecting reads with shared k-mer..." << endl;
    time_t beforeSearch = time(nullptr);

    // Checked here, before hours of work, because the alternative is finding out inside the
    // OpenMP region where ReadBuffer's exit turns into a segmentation fault with no diagnosis.
    // Both scans below -- the quantile pre-pass and the emit pass -- hold every split of every
    // thread open at once; see scanKmerRunsFdDemand for why lowering the thread count only in
    // this stage would drop whole partitions instead of fixing it.
    {
        const size_t demand = scanKmerRunsFdDemand(numOfSplits, par.threads);
        const size_t limit = getOpenFileLimit();
        const size_t held = openFdCount();
        // Room for the subGraph writers, the relations writers, stdio and the taxonomy handles.
        const size_t reserve = emitUnitCount(partitionCnt) + 64;
        if (demand + held + reserve > limit) {
            const size_t needed = demand + held + reserve;
            cerr << "Error: this run needs about " << needed << " open files but the soft limit is "
                 << limit << "." << endl;
            cerr << "       The k-mer scan holds every split open on every thread: "
                 << par.threads << " threads x " << numOfSplits << " splits x 2 = " << demand
                 << ", plus " << held << " already open and " << reserve << " in reserve." << endl;
            const size_t spare = (limit > held + reserve) ? (limit - held - reserve) : 0;
            const size_t fits = (numOfSplits > 0) ? spare / (numOfSplits * 2) : 0;
            cerr << "       Raise it with 'ulimit -n " << needed << "' (hard limit permitting)";
            if (fits >= 1) {
                cerr << ", or run with --threads " << fits
                     << " or fewer -- the scan's demand is proportional to it, and the grouping"
                     << " output does not depend on the thread count." << endl;
            } else {
                cerr << ". No thread count fits at this limit: even one thread needs "
                     << (numOfSplits * 2) << " descriptors for the scan." << endl;
            }
            exit(EXIT_FAILURE);
        }
    }

    // The per-read counts have to add up to the total the same filter reported, or the split
    // offset is wrong and every hub choice built on them is wrong with it. Cheap to check and
    // impossible to notice otherwise -- a mis-offset array still produces a plausible graph.
    {
        uint64_t perReadSum = 0;
        for (size_t i = 0; i < kmerCntPerRead.size(); ++i) {
            perReadSum += kmerCntPerRead[i];
        }
        cout << "[kmers] per-read counts: " << perReadSum << " summed over "
             << (kmerCntPerRead.empty() ? 0 : kmerCntPerRead.size() - 1) << " reads, "
             << totalFilteredKmers << " filtered in total ("
             << (perReadSum == totalFilteredKmers ? "match" : "MISMATCH") << "), "
             << humanBytes(kmerCntPerRead.size() * sizeof(uint32_t)) << " held" << endl;
        if (perReadSum != totalFilteredKmers) {
            cerr << "Error: per-read k-mer counts do not sum to the filtered k-mer total." << endl;
            cerr << "       The star edge mode picks its hub from these counts, so a wrong sum" << endl;
            cerr << "       means a wrong graph rather than a wrong log line." << endl;
            exit(EXIT_FAILURE);
        }
    }

    // Shared across threads, not per thread: see getRelationBudget. Each thread still keeps its
    // own map -- concurrent writes to one map would need a lock on the hot path -- but the
    // decision to flush is made against the running total, so the number of flushes depends on
    // how much data there is and not on how many workers are chewing through it.
    const size_t RELATION_BUDGET = getRelationBudget(par.ramUsage);
    std::atomic<size_t> bufferedPairs(0);

    // Emit buffers, one per output unit and shared by every thread. This is the whole point of
    // the 2026-08-27 change: with a map per thread, a flush emptied one thread's share and the
    // shared budget was reached again after another 1/N of it, so flushes -- and the files each
    // one writes -- scaled with --threads. Held per unit, a flush empties all of them at once,
    // so a file carries budget/units whatever the thread count is.
    const size_t emitUnitCnt = emitUnitCount(partitionCnt);
    const size_t emitRouteCnt = static_cast<size_t>(partitionCnt) * 2 + 1;
    const size_t emitRangeSize = getRouteRangeSize(processedReadCnt, partitionCnt);
    // A unit's buffer is striped over this many mutexes. The stripe of a pair is a function of
    // its key, so a unit's stripes hold disjoint keys and together hold exactly what one map
    // held: same entries, same weights, same count against the budget. Only the contention
    // changes -- see emitUnitShards and emitShardOf.
    const size_t emitShardCnt = emitUnitShards(emitUnitCnt, par.threads);
    const size_t emitBufferCnt = emitUnitCnt * emitShardCnt;
    std::vector<std::unique_ptr<EmitUnitBuffer>> emitUnits;
    emitUnits.reserve(emitBufferCnt);
    for (size_t b = 0; b < emitBufferCnt; ++b) {
        emitUnits.push_back(std::unique_ptr<EmitUnitBuffer>(new EmitUnitBuffer()));
        // The budget is shared across units, so a buffer's fair share is its slice of it. The
        // per-thread reserve this replaces divided by --threads instead, which is the same
        // arithmetic applied to the wrong denominator. Striping divides the same total further,
        // so the reserved capacity summed over buffers does not change.
        emitUnits.back()->pairs.reserve(std::max<size_t>(1, RELATION_BUDGET / emitBufferCnt));
    }
    // (route, shard) of each flat unit index, so the file name survives the flattening.
    std::vector<std::pair<size_t, size_t>> emitUnitKey(emitUnitCnt, {0, 0});
    for (size_t route = 0; route < emitRouteCnt; ++route) {
        const size_t shardCnt = shardsForRoute(route, partitionCnt);
        for (size_t shard = 0; shard < shardCnt; ++shard) {
            emitUnitKey[emitUnitIndex(route, shard, partitionCnt)] = {route, shard};
        }
    }
    // Serialises the harvest half of a flush round. Blocking, not try_lock: a thread that found a
    // round already running used to return and keep inserting, so the buffers grew past
    // RELATION_BUDGET for as long as the round took to write. Measured on CAMI2 strain-madness at
    // --partitions 16, the worst round carried 2.3x the budget at 4 threads, 5.0x at 16 and 38.6x
    // at 64 -- 7,721,398,786 pairs, which is 370 GB at this budget's own 48 B per entry against a
    // --max-ram of 128 GiB. The overrun grew with the worker count because that is how many
    // threads were racing past the closed door.
    std::mutex emitFlushMutex;
    // Rounds that had to wait for the harvest lock, and how long in total. This is the cost the
    // blocking lock adds; it is bounded by the harvest being emitUnitCnt pointer swaps.
    std::atomic<uint64_t> flushHarvestWaits(0);
    std::atomic<uint64_t> flushHarvestWaitMillis(0);
    // Time inside saveSubGraphToFile, summed over rounds and at its worst for one round. Rounds
    // overlap now, so the total can exceed the wall clock.
    std::atomic<uint64_t> flushWriteMillisTotal(0);
    std::atomic<uint64_t> flushWriteMillisMax(0);
    // Time spent releasing the harvested maps. A round holds up to RELATION_BUDGET nodes, each an
    // individual allocation, and freeing them is neither counted by `flush write` nor visible in
    // its build/sort/io split -- but it happens while the round still holds an in-flight slot, so
    // it keeps the cap closed exactly as the write does.
    std::atomic<uint64_t> flushFreeMillisTotal(0);
    std::atomic<uint64_t> flushFreeMillisMax(0);

    // Rounds harvested but not yet written. The harvest empties the buffers, so the budget is
    // honoured as soon as it runs -- but that data is still resident until the write finishes.
    // Without a cap the peak would again be bounded only by how fast the disk keeps up, which is
    // the same unbounded growth entered from the other side.
    const size_t maxFlushRoundsInFlight = concurrentFlushRounds();
    std::mutex flushInFlightMutex;
    std::condition_variable flushInFlightCv;
    size_t flushRoundsInFlight = 0;
    std::atomic<uint64_t> flushInFlightWaits(0);
    std::atomic<uint64_t> flushInFlightWaitMillis(0);
    // Budget-triggered calls that found the buffers already emptied and returned without
    // queueing for a slot. This is what the pre-queue check keeps out of the in-flight queue,
    // so it is the measure of what that check is worth.
    std::atomic<uint64_t> flushEarlyOuts(0);
    // Contention on the per-unit locks, in the only terms that matter here -- how often a
    // thread had to wait at all. Reported so a configuration with few units (--partitions 1
    // gives three) shows the cost instead of hiding it.
    std::atomic<uint64_t> unitLockWaits(0);
    std::atomic<uint64_t> unitLockTakes(0);
    // Thread-microseconds spent blocked on a buffer's mutex. The counts above say how often a
    // thread had to wait but not what the waiting cost, and without that the emit's thread-time
    // cannot be split between scanning and contention -- job 550906 left 36.9% of it unattributed.
    // Microseconds because a hold is a few thousand hash inserts, which is under a millisecond.
    std::atomic<uint64_t> unitLockWaitMicros(0);

    // The projected peak is printed next to the budget on purpose. Until 2026-09-04 the budget was
    // whatever a 200 M clamp allowed and the peak was nowhere in the log, so a run that declared
    // --max-ram 128 quietly used a fifth of it and nothing said so.
    cout << "Flush budget: " << RELATION_BUDGET << " pairs across all threads (memory budget "
         << humanBytes(getMemoryBudgetBytes(par.ramUsage)) << ", --max-ram "
         << par.ramUsage << " GiB; projected emit peak "
         << humanBytes(projectedEmitPeakBytes(par.ramUsage)) << " = "
         << peakBytesPerBufferedPair() << " B/pair x " << (concurrentFlushRounds() + 1)
         << " rounds)" << endl;

    // Disk counterpart of the memory budget above. A safety ceiling, not a target: a run that
    // fits under it behaves exactly as it does without one.
    const size_t tmpDiskBudget = getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB);
    if (par.maxTmpDiskMiB > 0) {
        cout << "Tmp disk budget: " << humanBytes(tmpDiskBudget)
             << " (--max-tmp-disk " << par.maxTmpDiskMiB << " MiB)" << endl;
    } else if (tmpDiskBudget == TMP_DISK_UNLIMITED) {
        cout << "[WARN] could not measure free space at " << outDir
             << "; running without a tmp disk cap. Pass --max-tmp-disk to set one." << endl;
    } else {
        cout << "Tmp disk budget: " << humanBytes(tmpDiskBudget) << " (auto: 80% of "
             << humanBytes(getFreeDiskBytes(outDir)) << " free at " << outDir << ")" << endl;
    }

    // Second fold trigger, independent of the budget: the merge can only open so many streams
    // per unit, so a prefix longer than the fan-in has to be folded regardless of disk. Note
    // that getMergeFanIn looks at the fd limit and the thread count only -- it knows nothing
    // about how much data is on disk, which is why the budget above is the primary trigger.
    // One merger: maybeFoldEmitted walks the units serially, so only one unit's streams are
    // open at any moment. Passing the thread count here shrank the fan-in for no reason.
    const size_t foldFanIn = getMergeFanIn(1);
    std::atomic<int> counter(0);

    // High-frequency k-mers are dropped by --max-kmer-freq-ratio. Record which ones so the
    // information loss stays auditable instead of silent.
    size_t skippedKmerCnt = 0;
    size_t skippedMaxM = 0;
    size_t skippedSumM = 0;
    double skippedPairEst = 0.0; // sum of C(m,2) over skipped k-mers; double because it overflows uint64

    // Edge-volume instrumentation: how much the C(m,2) clique actually produced, and the
    // largest m that survived the skip -- that m sets the per-k-mer worst case.
    size_t emittedEdgeCnt = 0; // Relation records written to subGraph_*
    size_t maxKeptM = 0;       // largest reads-per-k-mer that was NOT skipped
    // Pairs carried by each flush. A round empties every unit at once, so this is what one
    // subGraph_* generation holds and what the merge then has to work through. Reported as a
    // distribution because the mean hides the final partial rounds.
    std::vector<size_t> flushPairCounts;
    // Guards the two counters above. They used to be covered by emitFlushMutex, which is no
    // longer held while a round is written -- rounds now overlap, so they need a lock of their
    // own. Taken once per round, off the hot path.
    std::mutex emitStatsMutex;

    // A round has two halves with very different costs. Harvesting is a swap per unit and is
    // over in microseconds; writing is the whole round's data through zstd and is not. Only the
    // harvest is serialised: it decides which pairs belong to which flush index, and two threads
    // interleaving there would split a unit across indices. The write is left outside, so a
    // round in progress does not stop the next one from being harvested.
    //
    // Concurrent writes are safe because each round owns a distinct counter index and every file
    // name carries it (subGraphName(counter, route, shard)), so no two rounds touch the same
    // file. saveSubGraphToFile keeps its own state under subGraphSizeMutex and marks the index
    // complete only after the last file is closed, which is what emitWatermark walks -- a fold
    // can therefore never list a file that is still being written.
    //
    // Releases an in-flight slot however the round leaves. saveSubGraphToFile is fatal on error
    // rather than throwing, but a slot leaked on any path would stall every other thread, so the
    // release is not left to the control flow.
    struct FlushSlot {
        std::mutex & mutex;
        std::condition_variable & cv;
        size_t & count;
        ~FlushSlot() {
            {
                const std::lock_guard<std::mutex> guard(mutex);
                --count;
            }
            cv.notify_one();
        }
    };

    // `force` distinguishes the last round, after the parallel region joins, from the rounds the
    // budget triggers. A budget-triggered round is pointless once another thread has already
    // taken the buffers below the budget, and skipping it is what keeps the round count off the
    // thread count: without the re-check below, every thread that crossed the threshold queues on
    // the harvest lock and each spends a flush index on the crumbs the one before it left. On the
    // 20x fixture at --max-ram 1 that was 2 / 3 / 4 / 8 rounds for 1 / 4 / 16 / 64 threads --
    // exactly the dependency the unit buffers were introduced to remove.
    auto flushEmitUnits = [&](const bool force) {
        // The same re-check the harvest performs below, applied before this thread commits to the
        // in-flight queue. Without it every thread that crossed the budget queues for a slot,
        // wins it, finds the buffers already emptied and returns -- measured on CAMI2
        // strain-madness (job 550906) as 33,782 queue entries against 297 rounds, so 99.1% of
        // them did no work, and 547,187 thread-seconds spent in the queue against 1,071 seconds
        // of round writing to wait for. The threads were waiting on each other, not on the disk.
        //
        // The harvest subtracts from bufferedPairs before the write starts, precisely so that a
        // thread testing the budget during a write sees the emptied buffers; reserving a slot
        // first denied it the chance to look. Checking here narrows the window in which a herd
        // can form from a round's write to the harvest's pointer swaps.
        //
        // Cannot skip a needed round: the buffers are above the budget only if no thread has
        // harvested them, and every thread re-tests the budget after each k-mer, so one of them
        // arrives here and proceeds. `force` is the final flush after the parallel region joins
        // and has to run whatever the counters say.
        if (!force && bufferedPairs.load(std::memory_order_relaxed) < RELATION_BUDGET) {
            flushEarlyOuts.fetch_add(1, std::memory_order_relaxed);
            return;
        }
        // Reserve a slot before touching emitFlushMutex. Waiting with the harvest lock held would
        // stop every other thread from harvesting as well, which is the stall this cap exists to
        // bound rather than to create.
        {
            std::unique_lock<std::mutex> inFlightGuard(flushInFlightMutex);
            if (flushRoundsInFlight >= maxFlushRoundsInFlight) {
                const std::chrono::steady_clock::time_point waitStart =
                    std::chrono::steady_clock::now();
                flushInFlightWaits.fetch_add(1, std::memory_order_relaxed);
                flushInFlightCv.wait(inFlightGuard, [&]() {
                    return flushRoundsInFlight < maxFlushRoundsInFlight;
                });
                flushInFlightWaitMillis.fetch_add(
                    static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - waitStart).count()),
                    std::memory_order_relaxed);
            }
            ++flushRoundsInFlight;
        }
        const FlushSlot slot{flushInFlightMutex, flushInFlightCv, flushRoundsInFlight};

        std::vector<std::vector<std::unordered_map<uint64_t, uint16_t>>> taken(
                emitUnitCnt, std::vector<std::unordered_map<uint64_t, uint16_t>>(emitShardCnt));
        size_t movedPairs = 0;
        size_t counter_now = 0;
        {
            // The harvest, and nothing else. Held for emitUnitCnt pointer swaps. Blocking: a
            // thread that gives up here goes back to inserting into buffers already at the
            // budget, which is exactly the overrun this replaces.
            std::unique_lock<std::mutex> flushGuard(emitFlushMutex, std::try_to_lock);
            if (!flushGuard.owns_lock()) {
                const std::chrono::steady_clock::time_point waitStart =
                    std::chrono::steady_clock::now();
                flushHarvestWaits.fetch_add(1, std::memory_order_relaxed);
                flushGuard.lock();
                flushHarvestWaitMillis.fetch_add(
                    static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - waitStart).count()),
                    std::memory_order_relaxed);
            }
            // Whoever held the lock may have emptied the buffers already.
            if (!force
                && bufferedPairs.load(std::memory_order_relaxed) < RELATION_BUDGET) {
                return;
            }
            for (size_t u = 0; u < emitUnitCnt; ++u) {
                for (size_t s = 0; s < emitShardCnt; ++s) {
                    EmitUnitBuffer & buffer = *emitUnits[u * emitShardCnt + s];
                    const std::lock_guard<std::mutex> unitGuard(buffer.lock);
                    taken[u][s].swap(buffer.pairs);
                    movedPairs += taken[u][s].size();
                }
            }
            if (movedPairs == 0) {
                return; // do not spend a flush index on a round that would write only empty files
            }
            // Subtracted here rather than after the write, so a thread testing the budget while
            // this round is still being written sees the buffers as the emptied ones they are.
            bufferedPairs.fetch_sub(movedPairs, std::memory_order_relaxed);
            counter_now = static_cast<size_t>(counter.fetch_add(1, std::memory_order_relaxed));
        }

        // Milliseconds, not the project's usual time(nullptr) seconds: a round on the fixture is
        // under a second, and this is the half of a round that the harvest no longer waits on.
        const std::chrono::steady_clock::time_point writeStart = std::chrono::steady_clock::now();
        saveSubGraphToFile(taken, emitUnitKey, counter_now);
        const uint64_t writeMillis = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - writeStart).count());
        flushWriteMillisTotal.fetch_add(writeMillis, std::memory_order_relaxed);
        // Compare-exchange rather than a plain store: another round may finish while this one is
        // reading, and the maximum has to survive that.
        uint64_t prevMax = flushWriteMillisMax.load(std::memory_order_relaxed);
        while (prevMax < writeMillis
               && !flushWriteMillisMax.compare_exchange_weak(prevMax, writeMillis,
                                                             std::memory_order_relaxed)) {
        }

        // Released here rather than left to the scope exit, so the cost is visible. It is the same
        // work either way -- the maps die inside this call either at this line or a few lines
        // below -- but `flush write` above does not include it and the in-flight slot does, so a
        // round's hold on the cap was longer than any counter reported.
        {
            const std::chrono::steady_clock::time_point freeStart =
                std::chrono::steady_clock::now();
            taken.clear();
            taken.shrink_to_fit();
            const uint64_t freeMillis = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - freeStart).count());
            flushFreeMillisTotal.fetch_add(freeMillis, std::memory_order_relaxed);
            uint64_t prevFreeMax = flushFreeMillisMax.load(std::memory_order_relaxed);
            while (prevFreeMax < freeMillis
                   && !flushFreeMillisMax.compare_exchange_weak(prevFreeMax, freeMillis,
                                                                std::memory_order_relaxed)) {
            }
        }

        // Read only after the parallel region joins, but written here by rounds that now overlap.
        {
            const std::lock_guard<std::mutex> statsGuard(emitStatsMutex);
            emittedEdgeCnt += movedPairs;
            flushPairCounts.push_back(movedPairs);
        }

        checkDiskHeadroom();
        // Fold what is already complete instead of letting the whole emit pile up. Safe to call
        // with other rounds in flight: it folds only [foldedUpTo, emitWatermark), and
        // emitWatermark is the contiguous prefix of indices saveSubGraphToFile has finished.
        maybeFoldEmitted(tmpDiskBudget, foldFanIn);
    };
    // Sum of C(m,2) over kept k-mers = pair OCCURRENCES the clique produced. Records written
    // are fewer, because pair2weight collapses a pair seen across several k-mers into one
    // entry. emittedEdgeCnt / keptPairSum is that collapse factor, measured on this run.
    uint64_t keptPairSum = 0;
    // Sum of (m-1) over the same kept k-mers: what a center-star would have emitted where the
    // clique emitted C(m,2). Counted here rather than estimated because the ratio depends on the
    // whole m distribution, not on its mean -- one m = 255 k-mer is worth 127 m = 3 k-mers.
    uint64_t keptStarSum = 0;

    // m distribution over ALL k-mers, skipped or not. The [edges] summary only reports what
    // the chosen thresholds produced; this reports what any other threshold would have.
    std::vector<uint64_t> mHistKmers(M_HIST_BUCKETS, 0);
    std::vector<uint64_t> mHistPairs(M_HIST_BUCKETS, 0);
    std::vector<uint64_t> mHistStar(M_HIST_BUCKETS, 0);

    // Skip threshold, resolved once. An explicit --max-kmer-reads wins and skips the pre-pass
    // entirely; otherwise --max-kmer-quantile picks the cap from the reads-per-k-mer
    // distribution, which has to be measured before any edge is emitted. The cap stays an
    // absolute m, only its choice becomes a property of the data: m does not scale with the
    // size of the dataset but with per-genome coverage, so neither a hand-picked constant nor a
    // fraction of the read count transfers (--max-kmer-freq-ratio failed on exactly that).
    size_t absThr = (par.maxKmerReads > 0) ? static_cast<size_t>(par.maxKmerReads) : 0;
    bool capIsAuto = false;
    if (par.maxKmerReads <= 0 && par.maxKmerQuantile > 0.0f && par.maxKmerQuantile < 1.0f) {
        cout << "[skip] pre-pass: measuring the reads-per-k-mer distribution..." << endl;
        const time_t preStart = time(nullptr);
        std::vector<uint64_t> preHistKmers;
        std::vector<uint64_t> preHistPairs;
        scanMHistogram(processedReadCnt, preHistKmers, preHistPairs);
        // The pre-pass counted every k-mer value once already; the emit pass counts them again.
        // Leaving both in would double the denominator the disk projection divides by.
        kmerValuesTotal.store(0, std::memory_order_relaxed);

        absThr = capFromQuantile(preHistKmers, par.maxKmerQuantile);
        capIsAuto = (absThr > 0);

        uint64_t totalKmers = 0, totalPairs = 0, keptKmers = 0, keptPairs = 0;
        for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
            totalKmers += preHistKmers[b];
            totalPairs += preHistPairs[b];
            // Bucket b covers m in [2^b, 2^(b+1) - 1]; it survives when its top is under the cap.
            const bool kept = (absThr == 0)
                || (b + 1 < 63 && ((static_cast<size_t>(1) << (b + 1)) - 1) <= absThr);
            if (kept) {
                keptKmers += preHistKmers[b];
                keptPairs += preHistPairs[b];
            }
        }
        cout << "[skip] pre-pass done: " << double(time(nullptr) - preStart) << " s, "
             << totalKmers << " k-mers with m >= 2" << endl;
        cout << "[skip] auto cap " << (absThr > 0 ? to_string(absThr) : string("none"))
             << " (quantile " << par.maxKmerQuantile << "): keeps "
             << (totalKmers ? 100.0 * static_cast<double>(keptKmers)
                                    / static_cast<double>(totalKmers) : 0.0)
             << "% of k-mers, "
             << (totalPairs ? 100.0 * static_cast<double>(keptPairs)
                                    / static_cast<double>(totalPairs) : 0.0)
             << "% of pair occurrences" << endl;
    }
    cout << "[skip] threshold: --max-kmer-reads "
         << (absThr > 0 ? to_string(absThr) : string("off"))
         << (capIsAuto ? " (auto)" : "") << endl;

    #pragma omp parallel num_threads(par.threads)
    {
        int threadIdx = omp_get_thread_num();
        const string skippedPartName = outDir + "/skipped_kmers_" + to_string(threadIdx);
        ofstream skippedPart(skippedPartName);
        if (!skippedPart.is_open()) {
            cerr << "Error opening file: " << skippedPartName << endl;
        }
        size_t localSkippedCnt = 0;
        size_t localMaxM = 0;
        size_t localSumM = 0;
        double localPairEst = 0.0;
        size_t localMaxKeptM = 0;
        // Pair occurrences staged for each (unit, key shard), flushed into the shared buffer in
        // batches. The batch is what makes the per-buffer lock affordable: the alternative is
        // taking it inside the C(m,2) loop below, once per occurrence. The sub-batch is the
        // whole batch divided by the shard count, so what a thread stages in total does not
        // grow with the striping.
        const size_t stageBatch = std::max<size_t>(1, EMIT_STAGE_BATCH / emitShardCnt);
        std::vector<std::vector<uint64_t>> stage(emitBufferCnt);
        for (size_t b = 0; b < emitBufferCnt; ++b) {
            stage[b].reserve(stageBatch);
        }

        // Merge one buffer's staged occurrences into its shared map. Only newly created keys
        // count towards the budget: a repeat of a pair already buffered adds no footprint. The
        // count is exact rather than conservative because a key reaches one buffer only.
        // Contention counters kept per thread and added once after the scan: a drain happens every
        // few thousand pair occurrences, and touching an atomic that often would make the
        // instrumentation a source of the contention it measures.
        uint64_t localLockTakes = 0;
        uint64_t localLockWaits = 0;
        uint64_t localLockWaitMicros = 0;

        auto drainStage = [&](size_t b) {
            if (stage[b].empty()) {
                return;
            }
            size_t newKeys = 0;
            {
                std::unique_lock<std::mutex> guard(emitUnits[b]->lock, std::try_to_lock);
                if (!guard.owns_lock()) {
                    // Timed only on the contended path, so an uncontended drain does not pay for
                    // two clock reads.
                    const std::chrono::steady_clock::time_point waitStart =
                        std::chrono::steady_clock::now();
                    ++localLockWaits;
                    guard.lock();
                    localLockWaitMicros += static_cast<uint64_t>(
                        std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::steady_clock::now() - waitStart).count());
                }
                ++localLockTakes;
                std::unordered_map<uint64_t, uint16_t> & shared = emitUnits[b]->pairs;
                for (const uint64_t pairKey : stage[b]) {
                    const size_t before = shared.size();
                    addSat16(shared[pairKey], 1);
                    if (shared.size() != before) { newKeys++; }
                }
            }
            stage[b].clear();
            if (newKeys != 0) {
                bufferedPairs.fetch_add(newKeys, std::memory_order_relaxed);
            }
        };
        std::vector<uint64_t> localHistKmers(M_HIST_BUCKETS, 0);
        std::vector<uint64_t> localHistPairs(M_HIST_BUCKETS, 0);
        std::vector<uint64_t> localHistStar(M_HIST_BUCKETS, 0);
        uint64_t localKeptPairs = 0;
        uint64_t localKeptStar = 0;

        uint64_t localValuesDone = 0;

        scanKmerRuns(threadIdx, processedReadCnt, localValuesDone,
                     [&](uint64_t minKmer, const std::vector<uint32_t> & currentQueryIds) {
            const size_t m = currentQueryIds.size();

            // Histogram every k-mer, skipped or not, so the summary table can answer
            // "what would a different cap have cost".
            if (m >= 2) {
                const size_t bucket = mHistBucket(m);
                localHistKmers[bucket]++;
                localHistPairs[bucket] += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
                localHistStar[bucket] += static_cast<uint64_t>(m - 1);
            }

            if (absThr > 0 && m > absThr) {
                skippedPart << minKmer << "\t" << m << "\n";
                localSkippedCnt++;
                localSumM += m;
                if (m > localMaxM) { localMaxM = m; }
                localPairEst += static_cast<double>(m) * static_cast<double>(m - 1) * 0.5;
                return; // `continue` before the scan became a callback
            }
            // Every read pair sharing this k-mer gets an edge: the full C(m,2) clique.
            // Edge weight therefore counts shared k-mers directly, which is what the
            // Phase 1/2 thresholds and the knee detection assume. currentQueryIds is
            // sorted and deduplicated above, so id[i] < id[j] holds and pairKey needs no
            // min/max normalization. The quadratic term is bounded by the high-frequency
            // k-mer skip above, not by reducing the clique.
            if (m > localMaxKeptM) { localMaxKeptM = m; }
            if (m >= 2) {
                localKeptPairs += static_cast<uint64_t>(m) * static_cast<uint64_t>(m - 1) / 2;
                localKeptStar += static_cast<uint64_t>(m - 1);
            }
            if (edgeMode == EDGE_MODE_STAR) {
                // One hub, m-1 spokes. m < 2 leaves the hub alone with no edge to draw, which is
                // the same thing the clique loop below does with a one-read list.
                if (m >= 2) {
                    const uint32_t hub = pickStarHub(currentQueryIds, kmerCntPerRead);
                    for (size_t i = 0; i < m; ++i) {
                        const uint32_t other = currentQueryIds[i];
                        if (other == hub) {
                            continue;
                        }
                        // The hub is chosen by k-mer count, so it is not necessarily the smaller
                        // id. The clique path below can skip this because its list is sorted and
                        // it walks i < j; here the pair has to be ordered explicitly or the same
                        // pair would route to two different units and the merge would count it
                        // twice.
                        const uint32_t id1 = (hub < other) ? hub : other;
                        const uint32_t id2 = (hub < other) ? other : hub;
                        const size_t route = routeOf(id1, id2, partitionCnt, emitRangeSize);
                        const size_t u = emitUnitIndex(
                            route, shardOf(id1, shardsForRoute(route, partitionCnt)), partitionCnt);
                        const uint64_t pairKey = (static_cast<uint64_t>(id1) << 32)
                                                 | static_cast<uint64_t>(id2);
                        const size_t b = u * emitShardCnt
                                         + emitShardOf(pairKey, emitShardCnt);
                        stage[b].push_back(pairKey);
                        if (stage[b].size() >= stageBatch) {
                            drainStage(b);
                        }
                    }
                }
                return; // the clique loop below is the other branch of this choice
            }
            for (size_t i = 0; i + 1 < m; ++i) {
                const uint32_t id1 = currentQueryIds[i];
                const uint64_t idHi = static_cast<uint64_t>(id1) << 32;
                for (size_t j = i + 1; j < m; ++j) {
                    const uint32_t id2 = currentQueryIds[j];
                    // Route at insertion time, not at write time. The unit a pair belongs to
                    // is a function of (id1, id2) alone, so deciding it here changes nothing
                    // about where the pair ends up -- it only moves the decision earlier, which
                    // is what lets the buffer be held per unit.
                    const size_t route = routeOf(id1, id2, partitionCnt, emitRangeSize);
                    const size_t u = emitUnitIndex(
                        route, shardOf(id1, shardsForRoute(route, partitionCnt)), partitionCnt);
                    const uint64_t pairKey = idHi | static_cast<uint64_t>(id2);
                    const size_t b = u * emitShardCnt + emitShardOf(pairKey, emitShardCnt);
                    stage[b].push_back(pairKey);
                    if (stage[b].size() >= stageBatch) {
                        drainStage(b);
                    }
                }
            }
            // The budget is tested against the shared total, which is what every thread's
            // drains have contributed. Staged occurrences not yet drained are not counted:
            // they are bounded by EMIT_STAGE_BATCH per unit per thread and are not in a map.
            if (bufferedPairs.load(std::memory_order_relaxed) >= RELATION_BUDGET) {
                // Everything this thread is holding goes in before the round starts, so a
                // flush carries as much as it can rather than leaving a batch behind.
                for (size_t b = 0; b < emitBufferCnt; ++b) {
                    drainStage(b);
                }
                // Publish progress before the guards read it, so the projection matches the
                // bytes that are about to be written.
                kmerValuesDone.fetch_add(localValuesDone, std::memory_order_relaxed);
                localValuesDone = 0;
                // No barrier: one thread writes the round while the others keep filling the
                // now-empty buffers, and their data goes out in the next one.
                flushEmitUnits(false);
            }
        });

        // Nothing may stay staged past the end of this thread's scan.
        for (size_t b = 0; b < emitBufferCnt; ++b) {
            drainStage(b);
        }
        // After the last drain, so nothing this thread did is left out.
        unitLockTakes.fetch_add(localLockTakes, std::memory_order_relaxed);
        unitLockWaits.fetch_add(localLockWaits, std::memory_order_relaxed);
        unitLockWaitMicros.fetch_add(localLockWaitMicros, std::memory_order_relaxed);

        // No per-thread final write: what this thread built is in the shared unit buffers, and
        // the last round after the parallel region joins carries it out.
        skippedPart.close();

        #pragma omp critical
        {
            skippedKmerCnt  += localSkippedCnt;
            skippedSumM     += localSumM;
            skippedPairEst  += localPairEst;
            if (localMaxM > skippedMaxM) { skippedMaxM = localMaxM; }
            if (localMaxKeptM > maxKeptM) { maxKeptM = localMaxKeptM; }
            keptPairSum += localKeptPairs;
            keptStarSum += localKeptStar;
            for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
                mHistKmers[b] += localHistKmers[b];
                mHistPairs[b] += localHistPairs[b];
                mHistStar[b] += localHistStar[b];
            }
        }
    }
    // Final round. Every thread drained its staging before exiting, so whatever is left in the
    // unit buffers is the tail below the budget -- one more round writes it out. try_lock inside
    // flushEmitUnits cannot fail here: the parallel region has joined.
    flushEmitUnits(true);

    this->numOfGraph = counter.load(std::memory_order_relaxed);

    // Concatenate the per-thread skip logs. Thread t owns k-mer value range t (writeKmers
    // fixes kmerBoundaries once, so the ranges are disjoint and ordered), which makes thread
    // order equal to ascending k-mer order -- no merge sort needed. Unlike the other
    // intermediates this file is diagnostic output and is kept, not deleted.
    {
        const string skippedFileName = outDir + "/skipped_kmers";
        ofstream skippedAll(skippedFileName);
        if (!skippedAll.is_open()) {
            cerr << "Error opening file: " << skippedFileName << endl;
        } else {
            std::vector<std::string> skippedParts;
            for (int t = 0; t < par.threads; t++) {
                skippedParts.push_back(outDir + "/skipped_kmers_" + to_string(t));
            }
            for (const std::string & part : skippedParts) {
                ifstream in(part);
                // Streaming an empty rdbuf sets failbit on the destination, so guard it.
                if (in.is_open() && in.peek() != std::ifstream::traits_type::eof()) {
                    skippedAll << in.rdbuf();
                }
            }
            skippedAll.close();
            reportAndRemoveFiles(skippedParts, "skipped_kmers parts");
        }
    }

    if (skippedKmerCnt > 0) {
        cout << "[skip] " << skippedKmerCnt << " k-mers skipped (max m: " << skippedMaxM
             << ", sum m: " << skippedSumM
             << ", est. dropped read-pairs: " << skippedPairEst << ")" << endl;
        cout << "[skip] see " << outDir << "/skipped_kmers" << endl;
    } else if (absThr == 0) {
        cout << "[skip] 0 k-mers skipped (--max-kmer-reads off)" << endl;
    } else {
        cout << "[skip] 0 k-mers skipped (no k-mer exceeded the threshold)" << endl;
    }

    // k-mer files are fully consumed by the merge above; report + remove to free disk.
    {
        std::vector<std::string> deltaFiles, infoFiles;
        for (size_t i = 0; i < this->numOfSplits; i++) {
            for (int t = 0; t < par.threads; t++) {
                deltaFiles.push_back(outDir + "/kmer_delta_" + to_string(i) + "_" + to_string(t));
                infoFiles.push_back(outDir + "/kmer_info_"  + to_string(i) + "_" + to_string(t));
            }
        }
        reportAndRemoveFiles(deltaFiles, "kmer_delta");
        reportAndRemoveFiles(infoFiles, "kmer_info");
    }

    // m-distribution table. Each row is a candidate --max-kmer-reads value (the bucket's
    // upper bound) with the edge volume that cap would produce: the cumulative sum of
    // C(m,2) over every k-mer with m at or below it. Pick the row that fits the disk budget.
    // Counts cover all k-mers with m >= 2, whether or not the current run skipped them, so
    // the table stays valid for thresholds other than the one just used.
    {
        size_t topBucket = 0;
        uint64_t totalKmers = 0;
        uint64_t totalPairs = 0;
        for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
            totalKmers += mHistKmers[b];
            totalPairs += mHistPairs[b];
            if (mHistKmers[b] > 0) { topBucket = b; }
        }

        // Pair occurrences over-count what actually lands on disk, because pair2weight
        // collapses a pair seen across several k-mers into one record. Measure that factor
        // on this run and apply it, so the table predicts records rather than an upper bound.
        // Against the occurrences this run actually produced: star emitted (m-1) per k-mer, not
        // C(m,2), so dividing by the clique sum there would report a collapse that never happened.
        const uint64_t emittedOccurrences =
            (edgeMode == EDGE_MODE_STAR) ? keptStarSum : keptPairSum;
        const double dedup = (emittedOccurrences > 0)
            ? static_cast<double>(emittedEdgeCnt) / static_cast<double>(emittedOccurrences)
            : 1.0;

        cout << "[mhist] reads-per-k-mer distribution (m >= 2): " << totalKmers
             << " k-mers, " << totalPairs << " pair occurrences uncapped" << endl;
        cout << "[mhist] measured collapse factor: " << emittedEdgeCnt << " records / "
             << emittedOccurrences << " kept "
             << (edgeMode == EDGE_MODE_STAR ? "star" : "pair") << " occurrences = "
             << dedup << endl;

        // What a center-star would have emitted on exactly the k-mers this run kept. The clique
        // figure is the same sum with C(m,2) in place of (m-1), so the ratio is the volume this
        // pipeline would trade away for the graph precision a star gives up. Measured, not
        // estimated: it is set by the whole m distribution, and one m = 255 k-mer is worth 127
        // m = 3 ones.
        const double starRatio = (keptStarSum > 0)
            ? static_cast<double>(keptPairSum) / static_cast<double>(keptStarSum)
            : 0.0;
        cout << "[mhist] star vs clique on kept k-mers: " << keptPairSum
             << " clique pair occurrences / " << keptStarSum << " star (m-1) = "
             << starRatio << "x" << endl;
        cout << "[mhist] estimates below extrapolate that factor to other caps; it was"
             << " measured only on k-mers this run kept, so expect roughly +/-20%" << endl;
        if (totalKmers > 0) {
            cout << "[mhist] " << setw(16) << "cap (max m)"
                 << setw(14) << "k-mers"
                 << setw(18) << "pairs in range"
                 << setw(18) << "cum pairs"
                 << setw(18) << "cum star (m-1)"
                 << setw(16) << "est. records"
                 << setw(12) << "est. disk" << endl;
            uint64_t cumPairs = 0;
            uint64_t cumStar = 0;
            for (size_t b = 0; b <= topBucket; ++b) {
                cumPairs += mHistPairs[b];
                cumStar += mHistStar[b];
                if (mHistKmers[b] == 0) { continue; }
                const size_t lo = (size_t(1) << b);
                const size_t hi = (size_t(1) << (b + 1)) - 1;
                const uint64_t estRecords = static_cast<uint64_t>(static_cast<double>(cumPairs) * dedup);
                cout << "[mhist] " << setw(16) << (to_string(hi) + " (" + to_string(lo) + "-" + to_string(hi) + ")")
                     << setw(14) << mHistKmers[b]
                     << setw(18) << mHistPairs[b]
                     << setw(18) << cumPairs
                     << setw(18) << cumStar
                     << setw(16) << estRecords
                     << setw(12) << humanBytes(estRecords * sizeof(Relation)) << endl;
            }
            cout << "[mhist] read off the row whose est. disk fits the budget,"
                 << " then pass its cap as --max-kmer-reads" << endl;
        }
    }

    // Edge-volume summary. maxKeptM is the worst per-k-mer case that survived the skip:
    // that single k-mer contributed C(maxKeptM, 2) edges, so it is the knob to watch when
    // tuning --max-kmer-freq-ratio against a disk budget.
    const uint64_t rawEdgeBytes = static_cast<uint64_t>(emittedEdgeCnt) * sizeof(Relation);
    const uint64_t diskEdgeBytes = this->subGraphBytesOnDisk.load(std::memory_order_relaxed);
    cout << "[edges] emitted " << emittedEdgeCnt << " edge records into "
         << this->numOfGraph << " subgraphs ("
         << humanBytes(rawEdgeBytes) << " before merge, "
         << humanBytes(diskEdgeBytes) << " on disk, zstd level " << RELATION_ZSTD_LEVEL;
    if (diskEdgeBytes > 0) {
        cout << ", " << fixed << setprecision(2)
             << (static_cast<double>(rawEdgeBytes) / static_cast<double>(diskEdgeBytes)) << "x";
    }
    cout << ")" << endl;
    // Files, not flushes. Every flush writes one file per merge unit, and the unit count is
    // 3 x threads today (2*threads+1 routes, with the cross bucket sharded `threads` ways), so
    // the file count carries a factor of the thread count that the flush count alone hides.
    // Raising --threads shrinks each thread's share of the memory budget, which raises the flush
    // count too, and the product is what the merge then has to fold: CAMI2 strain-madness went
    // from 389 flushes x 48 units at 16 threads to 2,354 x 192 at 64, so 18,672 files became
    // 451,968 for the same 117.29 GB, and the run took 8h 08m instead of 3h 40m. Printed here
    // rather than only at merge time so the cost is visible before the fold starts.
    {
        const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;
        size_t unitCnt = 0;
        for (size_t r = 0; r < routeCnt; ++r) {
            unitCnt += shardsForRoute(r, partitionCnt);
        }
        cout << "[edges] " << this->numOfGraph << " flushes x " << unitCnt << " units ("
             << routeCnt << " routes, cross bucket sharded x"
             << shardsForRoute(static_cast<size_t>(partitionCnt), partitionCnt) << ") = "
             << (this->numOfGraph * unitCnt) << " subGraph files, at --threads "
             << par.threads << " --partitions " << partitionCnt << endl;

        // The two distributions the totals above cannot show. Whether raising --threads split
        // the same data into more, smaller pieces -- or actually moved more data -- is the
        // difference between these lines and the byte total, and it is the question the merge
        // time depends on. Both should be flat across --threads for a given input; anything
        // else is the flush unit tracking the thread count.
        std::vector<uint64_t> flushWidths;
        flushWidths.reserve(flushPairCounts.size());
        for (const size_t n : flushPairCounts) { flushWidths.push_back(n); }
        const OrderStats flushStats = orderStats(flushWidths);
        cout << "[edges] pairs per flush: min " << flushStats.min
             << ", median " << flushStats.median
             << ", max " << flushStats.max
             << " (budget " << RELATION_BUDGET << " shared across " << par.threads
             << " threads)" << endl;

        // What the budget actually bought. A ratio above 1 means rounds carried more than the
        // budget allows, and the excess is held in memory at RELATION_BUDGET's own 48 bytes per
        // entry -- so this ratio times --max-ram is what the run really asked the machine for.
        // Flat across --threads is the property that makes the tool portable; anything else means
        // the overrun tracks the worker count and the same input needs a bigger machine on a
        // bigger node.
        const double budgetRatio = RELATION_BUDGET > 0
            ? static_cast<double>(flushStats.max) / static_cast<double>(RELATION_BUDGET)
            : 0.0;
        cout << "[edges] flush budget: max/budget " << fixed << setprecision(3) << budgetRatio
             << "x over " << flushPairCounts.size() << " rounds, in-flight cap "
             << maxFlushRoundsInFlight << endl;

        // What enforcing the budget costs. Harvest waits are threads queued behind a swap loop;
        // cap waits are threads held back because two rounds were already unwritten. The second
        // is the one that trades wall time for a bounded peak, so it is reported separately.
        // Early-outs are the threads that crossed the budget but found the buffers already
        // emptied: they used to queue for a slot to learn that. Cap waits near the round count
        // mean the queue now holds only threads that go on to write something.
        cout << "[edges] flush waits: harvest " << flushHarvestWaits.load(std::memory_order_relaxed)
             << " (" << flushHarvestWaitMillis.load(std::memory_order_relaxed) << " ms), cap "
             << flushInFlightWaits.load(std::memory_order_relaxed)
             << " (" << flushInFlightWaitMillis.load(std::memory_order_relaxed) << " ms), "
             << flushEarlyOuts.load(std::memory_order_relaxed)
             << " early-outs before the queue" << endl;

        cout << "[edges] flush write: " << flushWriteMillisTotal.load(std::memory_order_relaxed)
             << " ms total, " << flushWriteMillisMax.load(std::memory_order_relaxed)
             << " ms worst round" << endl;

        // Reported next to the write because it is the other half of what keeps the in-flight cap
        // closed, and it is single-threaded where the write is not.
        cout << "[edges] flush free: " << flushFreeMillisTotal.load(std::memory_order_relaxed)
             << " ms total, " << flushFreeMillisMax.load(std::memory_order_relaxed)
             << " ms worst round" << endl;

        // The same time, split by what spends it. These are summed over every worker of every
        // round, so they exceed both the line above and the wall clock once the unit loop runs in
        // parallel -- the ratio between the four is what they are for, not their total.
        const uint64_t buildMs = emitBuildMicros.load(std::memory_order_relaxed) / 1000;
        const uint64_t sortMs = emitSortMicros.load(std::memory_order_relaxed) / 1000;
        const uint64_t writeMs = emitWriteMicros.load(std::memory_order_relaxed) / 1000;
        const uint64_t freeMs = emitFreeMicros.load(std::memory_order_relaxed) / 1000;
        cout << "[edges] flush write split: build " << buildMs << " ms, sort " << sortMs
             << " ms, compress+io " << writeMs << " ms, free " << freeMs << " ms" << endl;

        // What the split costs in memory. The maps of a round are already resident at 48 B per
        // pair; these vectors are the same pairs at sizeof(Relation), held only while a worker
        // sorts and writes one unit.
        cout << "[edges] flush writers: " << emitWriterCnt.load(std::memory_order_relaxed)
             << " per round, vector peak "
             << humanBytes(emitVectorPeakBytes.load(std::memory_order_relaxed)) << endl;

        // The ceiling on splitting the unit loop across workers: a round is never finished before
        // its largest unit is, so the largest unit's share of a round bounds the speed-up at
        // 1/share whatever the worker count.
        const uint64_t maxUnitRec = emitMaxUnitRecords.load(std::memory_order_relaxed);
        const uint64_t roundRec = emitRoundRecords.load(std::memory_order_relaxed);
        const double meanShare = roundRec ? static_cast<double>(maxUnitRec)
                                                / static_cast<double>(roundRec) : 0.0;
        const double worstShare =
            static_cast<double>(emitWorstShareScaled.load(std::memory_order_relaxed)) / 100000.0;
        cout << "[edges] unit skew: largest unit is " << fixed << setprecision(2)
             << (100.0 * meanShare) << "% of a round on average (worst " << (100.0 * worstShare)
             << "%), so a parallel unit loop caps at "
             << (meanShare > 0.0 ? 1.0 / meanShare : 0.0) << "x" << endl;

        OrderStats fileStats{0, 0, 0};
        size_t fileCnt = 0;
        {
            const std::lock_guard<std::mutex> guard(subGraphSizeMutex);
            fileCnt = subGraphFileBytes.size();
            fileStats = orderStats(subGraphFileBytes);
        }
        cout << "[edges] subGraph file size: median " << humanBytes(fileStats.median)
             << ", max " << humanBytes(fileStats.max)
             << " over " << fileCnt << " files, "
             << this->foldRounds << " folds" << endl;

        // Once per run, next to the other emit costs. The [disk] fold line prints once per fold
        // and carries no duration; this is what says whether folding is worth attacking.
        if (this->foldRounds > 0) {
            const uint64_t foldMs = foldMillisTotal.load(std::memory_order_relaxed);
            const uint64_t units = foldUnitsDone.load(std::memory_order_relaxed);
            // Wall clock per fold, not summed unit time: the unit loop runs several units at
            // once, so dividing the fold's duration by its units would understate what a unit
            // costs by the worker count. The worst single unit is the honest per-unit figure and
            // the ceiling on what parallelising the loop can win.
            cout << "[edges] fold cost: " << foldMs << " ms of wall clock over "
                 << this->foldRounds << " folds (mean " << (this->foldRounds ? foldMs / this->foldRounds : 0)
                 << " ms, worst " << foldMillisMax.load(std::memory_order_relaxed)
                 << " ms), worst single unit " << foldUnitMillisMax.load(std::memory_order_relaxed)
                 << " ms over " << units << " unit merges" << endl;
        }

        // Cost of holding the buffers per unit rather than per thread. The lock domain is a
        // (unit, key shard), so few units and many threads is the configuration to watch:
        // --partitions 1 leaves three units. Staging keeps the take count down to roughly
        // occurrences / EMIT_STAGE_BATCH, and a wait share near zero means the batch is doing
        // its job. Takes rise with the shard count because a batch is split that many ways;
        // the share is what to compare across configurations, not the count.
        const uint64_t takes = unitLockTakes.load(std::memory_order_relaxed);
        const uint64_t waits = unitLockWaits.load(std::memory_order_relaxed);
        cout << "[edges] unit lock: " << takes << " takes, " << waits << " contended ("
             << fixed << setprecision(2)
             << (takes > 0 ? 100.0 * static_cast<double>(waits) / static_cast<double>(takes) : 0.0)
             << "%), " << emitUnitCnt << " units x " << emitShardCnt << " key shards = "
             << emitBufferCnt << " lock domains for " << par.threads << " threads" << endl;
        // The share is what the contention costs, which the counts above cannot say: a wait is
        // reported whether it lasted a microsecond or a second. The denominator is the emit's
        // whole thread-time, so this is directly comparable with the cap wait on the line above.
        {
            const uint64_t lockMicros = unitLockWaitMicros.load(std::memory_order_relaxed);
            const double emitThreadSeconds = static_cast<double>(time(nullptr) - beforeSearch)
                                             * static_cast<double>(par.threads < 1 ? 1
                                                                                  : par.threads);
            cout << "[edges] unit lock wait: " << (lockMicros / 1000) << " ms of thread time";
            if (emitThreadSeconds > 0.0) {
                cout << " (" << fixed << setprecision(2)
                     << (100.0 * (static_cast<double>(lockMicros) / 1e6) / emitThreadSeconds)
                     << "% of emit thread time)";
            }
            cout << endl;
        }
    }
    // What the run actually held at once, which is the number a shared filesystem cares about.
    // "0 folds" means the data fit under the budget and nothing was folded early -- the path
    // this build takes when there is room, identical to the one before incremental folding.
    cout << "[disk] peak subGraph footprint: "
         << humanBytes(this->subGraphPeakBytes.load(std::memory_order_relaxed)) << " (budget "
         << (getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB) == TMP_DISK_UNLIMITED
             ? string("unlimited")
             : humanBytes(getTmpDiskBudgetBytes(outDir, par.maxTmpDiskMiB)))
         << ", " << this->foldRounds << " folds)" << endl;
    cout << "[edges] largest kept k-mer: m = " << maxKeptM << " -> "
         << (static_cast<double>(maxKeptM) * static_cast<double>(maxKeptM > 0 ? maxKeptM - 1 : 0) * 0.5)
         << " edges from that k-mer alone" << endl;

    cout << "Relations generated from files successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupGenerator::saveSubGraphToFile(
        std::vector<std::vector<std::unordered_map<uint64_t, uint16_t>>> & units,
        const std::vector<std::pair<size_t, size_t>> & unitKey,
        const size_t counter_now) {
    // The units arrive already partitioned: routeOf/shardOf ran at insertion time, so a pair
    // reached exactly one unit and cannot appear under two. That is what still makes the
    // per-route merges independent.
    // One entry per file, indexed by unit rather than appended: the units are written by several
    // workers below and the order they finish in is not fixed, but the log has to be.
    std::vector<uint64_t> flushFileBytes(units.size(), 0);

    // A unit owns its maps, its sort and its file, and two units never touch the same bytes, so the
    // loop below splits with nothing to synchronise. That the output cannot change with the split
    // rests on three facts, in order of how easy they are to get wrong:
    //   - a unit's file name is subGraphName(counter_now, route, shard) and no two units share one;
    //   - a file's contents are a function of its unit's maps alone. Relation::compare is a total
    //     order on (id1, id2), a map cannot hold the same key twice, and the unit's key shards hold
    //     disjoint keys (emitShardOf is a function of the key), so the sort has one answer;
    //   - everything written outside a worker's own slots is either a sum (compressedBytes, the
    //     three timing totals) or a maximum (the skew), both order-independent.
    // markEmitComplete stays after the join for the same reason it was after the loop: it is what
    // lets a fold list this index, and a fold may not see an index whose files are still open.
    // RelationWriter failure calls relationStreamFatal, which exits the process -- dying inside a
    // worker ends the run just as it did inside the loop, and this repository does not use
    // exceptions, so there is nothing to propagate across the join.
    const size_t writerCnt = subGraphWriters(units.size(), par.threads);
    std::vector<uint64_t> wBuild(writerCnt, 0), wSort(writerCnt, 0), wWrite(writerCnt, 0);
    std::vector<uint64_t> wFree(writerCnt, 0);
    std::vector<uint64_t> wRecords(writerCnt, 0), wMaxUnit(writerCnt, 0), wBytes(writerCnt, 0);

    const auto microsSince = [](const std::chrono::steady_clock::time_point & from) -> uint64_t {
        return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now() - from).count());
    };

    // Strided, not blocked: adjacent units are the shards of one route and carry similar loads, so
    // a stride spreads the heavy cross shards over every worker. Static either way, which keeps a
    // run reproducible.
    const auto writeUnits = [&](const size_t w) {
        for (size_t u = w; u < units.size(); u += writerCnt) {
            const std::chrono::steady_clock::time_point buildStart =
                std::chrono::steady_clock::now();
            // The unit's shards are concatenated, not merged: their key sets are disjoint, so no
            // pair can appear twice in the result and no weight needs re-summing.
            size_t unitRecords = 0;
            for (const auto & shard : units[u]) {
                unitRecords += shard.size();
            }
            std::vector<Relation> relations;
            relations.reserve(unitRecords);
            for (const auto & shard : units[u]) {
                for (const auto & entry : shard) {
                    const uint64_t pairKey = entry.first;
                    relations.push_back(makeRelation(static_cast<uint32_t>(pairKey >> 32),
                                                     static_cast<uint32_t>(pairKey & 0xFFFFFFFF),
                                                     entry.second));
                }
            }
            wBuild[w] += microsSince(buildStart);

            // The records are copied out, so the maps are dead weight from here on. Releasing them
            // now rather than at the caller's scope exit is the same work on a different thread:
            // this worker owns unit u alone, and no other worker or round reads it. It also takes
            // the unit's footprint down to the vector before the sort touches it.
            const std::chrono::steady_clock::time_point freeStart =
                std::chrono::steady_clock::now();
            for (auto & shard : units[u]) {
                std::unordered_map<uint64_t, uint16_t>().swap(shard);
            }
            wFree[w] += microsSince(freeStart);

            const std::chrono::steady_clock::time_point sortStart =
                std::chrono::steady_clock::now();
            sort(relations.begin(), relations.end(), Relation::compare);
            wSort[w] += microsSince(sortStart);

            wRecords[w] += relations.size();
            if (relations.size() > wMaxUnit[w]) { wMaxUnit[w] = relations.size(); }

            // Empty units still get an empty file: the merge's input list is then simply
            // subGraph_{0..numOfGraph-1}_{route}_{shard} with no existence checks.
            const string subGraphFileName = subGraphName(counter_now, unitKey[u].first,
                                                         unitKey[u].second);
            const std::chrono::steady_clock::time_point writeStart =
                std::chrono::steady_clock::now();
            RelationWriter out(subGraphFileName);
            for (size_t i = 0; i < relations.size(); ++i) {
                out.write(relations[i]);
            }
            out.finish();
            wWrite[w] += microsSince(writeStart);
            wBytes[w] += out.compressedBytes();
            flushFileBytes[u] = out.compressedBytes();
        }
    };

    if (writerCnt <= 1) {
        writeUnits(0);
    } else {
        // Worker 0 runs here rather than in a thread of its own: one fewer thread per round, and
        // --threads 1 then takes exactly the serial path above.
        std::vector<std::thread> writers;
        writers.reserve(writerCnt - 1);
        for (size_t w = 1; w < writerCnt; ++w) {
            writers.emplace_back(writeUnits, w);
        }
        writeUnits(0);
        for (size_t w = 0; w < writers.size(); ++w) {
            writers[w].join();
        }
    }

    uint64_t compressedBytes = 0;
    uint64_t buildMicros = 0, sortMicros = 0, writeMicros = 0, freeMicros = 0;
    uint64_t roundRecords = 0, maxUnitRecords = 0;
    for (size_t w = 0; w < writerCnt; ++w) {
        compressedBytes += wBytes[w];
        buildMicros += wBuild[w];
        sortMicros += wSort[w];
        writeMicros += wWrite[w];
        freeMicros += wFree[w];
        roundRecords += wRecords[w];
        if (wMaxUnit[w] > maxUnitRecords) { maxUnitRecords = wMaxUnit[w]; }
    }

    // Records held as vectors at once, which is what this change adds to the peak. Each worker
    // holds one unit at a time, so the bound is writers x the largest unit -- against the maps of
    // the same round, already resident at 48 B per pair to this 12.
    {
        const uint64_t vectorBound = static_cast<uint64_t>(writerCnt) * maxUnitRecords
                                     * sizeof(Relation);
        uint64_t peak = emitVectorPeakBytes.load(std::memory_order_relaxed);
        while (peak < vectorBound
               && !emitVectorPeakBytes.compare_exchange_weak(peak, vectorBound,
                                                             std::memory_order_relaxed)) {
        }
        emitWriterCnt.store(writerCnt, std::memory_order_relaxed);
    }

    emitBuildMicros.fetch_add(buildMicros, std::memory_order_relaxed);
    emitSortMicros.fetch_add(sortMicros, std::memory_order_relaxed);
    emitWriteMicros.fetch_add(writeMicros, std::memory_order_relaxed);
    emitFreeMicros.fetch_add(freeMicros, std::memory_order_relaxed);
    emitMaxUnitRecords.fetch_add(maxUnitRecords, std::memory_order_relaxed);
    emitRoundRecords.fetch_add(roundRecords, std::memory_order_relaxed);
    if (roundRecords > 0) {
        const uint64_t shareScaled = (maxUnitRecords * 100000) / roundRecords;
        uint64_t worst = emitWorstShareScaled.load(std::memory_order_relaxed);
        while (worst < shareScaled
               && !emitWorstShareScaled.compare_exchange_weak(worst, shareScaled,
                                                              std::memory_order_relaxed)) {
        }
    }

    {
        const std::lock_guard<std::mutex> guard(subGraphSizeMutex);
        subGraphFileBytes.insert(subGraphFileBytes.end(),
                                 flushFileBytes.begin(), flushFileBytes.end());
    }
    subGraphBytesOnDisk.fetch_add(compressedBytes, std::memory_order_relaxed);
    noteSubGraphBytes(static_cast<int64_t>(compressedBytes));
    // Only now, with every unit file of this flush closed, may the folder see this index.
    markEmitComplete(counter_now);
}

void GroupGenerator::noteSubGraphBytes(int64_t delta) {
    const int64_t live = subGraphLiveBytes.fetch_add(delta, std::memory_order_relaxed) + delta;
    if (live <= 0) {
        return;
    }
    const uint64_t liveBytes = static_cast<uint64_t>(live);
    uint64_t peak = subGraphPeakBytes.load(std::memory_order_relaxed);
    while (liveBytes > peak
           && !subGraphPeakBytes.compare_exchange_weak(peak, liveBytes,
                                                       std::memory_order_relaxed)) {
        // compare_exchange_weak refreshed `peak`; retry until it stops growing.
    }
}

void GroupGenerator::checkDiskHeadroom() {
    const size_t freeBytes = getFreeDiskBytes(outDir);
    if (freeBytes == 0) {
        return; // cannot measure; getTmpDiskBudgetBytes already warned about that
    }

    // Hard guard first. It needs no estimate, so it holds even when the projection is wrong.
    const size_t reserve = getDiskReserveBytes(outDir);
    if (freeBytes <= reserve) {
        cerr << "[ERROR] out of disk space at " << outDir << ": "
             << humanBytes(freeBytes) << " free, keeping " << humanBytes(reserve)
             << " in reserve.\n"
             << "        Intermediates hold " << humanBytes(subGraphBytesOnDisk.load())
             << " so far. Point the output directory at a larger filesystem, or lower\n"
             << "        --max-kmer-reads (currently " << par.maxKmerReads
             << ") to cut the edge volume -- see the [mhist] table for what each cap costs."
             << endl;
        exit(EXIT_FAILURE);
    }

    // Projection. Warn once, do not stop: folding decouples the footprint from the total.
    if (volumeWarned.load(std::memory_order_relaxed) != 0) {
        return;
    }
    const uint64_t total = kmerValuesTotal.load(std::memory_order_relaxed);
    const uint64_t done = kmerValuesDone.load(std::memory_order_relaxed);
    if (total == 0 || done == 0) {
        return;
    }
    const double progress = static_cast<double>(done) / static_cast<double>(total);
    if (progress < 0.02 || progress >= 1.0) {
        return; // too early to extrapolate, or already finished
    }
    const double projected = static_cast<double>(subGraphBytesOnDisk.load()) / progress;
    if (projected <= static_cast<double>(freeBytes)) {
        return;
    }
    int expected = 0;
    if (!volumeWarned.compare_exchange_strong(expected, 1)) {
        return; // another thread is printing it
    }
    cout << "[WARN] at " << fixed << setprecision(1) << (progress * 100.0)
         << "% of the emit, the intermediates project to "
         << humanBytes(static_cast<uint64_t>(projected)) << " in total but only "
         << humanBytes(freeBytes) << " is free at " << outDir << ".\n"
         << "       Folding keeps what is held at once well below that total, so the run may\n"
         << "       still finish. If it does not, lower --max-kmer-reads (currently "
         << par.maxKmerReads << ") -- the [mhist] table prints what each cap would cost --\n"
         << "       or use a larger filesystem. This early in the run the projection is rough\n"
         << "       and tends to overshoot; it is a heads-up, not a verdict." << endl;
}

void GroupGenerator::markEmitComplete(size_t flushIdx) {
    std::lock_guard<std::mutex> guard(foldMutex);
    if (emitDone.size() <= flushIdx) {
        emitDone.resize(flushIdx + 1, 0);
    }
    emitDone[flushIdx] = 1;
    while (emitWatermark < emitDone.size() && emitDone[emitWatermark] != 0) {
        ++emitWatermark;
    }
}

void GroupGenerator::maybeFoldEmitted(size_t tmpDiskBudget, size_t maxFanIn) {
    const int64_t live = subGraphLiveBytes.load(std::memory_order_relaxed);
    const bool overBudget = tmpDiskBudget != TMP_DISK_UNLIMITED
                            && live > 0
                            && static_cast<size_t>(live) >= tmpDiskBudget;

    // Cheap pre-check outside the lock. The authoritative test is repeated under it.
    if (!overBudget) {
        std::lock_guard<std::mutex> guard(foldMutex);
        if (emitWatermark - foldedUpTo < maxFanIn) {
            return;
        }
    }

    std::unique_lock<std::mutex> guard(foldMutex, std::try_to_lock);
    if (!guard.owns_lock()) {
        return; // another thread is folding; keep emitting instead of waiting on it
    }

    const size_t from = foldedUpTo;
    const size_t to = emitWatermark;
    // Nothing to gain from folding a single flush: it would rewrite the same records.
    if (to - from < 2) {
        return;
    }

    const size_t round = foldRounds++;
    foldedUpTo = to;
    if (foldedOutputs.empty()) {
        foldedOutputs.resize(emitUnitCount(partitionCnt));
    }

    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;

    // The units, flattened, so the loop below can be split across workers. Order is the order the
    // serial loop visited them in, and every worker writes only to its own units' slots, so the
    // result does not depend on how they are shared out.
    std::vector<std::pair<size_t, size_t>> foldUnits;   // (route, shard)
    foldUnits.reserve(emitUnitCount(partitionCnt));
    for (size_t route = 0; route < routeCnt; ++route) {
        const size_t shardCnt = shardsForRoute(route, partitionCnt);
        for (size_t shard = 0; shard < shardCnt; ++shard) {
            foldUnits.emplace_back(route, shard);
        }
    }

    // Units are independent -- a unit's inputs are its own subGraph_* files, its output name
    // carries its route and shard, and foldedOutputs is indexed by unit -- so the only thing that
    // kept this serial was descriptors. foldWorkers asks how many are actually free.
    //
    // Measured on CAMI2 plant-associated, where the single fold took 34,003 s of a 43,463 s emit:
    // 48 units at 708 s each, one at a time, while 127 cores idled.
    const size_t workerCnt = foldWorkers(foldUnits.size(), maxFanIn);
    // One extra stream per unit for the previous fold's output, which is re-folded below. The
    // merger count is what it now is: passing 1 while running several would size every worker's
    // read buffers as if it were alone and multiply the memory by the worker count.
    const size_t bufElems = getMergeBufferElems(to - from + 1, par.ramUsage, workerCnt);

    // Timed because nothing timed it before. A fold rewrites every unit's accumulated output, so
    // it is the single longest thing this stage does on a large input, and it holds an in-flight
    // slot while it runs -- yet the only trace it left was the write time of whichever round
    // happened to overlap it.
    std::vector<size_t> wFiles(workerCnt, 0);
    std::vector<uint64_t> wUnitMax(workerCnt, 0);
    const std::chrono::steady_clock::time_point foldStart = std::chrono::steady_clock::now();

    // Strided, and static either way, which keeps a run reproducible. mergeSubGraphBatch is fatal
    // on error rather than throwing, so a worker that dies ends the process just as it did inside
    // the serial loop; this repository does not use exceptions and there is nothing to propagate
    // across the join.
    const auto foldRange = [&](const size_t w) {
        for (size_t u = w; u < foldUnits.size(); u += workerCnt) {
            const size_t route = foldUnits[u].first;
            const size_t shard = foldUnits[u].second;
            const std::chrono::steady_clock::time_point unitStart =
                std::chrono::steady_clock::now();
            // Earlier folds' outputs go back in. Folding only the new flushes would let those
            // outputs pile up untouched, and the live footprint would climb past the budget one
            // fold at a time. Re-folding keeps exactly one file per unit, so the footprint
            // tracks the merged-so-far size -- the floor no amount of folding can go below.
            std::vector<std::string> & folded =
                foldedOutputs[emitUnitIndex(route, shard, partitionCnt)];
            std::vector<std::string> inputs = folded;
            inputs.reserve(inputs.size() + (to - from));
            for (size_t i = from; i < to; ++i) {
                inputs.push_back(subGraphName(i, route, shard));
            }
            // A distinct prefix from reduceSubGraphFanIn's subGraph_p*: that one numbers its
            // rounds independently during the final merge, and a shared name would let the two
            // overwrite each other.
            const std::string output = outDir + "/subGraph_f" + to_string(round)
                                     + "_" + to_string(route) + "_" + to_string(shard);
            mergeSubGraphBatch(inputs, output, bufElems);

            for (size_t i = 0; i < inputs.size(); ++i) {
                noteSubGraphBytes(-static_cast<int64_t>(FileUtil::getFileSize(inputs[i])));
                std::remove(inputs[i].c_str());
            }
            noteSubGraphBytes(static_cast<int64_t>(FileUtil::getFileSize(output)));
            folded.assign(1, output);
            ++wFiles[w];

            const uint64_t unitMillis = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - unitStart).count());
            if (unitMillis > wUnitMax[w]) { wUnitMax[w] = unitMillis; }
        }
    };

    if (workerCnt <= 1) {
        foldRange(0);
    } else {
        // Worker 0 runs here rather than in a thread of its own, so a single worker takes exactly
        // the serial path above.
        std::vector<std::thread> folders;
        folders.reserve(workerCnt - 1);
        for (size_t w = 1; w < workerCnt; ++w) {
            folders.emplace_back(foldRange, w);
        }
        foldRange(0);
        for (size_t w = 0; w < folders.size(); ++w) {
            folders[w].join();
        }
    }

    size_t foldedFiles = 0;
    uint64_t foldUnitMaxMillis = 0;
    for (size_t w = 0; w < workerCnt; ++w) {
        foldedFiles += wFiles[w];
        if (wUnitMax[w] > foldUnitMaxMillis) { foldUnitMaxMillis = wUnitMax[w]; }
    }

    // Accumulated once per fold rather than once per unit: a fold is rare and this keeps the
    // counters off the unit loop, which step 5 of this plan parallelises.
    const uint64_t foldMillis = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - foldStart).count());
    foldMillisTotal.fetch_add(foldMillis, std::memory_order_relaxed);
    foldUnitsDone.fetch_add(foldedFiles, std::memory_order_relaxed);
    {
        uint64_t prev = foldMillisMax.load(std::memory_order_relaxed);
        while (prev < foldMillis
               && !foldMillisMax.compare_exchange_weak(prev, foldMillis,
                                                       std::memory_order_relaxed)) {
        }
    }
    {
        uint64_t prev = foldUnitMillisMax.load(std::memory_order_relaxed);
        while (prev < foldUnitMaxMillis
               && !foldUnitMillisMax.compare_exchange_weak(prev, foldUnitMaxMillis,
                                                           std::memory_order_relaxed)) {
        }
    }

    // fan-in and peak fds in the same terms the final merge reports, so the two stages can be
    // compared. Peak is one unit's streams plus its output, times the workers folding at once.
    // Predicted, then what the process is actually holding, because a wrong prediction is exactly
    // how a run dies on a host with a smaller descriptor limit than the one it was tuned on --
    // and the worker count is derived from that same reading taken a moment earlier.
    cout << "[disk] fold " << round << ": flushes [" << from << ", " << to << ") -> "
         << foldedFiles << " files, live "
         << humanBytes(static_cast<uint64_t>(std::max<int64_t>(0, subGraphLiveBytes.load())))
         << " (" << workerCnt << " workers, fan-in " << maxFanIn << ", peak fds "
         << (mergeFdsPerUnit(maxFanIn) * workerCnt)
         << " predicted / " << openFdCount() << " open, soft limit " << getOpenFileLimit() << ")"
         << endl;
}

// No lock: this runs after makeSubGraph's parallel region has joined, so folding is over and
// foldedUpTo/foldedOutputs are stable.
std::vector<std::string> GroupGenerator::unitInputFiles(size_t route, size_t shard) const {
    std::vector<std::string> files;
    const size_t unit = emitUnitIndex(route, shard, partitionCnt);
    if (unit < foldedOutputs.size()) {
        files = foldedOutputs[unit];
    }
    files.reserve(files.size() + (this->numOfGraph - foldedUpTo));
    for (size_t i = foldedUpTo; i < this->numOfGraph; ++i) {
        files.push_back(subGraphName(i, route, shard));
    }
    return files;
}

size_t GroupGenerator::mergeSubGraphBatch(const std::vector<std::string> & inputs,
                                          const std::string & output,
                                          size_t bufElems) {
    std::vector<std::unique_ptr<RelationReader>> readers(inputs.size());
    for (size_t i = 0; i < inputs.size(); ++i) {
        readers[i].reset(new RelationReader(inputs[i], bufElems));
    }

    RelationWriter out(output);

    MergeHeap heap;
    for (size_t i = 0; i < readers.size(); ++i) {
        const Relation first = readers[i]->getNext();
        if (!(first == Relation())) {
            MergeHeapEntry entry;
            entry.rel = first;
            entry.stream = i;
            heap.push(entry);
        }
    }

    size_t written = 0;
    while (!heap.empty()) {
        const Relation minRelation = heap.top().rel;
        uint16_t totalWeight = 0;
        while (!heap.empty() && heap.top().rel == minRelation) {
            const size_t stream = heap.top().stream;
            addSat16(totalWeight, heap.top().rel.weight);
            heap.pop();
            const Relation next = readers[stream]->getNext();
            if (!(next == Relation())) {
                MergeHeapEntry entry;
                entry.rel = next;
                entry.stream = stream;
                heap.push(entry);
            }
        }
        out.write(makeRelation(minRelation.id1, minRelation.id2, totalWeight));
        ++written;
    }
    out.finish();

    return written;
}

void GroupGenerator::reduceSubGraphFanIn(std::vector<std::string> & files, size_t maxFanIn,
                                         size_t route, size_t shard) {
    if (files.size() <= maxFanIn) {
        return;
    }

    const size_t bufElems = getMergeBufferElems(maxFanIn, par.ramUsage, par.threads);
    size_t round = 0;
    while (files.size() > maxFanIn) {
        std::vector<std::string> next;
        size_t batchIdx = 0;
        for (size_t start = 0; start < files.size(); start += maxFanIn) {
            const size_t end = std::min(start + maxFanIn, files.size());
            const std::vector<std::string> batch(files.begin() + start, files.begin() + end);
            // Route AND shard are part of the name: units fold at the same time, and a shared
            // name would have two of them write the same file. Route alone is not enough --
            // the cross bucket's shards all carry the same route.
            const std::string out = outDir + "/subGraph_p" + to_string(route)
                                  + "_s" + to_string(shard)
                                  + "_r" + to_string(round) + "_" + to_string(batchIdx);
            mergeSubGraphBatch(batch, out, bufElems);
            // Delete this batch now, not after the whole round: peak disk then grows by one
            // batch's output instead of by a second copy of every subgraph.
            for (size_t i = 0; i < batch.size(); ++i) {
                std::remove(batch[i].c_str());
            }
            next.push_back(out);
            ++batchIdx;
        }
        files.swap(next);
        ++round;
    }
}

void GroupGenerator::mergeUnit(size_t route, size_t shard, const std::string & outPath,
                               const std::vector<std::string> & inputFiles,
                               size_t bufElems, size_t maxFanIn,
                               size_t & mergedOut, size_t & ceilingOut,
                               std::vector<uint64_t> & weightHistOut) {
    mergedOut = 0;
    ceilingOut = 0;
    weightHistOut.assign(M_HIST_BUCKETS, 0);

    std::vector<std::string> files = inputFiles;
    reduceSubGraphFanIn(files, maxFanIn, route, shard);

    const size_t streamCnt = files.size();
    std::vector<std::unique_ptr<RelationReader>> readers(streamCnt);
    for (size_t i = 0; i < streamCnt; ++i) {
        readers[i].reset(new RelationReader(files[i], bufElems));
    }

    // Output is relations_*, which makeGroups/mergeBySupport/makeGroupsPhase3 read with
    // ReadBuffer<Relation>. That format stays raw -- only the subGraph_* inputs above change.
    WriteBuffer<Relation> out(outPath, 1024 * 1024);

    // No routing decision here: every edge in these files already belongs to this unit.
    MergeHeap heap;
    for (size_t i = 0; i < streamCnt; ++i) {
        const Relation first = readers[i]->getNext();
        if (!(first == Relation())) {
            MergeHeapEntry entry;
            entry.rel = first;
            entry.stream = i;
            heap.push(entry);
        }
    }

    while (!heap.empty()) {
        const Relation minRelation = heap.top().rel;
        uint16_t totalWeight = 0;
        while (!heap.empty() && heap.top().rel == minRelation) {
            const size_t stream = heap.top().stream;
            addSat16(totalWeight, heap.top().rel.weight);
            heap.pop();
            const Relation next = readers[stream]->getNext();
            if (!(next == Relation())) {
                MergeHeapEntry entry;
                entry.rel = next;
                entry.stream = stream;
                heap.push(entry);
            }
        }
        if (totalWeight == UINT16_MAX) ceilingOut++;
        weightHistOut[mHistBucket(totalWeight)]++;
        mergedOut++;

        // Zeroed in place, not via a helper's return value -- see mergeGraph_one for why.
        Relation rel;
        std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
        rel.id1 = minRelation.id1;
        rel.id2 = minRelation.id2;
        rel.weight = totalWeight;
        out.write(&rel);
    }
    out.flush();

    readers.clear(); // close every input before unlinking it
    for (size_t i = 0; i < streamCnt; ++i) {
        std::remove(files[i].c_str());
    }
}

void GroupGenerator::mergeGraph() {
    cout << "Merging subgraphs" << endl;
    time_t before = time(nullptr);

    // Units merge independently: saveSubGraphToFile already split every flush by
    // (relations_* route, shard), and both are functions of the ids alone, so no pair spans
    // two units. Sharding exists because the routing rule funnels ~88% of all edges into the
    // cross bucket, which would otherwise cap the whole parallel merge at ~1.13x.
    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;
    std::vector<std::pair<size_t, size_t>> units; // (route, shard)
    for (size_t r = 0; r < routeCnt; ++r) {
        const size_t shardCnt = shardsForRoute(r, partitionCnt);
        for (size_t s = 0; s < shardCnt; ++s) {
            units.emplace_back(r, s);
        }
    }
    const size_t unitCnt = units.size();
    const size_t concurrentMergers = std::min(static_cast<size_t>(par.threads), unitCnt);
    const size_t maxFanIn = getMergeFanIn(concurrentMergers);
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->unitStreamCount(), maxFanIn),
                                                     par.ramUsage, concurrentMergers);
    const size_t streamsPerUnit = std::min(this->unitStreamCount(), maxFanIn);

    // A shard writes to a temp; the cross bucket's shards are concatenated afterwards.
    // Concatenation is legal because the relations_* consumers scan the file end to end and
    // never depend on ordering.
    auto unitOutPath = [&](size_t route, size_t shard) {
        return (shardsForRoute(route, partitionCnt) == 1)
             ? (outDir + "/relations_" + to_string(route) + ".bin")
             : (outDir + "/relations_" + to_string(route) + "_s" + to_string(shard) + ".tmp");
    };

    const size_t peakFds = concurrentMergers * mergeFdsPerUnit(std::min(this->unitStreamCount(), maxFanIn));
    cout << "Merge: " << unitCnt << " units (" << routeCnt << " routes, cross bucket sharded x"
         << shardsForRoute(static_cast<size_t>(partitionCnt), partitionCnt) << ") x "
         << this->numOfGraph << " subgraphs, " << concurrentMergers
         << " concurrent (fan-in " << maxFanIn << ", peak fds " << peakFds
         << " of soft limit " << getOpenFileLimit() << ")" << endl;
    if (peakFds >= getOpenFileLimit()) {
        cout << "Warning: peak descriptor use " << peakFds << " meets the soft limit "
             << getOpenFileLimit() << "; raise it or lower --threads." << endl;
    }
    cout << "Merge read buffers: " << streamsPerUnit << " streams x " << mergeBufElems
         << " Relations per unit, "
         << humanBytes(static_cast<uint64_t>(streamsPerUnit) * mergeBufElems
                       * sizeof(Relation) * concurrentMergers)
         << " total" << endl;
    if (this->numOfGraph > maxFanIn) {
        cout << "Folding each unit's " << this->numOfGraph << " subgraphs down to <= "
             << maxFanIn << endl;
    }
    if (mergeBufElems == MERGE_BUFFER_MIN_ELEMS) {
        cout << "Warning: read buffers are at the " << MERGE_BUFFER_MIN_ELEMS
             << "-element floor. Merging will be I/O bound." << endl;
    }

    // Each unit deletes its own inputs as it finishes, so total subgraph bytes are measured
    // up front to keep the [disk] report the single-pass merge used to print.
    uint64_t subGraphBytes = 0;
    for (size_t i = 0; i < this->numOfGraph; ++i) {
        for (size_t u = 0; u < unitCnt; ++u) {
            struct stat st;
            const string name = outDir + "/subGraph_" + to_string(i) + "_"
                              + to_string(units[u].first) + "_" + to_string(units[u].second);
            if (stat(name.c_str(), &st) == 0) {
                subGraphBytes += static_cast<uint64_t>(st.st_size);
            }
        }
    }

    std::vector<size_t> unitMerged(unitCnt, 0);
    std::vector<size_t> unitCeiling(unitCnt, 0);
    std::vector<double> unitSeconds(unitCnt, 0.0);
    std::vector<std::vector<uint64_t>> unitWeightHist(unitCnt);

    // schedule(dynamic): units still carry different loads even after sharding.
    #pragma omp parallel for schedule(dynamic) num_threads(concurrentMergers)
    for (size_t u = 0; u < unitCnt; ++u) {
        const time_t unitStart = time(nullptr);
        mergeUnit(units[u].first, units[u].second, unitOutPath(units[u].first, units[u].second),
                  unitInputFiles(units[u].first, units[u].second),
                  mergeBufElems, maxFanIn, unitMerged[u], unitCeiling[u], unitWeightHist[u]);
        unitSeconds[u] = double(time(nullptr) - unitStart);
    }

    // Fold the cross bucket's shard temps into its relations_*.bin.
    {
        const size_t crossRoute = static_cast<size_t>(partitionCnt);
        const size_t crossShards = shardsForRoute(crossRoute, partitionCnt);
        if (crossShards > 1) {
            const string finalPath = outDir + "/relations_" + to_string(crossRoute) + ".bin";
            ofstream dst(finalPath, std::ios::binary);
            if (!dst.is_open()) {
                cerr << "Error opening file: " << finalPath << endl;
            } else {
                for (size_t s = 0; s < crossShards; ++s) {
                    const string part = unitOutPath(crossRoute, s);
                    ifstream in(part, std::ios::binary);
                    // Streaming an empty rdbuf sets failbit on the destination, so guard it.
                    if (in.is_open() && in.peek() != std::ifstream::traits_type::eof()) {
                        dst << in.rdbuf();
                    }
                }
                dst.close();
                for (size_t s = 0; s < crossShards; ++s) {
                    std::remove(unitOutPath(crossRoute, s).c_str());
                }
            }
        }
    }

    size_t ceilingEdgeCnt = 0; // edges whose merged weight sits at the uint16 ceiling
    size_t mergedEdgeCnt = 0;  // distinct (id1,id2) pairs after merging
    size_t heaviest = 0;
    for (size_t u = 0; u < unitCnt; ++u) {
        mergedEdgeCnt += unitMerged[u];
        ceilingEdgeCnt += unitCeiling[u];
        if (unitMerged[u] > unitMerged[heaviest]) { heaviest = u; }
    }

    // Report the load split by EDGE COUNT, not just by time. On small inputs every unit
    // reports 0 s, which is how a 88%-of-everything unit stayed invisible before.
    {
        std::vector<size_t> order(unitCnt);
        for (size_t u = 0; u < unitCnt; ++u) { order[u] = u; }
        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b) { return unitMerged[a] > unitMerged[b]; });
        cout << "[edges] unit load (top 3 of " << unitCnt << "), max share "
             << (mergedEdgeCnt ? 100.0 * static_cast<double>(unitMerged[heaviest])
                                     / static_cast<double>(mergedEdgeCnt) : 0.0)
             << "%:" << endl;
        for (size_t i = 0; i < order.size() && i < 3; ++i) {
            const size_t u = order[i];
            cout << "[edges]   route " << units[u].first << " shard " << units[u].second
                 << ": " << unitMerged[u] << " edges, " << unitSeconds[u] << " s" << endl;
        }
    }
    cout << "[disk] subGraph: " << (this->numOfGraph * unitCnt) << " files, "
         << humanBytes(subGraphBytes) << " freed" << endl;

    cout << "[edges] merged into " << mergedEdgeCnt << " distinct edges ("
         << humanBytes(static_cast<uint64_t>(mergedEdgeCnt) * sizeof(Relation)) << ")" << endl;

    edgeWeightHist.assign(M_HIST_BUCKETS, 0);
    for (size_t u = 0; u < unitCnt; ++u) {
        if (unitWeightHist[u].size() != M_HIST_BUCKETS) { continue; }
        for (size_t b = 0; b < M_HIST_BUCKETS; ++b) {
            edgeWeightHist[b] += unitWeightHist[u][b];
        }
    }

    if (ceilingEdgeCnt > 0) {
        cout << "Warning: " << ceilingEdgeCnt << " edges reached the weight ceiling ("
             << UINT16_MAX << "); their weights are saturated." << endl;
    } else {
        cout << "Edges at weight ceiling (" << UINT16_MAX << "): 0" << endl;
    }

    cout << "Query relation graph merged successfully" << endl;
    cout << "Time spent: " << double(time(nullptr) - before) << " seconds." << endl;
    return;
}

void GroupGenerator::mergeGraph_one(size_t processedReadCnt) {
    cout << "Merging subgraphs" << endl;
    time_t before = time(nullptr);

    // Subgraphs are now split by route, so this path walks the routes in turn and appends
    // each one's merge to relations.bin. Records are sorted within a route but the file as a
    // whole is segmented by route -- this dump exists for weight-histogram analysis, which
    // does not depend on global ordering. The merge itself stays a linear scan.
    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;
    std::vector<std::pair<size_t, size_t>> units; // (route, shard)
    for (size_t r = 0; r < routeCnt; ++r) {
        const size_t shardCnt = shardsForRoute(r, partitionCnt);
        for (size_t s = 0; s < shardCnt; ++s) {
            units.emplace_back(r, s);
        }
    }
    const size_t maxFanIn = getMergeFanIn(1);
    const size_t mergeBufElems = getMergeBufferElems(std::min(this->unitStreamCount(), maxFanIn),
                                                     par.ramUsage, 1);
    cout << "Merge: " << units.size() << " units x " << this->numOfGraph
         << " subgraphs, sequential (fd soft limit " << getOpenFileLimit()
         << ", fan-in " << maxFanIn << ")" << endl;

    WriteBuffer<Relation> relationLog(outDir + "/relations.bin", 1024 * 1024);

    size_t ceilingEdgeCnt = 0; // edges whose merged weight sits at the uint16 ceiling

    for (size_t u = 0; u < units.size(); ++u) {
        const size_t route = units[u].first;
        const size_t shard = units[u].second;
        std::vector<std::string> subGraphFiles = unitInputFiles(route, shard);
        reduceSubGraphFanIn(subGraphFiles, maxFanIn, route, shard);

        const size_t streamCnt = subGraphFiles.size();
        std::vector<std::unique_ptr<RelationReader>> relationBuffers(streamCnt);
        std::vector<Relation> currentRelations(streamCnt);
        for (size_t i = 0; i < streamCnt; ++i) {
            relationBuffers[i].reset(new RelationReader(subGraphFiles[i], mergeBufElems));
            currentRelations[i] = relationBuffers[i]->getNext();
            if (currentRelations[i] == Relation()) {
                currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
            }
        }

        while (true) {
            Relation minRelation(UINT32_MAX, UINT32_MAX, 0);
            for (size_t i = 0; i < streamCnt; ++i) {
                if (currentRelations[i] < minRelation) {
                    minRelation = currentRelations[i];
                }
            }
            if (minRelation.id1 == UINT32_MAX) break;
            uint16_t totalWeight = 0;
            for (size_t i = 0; i < streamCnt; ++i) {
                if (currentRelations[i] == minRelation) {
                    addSat16(totalWeight, currentRelations[i].weight);
                    currentRelations[i] = relationBuffers[i]->getNext();
                    if (currentRelations[i] == Relation()) {
                        currentRelations[i] = Relation(UINT32_MAX, UINT32_MAX, UINT16_MAX);
                    }
                }
            }
            if (totalWeight == UINT16_MAX) ceilingEdgeCnt++;

            // Zero in place rather than copying from a helper's return value: a return-by-value
            // copy is free to leave the padding indeterminate, and this file is written byte
            // for byte.
            Relation rel;
            std::memset(static_cast<void *>(&rel), 0, sizeof(rel));
            rel.id1 = minRelation.id1;
            rel.id2 = minRelation.id2;
            rel.weight = totalWeight;
            relationLog.write(&rel);
        }

        relationBuffers.clear(); // close every input before unlinking it
        for (size_t i = 0; i < streamCnt; ++i) {
            std::remove(subGraphFiles[i].c_str());
        }
    }
    relationLog.flush();

    // This path exists to dump relations.bin for weight-histogram analysis, so a
    // silent saturation would quietly distort that histogram. Report it.
    if (ceilingEdgeCnt > 0) {
        cout << "Warning: " << ceilingEdgeCnt << " edges reached the weight ceiling ("
             << UINT16_MAX << "); their weights are saturated." << endl;
    } else {
        cout << "Edges at weight ceiling (" << UINT16_MAX << "): 0" << endl;
    }

    cout << "Query relation graph merged successfully" << endl;
    cout << "Time spent: " << double(time(nullptr) - before) << " seconds." << endl;
    return;
}

void GroupGenerator::makeGroups(int groupKmerThr,
                                size_t processedReadCnt,
                                unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo, 
                                vector<uint32_t> &queryGroupInfo) {
    cout << "Creating groups from relation file..." << endl;
    time_t beforeSearch = time(nullptr);

    DisjointSet ds(processedReadCnt);

    auto processFile = [&](const std::string& fname, DisjointSet& subDs) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            uint16_t weight = r.weight;
            if (static_cast<int>(weight) > groupKmerThr) {
                if (id1 == 0 || id2 == 0) continue;
                if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
                subDs.unionSets(id1, id2);
            }
        }
    };


    // Routes, not threads. There are 2*partitionCnt+1 of them and the count no longer follows
    // --threads, so a thread can no longer be identified with a route by index: the loops below
    // hand routes out with `omp for` instead. Union-find is order-independent and the relabelling
    // below keys on the minimum read id, so which worker took which route cannot reach the output.
    //
    // The two loops are kept separate, and the barrier between them is load-bearing. The second
    // one starts from a copy of the global `ds`, so it needs every route of the first group to
    // have been merged in already. Routes 0..partitionCnt-1 hold the same-residue edges and
    // partitionCnt..2*partitionCnt-1 the same-range ones; route 2*partitionCnt is read last,
    // single-threaded, as before.
    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;

    #pragma omp parallel num_threads(par.threads)
    {
        DisjointSet subDs(processedReadCnt);

        #pragma omp for schedule(dynamic)
        for (size_t route = 0; route < static_cast<size_t>(partitionCnt); ++route) {
            processFile(outDir + "/relations_" + std::to_string(route) + ".bin", subDs);
        }

        subDs.flatten();

        #pragma omp critical
        {
            ds += subDs;
        }
    }

    #pragma omp parallel num_threads(par.threads)
    {
        DisjointSet subDs(processedReadCnt);
        #pragma omp critical
        { subDs = ds;}

        #pragma omp for schedule(dynamic)
        for (size_t route = static_cast<size_t>(partitionCnt);
             route < routeCnt - 1; ++route) {
            processFile(outDir + "/relations_" + std::to_string(route) + ".bin", subDs);
        }

        subDs.flatten();

        #pragma omp critical
        { ds += subDs; }
    }

    {
        ReadBuffer<Relation> rb(outDir + "/relations_" + std::to_string(routeCnt - 1) + ".bin",
                                1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            if (static_cast<int>(r.weight) > groupKmerThr) {
                ds.unionSets(r.id1, r.id2);
            }
        }
    }

    // Pass 1: map each component root -> min node ID in that component (canonical label)
    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) {
            rootToMin[root] = i;
        }
    }
    // Pass 2: assign canonical group IDs
    for (uint32_t queryId = 1; queryId < ds.parent.size(); queryId++) {
        if (ds.grouped[queryId]){
            uint32_t groupId = rootToMin[ds.find(queryId)];
            groupInfo[groupId].insert(queryId);
            queryGroupInfo[queryId] = groupId;
        }
    }

    cout << "Query groups created successfully: " << groupInfo.size() << " groups." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

// Phase 2 -- see the declaration in GroupGenerator.h for the rationale.
//
// A "unit" is a Phase-1 core group (labelled by its minimum read id) or, for a read that Phase 1
// left alone, the read itself. Group ids are read ids of grouped reads, so a singleton's own id can
// never collide with a group id. Weak edges are those in (weakThr, coreThr] -- above Phase 3's
// floor, below the core threshold. Pairs where both sides are singletons are skipped: distinct edges
// are merged upstream, so two reads share exactly one edge and can never reach a support of 2. That
// also keeps the counting map bounded to pairs involving at least one multi-read group.
//
// Pass 1 counts weak EDGES per unit pair. With supportRatio == 0 that count, against the floor, is
// the whole rule. With supportRatio > 0 the requirement becomes a fraction of the smaller unit's
// read count and pass 2 recounts the qualifying pairs by DISTINCT READS on the smaller side. Both
// changes are needed together:
//   - the requirement has to scale, because chance links between units u and v grow with
//     |u| * |v|; a fixed count is met automatically once coverage is high.
//   - the count has to be distinct reads, because one repeat-bearing read linked to 50 reads of the
//     other unit yields 50 edges but explains only one read -- exactly the case being excluded.
// Pass 1's edge count is an upper bound on pass 2's distinct-read count, so it is a sound prefilter:
// a pair that fails on edges cannot pass on distinct reads.
//
// Unit sizes are the STATIC Phase-1 sizes, never the disjoint set's evolving component sizes. The
// support map is an unordered_map, so the order in which pairs are merged is unspecified; a
// size test against a changing component would make the output depend on that order.
void GroupGenerator::mergeBySupport(int coreThr,
                                   int weakThr,
                                   size_t processedReadCnt,
                                   unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                   vector<uint32_t> &queryGroupInfo) {
    // Both of these were parameters (--merge-support-ratio, --merge-max-unit-reads) and are now
    // fixed at the values every measured result used. The ratio form of the support requirement was
    // never measured against the plain count outside noise, and the size gate was off in every run
    // that produced a published number; keeping them as knobs invited sweeps over axes that had
    // already been decided. The branches they select are left in place so the reasoning stays
    // readable -- the compiler removes them.
    const float supportRatio = 0.0f;
    const size_t maxUnitReads = 0;
    const bool useRatio = (supportRatio > 0.0f);
    if (useRatio) {
        cout << "Phase 2: merging units whose smaller side has >= max(" << MERGE_SUPPORT_FLOOR
             << ", " << supportRatio << " x its read count) distinct reads carrying a weak link"
             << " (weight in (" << weakThr << ", " << coreThr << "])..." << endl;
    } else {
        cout << "Phase 2: merging units with >= " << MERGE_SUPPORT_FLOOR
             << " weak links (weight in (" << weakThr << ", " << coreThr << "])..." << endl;
    }
    time_t t0 = time(nullptr);

    const size_t groupsBefore = groupInfo.size();

    // unit[i]: Phase-1 group id, or the read's own id when Phase 1 left it alone.
    std::vector<uint32_t> unit(processedReadCnt + 1, 0);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        unit[i] = queryGroupInfo[i] ? queryGroupInfo[i] : i;
    }

    // The disjoint set starts from the Phase-1 components and the supported pairs are added to it
    // further down. It is built here, before the counting passes, because `compSize` is indexed by
    // its roots.
    DisjointSet ds(processedReadCnt);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        if (queryGroupInfo[i] != 0 && queryGroupInfo[i] != i) {
            ds.unionSets(i, queryGroupInfo[i]);
        }
    }

    // Reads per component, indexed by disjoint-set root. One array serves both roles: until the
    // merge loop runs, every root is still a Phase-1 root, so a lookup returns the static unit size
    // the ratio rule needs; from then on the loop keeps it current so the size gate can see what a
    // component has already grown to. Keeping two arrays instead would cost another 4 B per read,
    // which is ~666 MB on a 166 M-read run.
    // The order matters and is load-bearing: `requiredSupport` and `smallSideOf` must see static
    // sizes, or the support rule would depend on the merge order. Both are only called while
    // building `candidates`, which finishes before the merge loop starts.
    std::vector<uint32_t> compSize(processedReadCnt + 1, 0u);
    for (uint32_t i = 1; i <= processedReadCnt; ++i) { compSize[ds.find(i)]++; }
    const auto sizeOf = [&](uint32_t node) -> uint32_t { return compSize[ds.find(node)]; };

    // Support counting, exact, and independent of every knob that exists for memory or speed.
    //
    // This used to be a capped map. The cap made the answer depend on how the edges happened to be
    // divided up: relations_* is partitioned by (id1, id2) -- READ ids -- while support is keyed on
    // the UNIT pair, so one unit pair's weak edges are scattered across many partitions. Each
    // partition counted its own share against its own cap, so overflow dropped different pairs
    // depending on how many partitions there were and how many readers shared them. First that was
    // --threads (16 threads scored species purity 0.932 on strain-madness, 64 scored 0.714 on a
    // byte-identical k-mer set); sizing the caps per route only moved the dependency to
    // --partitions. A device that exists to bound memory must not change the result, so the cap is
    // gone.
    //
    // What replaces it: shard by the support key itself. hash(key) sends every occurrence of a
    // pair to exactly one shard, so a shard can be counted to completion in memory with no cap at
    // all, and the shard count changes only how much memory and disk the counting uses -- never
    // which pairs are counted or how far. Shards are 8 bytes per weak edge on disk, against the
    // 117 GB of relations_* already there.
    std::unordered_map<uint64_t, uint32_t> support;
    size_t weakEdgeCnt = 0;

    auto packPair = [](uint32_t a, uint32_t b) -> uint64_t {
        return (a < b) ? ((uint64_t)a << 32 | b) : ((uint64_t)b << 32 | a);
    };

    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;
    const auto routePath = [&](size_t route) {
        return outDir + "/relations_" + std::to_string(route) + ".bin";
    };

    // How many shards. relations_* records are a fixed 12 bytes, so the file sizes give an exact
    // upper bound on the weak-edge count without reading anything, and distinct keys can never
    // outnumber edges. One shard whenever the whole thing fits -- the common case, and it skips
    // the shard files entirely.
    size_t upperBoundEdges = 0;
    for (size_t r = 0; r < routeCnt; ++r) {
        upperBoundEdges += static_cast<size_t>(FileUtil::getFileSize(routePath(r)))
                           / sizeof(Relation);
    }
    const size_t entryBudget = getRelationBudget(par.ramUsage);
    size_t shardCnt = (upperBoundEdges + entryBudget - 1) / std::max<size_t>(1, entryBudget);
    if (shardCnt < 1) { shardCnt = 1; }
    // The bound above counts every record in relations_*, but only the ones inside the weak band
    // are keyed here -- on strain-madness that was 450 M of 10.49 G, so the estimate is roughly
    // 20x conservative and the shard count can be trimmed hard without a shard ever coming close
    // to the budget. Trimmed by what phase A can hold open: every writer keeps one file per shard.
    const size_t fdBudget = getOpenFileLimit() / 2;
    const size_t writers = static_cast<size_t>(std::max(1, par.threads));
    const size_t shardFdCap = std::max<size_t>(1, fdBudget / writers);
    if (shardCnt > shardFdCap) {
        cout << "Phase 2: " << shardCnt << " key shards would need " << (shardCnt * writers)
             << " open files; trimming to " << shardFdCap << " (limit "
             << getOpenFileLimit() << "). The estimate counts every relations_* record, not just"
             << " the weak band, so it is far above what the shards actually hold." << endl;
        shardCnt = shardFdCap;
    }
    // Write buffers are charged against a fixed total per writer rather than a fixed size each,
    // so raising the shard count cannot multiply the memory it uses.
    const size_t SHARD_WRITE_BYTES_PER_WRITER = 16u << 20;
    const size_t shardBufElems = std::max<size_t>(
        4096, SHARD_WRITE_BYTES_PER_WRITER / (shardCnt * sizeof(uint64_t)));

    // Per-route weak-edge counts, for the load-skew report. Instrumentation only.
    std::vector<size_t> routeWeak(routeCnt, 0);

    // The filter every pass applies to a Relation, in one place so the two passes cannot drift.
    // Returns false when the edge carries no support for any pair.
    const auto weakPairOf = [&](const Relation & r, uint64_t & keyOut) -> bool {
        const uint32_t id1 = r.id1, id2 = r.id2;
        if (id1 == 0 || id2 == 0) { return false; }
        if (id1 > processedReadCnt || id2 > processedReadCnt) { return false; }
        const int w = static_cast<int>(r.weight);
        if (w <= weakThr || w > coreThr) { return false; }
        if (queryGroupInfo[id1] == 0 && queryGroupInfo[id2] == 0) { return false; }
        const uint32_t u = unit[id1], v = unit[id2];
        if (u == v) { return false; }
        keyOut = packPair(u, v);
        return true;
    };

    // Splitmix64. Any mixer would do; what matters is that it depends on the key alone, so a
    // pair's shard is the same on every run and in every configuration.
    const auto shardOfKey = [shardCnt](uint64_t key) -> size_t {
        if (shardCnt == 1) { return 0; }
        uint64_t x = key + 0x9E3779B97F4A7C15ull;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
        x ^= x >> 31;
        return static_cast<size_t>(x % shardCnt);
    };

    if (shardCnt == 1) {
        // Everything fits. Threads count into private maps and the merge sums them, which is
        // exact because nothing is discarded: `support[key] += n` recovers the total however the
        // occurrences were spread across readers. A key can appear under several routes, so the
        // floor cannot be applied per thread -- it is applied to the merged counts below.
        #pragma omp parallel for schedule(dynamic) num_threads(par.threads)
        for (size_t route = 0; route < routeCnt; ++route) {
            std::unordered_map<uint64_t, uint32_t> local;
            size_t localWeak = 0;
            uint64_t key = 0;
            ReadBuffer<Relation> rb(routePath(route), 1024 * 1024);
            for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
                if (!weakPairOf(r, key)) { continue; }
                localWeak++;
                local[key]++;
            }
            routeWeak[route] = localWeak;
            #pragma omp critical
            {
                for (const auto& kv : local) { support[kv.first] += kv.second; }
                weakEdgeCnt += localWeak;
            }
        }
        // Same filter the sharded path applies, so `support` means the same thing either way:
        // pairs that could still be merged. requiredSupport never returns less than the floor,
        // so nothing droppable here could have qualified.
        for (auto it = support.begin(); it != support.end(); ) {
            if (it->second < MERGE_SUPPORT_FLOOR) { it = support.erase(it); } else { ++it; }
        }
    } else {
        // Two passes. Scatter keys to shard files, then count one shard at a time. Each writer
        // owns its own file so no lock is needed on the hot path; the reader takes every writer's
        // file for its shard.
        const auto shardPath = [&](size_t shard, size_t writer) {
            return outDir + "/support_" + std::to_string(shard) + "_" + std::to_string(writer);
        };
        std::vector<std::string> shardFiles;

        #pragma omp parallel num_threads(par.threads)
        {
            const int writer = omp_get_thread_num();
            std::vector<std::unique_ptr<WriteBuffer<uint64_t>>> out(shardCnt);
            for (size_t sh = 0; sh < shardCnt; ++sh) {
                out[sh].reset(new WriteBuffer<uint64_t>(shardPath(sh, writer), shardBufElems));
            }
            size_t localWeak = 0;

            #pragma omp for schedule(dynamic)
            for (size_t route = 0; route < routeCnt; ++route) {
                size_t routeLocalWeak = 0;
                uint64_t key = 0;
                ReadBuffer<Relation> rb(routePath(route), 1024 * 1024);
                for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
                    if (!weakPairOf(r, key)) { continue; }
                    routeLocalWeak++;
                    out[shardOfKey(key)]->write(&key, 1);
                }
                routeWeak[route] = routeLocalWeak;
                localWeak += routeLocalWeak;
            }

            for (size_t sh = 0; sh < shardCnt; ++sh) {
                out[sh]->flush();
                out[sh].reset();
            }
            #pragma omp critical
            {
                weakEdgeCnt += localWeak;
                for (size_t sh = 0; sh < shardCnt; ++sh) {
                    shardFiles.push_back(shardPath(sh, writer));
                }
            }
        }

        cout << "Phase 2: " << weakEdgeCnt << " weak edges scattered over " << shardCnt
             << " key shards x " << writers << " writers ("
             << humanBytes(weakEdgeCnt * sizeof(uint64_t)) << " on disk, "
             << shardBufElems << " elems per write buffer)" << endl;

        // One shard at a time so peak memory is one shard's distinct keys, whatever shardCnt is.
        // Only pairs that could still matter are carried out of the shard: a pair below the floor
        // can never be merged, so keeping it would only grow `support` for nothing.
        for (size_t sh = 0; sh < shardCnt; ++sh) {
            std::unordered_map<uint64_t, uint32_t> shardCount;
            for (int writer = 0; writer < par.threads; ++writer) {
                ReadBuffer<uint64_t> rb(shardPath(sh, writer), 1024 * 1024);
                for (uint64_t key = rb.getNext(); key != 0; key = rb.getNext()) {
                    shardCount[key]++;
                }
            }
            for (const auto& kv : shardCount) {
                if (kv.second >= MERGE_SUPPORT_FLOOR) { support.emplace(kv.first, kv.second); }
            }
        }
        reportAndRemoveFiles(shardFiles, "support shards");
    }

    // Where the weak edges are, what each route was allowed to hold, and whether that was
    // enough. The totals below say how much was dropped but not that one route holds most of
    // it, which is the fact that used to make the cap thread-dependent.
    {
        std::vector<size_t> order(routeCnt);
        for (size_t r = 0; r < routeCnt; ++r) { order[r] = r; }
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (routeWeak[a] != routeWeak[b]) { return routeWeak[a] > routeWeak[b]; }
            return a < b;
        });
        const size_t shown = std::min<size_t>(3, routeCnt);
        cout << "Phase 2: weak edges by route, top " << shown << " of " << routeCnt;
        for (size_t i = 0; i < shown; ++i) {
            const size_t r = order[i];
            cout << " | route " << r << ": " << routeWeak[r];
            if (weakEdgeCnt > 0) {
                cout << " (" << fixed << setprecision(1)
                     << (100.0 * static_cast<double>(routeWeak[r])
                             / static_cast<double>(weakEdgeCnt))
                     << "%)";
            }
        }
        cout << endl;
        // Load skew only. There is nothing to drop any more, so this is about where the work is,
        // not about what was lost -- the cross bucket carrying ~88% is why the merge shards it.
    }

    // Support required of a pair. Without the ratio this is the bare floor, which is the rule the
    // benchmark operating point was measured with.
    auto requiredSupport = [&](uint32_t u, uint32_t v) -> uint32_t {
        if (!useRatio) { return MERGE_SUPPORT_FLOOR; }
        const uint32_t smaller = std::min(sizeOf(u), sizeOf(v));
        const double need = std::ceil(static_cast<double>(supportRatio) * static_cast<double>(smaller));
        const uint32_t needed = (need >= static_cast<double>(UINT32_MAX))
            ? UINT32_MAX : static_cast<uint32_t>(need);
        return std::max(MERGE_SUPPORT_FLOOR, needed);
    };

    // Which side of a pair the distinct reads are counted on. Sizes decide it; equal sizes fall back
    // to the smaller unit id so the choice does not depend on map order.
    auto smallSideOf = [&](uint32_t u, uint32_t v) -> uint32_t {
        const uint32_t su = sizeOf(u), sv = sizeOf(v);
        if (su != sv) { return (su < sv) ? u : v; }
        return std::min(u, v);
    };

    // Pairs that cleared pass 1. With the ratio on, this is only a prefilter -- pass 2 recounts them
    // by distinct reads, which can only be smaller.
    // The size gate is applied here as well as in the merge loop. Here it is the cheap half: a pair
    // whose Phase-1 units are already too big can never qualify, whatever happens later, so dropping
    // it now also keeps it out of pass 2's distinct-read counting.
    std::unordered_map<uint64_t, uint32_t> candidates; // pair -> required support
    size_t gatedStatic = 0, gatedStaticReads = 0;
    for (const auto& kv : support) {
        const uint32_t u = static_cast<uint32_t>(kv.first >> 32);
        const uint32_t v = static_cast<uint32_t>(kv.first & 0xFFFFFFFFull);
        const uint32_t need = requiredSupport(u, v);
        if (kv.second < need) { continue; }
        if (maxUnitReads > 0
            && (sizeOf(u) > maxUnitReads || sizeOf(v) > maxUnitReads)) {
            gatedStatic++;
            gatedStaticReads += static_cast<size_t>(sizeOf(u)) + sizeOf(v);
            continue;
        }
        candidates.emplace(kv.first, need);
    }
    if (maxUnitReads > 0) {
        cout << "Phase 2: size gate " << maxUnitReads << " reads/unit dropped " << gatedStatic
             << " supported pairs spanning " << gatedStaticReads << " unit-reads" << endl;
    }

    // Pass 2: distinct reads on the smaller side, for candidate pairs only.
    std::unordered_set<uint64_t> satisfied;
    if (useRatio && !candidates.empty()) {
        // No cap. This pass counts distinct reads per candidate pair, and a cap here would put
        // the same dependency back that the support pass just lost: which reads get dropped would
        // be decided by how the edges were partitioned. `candidates` is already the pairs that
        // cleared the support floor, so the read sets are bounded by that -- far smaller than the
        // edge stream. If this ever does not fit, shard it by key as the support pass does; do not
        // reintroduce a cap.
        // A pair is satisfied as soon as `required` distinct reads are seen, so a reader that
        // stops inserting at that point cannot make any pair's union fall short: no reader stops
        // earlier.
        auto distinctRoute = [&](size_t route,
                                 std::unordered_map<uint64_t, std::unordered_set<uint32_t>>& local,
                                 std::unordered_set<uint64_t>& localSatisfied,
                                 size_t& held) {
            ReadBuffer<Relation> rb(routePath(route), 1024 * 1024);
            for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
                const uint32_t id1 = r.id1, id2 = r.id2;
                if (id1 == 0 || id2 == 0) continue;
                if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
                const int w = static_cast<int>(r.weight);
                if (w <= weakThr || w > coreThr) continue;
                if (queryGroupInfo[id1] == 0 && queryGroupInfo[id2] == 0) continue;
                const uint32_t u = unit[id1], v = unit[id2];
                if (u == v) continue;
                const uint64_t key = packPair(u, v);
                const auto cand = candidates.find(key);
                if (cand == candidates.end()) continue;
                if (localSatisfied.count(key)) continue;
                const uint32_t small = smallSideOf(u, v);
                const uint32_t readId = (unit[id1] == small) ? id1 : id2;
                auto& reads = local[key];
                if (reads.insert(readId).second) { held++; }
                if (reads.size() >= cand->second) {
                    localSatisfied.insert(key);
                    held -= reads.size();
                    local.erase(key);
                }
            }
        };

        size_t heldReads = 0;
        std::unordered_map<uint64_t, std::unordered_set<uint32_t>> partial;
        // Every route in one loop, the cross bucket included. It used to be read separately and
        // last, against the whole PASS2_READ_CAP, while the others shared it by thread; both of
        // those made the outcome depend on --threads.
        #pragma omp parallel for schedule(dynamic) num_threads(par.threads)
        for (size_t route = 0; route < routeCnt; ++route) {
            std::unordered_map<uint64_t, std::unordered_set<uint32_t>> local;
            std::unordered_set<uint64_t> localSatisfied;
            size_t localHeld = 0;

            distinctRoute(route, local, localSatisfied, localHeld);

            #pragma omp critical
            {
                satisfied.insert(localSatisfied.begin(), localSatisfied.end());
                for (const auto& kv : local) {
                    partial[kv.first].insert(kv.second.begin(), kv.second.end());
                }
                heldReads += localHeld;
            }
        }
        // Unions can have crossed the requirement even when no single route's share did.
        for (auto it = partial.begin(); it != partial.end(); ) {
            if (it->second.size() >= candidates[it->first]) {
                satisfied.insert(it->first);
                it = partial.erase(it);
            } else {
                ++it;
            }
        }

        cout << "Phase 2: " << candidates.size() << " candidate pairs, " << satisfied.size()
             << " reached the distinct-read requirement (" << heldReads
             << " reads still held for the rest)" << endl;
    }

    // Merge the supported pairs into the Phase-1 components built above.
    //
    // The order is fixed on purpose. The gate below asks how big a component already is, so with
    // the plain `unordered_map` walk the outcome would depend on the map's bucket layout. Sorting
    // by support puts the best-evidenced pairs first, which is also where the size budget should go.
    const size_t LOG_BUCKETS = 8;
    std::vector<std::pair<uint64_t, uint32_t>> ordered; // (pair key, support observed in pass 1)
    ordered.reserve(candidates.size());
    for (const auto& kv : candidates) {
        if (useRatio && satisfied.count(kv.first) == 0) continue;
        const auto found = support.find(kv.first);
        ordered.emplace_back(kv.first, (found == support.end()) ? 0u : found->second);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const std::pair<uint64_t, uint32_t>& a, const std::pair<uint64_t, uint32_t>& b) {
                  if (a.second != b.second) { return a.second > b.second; }
                  return a.first < b.first;
              });

    // Units touched by an accepted merge. A bitmap, not a set of ids: the id space is the read
    // space, so this is one bit per read where a hash set would be tens of bytes per merged pair.
    std::vector<bool> unitMerged(processedReadCnt + 1, false);
    std::vector<size_t> mergeSizeHist(LOG_BUCKETS, 0); // sizes the gate actually saw, both sides
    size_t mergedPairs = 0, gatedDynamic = 0;

    const auto bucketOf = [LOG_BUCKETS](uint32_t sz) -> size_t {
        // bucket 0 is the singleton component; from 2 on it is 1 + floor(log10(sz)), capped.
        if (sz < 2) { return 0; }
        size_t bucket = 1;
        for (uint32_t bound = 10; bucket + 1 < LOG_BUCKETS && sz >= bound; ++bucket) {
            bound *= 10;
        }
        return bucket;
    };

    for (const auto& entry : ordered) {
        const uint32_t u = static_cast<uint32_t>(entry.first >> 32);
        const uint32_t v = static_cast<uint32_t>(entry.first & 0xFFFFFFFFull);
        const uint32_t ru = ds.find(u), rv = ds.find(v);
        if (ru == rv) { continue; } // an earlier merge already joined them
        // Re-checked here, not just on the Phase-1 sizes: a component that grew past the bound has
        // to stop growing, and that is what makes the bound hold. Both sides are at most
        // maxUnitReads, so no component can end up with 2 * maxUnitReads reads or more, and one
        // that reaches the bound is out of the running for every remaining pair.
        if (maxUnitReads > 0
            && (compSize[ru] > maxUnitReads || compSize[rv] > maxUnitReads)) {
            gatedDynamic++;
            continue;
        }
        mergeSizeHist[bucketOf(compSize[ru])]++;
        mergeSizeHist[bucketOf(compSize[rv])]++;
        const uint64_t merged = static_cast<uint64_t>(compSize[ru]) + compSize[rv];
        ds.unionSets(u, v);
        compSize[ds.find(u)] = static_cast<uint32_t>(merged); // read ids are uint32, so is any sum
        unitMerged[u] = true;
        unitMerged[v] = true;
        mergedPairs++;
    }

    // A Phase-1 group of exactly one read (its own label) never entered the disjoint set above, so
    // re-mark those reads as grouped to keep them out of Phase 3's singleton pass.
    for (uint32_t i = 1; i <= processedReadCnt; ++i) {
        if (queryGroupInfo[i] != 0) ds.grouped[i] = true;
    }

    // Relabel every component by its minimum read id, as Phase 1 does.
    groupInfo.clear();
    std::fill(queryGroupInfo.begin(), queryGroupInfo.end(), 0u);
    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        const uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) rootToMin[root] = i;
    }
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        const uint32_t groupId = rootToMin[ds.find(i)];
        groupInfo[groupId].insert(i);
        queryGroupInfo[i] = groupId;
    }

    cout << "Phase 2: " << weakEdgeCnt << " weak edges over " << support.size()
         << " unit pairs at or above the support floor, " << mergedPairs << " pairs merged"
         << "; groups " << groupsBefore << " -> " << groupInfo.size()
         << " [counted exactly, no cap; independent of --threads and --partitions]" << endl;

    // Read mass this phase moved, and the sizes it moved it from. Pair counts alone hide the damage:
    // a few merges between large units cost more purity than many merges between small ones, because
    // every read of a wrongly merged unit becomes impure. The share below is the worst case that
    // purity can lose here, so it is readable before any grading run.
    {
        // Reads sitting in a component that this phase enlarged. Counted per component, from the
        // final sizes, so a component merged several times is counted once.
        std::unordered_set<uint32_t> mergedRoots;
        size_t mergedUnits = 0;
        for (uint32_t i = 1; i <= processedReadCnt; ++i) {
            if (!unitMerged[i]) continue;
            mergedUnits++;
            mergedRoots.insert(ds.find(i));
        }
        size_t mergedReads = 0, maxMergedComp = 0;
        for (const uint32_t root : mergedRoots) {
            mergedReads += compSize[root];
            if (compSize[root] > maxMergedComp) { maxMergedComp = compSize[root]; }
        }

        size_t groupedReads = 0, maxGroupSize = 0;
        std::vector<size_t> groupSizes;
        groupSizes.reserve(groupInfo.size());
        for (const auto& g : groupInfo) {
            groupSizes.push_back(g.second.size());
            groupedReads += g.second.size();
            if (g.second.size() > maxGroupSize) { maxGroupSize = g.second.size(); }
        }
        const size_t topN = std::min<size_t>(5, groupSizes.size());
        std::partial_sort(groupSizes.begin(), groupSizes.begin() + topN, groupSizes.end(),
                          std::greater<size_t>());

        const double massShare = groupedReads
            ? (100.0 * static_cast<double>(mergedReads) / static_cast<double>(groupedReads)) : 0.0;
        cout << "Phase 2: " << mergedUnits << " units merged, holding " << mergedReads
             << " reads (" << massShare << "% of grouped reads -- upper bound on the purity this"
             << " phase can cost)" << endl;

        // Sizes as the gate saw them: one entry per side of every accepted merge, taken from the
        // component size at that moment. This is the distribution --merge-max-unit-reads cuts, so
        // it is the table to read when picking a value for it.
        cout << "Phase 2: merged component sizes";
        static const char* BUCKET_NAME[] = {
            "1", "2-9", "10-99", "100-999", "1e3-1e4", "1e4-1e5", "1e5-1e6", "1e6+"
        };
        for (size_t b = 0; b < LOG_BUCKETS; ++b) {
            cout << " | " << BUCKET_NAME[b] << ": " << mergeSizeHist[b];
        }
        cout << endl;
        if (maxUnitReads > 0) {
            cout << "Phase 2: size gate also blocked " << gatedDynamic
                 << " pairs at merge time (component had grown past the bound)" << endl;
        }

        // Two different numbers on purpose. The first covers every group, including Phase-1 groups
        // this phase never touched -- those can be larger than the bound and are left alone. The
        // second is the one the bound governs: no component this phase enlarged reaches
        // 2 * maxUnitReads reads.
        cout << "Phase 2: max group size " << maxGroupSize << ", top " << topN << ":";
        for (size_t i = 0; i < topN; ++i) { cout << " " << groupSizes[i]; }
        cout << endl;
        cout << "Phase 2: largest component this phase enlarged: " << maxMergedComp;
        if (maxUnitReads > 0) {
            cout << " (bound 2 x " << maxUnitReads << " = " << (2 * maxUnitReads) << ")";
        }
        cout << endl;
    }

    cout << "Phase 2 done: " << double(time(nullptr) - t0) << " s" << endl;
}

void GroupGenerator::makeGroupsPhase3(
    int groupKmerThr,
    size_t processedReadCnt,
    const std::vector<bool>& isSingleton,
    unordered_map<uint32_t, unordered_set<uint32_t>>& groupInfo,
    vector<uint32_t>& queryGroupInfo)
{
    cout << "Phase 3: linking singletons..." << endl;
    time_t t0 = time(nullptr);

    DisjointSet ds(processedReadCnt);

    auto processFile = [&](const std::string& fname, DisjointSet& subDs) {
        ReadBuffer<Relation> rb(fname, 1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            uint32_t id1 = r.id1, id2 = r.id2;
            if (id1 == 0 || id2 == 0) continue;
            if (id1 > processedReadCnt || id2 > processedReadCnt) continue;
            if (!isSingleton[id1] || !isSingleton[id2]) continue;
            if (static_cast<int>(r.weight) > groupKmerThr) {
                subDs.unionSets(id1, id2);
            }
        }
    };

    // Routes, handed out by `omp for`, for the same reason as in makeGroups: the route count is
    // 2*partitionCnt+1 and no longer equals the thread count, so a thread index is not a route
    // index. The barrier between the two loops stays -- the second starts from a copy of `ds`.
    const size_t routeCnt = static_cast<size_t>(partitionCnt) * 2 + 1;

    #pragma omp parallel num_threads(par.threads)
    {
        DisjointSet subDs(processedReadCnt);
        #pragma omp for schedule(dynamic)
        for (size_t route = 0; route < static_cast<size_t>(partitionCnt); ++route) {
            processFile(outDir + "/relations_" + std::to_string(route) + ".bin", subDs);
        }
        subDs.flatten();
        #pragma omp critical
        { ds += subDs; }
    }

    #pragma omp parallel num_threads(par.threads)
    {
        DisjointSet subDs(processedReadCnt);
        #pragma omp critical
        { subDs = ds; }
        #pragma omp for schedule(dynamic)
        for (size_t route = static_cast<size_t>(partitionCnt);
             route < routeCnt - 1; ++route) {
            processFile(outDir + "/relations_" + std::to_string(route) + ".bin", subDs);
        }
        subDs.flatten();
        #pragma omp critical
        { ds += subDs; }
    }

    {
        ReadBuffer<Relation> rb(outDir + "/relations_" + std::to_string(routeCnt - 1) + ".bin",
                                1024 * 1024);
        for (Relation r = rb.getNext(); !(r == Relation()); r = rb.getNext()) {
            if (r.id1 == 0 || r.id2 == 0) continue;
            if (r.id1 > processedReadCnt || r.id2 > processedReadCnt) continue;
            if (!isSingleton[r.id1] || !isSingleton[r.id2]) continue;
            if (static_cast<int>(r.weight) > groupKmerThr) {
                ds.unionSets(r.id1, r.id2);
            }
        }
    }

    std::unordered_map<uint32_t, uint32_t> rootToMin;
    for (uint32_t i = 1; i < ds.parent.size(); ++i) {
        if (!ds.grouped[i]) continue;
        uint32_t root = ds.find(i);
        auto it = rootToMin.find(root);
        if (it == rootToMin.end() || i < it->second) rootToMin[root] = i;
    }
    for (uint32_t queryId = 1; queryId < ds.parent.size(); queryId++) {
        if (ds.grouped[queryId]) {
            uint32_t groupId = rootToMin[ds.find(queryId)];
            groupInfo[groupId].insert(queryId);
            queryGroupInfo[queryId] = groupId;
        }
    }

    cout << "Phase 3 done: " << double(time(nullptr) - t0) << " s" << endl;
}

void GroupGenerator::saveGroupsToFile(const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                      const vector<uint32_t> &queryGroupInfo) {
    // save group in txt file
    const string& groupInfoFileName = outDir + "/groups";
    ofstream outFile1(groupInfoFileName);
    if (!outFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    for (const auto& [groupId, queryIds] : groupInfo) {
        outFile1 << groupId << "\t";
        for (const auto& queryId : queryIds) {
            outFile1 << queryId << "\t";
        }
        outFile1 << endl;
    }
    outFile1.close();
    cout << "Query group saved to " << groupInfoFileName << " successfully." << endl;
    

    const string& queryGroupInfoFileName = outDir + "/groupMap";
    ofstream outFile2(queryGroupInfoFileName);
    if (!outFile2.is_open()) {
        cerr << "Error opening file: " << queryGroupInfoFileName << endl;
        return;
    }

    for (size_t i = 1; i < queryGroupInfo.size(); ++i) {
        outFile2 << i << "\t" << queryGroupInfo[i] << "\n";
    }
    outFile2.close();
    cout << "Query group saved to " << queryGroupInfoFileName << " successfully." << endl;

    return;
}