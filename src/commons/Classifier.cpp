#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"

#include <variant>

Classifier::Classifier(LocalParameters & par) : par(par) {
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);
    kmerFormat = par.kmerFormat;
    if (par.dbTotalLength == 0) {
        par.dbTotalLength = readDbSize(par.filenames[1 + (par.seqMode == 2)]);
    }

    cout << "Database name : " << par.dbName << endl;
    cout << "Creation date : " << par.dbDate << endl;
    
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    if (par.precisionMode != 0) {
        preciseModePreset(par);
    }

    if (kmerFormat == 1) {
        metamerPattern = new LegacyPattern(std::make_unique<RegularGeneticCode>(), 8);
        cerr << "Warning: the specified database uses an old k-mer format." << endl;
        cerr << "   E-value calculation is not supported." << endl;
        par.maxEValue = 0;
    } else if (!par.customMetamer.empty()) {
        int codeNum = getCodeNum(par.customMetamer);
        if (codeNum == 1) {
            if (par.spaceMask.empty()) {
                cout << "Using SingleCodePattern with custom metamer." << endl;
                metamerPattern = new SingleCodePattern(par.customMetamer);
            } else {
                cout << "Using SpacedPattern with custom metamer." << endl;
                uint32_t mask = parseMask(par.spaceMask.c_str());
                metamerPattern = new SpacedPattern(par.customMetamer, mask);
            }
        } else if (codeNum > 1) {
            metamerPattern = new MultiCodePattern(par.customMetamer);
        }
    } else if (kmerFormat == 2) {
        if (par.spaceMask.empty()) {
            cout << "Using SingleCodePattern with RegularGeneticCode" << endl;
            metamerPattern = new SingleCodePattern(std::make_unique<RegularGeneticCode>(), 8);
        } else {
            cout << "Using SpacedPattern with RegularGeneticCode" << endl;
            uint32_t mask = parseMask(par.spaceMask.c_str());
            metamerPattern = new SpacedPattern(std::make_unique<RegularGeneticCode>(), __builtin_popcount(mask), mask);
        }
    }
    if (par.storeKmerPos) {
        C_LOG2_C.resize(256);
        for (int i = 1; i < 256; ++i) {
            C_LOG2_C[i] = i * std::log2(static_cast<double>(i));
        }
    }
    kmerExtractor = new KmerExtractor(par, metamerPattern);
    queryIndexer = new QueryIndexer(par);
    kmerMatcher = new KmerMatcher(par, taxonomy, metamerPattern);
    reporter = new Reporter(par, taxonomy);
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete reporter;
    delete geneticCode;
    delete metamerPattern;
    if (mappingResList) delete[] mappingResList;
}

void Classifier::preciseModePreset(LocalParameters & par) {
    uint32_t mask = parseMask(par.spaceMask.c_str());
    size_t windowSize = par.spaceMask.length();
    size_t kmerLen = __builtin_popcount(mask);

    float minScoreCp = par.minScore;
    float minSpScoreCp = par.minSpScore;
    double maxEValueCp = par.maxEValue;

    if (par.precisionMode == 1) { // short-read preset
        if (windowSize == kmerLen || kmerLen == 0) {
            std::cout << "Using short-read presets for contiguous k-mer search: " << std::endl;
            std::cout << "   --min-score 0.15 --min-sp-score 0.5 -e 0.001" << std::endl;
            par.minScore = 0.15;
            par.minSpScore = 0.5;
            par.maxEValue = 0.001;
        } else {
            std::cout << "Using short-read presets for spaced k-mer search: " << std::endl;
            std::cout << "   --min-score 0.2 --min-sp-score 0.6 -e 0.001" << std::endl;
            par.minScore = 0.2;
            par.minSpScore = 0.6;
            par.maxEValue = 0.001;
        }
    } else if (par.precisionMode == 2) { // HiFi long-read preset
        std::cout << "Using HiFi long-read presets: " << std::endl;
        std::cout << "   --min-score 0.07 --min-sp-score 0.3 -e 0.001" << std::endl;
        par.minScore = 0.07;
        par.minSpScore = 0.3;
        par.maxEValue = 0.001;
    }

    if (minScoreCp != 0 && minScoreCp != par.minScore) {
        std::cout << "Overriding preset --min-score " << par.minScore << " with user specified value " << minScoreCp << std::endl;
        par.minScore = minScoreCp;
    }
    if (minSpScoreCp != 0 && minSpScoreCp != par.minSpScore) {
        std::cout << "Overriding preset --min-sp-score " << par.minSpScore << " with user specified value " << minSpScoreCp << std::endl;
        par.minSpScore = minSpScoreCp;
    }
    if (maxEValueCp != 1 && maxEValueCp != par.maxEValue) {
        std::cout << "Overriding preset --e-value " << par.maxEValue << " with user specified value " << maxEValueCp << std::endl;
        par.maxEValue = maxEValueCp;
    }

}

uint64_t Classifier::calculateBufferSize(
    uint64_t queryListSize,
    uint64_t matchPerKmer) 
{
    size_t totalBytes = (size_t) par.ramUsage * 1024 * 1024 * 1024;

    size_t bytesPerThread = 
        1024 * 1024 * sizeof(Match) + // local match buffer 32 MB
        1024 * 1024 * 4 * sizeof(Kmer) +  // DeltaIdxReader::valueBuffer 64 MB
        1024 * 1024 * 4 * sizeof(uint16_t) +   // DeltaIdxReader::deltaBuffer 2 MB
        1024 * 1024 * 4 * sizeof(uint32_t);    // DeltaIdxReader::infoBuffer 4 MB

    size_t overhead = 128 * 1024 * 1024; // 128MB
    size_t queryListBytes = queryListSize * (sizeof(Query) + 150); //  104,857,600
    size_t availableBytes = totalBytes - (par.threads * bytesPerThread) - overhead - queryListBytes; 
    size_t bytePerKmer = sizeof(Kmer) + matchPerKmer * sizeof(Match);
    size_t totalSize = availableBytes / bytePerKmer;

    return totalSize;
}

void Classifier::classifyReads() {
    Buffer<Kmer> queryKmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;
    queryList.reserve(512 * 1024);

    reporter->openReadClassificationFile();

    bool complete = false;
    SeqEntry * savedSeq_1 = new SeqEntry();
    SeqEntry * savedSeq_2 = new SeqEntry(); 
    uint64_t processedReadCnt = 0;

    std::cout << "--------------------" << std::endl;
    while (!complete) {
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = par.seqMode == 2 ? KSeqFactory(par.filenames[1].c_str()) : nullptr;

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        size_t queryKmerBufferSize = calculateBufferSize(queryList.size(), matchPerKmer);
        queryKmerBuffer.reallocateMemory(queryKmerBufferSize);
        matchBuffer.reallocateMemory(queryKmerBufferSize * matchPerKmer);

        bool moreReads = true;
        while (moreReads) {
            queryList.clear(); // TODO: use reserve and push_back instead
            queryKmerBuffer.init();
            matchBuffer.startIndexOfReserve = 0;
            // matchBuffer.init();

            // 1) Extract query k-mers
            time_t start = time(nullptr);
            std::cout << "Query k-mer extraction : " << std::flush;
            uint64_t seqCnt = 0;
            moreReads = kmerExtractor->extractQueryKmers(
                    queryKmerBuffer,
                    queryList,
                    seqCnt,
                    savedSeq_1,
                    savedSeq_2,
                    kseq1,
                    kseq2);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // 2) Sort k-mers
            start = time(nullptr);
            std::cout << "Query k-mer sorting    : " << std::flush;
            SORT_PARALLEL(queryKmerBuffer.buffer, queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, Kmer::compareQueryKmer);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // for (size_t i = 0; i < queryKmerBuffer.startIndexOfReserve; i++) {
            //     if (queryKmerBuffer.buffer[i].qInfo.sequenceID == 1003) {
            //         cout << "Sorted Kmer: ";
            //         metamerPattern->printDNA(queryKmerBuffer.buffer[i].value); cout << endl;
            //     }
            // }
            // 3) Match k-mers
            start = time(nullptr);
            bool searchComplete = kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer, dbDir);
            if (searchComplete) {
                std::cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << std::endl;

                // 4) Sort matches
                kmerMatcher->sortMatches(&matchBuffer);

                // 5) Assign taxonomy
                assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);

                // 6) Write classification results
                start = time(nullptr);
                std::cout << "Writing results        : " << std::flush;
                reporter->writeReadClassification(queryList);
                processedReadCnt += seqCnt;                    
                std::cout << difftime(time(nullptr), start) << " s" << std::endl;
                std::cout << "Processed read count   : " << processedReadCnt << std::endl;
            } else {
                matchPerKmer *= 2;
                moreReads = true;
                std::cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << std::endl;
                break;
            }
            std::cout << "--------------------" << std::endl;

        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (!moreReads) {
            complete = true;
        }
        std::cout << "--------------------" << std::endl;
    }
    delete savedSeq_1;
    delete savedSeq_2;
    reporter->closeReadClassificationFile();
    std::cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << std::endl;    


    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCounts, parentToChildren);
    for (const auto & it : sp2scoreSum_global) {
        sp2scoreSum_global[it.first] = it.second / cladeCounts[it.first].cladeCount;
    }
    rollUpScore(sp2scoreSum_global, parentToChildren, 1);

    // Finalize original classification results
    reporter->writeReportFile(processedReadCnt, taxCounts, ReportType::Default);
    std::cout << "Taxonomic classification completed." << std::endl;

}

void Classifier::startClassify(const LocalParameters &par) {
    Buffer<Kmer> queryKmerBuffer;
    Buffer<Match> matchBuffer;
    Buffer<MatchWithPos> matchBufferWithPos;
    vector<Query> queryList;
    size_t numOfTatalQueryKmerCnt = 0;
    reporter->openReadClassificationFile();

    bool complete = false;
    size_t processedReadCnt = 0;
    size_t tries = 0;
    size_t totalSeqCnt = 0;
    
    // Extract k-mers from query sequences and compare them to target k-mer DB
    while (!complete) {
        tries++;
        queryIndexer->setBytesPerKmer(matchPerKmer);
        queryIndexer->indexQueryFile(processedReadCnt);
        const vector<QuerySplit> & queryReadSplit = queryIndexer->getQuerySplits();
        if (tries == 1) {
            totalSeqCnt = queryIndexer->getReadNum_1();
            cout << "--------------------" << endl;
            cout << "Total read count : " << queryIndexer->getReadNum_1() << endl;
            cout << "Total read length: " << queryIndexer->getTotalReadLength() <<  "nt" << endl;
            cout << "--------------------" << endl;
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
            // memset(queryKmerBuffer.buffer, 0, queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer));

            // Allocate memory for match buffer
            if (par.storeKmerPos) {
                if (queryReadSplit.size() == 1) {
                    size_t remain = queryIndexer->getAvailableRam()
                                    - (queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer))
                                    - (queryIndexer->getReadNum_1() * 200); // TODO: check it later
                    matchBufferWithPos.reallocateMemory(remain / sizeof(MatchWithPos));
                } else {
                    matchBufferWithPos.reallocateMemory(queryReadSplit[splitIdx].kmerCnt * matchPerKmer);
                }
                matchBufferWithPos.startIndexOfReserve = 0;
            } else {
                if (queryReadSplit.size() == 1) {
                    size_t remain = queryIndexer->getAvailableRam()
                                    - (queryReadSplit[splitIdx].kmerCnt * sizeof(Kmer))
                                    - (queryIndexer->getReadNum_1() * 200); // TODO: check it later
                    matchBuffer.reallocateMemory(remain / sizeof(Match));
                } else {
                    matchBuffer.reallocateMemory(queryReadSplit[splitIdx].kmerCnt * matchPerKmer);
                }
                matchBuffer.startIndexOfReserve = 0;
            }

            // Extract query k-mers
            kmerExtractor->extractQueryKmers(queryKmerBuffer,
                                             queryList,
                                             queryReadSplit[splitIdx],
                                             par,
                                             kseq1,
                                             kseq2);
            
            // Search matches between query and target k-mers
            bool searchComplete = false;
            if (par.storeKmerPos) {
                searchComplete = kmerMatcher->matchKmersWithPos(&queryKmerBuffer, &matchBufferWithPos, dbDir);
                if (searchComplete) {
                    cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << endl;
                    SORT_PARALLEL(matchBufferWithPos.buffer,
                                  matchBufferWithPos.buffer + matchBufferWithPos.startIndexOfReserve,
                                  MatchWithPos::compare);
                    assignTaxonomy<MatchWithPos>(matchBufferWithPos.buffer, matchBufferWithPos.startIndexOfReserve, queryList, par);
                }
            } else {
                searchComplete = kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer, dbDir);
                if (searchComplete) {
                    cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << endl;
                    kmerMatcher->sortMatches(&matchBuffer);
                    assignTaxonomy<Match>(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                }
            }

            if (searchComplete) {
                reporter->writeReadClassification(queryList);
                if (par.em) {
                    reporter->writeMappings(queryList, processedReadCnt);
                    getTopSpecies(queryList);
                }
                processedReadCnt += queryReadSplit[splitIdx].readCnt;
                cout << "Processed read count   : " << processedReadCnt << " (" << (double) processedReadCnt / (double) totalSeqCnt << ")" << endl;
            } else { // search was incomplete
                matchPerKmer += 4;
                cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << endl;
                break;
            }
            cout << "--------------------" << endl;
        }
         
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (processedReadCnt == totalSeqCnt) {
            complete = true;
        }

        cout << "--------------------" << endl;
    }

    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCounts, parentToChildren);

    // Finalize original classification results
    cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << endl;    
    reporter->closeReadClassificationFile();
    reporter->writeReportFile(totalSeqCnt, taxCounts, ReportType::Default);
    cout << "Taxonomic classification completed." << endl;

    // Perform EM algorithm    
    if (par.em) {
        cout << "Running EM algorithm for taxonomic re-assignment..." << endl;
        reporter->freeMappingWriteBuffer();        
        loadOriginalResults(reporter->getClassificationFileName(), totalSeqCnt);   
        em(totalSeqCnt);
        reporter->writeReclassifyResults(emResults);
        unordered_map<TaxID, TaxonCounts> emCladeCounts = taxonomy->getCladeCounts(emTaxCounts, parentToChildren);
        reporter->writeReportFile(totalSeqCnt, emTaxCounts, ReportType::EM);

        unordered_map<TaxID, TaxonCounts> reclassifyCladeCounts = taxonomy->getCladeCounts(reclassifyTaxCounts, parentToChildren);
        reporter->writeReportFile(totalSeqCnt, reclassifyTaxCounts, ReportType::EM_RECLASSIFY);
    }
    
    return;
}

template <typename MatchType>
void Classifier::assignTaxonomy(const MatchType *matchList,
                               size_t numOfMatches,
                               std::vector<Query> &queryList,
                               const LocalParameters &par) {
    time_t beforeAnalyze = time(nullptr);
    std::cout << "K-mer match analysis   : " << std::flush;
    // Divide matches into blocks for multi threading
    size_t seqNum = queryList.size();
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    if (par.printLog==1) {
#ifdef OPENMP
        omp_set_num_threads(1);
#endif
    }
    while (matchIdx < numOfMatches) {
        currentQuery = matchList[matchIdx].qKmer.qInfo.sequenceID;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList[matchIdx].qKmer.qInfo.sequenceID) && (matchIdx < numOfMatches)) ++matchIdx;
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }
    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx, par)
    {
        Taxonomer<MatchType> taxonomer(par, taxonomy, metamerPattern);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            taxonomer.chooseBestTaxon(
                            matchBlocks[i].id - 1,
                            matchBlocks[i].start,
                            matchBlocks[i].end,
                            matchList,
                            queryList);
        }

        #pragma omp critical
        {
            for (const auto & it : taxonomer.sp2scoreSum) {
                sp2scoreSum_global[it.first] += it.second;
            }   

            if (par.storeKmerPos) {
                for (const auto &it : taxonomer.sp2totalReadLength) {
                    sp2totalReadLength[it.first] += it.second;
                }

                for (auto& [spId, localArray] : taxonomer.sp2coverage) {
                    const vector<uint8_t> &localBins = localArray;
                    if (sp2coverage_global.find(spId) == sp2coverage_global.end()) {
                        sp2coverage_global[spId].resize(65536, 0);
                    }

                    uint8_t *globalBins = sp2coverage_global[spId].data();
                    for (int j = 0; j < 65536; ++j) {
                        int sum = globalBins[j] + localBins[j];
                        globalBins[j] = (sum > 255) ? 255 : static_cast<uint8_t>(sum);
                    }
                }
            }
        }    
    }

    if (par.printLog) {
#ifdef OPENMP
        omp_set_num_threads(par.threads);
#endif
    }

    for (size_t i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    
    delete[] matchBlocks;
    cout << double(time(nullptr) - beforeAnalyze) << " s" << endl;

}

// EM Algorithm Constants and Helper Functions
namespace {
    static const double ABUNDANCE_EXP_TAU = 3.0;
    static const double ABUNDANCE_SMOOTH_EPS = 1e-8;
    static const double COVERAGE_PRIOR_WEIGHT = 0.1;
    static const double SQUAREM_STEP_MIN = -1.0;
    static const double SQUAREM_STEP_MAX = 1.0;
    
    static inline double clamp01(double value) {
        return std::max(0.0, std::min(1.0, value));
    }
    
    static inline std::vector<double> projectSimplex(const std::vector<double> &x) {
        std::vector<double> projected = x;
        for (size_t i = 0; i < projected.size(); ++i) {
            if (projected[i] < 0.0) {
                projected[i] = 0.0;
            }
        }
        double sum = 0.0;
        for (size_t i = 0; i < projected.size(); ++i) {
            sum += projected[i];
        }
        if (sum > 0.0) {
            for (size_t i = 0; i < projected.size(); ++i) {
                projected[i] /= sum;
            }
        }
        return projected;
    }
    
    static std::vector<double> abundanceVectorFromMap(
        const std::unordered_map<TaxID, double> &abundanceMap,
        const std::vector<TaxID> &targetList) {
        std::vector<double> out;
        out.reserve(targetList.size());
        for (size_t i = 0; i < targetList.size(); ++i) {
            auto it = abundanceMap.find(targetList[i]);
            out.push_back(it != abundanceMap.end() ? it->second : 0.0);
        }
        return out;
    }
    
    static void setAbundanceFromVector(
        std::unordered_map<TaxID, double> &abundanceMap,
        const std::vector<TaxID> &targetList,
        const std::vector<double> &abundanceVector) {
        for (size_t i = 0; i < targetList.size(); ++i) {
            abundanceMap[targetList[i]] = abundanceVector[i];
        }
    }
}

void Classifier::em(size_t totalQueryCnt) {
    loadMappings(reporter->getMappingFileName());
    int maxTaxId = taxonomy->getMaxTaxID();
    vector<uint32_t> sp2uniqKmerCnt(maxTaxId + 1, 0);
    countUniqueKmerPerSpecies(sp2uniqKmerCnt);
    vector<double> sp2lengthFactor(maxTaxId + 1, 0.0);
    for (size_t i = 0; i < sp2uniqKmerCnt.size(); i++) {
        if (sp2uniqKmerCnt[i] > 0) {
            sp2lengthFactor[i] = 1.0 / log(sp2uniqKmerCnt[i]);
        } else {
            sp2lengthFactor[i] = 0.0;
        }
    }
    
    time_t st = time(nullptr);
    cout << "EM algorithm for taxonomic assignment (with SQUAREM acceleration): " << std::endl;
    size_t idx = 0;
    std::vector<std::pair<size_t, size_t>> queryRanges;
    std::unordered_set<TaxID> mappedSpeciesSet;
    while (idx < mappingResListSize) {
        TaxID currentQuery = mappingResList[idx].queryId;
        size_t startIdx = idx;
        while (idx < mappingResListSize && mappingResList[idx].queryId == currentQuery) {
            mappedSpeciesSet.insert(mappingResList[idx].speciesId);
            idx++;
        }
        queryRanges.emplace_back(startIdx, idx);
    }

    cout << "Species count: " << mappedSpeciesSet.size() << std::endl;
    std::vector<TaxID> spList;
    std::unordered_map<TaxID, double> abundance;
    if (mappedSpeciesSet.empty()) {
        cout << "No species found in mapping results. Skipping EM." << endl;
        return;
    }
    const double initProb = 1.0 / static_cast<double>(mappedSpeciesSet.size());
    for (const auto & sp : mappedSpeciesSet) {
        taxProbs[sp] = initProb;
        abundance[sp] = taxProbs[sp];
        spList.push_back(sp);
    }

    // Build coverage metrics from the original coverage stack when available.
    // This keeps the EM prior aligned with the original coverage report pipeline.
    if (sp2covMetric.empty() && !sp2coverage_global.empty()) {
        if (sp2genomeSize.empty()) {
            parseSp2GenomeSize();
        }
        const std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
        const unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCounts, parentToChildren);
        for (const auto &it : sp2coverage_global) {
            const TaxID spId = it.first;
            const auto countIt = cladeCounts.find(spId);
            if (countIt == cladeCounts.end() || countIt->second.cladeCount == 0) {
                continue;
            }
            const uint64_t genomeSize = (sp2genomeSize.find(spId) != sp2genomeSize.end()) ? sp2genomeSize[spId] : 65536;
            const uint64_t totalReadLen = (sp2totalReadLength.find(spId) != sp2totalReadLength.end()) ? sp2totalReadLength[spId] : 0;
            sp2covMetric[spId] = calCovMetrics(it.second, countIt->second.cladeCount, totalReadLen, genomeSize);
        }
        rollUpCoverageMetrics(parentToChildren, cladeCounts, sp2covMetric, 1);
    }

    // Species-level coverage prior used in M-step smoothing.
    std::unordered_map<TaxID, double> sp2coveragePrior;
    sp2coveragePrior.reserve(spList.size());
    for (const TaxID sp : spList) {
        auto covIt = sp2covMetric.find(sp);
        sp2coveragePrior[sp] = (covIt != sp2covMetric.end()) ? covIt->second.coverage : 0.0;
    }

    const double alphaMax = 1.0;
    const int maxIter = 1000;
    const double convergenceTol = 1e-6;
    
    // SQUAREM acceleration: track multiple iterations
    std::vector<double> x0 = abundanceVectorFromMap(abundance, spList);
    for (int iter = 0; iter < maxIter; ++iter) {
        // Annealed exponent increases from 0 to alphaMax over iterations
        const double abundanceExponent = alphaMax * (1.0 - std::exp(-static_cast<double>(iter + 1) / ABUNDANCE_EXP_TAU));

        // M-step with coverage-prior smoothing:
        // theta_tilde = responsibilities + kappa * coverage + eps
        std::unordered_map<TaxID, double> abundanceUpdate;
        abundanceUpdate.reserve(spList.size());
        for (const TaxID sp : spList) {
            abundanceUpdate[sp] = COVERAGE_PRIOR_WEIGHT * sp2coveragePrior[sp] + ABUNDANCE_SMOOTH_EPS;
        }

        #pragma omp parallel default(none) shared(queryRanges, mappingResList, abundance, sp2lengthFactor, abundanceExponent, abundanceUpdate)
        {
            std::unordered_map<TaxID, double> localResponsibilities;

            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < queryRanges.size(); ++i) {
                double denom = 0.0;
                std::vector<double> posteriors;

                for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                    TaxID sp = mappingResList[j].speciesId;
                    double abundSmoothed = std::max(abundance[sp], ABUNDANCE_SMOOTH_EPS);
                    double weight = mappingResList[j].score *
                                   std::pow(abundSmoothed, abundanceExponent) *
                                   sp2lengthFactor[sp];
                    posteriors.push_back(weight);
                    denom += weight;
                }

                if (denom <= 0.0) { continue; }

                for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                    TaxID sp = mappingResList[j].speciesId;
                    const double posterior = posteriors[j - queryRanges[i].first] / denom;
                    localResponsibilities[sp] += posterior;
                }
            }

            #pragma omp critical
            {
                for (const auto &spAbund : localResponsibilities) {
                    abundanceUpdate[spAbund.first] += spAbund.second;
                }
            }
        }

        abundance.swap(abundanceUpdate);
        
        // Normalize abundance to sum to 1
        double totalAbundance = 0.0;
        for (const auto & sp : spList) {
            totalAbundance += abundance[sp];
        }
        if (totalAbundance > 0.0) {
            for (auto & sp : spList) {
                abundance[sp] /= totalAbundance;
            }
        }
        
        // Compute log-likelihood for convergence monitoring
        std::vector<double> x1 = abundanceVectorFromMap(abundance, spList);
        
        // Perform second EM update for SQUAREM acceleration
        std::unordered_map<TaxID, double> abundanceX2 = abundance;
        std::unordered_map<TaxID, double> abundanceX2Update;
        abundanceX2Update.reserve(spList.size());
        for (const TaxID sp : spList) {
            abundanceX2Update[sp] = COVERAGE_PRIOR_WEIGHT * sp2coveragePrior[sp] + ABUNDANCE_SMOOTH_EPS;
        }

        #pragma omp parallel default(none) shared(queryRanges, mappingResList, abundanceX2, sp2lengthFactor, abundanceExponent, ABUNDANCE_SMOOTH_EPS, abundanceX2Update)
        {
            std::unordered_map<TaxID, double> localResponsibilities;
            
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < queryRanges.size(); ++i) {
                double denom = 0.0;
                std::vector<double> posteriors;
                
                for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                    TaxID sp = mappingResList[j].speciesId;
                    double abundSmoothed = std::max(abundanceX2[sp], ABUNDANCE_SMOOTH_EPS);
                    double weight = mappingResList[j].score * 
                                   std::pow(abundSmoothed, abundanceExponent) * 
                                   sp2lengthFactor[sp];
                    posteriors.push_back(weight);
                    denom += weight;
                }
                
                if (denom <= 0.0) { continue; }
                
                for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                    TaxID sp = mappingResList[j].speciesId;
                    double posterior = posteriors[j - queryRanges[i].first] / denom;
                    localResponsibilities[sp] += posterior;
                }
            }
            
            #pragma omp critical
            {
                for (const auto &spAbund : localResponsibilities) {
                    abundanceX2Update[spAbund.first] += spAbund.second;
                }
            }
        }
        abundanceX2.swap(abundanceX2Update);
        
        // Normalize abundanceX2
        totalAbundance = 0.0;
        for (const auto & sp : spList) {
            totalAbundance += abundanceX2[sp];
        }
        if (totalAbundance > 0.0) {
            for (auto & sp : spList) {
                abundanceX2[sp] /= totalAbundance;
            }
        }
        
        std::vector<double> x2 = abundanceVectorFromMap(abundanceX2, spList);
        
        // SQUAREM extrapolation
        std::vector<double> r(x0.size()), v(x0.size());
        double normR = 0.0, normV = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            r[i] = x1[i] - x0[i];
            v[i] = x2[i] - x1[i] - r[i];
            normR += r[i] * r[i];
            normV += v[i] * v[i];
        }
        normR = std::sqrt(normR);
        normV = std::sqrt(normV);
        
        double accel = (normV == 0.0) ? -1.0 : -(normR / normV);
        accel = std::max(SQUAREM_STEP_MIN, std::min(SQUAREM_STEP_MAX, accel));
        
        std::vector<double> xNew(x0.size());
        for (size_t i = 0; i < x0.size(); ++i) {
            xNew[i] = x0[i] - 2.0 * accel * r[i] + accel * accel * v[i];
        }
        xNew = projectSimplex(xNew);
        setAbundanceFromVector(abundance, spList, xNew);
        
        // Check convergence
        double parameterChange = 0.0;
        for (size_t i = 0; i < x0.size(); ++i) {
            parameterChange = std::max(parameterChange, std::fabs(xNew[i] - x0[i]));
        }
        
        cout << "Iteration " << iter + 1 << ": change=" << parameterChange << endl;
        
        if (parameterChange < convergenceTol && iter > 5) {
            cout << "Converged after " << (iter + 1) << " iterations." << endl;
            break;
        }
        
        x0 = xNew;
    }

    // Final update of tax probabilities
    for (const auto & sp : spList) {
        taxProbs[sp] = abundance[sp];
    }
    
    size_t explainedQueryCount = 0;
    for (size_t i = 0; i < spList.size(); i++) {
        emTaxCounts[spList[i]] = taxProbs[spList[i]] * queryRanges.size();
        explainedQueryCount += emTaxCounts[spList[i]];
    }
    emTaxCounts[0] = totalQueryCnt - explainedQueryCount;
    cout << "EM completed in " << double(time(nullptr) - st) << " s" << endl;
    reclassify(queryRanges, mappingResList, sp2lengthFactor, totalQueryCnt);
}


void Classifier::reclassify(
    const std::vector<std::pair<size_t, size_t>> & queryRanges,
    const MappingRes * mappingResList,
    const vector<double> & sp2lengthFactor,
    size_t totalQueryCnt) 
{
    #pragma omp parallel default(none) shared(queryRanges, mappingResList, sp2lengthFactor, std::cout)
    {
        std::vector<std::pair<TaxID, double>> sp2prob;
        std::vector<TaxID> candidateSpecies;
        std::unordered_map<TaxID, unsigned int> localTaxCounts;

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < queryRanges.size(); ++i) {
            double denom = 0.0;
            candidateSpecies.clear();
            sp2prob.clear();

            TaxID currentQuery = mappingResList[queryRanges[i].first].queryId;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                const TaxID sp = mappingResList[j].speciesId;
                const auto probIt = taxProbs.find(sp);
                const double prob = (probIt == taxProbs.end()) ? 0.0 : probIt->second;
                double score = prob * mappingResList[j].score * sp2lengthFactor[sp];
                denom += score;
                sp2prob.emplace_back(sp, score);
            }
            if (denom == 0.0) {
                emResults[currentQuery].taxId = 0;
                emResults[currentQuery].score = 0.0;
                continue;
            }

            for (auto & sp : sp2prob) {
                sp.second /= denom; // Normalize scores
            }

            std::sort(sp2prob.begin(), sp2prob.end(),
                      [](const std::pair<TaxID, double> &a, const std::pair<TaxID, double> &b) {
                          return a.second > b.second;
                      });

            double sum = 0.0;
            for (size_t j = 0; j < sp2prob.size() && sum < 0.5; ++j) {
                sum += sp2prob[j].second;
                candidateSpecies.push_back(sp2prob[j].first);
            }
            emResults[currentQuery].taxId = taxonomy->LCA(candidateSpecies)->taxId;
            emResults[currentQuery].score = sum;
            localTaxCounts[emResults[currentQuery].taxId] += 1;
        }
        #pragma omp critical
        {
            for (const auto & taxCount : localTaxCounts) {
                reclassifyTaxCounts[taxCount.first] += taxCount.second;
            }
        }
    }
    size_t explainedQueryCnt = 0;
    for (const auto & taxCount : reclassifyTaxCounts) {
        if (taxCount.first != 0) { // Skip unclassified
            explainedQueryCnt += taxCount.second;
        }
    }
//    std::cout << double(time(nullptr) - st) << " s" << endl;
}

void Classifier::countUniqueKmerPerSpecies(
    vector<uint32_t> & sp2uniqKmerCnt) 
{
    string sp2uniqKmerCntFileName = dbDir + "/sp2uniqKmerCnt";
    if (FileUtil::fileExists(sp2uniqKmerCntFileName.c_str())) {
        cout << "Loading unique k-mer count per species from file: " << sp2uniqKmerCntFileName << endl;
        ifstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << endl;
            return;
        }
        TaxID taxId;
        uint32_t count;
        while (sp2uniqKmerCntFile >> taxId >> count) {
            if (taxId < sp2uniqKmerCnt.size()) {
                sp2uniqKmerCnt[taxId] = count;
            } else {
                cerr << "Warning: TaxID " << taxId << " exceeds the size of sp2uniqKmerCnt vector." << endl;
            }
        }
        sp2uniqKmerCntFile.close();
        return;
    } else {
        string infoFileName = dbDir + "/info";
        ReadBuffer<TaxID> readBuffer(infoFileName, 1000000);
        TaxID taxId;
        time_t st = time(nullptr);
        unordered_map<TaxID, TaxID> & taxid2speciesId = kmerMatcher->getTaxId2SpeciesId();
        cout << "Count unique k-mers per species: " << std::flush;
        while((taxId = readBuffer.getNext()) != 0) {
            TaxID speciesTaxId = taxid2speciesId[taxId];
            if (speciesTaxId) {
                ++sp2uniqKmerCnt[speciesTaxId];
            }
        }
        cout << double(time(nullptr) - st) << " s" << endl;
    
        // Save the unique k-mer count per species to a file
        ofstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << " for writing." << endl;
            return;
        }
        for (size_t i = 0; i < sp2uniqKmerCnt.size(); ++i) {
            if (sp2uniqKmerCnt[i] > 0) {
                sp2uniqKmerCntFile << i << " " << sp2uniqKmerCnt[i] << "\n";
            }
        }
        sp2uniqKmerCntFile.close();
    }
}

void Classifier::loadMappings(const string & mappingResFileName) {
    size_t fileSize = FileUtil::getFileSize(mappingResFileName);
    mappingResListSize = fileSize / sizeof(MappingRes);
    mappingResList = new MappingRes[mappingResListSize];
    FILE * mappingResFile = fopen(mappingResFileName.c_str(), "rb");
    if (mappingResFile == nullptr) {
        cerr << "Error: Could not open mapping results file " << mappingResFileName << endl;
        return;
    }
    size_t readCount = fread(mappingResList, sizeof(MappingRes), mappingResListSize, mappingResFile);
    if (readCount != mappingResListSize) {
        cerr << "Error: Could not read all mapping results from file " << mappingResFileName << endl;
        fclose(mappingResFile);
        return;
    }
    fclose(mappingResFile);
}

void Classifier::loadOriginalResults(
    const string & classificationFileName,
    size_t seqNum) 
{
    ifstream classificationFile(classificationFileName);
    if (!classificationFile.is_open()) {
        cerr << "Error: Could not open classification results file " << classificationFileName << endl;
        return;
    }
    emResults.reserve(seqNum);
    string line;
    while (getline(classificationFile, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        vector<string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 4);
        if (columns.size() < 4) {
            cerr << "Error: Invalid line format in classification results file: " << line << endl;
            continue;
        }
        string queryName = columns[1];
        int length = atoi(columns[3].c_str());
        emResults.emplace_back(queryName, length);
    }
}

CovMetric Classifier::calCovMetrics(
    const std::vector<uint8_t>& bins,
    int readCnt,
    uint64_t totalReadLength,
    uint64_t genomeSize)
{
    uint64_t totalCount = 0;
    double sum_c_log_c = 0.0;
    int occupiedBins = 0;

    int occupiedMacroBins = 0;
    bool macro_seen[256] = {false};

    size_t binLimit = std::min(genomeSize + 1, static_cast<uint64_t>(65536));
    for (size_t i = 1; i < binLimit; ++i) {
        uint8_t count = bins[i];
        if (count == 0) continue;

        totalCount += count;
        occupiedBins++;
        sum_c_log_c += C_LOG2_C[count];

        size_t macro_bin_idx = i >> 8;
        if (!macro_seen[macro_bin_idx]) {
            macro_seen[macro_bin_idx] = true;
            occupiedMacroBins++;
        }
    }

    if (totalCount == 0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }

    double effective_bins = std::min(65535.0, static_cast<double>(genomeSize));
    if (effective_bins <= 1.0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }

    double max_macro_bins = std::ceil(effective_bins / 256.0);
    max_macro_bins = std::min(max_macro_bins, 256.0);
    double macro_coverage = static_cast<double>(occupiedMacroBins) / max_macro_bins;
    macro_coverage = std::min(1.0, macro_coverage);

    double coverage = static_cast<double>(occupiedBins) / effective_bins;
    coverage = std::min(1.0, coverage);

    double H_obs = std::log2(static_cast<double>(totalCount)) - (sum_c_log_c / totalCount);
    H_obs = std::max(0.0, H_obs);

    double max_H_standard = std::log2(effective_bins);
    double evenness = (max_H_standard > 0.0001) ? std::min(1.0, H_obs / max_H_standard) : 0.0;

    double read_length = static_cast<double>(totalReadLength) / readCnt;
    double bin_size_bp = static_cast<double>(genomeSize) / effective_bins;
    double bins_per_read = 1.0 + (read_length / bin_size_bp);
    double expected_occupied = effective_bins * (1.0 - std::exp(-(readCnt * bins_per_read) / effective_bins));
    double expected_max_H = std::log2(std::max(1.0, expected_occupied));

    double adjustedEvenness = 0.0;
    if (expected_max_H >= 0.0001) {
        adjustedEvenness = std::min(1.0, H_obs / expected_max_H);
    }

    double unified_score = std::pow(2.0, H_obs) / effective_bins;
    return {evenness, coverage, adjustedEvenness, unified_score, macro_coverage};
}

void Classifier::rollUpCoverageMetrics(
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren,
    const std::unordered_map<TaxID, TaxonCounts>& cladeCounts,
    std::unordered_map<TaxID, CovMetric>& allMetrics,
    TaxID currentTaxID)
{
    if (allMetrics.find(currentTaxID) != allMetrics.end()) {
        return;
    }

    auto childrenIt = parentToChildren.find(currentTaxID);
    if (childrenIt == parentToChildren.end() || childrenIt->second.empty()) {
        return;
    }

    double weightedEvenness = 0.0;
    double weightedCoverage = 0.0;
    double weightedAdjEvenness = 0.0;
    double weightedUnifiedScore = 0.0;
    double weightedAvgScore = 0.0;
    double weightedMacroCoverage = 0.0;
    uint64_t totalChildReads = 0;

    for (TaxID childID : childrenIt->second) {
        rollUpCoverageMetrics(parentToChildren, cladeCounts, allMetrics, childID);

        auto countIt = cladeCounts.find(childID);
        auto metricIt = allMetrics.find(childID);
        if (countIt != cladeCounts.end() && metricIt != allMetrics.end()) {
            uint64_t reads = countIt->second.cladeCount;
            if (reads > 0) {
                weightedEvenness += metricIt->second.evenness * reads;
                weightedCoverage += metricIt->second.coverage * reads;
                weightedAdjEvenness += metricIt->second.adjustedEvenness * reads;
                weightedUnifiedScore += metricIt->second.unifiedScore * reads;
                weightedAvgScore += metricIt->second.avgScore * reads;
                weightedMacroCoverage += metricIt->second.macroCoverage * reads;
                totalChildReads += reads;
            }
        }
    }

    if (totalChildReads > 0) {
        CovMetric parentMetric;
        parentMetric.evenness = weightedEvenness / static_cast<double>(totalChildReads);
        parentMetric.coverage = weightedCoverage / static_cast<double>(totalChildReads);
        parentMetric.adjustedEvenness = weightedAdjEvenness / static_cast<double>(totalChildReads);
        parentMetric.unifiedScore = weightedUnifiedScore / static_cast<double>(totalChildReads);
        parentMetric.avgScore = weightedAvgScore / static_cast<double>(totalChildReads);
        parentMetric.macroCoverage = weightedMacroCoverage / static_cast<double>(totalChildReads);
        allMetrics[currentTaxID] = parentMetric;
    }
}



void Classifier::rollUpScore(
    std::unordered_map<TaxID, double> & sp2score_global,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren,
    TaxID currentTaxID) 
{
    auto childrenIt = parentToChildren.find(currentTaxID);
    if (childrenIt == parentToChildren.end() || childrenIt->second.empty()) {
        return; 
    }

    double childrenSum = 0.0;
    // double weightedAvgScore = 0.0;
    // uint64_t totalChildReads = 0;

    // 1. Recursively process all children FIRST (Post-order traversal)
    for (TaxID childID : childrenIt->second) {
        rollUpScore(sp2score_global, parentToChildren, childID);
        
        auto scoreIt = sp2score_global.find(childID);
        if (scoreIt != sp2score_global.end()) {
            childrenSum += scoreIt->second;
        }
    }

    sp2score_global[currentTaxID] += childrenSum;
}
