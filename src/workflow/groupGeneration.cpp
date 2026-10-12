#include "GroupGenerator.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include <fstream>
#include <string>

void setGroupGenerationDefaults(LocalParameters & par){
    // The reads-per-k-mer cap, chosen from the sample rather than fixed. A k-mer shared by m reads
    // contributes C(m,2) edges, so the frequency tail decides the edge volume, and the tail is not
    // the same shape on two datasets: at this quantile the cap resolves to 63 on species-inclusion,
    // 127 on species-exclusion and 255 on CAMI2 strain-madness. An absolute default cannot be
    // right on all three. maxKmerReads is the absolute form and stays available, but it is off by
    // default: an explicit cap wins over the quantile, so leaving it set would make the quantile dead.
    par.maxKmerReads = 0;
    par.maxKmerQuantile = 0.995f;

    // Phase 1 core threshold as a fraction of k-mers per read. 0.5 was chosen for purity on CAMI2
    // strain-madness (purity 0.953 at 0.5 against 0.670 at 0.3); species-exclusion is 0.999 pure
    // at every value, so the default follows the harder case.
    par.minOverlapRatio = 0.5f;
    // Weak-band lower bound as a fraction of the core threshold, also Phase 3's floor.
    // 5/15 = 0.3333 reproduces the species-inclusion operating point, where the band was (5, 15].
    par.weakBandRatio = 0.3333f;
    // Intermediate partitions. 16 is the value the published numbers were produced at, and it is
    // a constant on purpose: 0 makes the routing follow --threads, so the same input on a machine
    // with a different core count would be split differently.
    par.partitions = 16;
    // K-mers discarded either side of a common-k-mer hit. 0 discards the hit alone; on CAMI2
    // strain-madness, 1 -> 0 took species Recall*c from 0.164 to 0.245 for 0.010 of purity.
    par.commonKmerSpan = 0;
    // Peak disk the subGraph_* intermediates may hold at once. 0 derives it from the free space
    // at the output directory (80%).
    par.maxTmpDiskMiB = 0;
    par.edgeMode = 0;
    par.syncmer = 1;
    par.smerLen = 5;
    par.seqMode = 2;
    par.verbosity = 3;
    par.ramUsage = LocalParameters::defaultRamUsage();
    par.printLog = 0;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4;
}

bool prepareGroupGeneration(LocalParameters & par, const std::string & commonKmerDb, std::string & error,
                            bool adoptDbKmerSettings) {
    par.kmerFormat = (par.syncmer == 0) ? 3 : 5;

    // Both algorithm thresholds are ratios and there is no absolute fallback for either, so a
    // value outside its usable range has to stop the run rather than silently pick something else.
    if (par.minOverlapRatio <= 0.0f) {
        error = "--min-overlap-ratio must be > 0 (given " + std::to_string(par.minOverlapRatio) + ").\n"
                "       The Phase 1 core threshold is derived from it as ratio x k-mers per read;\n"
                "       there is no absolute threshold to fall back to.";
        return false;
    }
    if (par.weakBandRatio <= 0.0f || par.weakBandRatio >= 1.0f) {
        error = "--weak-band-ratio must be in (0, 1) (given " + std::to_string(par.weakBandRatio) + ").\n"
                "       It is the weak band's lower bound as a fraction of the core threshold.";
        return false;
    }
    if (par.partitions < 0) {
        error = "--partitions must be >= 0 (given " + std::to_string(par.partitions) + "). 0 means follow --threads.";
        return false;
    }
    if (par.maxKmerQuantile < 0.0f || par.maxKmerQuantile > 1.0f) {
        error = "--max-kmer-quantile must be in [0, 1] (given " + std::to_string(par.maxKmerQuantile) + ").";
        return false;
    }
    if (par.edgeMode != 0 && par.edgeMode != 1) {
        error = "--edge-mode must be 0 (clique) or 1 (star).";
        return false;
    }
    if (!FileUtil::directoryExists(commonKmerDb.c_str())) {
        error = "common k-mer DB " + commonKmerDb + " is not a directory.";
        return false;
    }
    for (const char *name : {"/diffIdx", "/split"}) {
        if (!FileUtil::fileExists((commonKmerDb + name).c_str())) {
            error = "common k-mer DB " + commonKmerDb + " has no " + std::string(name + 1)
                    + " file. Build it with create-common-kmer-list.";
            return false;
        }
    }

    // The query k-mers have to be built the same way as the common k-mer DB, or the two k-mer
    // spaces do not line up and the common-k-mer filter silently removes nothing. DBs built
    // before this record existed cannot be checked and pass silently.
    const std::string paramFile = commonKmerDb + "/kmer_params";
    if (FileUtil::fileExists(paramFile.c_str())) {
        int dbSyncmer = -1, dbSmerLen = -1, dbKmerFormat = -1;
        std::ifstream in(paramFile);
        std::string key;
        int value = 0;
        while (in >> key >> value) {
            if (key == "syncmer") { dbSyncmer = value; }
            else if (key == "smer_len") { dbSmerLen = value; }
            else if (key == "kmer_format") { dbKmerFormat = value; }
        }
        if (adoptDbKmerSettings) {
            // The DB decides how the query k-mers are built.
            if (dbSyncmer >= 0) { par.syncmer = dbSyncmer; }
            if (dbSmerLen > 0) { par.smerLen = dbSmerLen; }
            par.kmerFormat = (par.syncmer == 0) ? 3 : 5;
            if (dbKmerFormat >= 0 && dbKmerFormat != par.kmerFormat) {
                error = "common k-mer DB " + commonKmerDb + " records k-mer format " + std::to_string(dbKmerFormat)
                        + ", which the read grouping cannot use (it needs 3 or 5).";
                return false;
            }
        } else if (dbSyncmer != par.syncmer || dbSmerLen != par.smerLen || dbKmerFormat != par.kmerFormat) {
            error = "k-mer settings do not match the common k-mer DB at " + commonKmerDb + ".\n"
                    "       DB:      --syncmer " + std::to_string(dbSyncmer) + " --smer-len " + std::to_string(dbSmerLen)
                    + " (k-mer format " + std::to_string(dbKmerFormat) + ")\n"
                    "       Request: --syncmer " + std::to_string(par.syncmer) + " --smer-len " + std::to_string(par.smerLen)
                    + " (k-mer format " + std::to_string(par.kmerFormat) + ")\n"
                    "       Rebuild the DB with these settings, or pass the DB's settings.";
            return false;
        }
    }
    return true;
}

bool runReadGrouping(LocalParameters & par,
                     const std::vector<std::string> & queryFiles,
                     const std::string & commonKmerDb,
                     const std::string & outDir)
{
    if (queryFiles.empty() || queryFiles.size() > 2) {
        std::cerr << "Error: read grouping needs one read file, or two for paired-end (given "
                  << queryFiles.size() << ")." << std::endl;
        return false;
    }
    if ((par.seqMode == 2) != (queryFiles.size() == 2)) {
        std::cerr << "Error: --seq-mode " << par.seqMode << " does not match " << queryFiles.size()
                  << " read file(s)." << std::endl;
        return false;
    }
    if (!FileUtil::directoryExists(outDir.c_str())) {
        FileUtil::makeDir(outDir.c_str());
    }

    // GroupGenerator and QueryIndexer take their inputs from par.filenames laid out as the
    // grouping command does: <reads...> <common k-mer DB> <output directory>.
    const std::vector<std::string> savedFilenames = par.filenames;
    par.filenames = queryFiles;
    par.filenames.push_back(commonKmerDb);
    par.filenames.push_back(outDir);

    GroupGenerator * groupGenerator = new GroupGenerator(par);
    groupGenerator->startGroupGeneration(par);
    delete groupGenerator;

    par.filenames = savedFilenames;
    return FileUtil::fileExists((outDir + "/groupMap").c_str());
}

// Everything the two grouping commands share. They differ only in the defaults their caller set,
// so the body lives here once rather than being copied into each entry point.
static int runGroupGeneration(int argc, const char **argv, const Command& command,
                              LocalParameters & par)
{
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    // Positional layout:
    //   --seq-mode 2   : <read_1> <read_2> <common k-mer DB> <output directory>
    //   --seq-mode 1/3 : <read>            <common k-mer DB> <output directory>
    const size_t readFileCnt = (par.seqMode == 2) ? 2 : 1;
    if (par.filenames.size() != readFileCnt + 2) {
        cerr << "Error: expected " << readFileCnt << " read file(s), a common k-mer DB and an output directory"
             << " for --seq-mode " << par.seqMode << ", got " << par.filenames.size() << " arguments." << endl;
        return 1;
    }
    if (par.seqMode == 2 && FileUtil::directoryExists(par.filenames[1].c_str())) {
        cerr << "Error: " << par.filenames[1] << " is a directory. Please specify a query file name." << endl;
        cerr << "       For '--seq-mode 2', please provide two query files." << endl;
        return 1;
    }
    if (par.seqMode != 2 && FileUtil::fileExists(par.filenames[1].c_str())
        && !FileUtil::directoryExists(par.filenames[1].c_str())) {
        cerr << "Error: " << par.filenames[1] << " is a file. Please specify a database directory." << endl;
        cerr << "       For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
        return 1;
    }

    const std::vector<std::string> queryFiles(par.filenames.begin(), par.filenames.begin() + readFileCnt);
    const std::string commonKmerDb = par.filenames[readFileCnt];
    const std::string outDir = par.filenames[readFileCnt + 1];

    std::string error;
    if (!prepareGroupGeneration(par, commonKmerDb, error)) {
        cerr << "Error: " << error << endl;
        return 1;
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    return runReadGrouping(par, queryFiles, commonKmerDb, outDir) ? 0 : 1;
}

int groupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupGenerationDefaults(par);
    return runGroupGeneration(argc, argv, command, par);
}

// The fast, light variant. Only the edge set differs: one shared k-mer produces a hub and its
// m-1 spokes instead of all C(m,2) pairs, so the edge volume is linear in m rather than quadratic.
// Every other default is grouping's, on purpose.
void setEasyGroupGenerationDefaults(LocalParameters & par){
    setGroupGenerationDefaults(par);
    par.edgeMode = 1; // EdgeMode::EDGE_MODE_STAR
}

int easyGroupGeneration(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setEasyGroupGenerationDefaults(par);
    return runGroupGeneration(argc, argv, command, par);
}
