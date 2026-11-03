#ifndef METABULI_CREATE_COMMON_KMER_LIST_H
#define METABULI_CREATE_COMMON_KMER_LIST_H

#include "editNames.h"
#include "accession2taxid.h"
#include "IndexCreator.h"


void setDefaults(LocalParameters & par){
    par.kmerFormat = 3;
    par.syncmer = 0;
    par.smerLen = 5;
    par.gtdb = 0;
    par.makeLibrary = 0;
    par.ramUsage = LocalParameters::defaultRamUsage();
    par.validateInput = 0;
    par.validateDb = 0;
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    // Six-frame extraction. The ORF path (fillTargetKmerBuffer) needs a MetamerPattern, which the
    // format-aware IndexCreator used below does not have.
    par.cdsInfo = "x";
    par.maskProb = 0.9;
    par.maskMode = 0;
    par.accessionLevel = 0;
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    // Not rand(). The name went into db.parameters, so the same inputs wrote a different
    // file every run, and nothing reads the field back except to round-trip it.
    par.dbName = "common-kmer-list";
}

int create_common_kmer_list(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string & dbDir = par.filenames[0];
    if (par.cdsInfo != "x") {
        cerr << "Error: create-common-kmer-list extracts k-mers from six-frame translations only;" << endl;
        cerr << "       --cds-info is not supported here." << endl;
        return 1;
    }
    const string taxonomyDir = par.filenames[3];

    #ifdef OPENMP
        omp_set_num_threads(par.threads);
    #endif

    if (!FileUtil::directoryExists(dbDir.c_str())) {
        FileUtil::makeDir(dbDir.c_str());
    }
    
    if (par.validateInput == 1) {
        const string & fnaListFileName = par.filenames[1];
        // Read the file line by line
        ifstream file(fnaListFileName);
        if (!file.is_open()) {
            cout << "Error: Unable to open file " << fnaListFileName << endl;
            return 1;
        }
        string singleFastaName;
        while (getline(file, singleFastaName)) {
            if (!LocalUtil::isFasta(singleFastaName)) {
                cout << "Error: " << singleFastaName << " has a unsupported extension." << endl;
                cout << "       Allowed extensions are .fna, .fasta, .fa, and their gzip versions (e.g., .fna.gz)" << endl;
                file.close();
                return 1;
            }
            if (!FileUtil::fileExists(singleFastaName.c_str())) {
                cout << "Error: " << singleFastaName << " does not exist." << endl;
                file.close();
                return 1;
            }
            int validateRes = validate_fasta_file(singleFastaName.c_str(), 1);
            if (validateRes != 0) {
                cout << "Error: " << singleFastaName << " is not a valid FASTA file." << endl;
                cout << "       Please check the entries in the file" << endl;
                file.close();
                return 1;
            }
        }
    }

    if (par.gtdb == 1) {
        editNames(taxonomyDir + "/names.dmp", par.filenames[2]);
        accession2taxid(par.filenames[1], par.filenames[2]);
        par.filenames[2] = par.filenames[2].substr(0, par.filenames[2].find_last_of('.')) + ".accession2taxid";
    }

    TaxonomyWrapper * taxonomy =  
        new TaxonomyWrapper(taxonomyDir + "/names.dmp",
                            taxonomyDir + "/nodes.dmp",
                            taxonomyDir + "/merged.dmp",
                            true);

    // Not the taxonomy-taking constructor: it builds its KmerExtractor from a MetamerPattern
    // fixed at 8 amino acids and ignores kmerFormat, so --kmer-format 5 produced 8-mer
    // metamers that grouping's 12-mer syncmers never matched.
    IndexCreator idxCre(par, par.kmerFormat);
    idxCre.setTaxonomy(taxonomy);
    idxCre.createCommonKmerIndex();
    const int kmerLen = idxCre.getKmerLen();

    // Record how the k-mers were built. The grouping workflow compares its own --syncmer /
    // --smer-len against this and refuses to run on a mismatch: the two k-mer spaces would not
    // line up and the common-k-mer filter would remove nothing, without any error. kmer_len is
    // what the extractor reported, not what was asked for.
    {
        const string paramFile = dbDir + "/kmer_params";
        ofstream out(paramFile);
        if (!out.is_open()) {
            cerr << "Error opening file: " << paramFile << endl;
            return 1;
        }
        out << "syncmer " << par.syncmer << "\n"
            << "smer_len " << par.smerLen << "\n"
            << "kmer_format " << par.kmerFormat << "\n"
            << "kmer_len " << kmerLen << "\n";
        out.close();
    }
    delete taxonomy;
    return 0;
}

#endif // METABULI_CREATE_COMMON_KMER_LIST_H