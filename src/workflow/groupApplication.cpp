#include "GroupApplier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"

void setGroupApplicationDefaults(LocalParameters & par){    
    par.seqMode = 2;    
    par.ramUsage = 128;
    par.scoreCol = 5;
    par.readIdCol = 2;
    par.taxidCol = 3;
    par.weightMode = 1; // 0: uniform, 1: score, 2: score^2
    par.minVoteScr = 0.15;
}

int groupApplication(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setGroupApplicationDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.weightMode != 0) {
        cout << "Warning: --weight-mode " << par.weightMode << " requires classification scores." << endl;
        cout << "         Make sure that score column is correctly set using --score-col." << endl;
    }

    // apply-group takes exactly five arguments -- group result, group mapping result, taxonomy
    // directory, read-by-read result, output directory -- and GroupApplier's constructor reads
    // filenames[0..4]. The output directory is therefore always index 4.
    //
    // This used to branch on --seq-mode and, for the default of 2, index filenames[5] and check
    // filenames[1] for being a directory. Both were wrong: there is no sixth argument, so
    // filenames[5] was an unchecked std::vector::operator[] past the end -- undefined behaviour
    // that happened not to crash -- and filenames[1] is the group mapping FILE, so the check
    // rejected the correct input rather than a wrong one. The branch was copied from
    // groupGeneration, where --seq-mode does shift the argument positions. apply-group never
    // reads par.seqMode at all (nothing in GroupApplier references it).
    if (!FileUtil::directoryExists(par.filenames[4].c_str())) {
        FileUtil::makeDir(par.filenames[4].c_str());
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif    
    GroupApplier * groupApplier = new GroupApplier(par);
    groupApplier->startGroupApplication(par);
    delete groupApplier;
    return 0;
}