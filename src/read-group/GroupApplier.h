#ifndef GROUP_APPLIER_H
#define GROUP_APPLIER_H

// apply-group: give every read of a group the group's representative taxonomy label.
// Ported from the Metabuli_work fork; the output is written here directly instead of
// through Reporter, whose grouped writer only exists in that fork.

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstdint>
#include "LocalParameters.h"
#include "TaxonomyWrapper.h"

struct OrgResult {
    int label;        // classification as written in the input (original taxonomy id)
    float score;
    int queryLength;
    std::string name;
};

class GroupApplier {
protected:
    const LocalParameters & par;
    std::string groupFileDir;
    std::string groupmapFileDir;
    std::string taxDbDir;
    std::string orgRes;
    std::string outDir;

    TaxonomyWrapper *taxonomy = nullptr;

    // Output
    std::string updatedResultFileName;

public:
    explicit GroupApplier(LocalParameters & par);

    ~GroupApplier();

    void startGroupApplication(const LocalParameters & par);

    void loadOrgResult(std::vector<OrgResult>& orgResults, size_t& processedReadCnt);

    void loadGroupsFromFile(std::unordered_map<uint32_t, std::unordered_set<uint32_t>> &groupInfo,
                            std::vector<uint32_t> &groupMappingInfo,
                            const std::string &groupInfoFileName,
                            const std::string &groupMappingInfoFileName);

    void getRepLabel(std::vector<OrgResult>& orgResults,
                     const std::unordered_map<uint32_t, std::unordered_set<uint32_t>>& groupInfo,
                     std::unordered_map<uint32_t, uint32_t>& repLabel,
                     std::unordered_map<int, int>& external2internalTaxId);

    void applyRepLabel(const std::vector<OrgResult>& orgResults,
                       const std::vector<uint32_t>& queryGroupInfo,
                       const std::unordered_map<uint32_t, uint32_t>& repLabel,
                       std::unordered_map<int, int>& external2internalTaxId);
};

#endif // GROUP_APPLIER_H
