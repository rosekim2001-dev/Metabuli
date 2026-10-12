#include "GroupApplier.h"
#include "FileUtil.h"
#include "common.h"

#include <fstream>
#include <iostream>
#include <sstream>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
using std::unordered_set;

GroupApplier::GroupApplier(LocalParameters & par) : par(par) {
    groupFileDir    = par.filenames[0];
    groupmapFileDir = par.filenames[1];
    taxDbDir        = par.filenames[2];
    orgRes          = par.filenames[3];
    outDir          = par.filenames[4];

    taxonomy = new TaxonomyWrapper(
                    taxDbDir + "/names.dmp",
                    taxDbDir + "/nodes.dmp",
                    taxDbDir + "/merged.dmp",
                    true);

    updatedResultFileName = outDir + "/updated_classifications.tsv";
}

GroupApplier::~GroupApplier() {
    delete taxonomy;
}

void GroupApplier::startGroupApplication(const LocalParameters &par) {
    size_t processedReadCnt = 0;

    vector<OrgResult> orgResult;
    loadOrgResult(orgResult, processedReadCnt);

    unordered_map<uint32_t, unordered_set<uint32_t>> groupInfo;
    vector<uint32_t> groupMappingInfo;
    loadGroupsFromFile(groupInfo, groupMappingInfo, groupFileDir, groupmapFileDir);

    // applyRepLabel walks orgResults and indexes groupMappingInfo by the same counter, so the two
    // have to describe the same reads in the same order. A shorter map is an unchecked read past
    // the end of a vector; a longer one means the two files are from different runs and every
    // label would be applied to the wrong read. Neither can be recovered from here.
    if (groupMappingInfo.size() != processedReadCnt) {
        cerr << "Error: " << groupmapFileDir << " has " << groupMappingInfo.size()
             << " read-to-group rows but " << orgRes << " has " << processedReadCnt
             << " classified reads." << endl;
        cerr << "       These must be the same reads in the same order -- the grouping run and"
             << " the classification run have to share their input." << endl;
        exit(EXIT_FAILURE);
    }

    std::unordered_map<int, int> external2internalTaxId;
    taxonomy->getExternal2internalTaxID(external2internalTaxId);
    unordered_map<uint32_t, uint32_t> repLabel;
    getRepLabel(orgResult, groupInfo, repLabel, external2internalTaxId);
    applyRepLabel(orgResult, groupMappingInfo, repLabel, external2internalTaxId);
    (void) par;
}


void GroupApplier::loadOrgResult(vector<OrgResult>& orgResults, size_t& processedReadCnt) {
    std::ifstream inFile(orgRes);
    if (!inFile.is_open()) {
        cerr << "Error opening file: " << orgRes << endl;
        return;
    }

    const int classificationCol = par.taxidCol - 1;
    const int scoreCol = par.scoreCol - 1;
    const int readNameCol = par.readIdCol - 1;
    // Metabuli's read-by-read result keeps the query length in the 4th column.
    const int lengthCol = 3;
    const bool useScore = !(par.weightMode == 0 || scoreCol < 0);
    string line;
    while (getline(inFile, line)) {
        if (line.empty()) continue;
        if (line.front() == '#') continue;
        std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 20);
        OrgResult result;
        result.label = stoi(columns[classificationCol]);
        result.score = useScore ? stof(columns[scoreCol]) : 1.0f;
        result.queryLength = (static_cast<int>(columns.size()) > lengthCol) ? atoi(columns[lengthCol].c_str()) : 0;
        result.name = columns[readNameCol];
        orgResults.push_back(std::move(result));
        processedReadCnt++;
    }
    inFile.close();
    cout << "Original Metabuli result loaded from " << orgRes << " successfully." << endl;
    cout << "Number of query result: " << orgResults.size() << endl;
}

void GroupApplier::loadGroupsFromFile(unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                                      vector<uint32_t> &groupMappingInfo,
                                      const string &groupInfoFileName,
                                      const string &groupMappingInfoFileName) {
    // 1. Load groupInfo
    std::ifstream inFile1(groupInfoFileName);
    if (!inFile1.is_open()) {
        cerr << "Error opening file: " << groupInfoFileName << endl;
        return;
    }

    string line;
    while (getline(inFile1, line)) {
        std::istringstream ss(line);
        uint32_t groupId;
        ss >> groupId;

        uint32_t queryId;
        while (ss >> queryId) {
            groupInfo[groupId].insert(queryId - 1); // convert to 0-based index
        }
    }
    inFile1.close();
    cout << "Group info loaded from " << groupInfoFileName << " successfully." << endl;

    // 2. Load groupMappingInfo
    std::ifstream inFile2(groupMappingInfoFileName);
    if (!inFile2.is_open()) {
        cerr << "Error opening file: " << groupMappingInfoFileName << endl;
        return;
    }

    groupMappingInfo.clear();
    int groupId;
    string readId;
    while (inFile2 >> readId >> groupId) {
        groupMappingInfo.emplace_back(groupId);
    }
    inFile2.close();
    cout << "Query-to-group map loaded from " << groupMappingInfoFileName << " successfully." << endl;
}

void GroupApplier::getRepLabel(vector<OrgResult> &orgResults,
                               const unordered_map<uint32_t, unordered_set<uint32_t>> &groupInfo,
                               unordered_map<uint32_t, uint32_t> &repLabel,
                               std::unordered_map<int, int>& external2internalTaxId) {
    cout << "Find query group representative labels..." << endl;
    time_t beforeSearch = time(nullptr);

    for (const auto& group : groupInfo) {
        uint32_t groupId = group.first;
        const unordered_set<uint32_t>& queryIds = group.second;

        vector<WeightedTaxHit> setTaxa;

        for (const auto& queryId : queryIds) {
            int query_label = external2internalTaxId[orgResults[queryId].label];
            if (par.weightMode == 0) {
                float score = 1;
                if (query_label != 0) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 1) {
                float score = orgResults[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score, 2);
                }
            } else if (par.weightMode == 2) {
                float score = orgResults[queryId].score;
                if (query_label != 0 && score >= par.minVoteScr) {
                    setTaxa.emplace_back(query_label, score * score, 2);
                }
            }
        }

        WeightedTaxResult result = taxonomy->weightedMajorityLCA(setTaxa, 0.5);

        if (result.taxon != 0 && result.taxon != 1) {
            repLabel[groupId] = result.taxon;
        }
        else{
            repLabel[groupId] = 0;
        }
    }

    const string& groupRepFileName = outDir + "/groupRep";
    std::ofstream outFile(groupRepFileName);
    if (!outFile.is_open()) {
        cerr << "Error opening file: " << groupRepFileName << endl;
        return;
    }

    for (const auto& [groupId, groupRep] : repLabel) {
        outFile << groupId << "\t" << taxonomy->getOriginalTaxID(groupRep) << "\n";
    }

    outFile.close();

    cout << "Query group representative label saved to " << groupRepFileName << " successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}

void GroupApplier::applyRepLabel(const vector<OrgResult>& orgResults,
                                 const vector<uint32_t>& groupMappingInfo,
                                 const unordered_map<uint32_t, uint32_t>& repLabel,
                                 std::unordered_map<int, int>& external2internalTaxId) {
    cout << "Apply query group representative labels..." << endl;
    time_t beforeSearch = time(nullptr);

    std::ofstream out(updatedResultFileName);
    if (!out.is_open()) {
        cerr << "Error opening file: " << updatedResultFileName << endl;
        return;
    }
    // Same columns as Metabuli's read-by-read result plus the group id ("-" = no group).
    out << "#is_classified\tname\ttaxID\tquery_length\tscore\trank\tgroup";
    if (par.printLineage) {
        out << "\tlineage";
    }
    out << "\ttaxID:match_count\n";

    size_t relabeled = 0;
    for (size_t queryIdx = 0; queryIdx < orgResults.size(); queryIdx++) {
        const OrgResult &org = orgResults[queryIdx];
        const uint32_t groupId = groupMappingInfo[queryIdx];
        auto repLabelIt = repLabel.find(groupId);
        const bool hasRepLabel = (repLabelIt != repLabel.end() && repLabelIt->second != 0);

        // Internal taxonomy id of the label this read ends up with; 0 = unclassified.
        int label = 0;
        if (hasRepLabel) {
            const bool applyRep = (par.weightMode == 0)
                || (org.label == 0 || org.score < par.minVoteScr);
            if (applyRep) {
                label = static_cast<int>(repLabelIt->second);
                ++relabeled;
            } else {
                label = external2internalTaxId[org.label];
            }
        } else if (org.label != 0) {
            label = external2internalTaxId[org.label];
        }

        const string groupText = (groupId == 0) ? string("-") : std::to_string(groupId);
        if (label != 0) {
            const TaxonNode *node = taxonomy->taxonNode(label);
            out << "1\t" << org.name << '\t' << taxonomy->getOriginalTaxID(label) << '\t'
                << org.queryLength << '\t' << org.score << '\t'
                << taxonomy->getString(node->rankIdx) << '\t' << groupText;
            if (par.printLineage) {
                out << '\t' << taxonomy->taxLineage2(node);
            }
            out << "\t\n";
        } else {
            out << "0\t" << org.name << "\t0\t" << org.queryLength << '\t' << org.score
                << "\t-\t" << groupText;
            if (par.printLineage) {
                out << "\t-";
            }
            out << "\t-\n";
        }
    }
    out.close();

    cout << "Reads relabeled with a group label: " << relabeled << "/" << orgResults.size() << endl;
    cout << "Result saved to " << updatedResultFileName << " successfully." << endl;
    cout << "Time spent: " << double(time(nullptr) - beforeSearch) << " seconds." << endl;
}
