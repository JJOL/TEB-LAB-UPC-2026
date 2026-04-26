#include "aligner.hpp"
#include <iostream>
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void printDPTable(int *dpCost, int m, int n, const std::string &read, const std::string &reference) {
    std::cout << "DP Table:" << std::endl;
    std::cout << "\t\t";
    for (int j = 0; j < m; j++) {
        std::cout << reference[j] << "\t";
    }
    std::cout << std::endl;
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            std::cout << read[i-1] << "\t";
        } else {
            std::cout << "\t";
        }
        for (int j = 0; j < m; j++) {
            std::cout << dpCost[i * m + j] << "\t";
        }
        std::cout << std::endl;
    }
}

int SUBSTITUION_MATRIX[4][4] = {
    {0, 1, 1, 1}, // A
    {1, 0, 1, 1}, // C
    {1, 1, 0, 1}, // G
    {1, 1, 1, 0}  // T
};

std::string compressCigar(std::vector<char> &cigarOps) {
    std::string compressedCigar;
    if (cigarOps.empty()) {
        return compressedCigar;
    }

    char currentOp = cigarOps[0];
    int count = 1;

    for (size_t i = 1; i < cigarOps.size(); ++i) {
        if (cigarOps[i] == currentOp) {
            count++;
        } else {
            compressedCigar += std::to_string(count) + currentOp;
            currentOp = cigarOps[i];
            count = 1;
        }
    }
    compressedCigar += std::to_string(count) + currentOp; // Add the last operation

    return compressedCigar;
}

std::string compressCigar(const std::string &cigarString) {
    std::vector<char> cigarOps(cigarString.begin(), cigarString.end());
    return compressCigar(cigarOps);
}

ReadAlignment alignReadToReference(const std::string &read, const std::string &reference) {
    ReadAlignment alignment;
    alignment.editDistance = 0; // Example value, replace with actual edit distance
    alignment.readStart = 0; // Example value, replace with actual read start position
    alignment.refStart = 0; // Example value, replace with actual reference start position
    alignment.alignmentLength = 0;
    alignment.cigarString = "";


    int m = reference.size();
    int n = read.size();

    // std::cout << "Aligning read of length " << n << " to reference of length " << m << std::endl;

    // int dpCostPtr[n + 1][m + 1];
    int *dpCostPtr = new int[(n + 1) * (m + 1)];
    // int dp_parent[n + 1][m + 1];
    int *dpParent = new int[(n + 1) * (m + 1)];
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            dpCostPtr[i * (n + 1) + j] = 0;
            dpParent[i * (n + 1) + j] = -1;
        }
    }

    // Semi-global aligment init.
    for (int i = 0; i <= m; i++) {
        dpCostPtr[0 * (m + 1) + 0] = 0; // No cost for leading gaps in the reference
    }
    for (int j = 0; j <= n; j++) {
        dpCostPtr[j * (m + 1) + 0] = j; // Cost of inserting all characters of the read
    }

    // ptr to dp_cost

    // std::cout << "Initial DP Table of size: " << (n + 1) << " x " << (m + 1) << std::endl;
    // printDPTable(dpCostPtr, m+1, n+1, read, reference);

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            // First Align Priority: Match/Mismatch
            int matchCost = dpCostPtr[(i - 1) * (m + 1) + (j - 1)] + (read[i - 1] == reference[j - 1] ? 0 : 1);
            int insertCost = dpCostPtr[(i - 1) * (m + 1) + j] + 1; // Insertion cost
            int deleteCost = dpCostPtr[i * (m + 1) + (j - 1)] + 1; // Deletion cost

            dpCostPtr[i * (m + 1) + j] = MIN(MIN(matchCost, insertCost), deleteCost);
            dpParent[i * (m + 1) + j] = (dpCostPtr[i * (m + 1) + j] == matchCost) ? 0 : 
                                        (dpCostPtr[i * (m + 1) + j] == insertCost) ? 1 : 2; // 0: match/mismatch, 1: insertion, 2: deletion
        }
    }
    // printDPTable(dpCostPtr, m+1, n+1, read, reference);

    int minAlignmentEndPos = 0;
    for (int j = 0; j <= m; j++) {
        if (dpCostPtr[n * (m + 1) + j] < dpCostPtr[n * (m + 1) + minAlignmentEndPos]) {
            minAlignmentEndPos = j;
        }
    }

    // recontruct path
    int i = n;
    int j = minAlignmentEndPos;
    std::vector<char> cigarOps;
    while (i > 0 && j > 0) {
        int parent = dpParent[i * (m + 1) + j];
        if (parent == 0) { // match/mismatch
            i--;
            j--;
            cigarOps.push_back((read[i] == reference[j]) ? 'M' : 'X'); // M for match, X for mismatch
        } else if (parent == 1) { // insertion
            i--;
            cigarOps.push_back('I'); // I for insertion
        } else { // deletion
            j--;
            cigarOps.push_back('D'); // D for deletion
        }
    }

    alignment.editDistance = dpCostPtr[n * (m + 1) + minAlignmentEndPos];
    alignment.readStart = i;
    alignment.refStart = j;
    alignment.alignmentLength = n - i; // Assuming the alignment length is the length of the read aligned
    alignment.cigarString = ""; // CIGAR string reconstruction can be implemented
    for (auto it = cigarOps.rbegin(); it != cigarOps.rend(); ++it) {
        alignment.cigarString += *it;
    }
    alignment.cigarString = compressCigar(alignment.cigarString);


    delete[] dpCostPtr;
    delete[] dpParent;

    return alignment;
}

ReadAlignment dpSimpleSemiGlobalEditDistance(const std::string &readSequence, const std::string &refSequence) {
    return alignReadToReference(readSequence, refSequence);
}