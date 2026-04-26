#include "writer.hpp"
#include "types.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

void writeAlignmentsToSAM(
    const std::vector<std::vector<ReadRecord>> &readChunks,
    const std::unordered_map<std::string, ReadBestMatches> &readBestMatchesMap,
    const std::string &samOutputPath
) {
    std::ofstream samFile(samOutputPath);
    if (!samFile.is_open()) {
        std::cerr << "Error: Could not open SAM output file: " << samOutputPath << std::endl;
        return;
    }

    // SAM consists of
    // readSequenceName \t refSeq1Name \t refSeq1Offset \t cigar1 \t readSequence \t readQuality <\t ALT:refSeq2Name,refSeq2Offset,cigar2 >
    // if no match is found, refSeq1Name and cigar1 are * and refSeq1Offset is 0. If none or only one match is found, the ALT field is omitted.
    // quality is set for now to F for all read bases.

    for (const auto &chunk : readChunks) {
        for (const auto &read : chunk) {
            const std::string &readName = std::get<0>(read);
            const std::string &readSequence = std::get<1>(read);
            const std::string readQuality(readSequence.size(), 'F'); // Placeholder quality string

            if (readBestMatchesMap.find(readName) != readBestMatchesMap.end()) {
                const ReadBestMatches &bestMatches = readBestMatchesMap.at(readName);
                const ReadAlignment &bestAlignment = bestMatches.bestAlignment;
                const ReadAlignment &secondBestAlignment = bestMatches.secondBestAlignment;

                std::string refSeq1Name = "*";
                int refSeq1Offset = 0;
                std::string cigar1 = "*";

                if (bestAlignment.editDistance != -1) {
                    refSeq1Name = bestAlignment.refName;
                    refSeq1Offset = bestAlignment.refStart;
                    cigar1 = bestAlignment.cigarString;
                }

                samFile << readName << "\t" << refSeq1Name << "\t" << refSeq1Offset << "\t" << cigar1 << "\t" << readSequence << "\t" << readQuality;

                if (secondBestAlignment.editDistance != -1) {
                    samFile << "\tALT:" << secondBestAlignment.refName << "," << secondBestAlignment.refStart << "," << secondBestAlignment.cigarString;
                }

                samFile << "\n";
            } else {
                // No matches found, write the read with default values
                samFile << readName << "\t*\t0\t*\t" << readSequence << "\t" << readQuality << "\n";
            }
        }
    }
    
    samFile.close();
}