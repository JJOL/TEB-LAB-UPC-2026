#ifndef TYPES_HPP
#define TYPES_HPP
#include <string>
#include <vector>
#include <tuple>

// readQuality is only used to write to SAM output. It is not used in the mapping process, so it can be ignored in the FASTA reader.
typedef std::tuple<std::string, std::string, std::string> ReadRecord; // (readName, readSequence, readQuality)
typedef std::tuple<int, int, int> ReadSeedHit; // (readId, readPosition, targetPosition)

struct ChromosomeInfo {
    std::string name;
    size_t length;
    size_t offset;
    ChromosomeInfo(const std::string &name, size_t length, size_t offset) : name(name), length(length), offset(offset) {}
};

struct ReadAlignment {
    int editDistance;
    int readStart;
    int refStart;
    int alignmentLength;
    std::string cigarString;
    std::string readName;
    std::string refName;

    ReadAlignment() : editDistance(-1), readStart(-1), refStart(-1), alignmentLength(-1), cigarString(""), readName(""), refName("") {}
};

struct ReadBestMatches {
    ReadAlignment bestAlignment;
    ReadAlignment secondBestAlignment;
    
    void updateBestAlignmentAgainst(const ReadAlignment &newAlignment) {
        if (bestAlignment.editDistance == -1 || newAlignment.editDistance < bestAlignment.editDistance) {
            secondBestAlignment = bestAlignment;
            bestAlignment = newAlignment;
        } else if (secondBestAlignment.editDistance == -1 || newAlignment.editDistance < secondBestAlignment.editDistance) {
            secondBestAlignment = newAlignment;
        }
    }
};

#endif