#ifndef READER_HPP
#define READER_HPP

#include <string>
#include <fstream>
#include <vector>
#include "types.hpp"

class FASTAReader {
public:
    FASTAReader(const std::string &filePath);
    ~FASTAReader();
    bool hasNextSequence() const;
    bool hasNextChunk() const;
    std::string readChunk(size_t chunkSize);
    std::string readSequenceName();
    int getCurrentLineNumber() const { return currentLineNumber; }
    int getCurrentLineOffset() const { return currentLineOffset; }
private:
    int nSeqs;
    int basesReadInChunk;
    int currentChunk;
    int totalBasesRead;
    int currentLineNumber; // Line number of the current line being read
    int currentLineOffset; // Offset in the current line of the next base to read
    std::string currentLine;
    std::ifstream file;

    bool pureAlphabetMode;
};

class FASTQReader {
public:
    FASTQReader(const std::string &filePath);
    ~FASTQReader();
    bool hasNext() const;
    std::vector<ReadRecord> readChunk(size_t chunkSize);
    int getNumReads() const;
    int getNumChunks() const;
private:
    int nReads;
    int nChunks;
    std::string currentLine;
    std::ifstream file;
};

std::vector<ChromosomeInfo> findChromosomesInFASTA(const std::string &filePath);

#endif // READER_HPP