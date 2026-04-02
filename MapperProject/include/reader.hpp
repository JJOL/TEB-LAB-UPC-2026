#ifndef READER_HPP
#define READER_HPP

#include <string>
#include <fstream>

class FASTAReader {
public:
    FASTAReader(const std::string &filePath);
    bool hasNextSequence() const;
    bool hasNextChunk() const;
    std::string readChunk(size_t chunkSize);
    std::string readSequenceName();
private:
    int nSeqs;
    int basesReadInChunk;
    int currentChunk;
    int totalBasesRead;
    std::string currentLine;
    std::ifstream file;
};

#endif // READER_HPP