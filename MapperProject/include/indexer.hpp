#ifndef INDEXER_HPP
#define INDEXER_HPP

#include <unordered_map>
#include <string>
#include <vector>

void createIndex(const std::unordered_map<std::string, std::string> &argsMap);

class SimpleKmerIndexer {
public:
    SimpleKmerIndexer(int kmerSize, const std::string &sequenceName, int chunkNumber, int fileChunkLine, int fileChunkLineOffset, int totalBasesRead)
        : kmerSize(kmerSize), sequenceName(sequenceName), chunkNumber(chunkNumber), fileChunkLine(fileChunkLine), fileChunkLineOffset(fileChunkLineOffset), totalBasesRead(totalBasesRead) {}
    void indexSequenceChunk(const std::string &sequenceChunk);
    void writeIndexToFile(std::ofstream &indexFile);
    void readIndexFromFile(const std::string &indexFilePath);
    std::string getSequenceName() const { return sequenceName; }
    int getChunkNumber() const { return chunkNumber; }
    int getKmerSize() const { return kmerSize; }
public:
    std::unordered_map<std::string, std::vector<int>> kmerIndex;
private:
    void _writeEfficientLaterLoadIndexFile(std::ofstream &indexFile);
    void _writeHumanReadableIndexFile(std::ostream &indexFile);
    void _readEfficientLaterLoadIndexFile(std::ifstream &indexFile);
private:
    int kmerSize;
    std::string sequenceName;
    int chunkNumber;
public:
    int fileChunkLine;
    int fileChunkLineOffset;
    int totalBasesRead;
};

#endif // INDEXER_HPP