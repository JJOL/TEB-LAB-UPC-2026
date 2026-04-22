#ifndef INDEXER_HPP
#define INDEXER_HPP

#include <unordered_map>
#include <string>
#include <vector>

void createIndex(const std::unordered_map<std::string, std::string> &argsMap);

class SimpleKmerIndexer {
public:
    SimpleKmerIndexer(int kmerSize, const std::string &sequenceName, int chunkNumber)
        : kmerSize(kmerSize), sequenceName(sequenceName), chunkNumber(chunkNumber) {}
    void indexSequenceChunk(const std::string &sequenceChunk);
    void writeIndexToFile(std::ofstream &indexFile);
private:
    int kmerSize;
    std::string sequenceName;
    int chunkNumber;
    std::unordered_map<std::string, std::vector<int>> kmerIndex;
};

#endif // INDEXER_HPP