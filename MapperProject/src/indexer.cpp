#include "indexer.hpp"
#include "reader.hpp"
#include <iostream>
#include <fstream>

void createIndex(const std::unordered_map<std::string, std::string> &argsMap) {
    // Placeholder for index creation logic
    // You can access the arguments from argsMap, for example:
    std::string referencePath = argsMap.at("reference");
    std::string indexPath = argsMap.at("index_path");
    int kmerSize = std::stoi(argsMap.at("kmer_size"));
    int chunkSize = std::stoi(argsMap.at("chunk_size"));

    // Implement the indexing logic here using the provided arguments
    std::cout << "Creating index with the following parameters:" << std::endl;
    std::cout << "Reference Path: " << referencePath << std::endl;
    std::cout << "Index Path: " << indexPath << std::endl;
    std::cout << "K-mer Size: " << kmerSize << std::endl;
    std::cout << "Chunk Size: " << chunkSize << std::endl;

    std::cout << "Reading reference genome and creating index..." << std::endl;
    auto reader = FASTAReader(referencePath);
    std::cout << "Handle obtained." << std::endl;

    // NOTE: For now we now that we only have 1 sequence in the reference for chromosome 1.
    while (reader.hasNextSequence()) {
        std::string baseSequenceName = reader.readSequenceName(); // Get the sequence name from the reader
        std::cout << "Sequence name: " << baseSequenceName << std::endl;
        std::cout << "Reading sequence..." << std::endl;
        
        int chunkNumber = 0;
        int totalBasesRead = 0;
        while (reader.hasNextChunk()) {
            std::cout << "Reading chunk " << chunkNumber << " of sequence " << baseSequenceName << std::endl;
            std::string sequenceChunk = reader.readChunk(chunkSize);
            std::cout << "Processing chunk with size " << sequenceChunk.size() << std::endl;
            totalBasesRead += sequenceChunk.size();

            SimpleKmerIndexer kmerIndexer(kmerSize, baseSequenceName, chunkNumber);
            kmerIndexer.indexSequenceChunk(sequenceChunk);
            // Process the sequence chunk and create the index
            // You can use the baseSequenceName to associate the chunks with the original sequence
            std::cout << "Processing done." << std::endl;
            std::cout << "Writing to file..." << std::endl;

            // Create the index file with nomenclature indexPath/baseSequenceName.chunkNumber.idx
            std::ofstream indexFile(indexPath + "/" + baseSequenceName + ".chunk" + std::to_string(chunkNumber) + ".idx");
            if (!indexFile.is_open()) {
                std::cerr << "Error: Could not create index file for chunk " << chunkNumber << " of sequence " << baseSequenceName << std::endl;
                continue; // Skip to the next chunk
            }
            kmerIndexer.writeIndexToFile(indexFile);
            indexFile.close();
            chunkNumber++;
        }
    }
}

void SimpleKmerIndexer::indexSequenceChunk(const std::string &sequenceChunk) {
    std::cout << "Indexing a sequence chunk of size " << sequenceChunk.size() << " with k-mer size " << kmerSize << std::endl;
    int timesN = 0, localN = 0;
    int totalKmers = sequenceChunk.size() - kmerSize + 1;
    for (int i = 0; i < totalKmers; i++) { // goes from 0 to 10M - 10 + 1 = 9,999,991
        if (i % 500000 == 0) {
            std::cout << i << "/" << sequenceChunk.size() << "..." << std::endl;
        }
        kmerIndex[sequenceChunk.substr(i, kmerSize)].push_back(i);
    }
}

void _writeHumanReadableIndexFile(std::ofstream &indexFile, const std::unordered_map<std::string, std::vector<int>> &kmerIndex) {
    indexFile << "#Kmers:" << kmerIndex.size() << std::endl;
    for (auto pair : kmerIndex) {
        indexFile << pair.first << ":" << pair.second.size() << std::endl;
    }
    std::cout << "Wrote to file" << std::endl;
}

void _writeEfficientLaterLoadIndexFile(std::ofstream &indexFile, const std::unordered_map<std::string, std::vector<int>> &kmerIndex) {
    // Placeholder for writing the index in a more efficient format for later loading
    // You can implement a binary format or a more compact representation here

}

void SimpleKmerIndexer::writeIndexToFile(std::ofstream &indexFile) {
    _writeHumanReadableIndexFile(indexFile, kmerIndex);
    // _writeEfficientLaterLoadIndexFile(indexFile, kmerIndex);
    std::cout << "Wrote to file" << std::endl;
}