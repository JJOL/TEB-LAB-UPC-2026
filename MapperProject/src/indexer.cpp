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
        std::cout << "Reading sequence..." << std::endl;
        std::string baseSequenceName = reader.readSequenceName(); // Get the sequence name from the reader
        std::cout << "Sequence name: " << baseSequenceName << std::endl;
        
        int chunkNumber = 0;
        while (reader.hasNextChunk()) {
            std::cout << "Reading chunk " << chunkNumber << " of sequence " << baseSequenceName << std::endl;
            std::string sequenceChunk = reader.readChunk(chunkSize);
            std::cout << "Processing chunk with size " << sequenceChunk.size() << std::endl;
            // Process the sequence chunk and create the index
            // You can use the baseSequenceName to associate the chunks with the original sequence
    
            // Create the index file with nomenclature indexPath/baseSequenceName.chunkNumber.idx
            std::ofstream indexFile(indexPath + "/" + baseSequenceName + ".chunk" + std::to_string(chunkNumber) + ".idx");
            chunkNumber++;
        }
    }

}