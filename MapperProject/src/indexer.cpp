#include "indexer.hpp"
#include "reader.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>

void createIndex(const std::unordered_map<std::string, std::string> &argsMap) {
    std::string referencePath = argsMap.at("reference");
    std::string indexPath = argsMap.at("index_path");
    std::string sequencesWhitelistString = argsMap.find("sequences_whitelist") != argsMap.end() ? argsMap.at("sequences_whitelist") : "";
    int kmerSize = std::stoi(argsMap.at("kmer_size"));
    int chunkSize = std::stoi(argsMap.at("chunk_size"));

    std::cout << "Creating index with the following parameters:" << std::endl;
    std::cout << "Reference Path: " << referencePath << std::endl;
    std::cout << "Index Path: " << indexPath << std::endl;
    std::cout << "K-mer Size: " << kmerSize << std::endl;
    std::cout << "Chunk Size: " << chunkSize << std::endl;
    std::cout << "Sequences Whitelist: " << sequencesWhitelistString << std::endl;
    std::cout << "Reading reference genome and creating index..." << std::endl;

    std::vector<std::string> whitelistSequences;
    if (!sequencesWhitelistString.empty()) {
        whitelistSequences = splitString(sequencesWhitelistString, ';');
    }

    auto reader = FASTAReader(referencePath);
    std::cout << "Handle obtained." << std::endl;

    while (reader.hasNextSequence()) {
        std::string baseSequenceName = reader.readSequenceName(); // Get the sequence name from the reader
        std::cout << "Sequence name: " << baseSequenceName << std::endl;
        if (!whitelistSequences.empty() && std::find(whitelistSequences.begin(), whitelistSequences.end(), baseSequenceName) == whitelistSequences.end()) {
            std::cout << "Skipping sequence " << baseSequenceName << " as it is not in the whitelist." << std::endl;
            while (reader.hasNextChunk()) {
                reader.readChunk(chunkSize); // Read and discard the chunks of the sequence to move the reader to the next sequence
            }
            continue;
        }
        std::cout << "Reading sequence..." << std::endl;
        
        int chunkNumber = 0;
        int totalBasesRead = 0;
        while (reader.hasNextChunk()) {
            std::cout << "Reading chunk " << chunkNumber << " of sequence " << baseSequenceName << std::endl;
            int fileChunkLine = reader.getCurrentLineNumber();
            int fileChunkLineOffset = reader.getCurrentLineOffset();
            std::string sequenceChunk = reader.readChunk(chunkSize);
            std::cout << "Processing chunk with size " << sequenceChunk.size() << std::endl;
            totalBasesRead += sequenceChunk.size();

            SimpleKmerIndexer kmerIndexer(kmerSize, baseSequenceName, chunkNumber, fileChunkLine, fileChunkLineOffset, sequenceChunk.size());
            kmerIndexer.indexSequenceChunk(sequenceChunk);
            // Process the sequence chunk and create the index
            // You can use the baseSequenceName to associate the chunks with the original sequence
            std::cout << "Processing done." << std::endl;
            std::cout << "Writing to file..." << std::endl;

            // Create the index file with nomenclature indexPath/baseSequenceName.chunkNumber.idx. It will be a binary file that contains the kmer-interger positions map
            std::ofstream indexFile(indexPath + "_" + baseSequenceName + "." + std::to_string(chunkNumber) + ".idx");
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

void SimpleKmerIndexer::_writeHumanReadableIndexFile(std::ostream &indexFile) {
    indexFile << "#Kmers:" << kmerIndex.size() << std::endl;
    int maxPositions = 0;
    unsigned long long totalPositions = 0;
    for (auto pair : kmerIndex) {
        // indexFile << pair.first << ":" << pair.second.size() << std::endl;
        if (pair.second.size() > maxPositions) {
            maxPositions = pair.second.size();
        }
        totalPositions += pair.second.size();
    }
    float avgPositions = (float)totalPositions / kmerIndex.size();
    std::cout << "Max positions for a kmer: " << maxPositions << std::endl;
    std::cout << "Average positions for a kmer: " << avgPositions << std::endl;
    std::cout << "Wrote to file" << std::endl;
}

void SimpleKmerIndexer::_writeEfficientLaterLoadIndexFile(std::ofstream &indexFile) {
    // Placeholder for writing the index in a more efficient format for later loading
    // You can implement a binary format or a more compact representation here
    // for each kmer, write the kmer string, followed by number of positions, followed by the positions as binary data.
    // First we write the reference sequence name and chunk number as a header, then we write the kmer index.
    // write sequenceName with trailing null character to ensure it can be read correctly later
    // write sequence retrieval information as a header. Reference sequence name (0 ending string), chunk number (4 bytes), file chunk line (4 bytes), file chunk line offset (4 bytes), total bases read (4 bytes)
    indexFile.write(sequenceName.c_str(), sequenceName.size() + 1); // Write the reference sequence name as a header
    indexFile.write(reinterpret_cast<const char*>(&chunkNumber), sizeof(chunkNumber)); // Write the chunk number as a header
    indexFile.write(reinterpret_cast<const char*>(&fileChunkLine), sizeof(fileChunkLine)); // Write the file chunk line as a header
    indexFile.write(reinterpret_cast<const char*>(&fileChunkLineOffset), sizeof(fileChunkLineOffset)); // Write the file chunk line offset as a header
    indexFile.write(reinterpret_cast<const char*>(&totalBasesRead), sizeof(totalBasesRead)); // Write the total bases read as a header
     std::cout << "Read header from file. Sequence name: " << sequenceName << ", chunk number: " << chunkNumber << ", file chunk line: " << fileChunkLine << ", file chunk line offset: " << fileChunkLineOffset << ", total bases read: " << totalBasesRead << std::endl;

    int n = kmerIndex.size();
    int i = 0;
    for (auto pair : kmerIndex) {
        if (i == n) break;
        i++;
        indexFile.write(pair.first.c_str(), pair.first.size()); // Write the kmer string
        int count = pair.second.size();
        indexFile.write(reinterpret_cast<const char*>(&count), sizeof(count)); // Write the number of positions
        indexFile.write(reinterpret_cast<const char*>(pair.second.data()), count * sizeof(int)); // Write the positions as binary data
    }
    std::cout << "Wrote  " << i << " kmers to file" << std::endl;
}

void SimpleKmerIndexer::writeIndexToFile(std::ofstream &indexFile) {
    _writeHumanReadableIndexFile(std::cout);
    _writeEfficientLaterLoadIndexFile(indexFile);
    std::cout << "Wrote to file" << std::endl;
}

void SimpleKmerIndexer::_readEfficientLaterLoadIndexFile(std::ifstream &indexFile) {
    // read sequence name and chunk number from the header
    char sequenceNameBuffer[1000]; // Assuming the sequence name will not exceed 1000 characters, you can adjust this size as needed
    for (int i = 0; i < 1000; i++) sequenceNameBuffer[i] = 0; // initialize the buffer with zeros to ensure null termination
    // read the sequence name as a null-terminated string
    indexFile.getline(sequenceNameBuffer, 1000, '\0'); // Read the reference sequence name from the header
    sequenceName = std::string(sequenceNameBuffer);
    indexFile.read(reinterpret_cast<char*>(&chunkNumber), sizeof(chunkNumber));
    indexFile.read(reinterpret_cast<char*>(&fileChunkLine), sizeof(fileChunkLine));
    indexFile.read(reinterpret_cast<char*>(&fileChunkLineOffset), sizeof(fileChunkLineOffset));
    indexFile.read(reinterpret_cast<char*>(&totalBasesRead), sizeof(totalBasesRead));
    // we know first kmer bytes are the kmer string, then we have 4 bytes for the number of positions, then we have the positions as binary data.
    // we will continue reading until we reach the end of the file.
    // probably we can accelerate by reading a big chunk of the file into memory and then parsing it, instead of reading kmer by kmer from the file.
    // first let's implement the simple version that reads kmer by kmer from the file.
    const int KMER_MAX_SIZE = 100; // we know that the kmer size is at most 100, so we can use a fixed size buffer to read the kmer strings.
    char kmerString[KMER_MAX_SIZE + 1]; // +1 for null terminator
    for (int i = 0; i < KMER_MAX_SIZE + 1; i++) kmerString[i] = 0; // initialize the buffer with zeros to ensure null termination
    
    int count;
    int kmersRead = 0;
    while (indexFile.read(kmerString, kmerSize)) {
        if (kmersRead % 50000 == 0) {
            std::cout << "Read " << kmersRead << " kmers from file..." << std::endl;
        }
        indexFile.read(reinterpret_cast<char*>(&count), sizeof(count)); // Read the number of positions
        std::vector<int> positions(count);
        indexFile.read(reinterpret_cast<char*>(positions.data()), count * sizeof(int)); // Read the positions as binary data
        kmerIndex[std::string(kmerString)] = positions;
        kmersRead++;
    }
}

void SimpleKmerIndexer::readIndexFromFile(const std::string &indexFilePath) {
    // Placeholder for reading the index from a file
    // You can implement the logic to read the index file and populate the kmerIndex map here
    std::cout << "Reading index from file: " << indexFilePath << std::endl;
    std::ifstream indexFile(indexFilePath, std::ios::binary);
    if (!indexFile.is_open()) {
        std::cerr << "Error: Could not open index file " << indexFilePath << std::endl;
        return;
    }
    _readEfficientLaterLoadIndexFile(indexFile);
    std::cout << "Finished reading index from file. Total kmers read: " << kmerIndex.size() << std::endl;
}