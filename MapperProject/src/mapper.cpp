#include "mapper.hpp"
#include "reader.hpp"
#include "types.hpp"
#include "indexer.hpp"
#include "aligner.hpp"
#include "writer.hpp"
#include "utils.hpp"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <limits>
#include <algorithm>
#include <limits>

std::vector<std::string> collectIndexFilePaths(const std::string &indexPath) {
    std::vector<std::string> indexFiles;

    std::string parentDirectory = indexPath.substr(0, indexPath.find_last_of('/'));
    std::cout << "Parent Directory: " << parentDirectory << std::endl;
    std::string fileNamePrefix = indexPath.substr(indexPath.find_last_of('/') + 1);

    for (const auto &entry : std::filesystem::directory_iterator(parentDirectory)) {
        if (entry.is_regular_file()) {
            std::string filePath = entry.path().string();
            std::string fileName = entry.path().filename().string();
            if (fileName.rfind(fileNamePrefix, 0) == 0) { // Check if fileName starts with fileNamePrefix
                indexFiles.push_back(filePath);
            }
        }
    }

    // sort the index files lexicographically to ensure chunk0 is always first to test the mapping.
    std::sort(indexFiles.begin(), indexFiles.end());

    return indexFiles;
}

std::vector<std::string> filterIndexFilesByWhitelist(const std::vector<std::string> &indexFiles, const std::vector<std::string> &whitelistSequences) {
    if (whitelistSequences.empty()) {
        return indexFiles; // No whitelist provided, return all index files
    }

    std::vector<std::string> filteredIndexFiles;
    for (const auto &file : indexFiles) {
        for (const auto &whitelistSeq : whitelistSequences) {
            if (file.find(whitelistSeq) != std::string::npos) {
                filteredIndexFiles.push_back(file);
                break; // Move to the next file after a match is found
            }
        }
    }

    return filteredIndexFiles;
}

std::string loadReferenceSequence(const std::string &referencePath, const std::string &refSequenceName, int chunkSize = 10000000) {
    std::string refSequence;
    FASTAReader fastaReader(referencePath);
    std::cout << "Loading reference sequence: " << refSequenceName << " from file: " << referencePath << std::endl;
    while (fastaReader.hasNextSequence()) {
        std::string sequenceName = fastaReader.readSequenceName();
        if (sequenceName.compare(refSequenceName) == 0) {
            std::cout << "Found matching reference sequence: " << sequenceName << std::endl;
            refSequence = "";
            // std::vector<std::string> chunks;
            while (fastaReader.hasNextChunk()) {
                std::string chunk = fastaReader.readChunk(chunkSize);
                // chunks.push_back(chunk);
                refSequence += chunk;
            }

            // refSequence = std::accumulate(chunks.begin(), chunks.end(), std::string());
            // std::cout << "Loaded reference sequence with length: " << refSequence.size() << std::endl;
            break;
        } else {
            while (fastaReader.hasNextChunk()) {
                fastaReader.readChunk(chunkSize); // Skip chunks of non-matching sequences
            }
        }
    }
    return refSequence;
}

std::string loadReferenceSequenceChunk(const std::string &referencePath, const std::string &refSequenceName, int chunkId, int chunkSize) {
    std::string refSequence;
    FASTAReader fastaReader(referencePath);

    std::cout << "Loading reference sequence chunk: " << refSequenceName << ", chunk ID: " << chunkId << " from file: " << referencePath << std::endl;
    while (fastaReader.hasNextSequence()) {
        std::string sequenceName = fastaReader.readSequenceName();
        if (sequenceName.compare(refSequenceName) == 0) {
            std::cout << "Found matching reference sequence: " << sequenceName << std::endl;
            int currentChunkId = 0;
            while (fastaReader.hasNextChunk()) {
                std::string chunk = fastaReader.readChunk(chunkSize); // Use the same
                if (currentChunkId == chunkId) {
                    refSequence = chunk;
                    break;
                }
                currentChunkId++;
            }
            break;
        } else {
            while (fastaReader.hasNextChunk()) {
                fastaReader.readChunk(chunkSize); // Skip chunks of non-matching sequences
            }
        }
    }
    return refSequence;
}

unsigned long long mismatchCount = 0;
unsigned long long processedSeeds = 0;
void dpExpand(const ReadSeedHit &seedHit, const std::string &refSequence, const std::string &refSequenceName) {
    int readChunkId = std::get<0>(seedHit);
    int readPosition = std::get<1>(seedHit);
    int targetPosition = std::get<2>(seedHit);

    processedSeeds++;
    if (processedSeeds == ULLONG_MAX) {
        std::cout << "Close to overflow: Processed " << processedSeeds << " seeds, mismatch count: " << mismatchCount << std::endl;
        processedSeeds = 0;
    }

    // std::cout << "DP Expand - Read Chunk ID: " << readChunkId << ", Read Position: " << readPosition << ", Target Position: " << targetPosition << std::endl;
    // std::cout << "Reference Sequence Name: " << refSequenceName << std::endl;
    // std::cout << "Reference Sequence Length: " << refSequence.size() << std::endl;
}

void mapReads(const std::unordered_map<std::string, std::string> &argsMap) {

    // std::string reference = "ATGATGAAGTGTGTAGCTCGCGGCCGATCGACTGCACGTACGTAGCCGCGACGATCTAGCTATATATAGCTAGTCGATCGCGACTGCATGCATCGCCGTCGCTCCTCCCGAATAACTAGCTACAGATAGAGAGAGAGATCGACTAGCTACGATCGCACTGTTTGACCACTGCAG";
    // std::string read = "AGCTAGTCGATCGCGACTGCATGCATCGCCGTCGCTCCTCCCGAATAACTAGCTACAGATAGAGAGAGAGATCGACTAGCTACGATCGCACTGTTTGACCACTGCAG";
    // reference = "AGCTAGTCGATCGCGACTGCATGCATCGCCGTCGCTCCTCCCGAATAACTAGCTACAGATAGAGAGAGAGATCGACTAGCTACGATCGCACTGTTTGACCACTGCAG";
    // // read                                               CGAAGACTAGCTC
    // // cigar                                              MMMMDXMMMMMMMX
    // read = "CGAAGACTAGCTC";
    // ReadAlignment alignment = alignReadToReference(read, reference);
    // alignment.readName = "read1";
    // alignment.refName = "ref1";
    // std::cout << "Alignment Results:" << std::endl;
    // std::cout << "Read Name: " << alignment.readName << std::endl;
    // std::cout << "Reference Name: " << alignment.refName << std::endl;
    // std::cout << "Edit Distance: " << alignment.editDistance << std::endl;
    // std::cout << "Read Start: " << alignment.readStart << std::endl;
    // std::cout << "Reference Start: " << alignment.refStart << std::endl;
    // std::cout << "Alignment Length: " << alignment.alignmentLength << std::endl;
    // std::cout << "CIGAR String: " << alignment.cigarString << std::endl;
    // return;

    // Placeholder for mapping logic
    // You can access the arguments from argsMap, for example:
    std::string referencePath = argsMap.at("reference");
    std::string sequencesWhitelistString = argsMap.find("sequences_whitelist") != argsMap.end() ? argsMap.at("sequences_whitelist") : "";
    std::string indexPath = argsMap.at("index_path");

    int kmerSize = std::stoi(argsMap.at("kmer_size"));
    std::string readsPath = argsMap.at("reads");
    std::string samOutputPath = argsMap.at("sam_output_path");
    // Implement the mapping logic here using the provided arguments
    std::cout << "Mapping reads with the following parameters:" << std::endl;
    std::cout << "Reference Path: " << referencePath << std::endl;
    std::cout << "Index Path: " << indexPath << std::endl;
    std::cout << "K-mer Size: " << kmerSize << std::endl;
    std::cout << "Reads Path: " << readsPath << std::endl;
    std::cout << "SAM Output Path: " << samOutputPath << std::endl;

    std::vector<std::string> whitelistSequences;
    if (!sequencesWhitelistString.empty()) {
        whitelistSequences = splitString(sequencesWhitelistString, ';');
    }

    // ------------------ Collect Index Files ----------------- //
    std::vector<std::string> indexFiles = collectIndexFilePaths(indexPath);
    indexFiles = filterIndexFilesByWhitelist(indexFiles, whitelistSequences);
    std::cout << "Collected " << indexFiles.size() << " index files." << std::endl;
    for (const auto &file : indexFiles) {
        std::cout << "- Index File: " << file << std::endl;
    }


    // ------------------ Read FASTQ File in Chunks ----------------- //
    const int READS_CHUNK_SIZE = 10000; // This is because 10,000 x 101 (read length) = 10,000,000 bases. Similar to the chunk size used for the reference sequences.
    std::cout << "Starting to read FASTQ file: " << readsPath << std::endl;
    std::vector<std::vector<ReadRecord>> readChunks;
    {
        FASTQReader qReader(readsPath);
        while (qReader.hasNext()) {
            readChunks.push_back(qReader.readChunk(READS_CHUNK_SIZE));
        }
        std::cout << "Finished reading FASTQ file." << std::endl;
        std::cout << "Read " << qReader.getNumReads() << " reads divided into " << qReader.getNumChunks() << " chunks." << std::endl;
    }

    std::string refSequenceName;
    std::string refSequence;

    std::unordered_map<std::string, ReadBestMatches> readBestMatchesMap; // Map to store the best and second-best alignments for each read
    std::vector<ReadSeedHit> seedHits;
    int readsProcessed = 0;
    for (int i = 0; i < indexFiles.size(); ++i) {
        std::cout << "-- IndexChunk " << i << "/" << indexFiles.size()-1 << std::endl;
        SimpleKmerIndexer indexer(kmerSize, "", i, -1, -1, -1);
        indexer.readIndexFromFile(indexFiles[i]);
        refSequenceName = indexer.getSequenceName();
        refSequence = loadReferenceSequenceChunk(referencePath, refSequenceName, indexer.getChunkNumber(), 320);

        // if (refSequence.empty() || refSequenceName.compare(indexer.getSequenceName()) != 0) {
        //     refSequenceName = indexer.getSequenceName();
        //     std::cout << "Reference sequence name: " << refSequenceName << std::endl;
        //     refSequence = loadReferenceSequence(referencePath, refSequenceName, 320);
        //     std::cout << "Loaded reference sequence with length: " << refSequence.size() << std::endl;
        // }
        
        for (int j = 0; j < readChunks.size(); ++j) {
            std::cout << "Processing READ chunk " << j << "/" << readChunks.size()-1 << std::endl;
            // SEED STAGE: Generate seeds for each read in the chunk and find candidate positions in the reference using the index
            
            readsProcessed = 0;
            seedHits.clear();
            for (const auto &read : readChunks[j]) {
                const std::string &readName = std::get<0>(read);
                const std::string &readSequence = std::get<1>(read);

                // Seeding
                for (int k = 0; k <= readSequence.size() - kmerSize; ++k) {
                    std::string kmer = readSequence.substr(k, kmerSize);
                    for (int pos : indexer.kmerIndex[kmer]) {
                        seedHits.emplace_back(j, k, pos); // (readChunkId, readPosition, targetPosition)
                    }
                }
                
                // std::cout << "Aligning with DP read '" << readName << "' with sequence: " << readSequence << std::endl;
                auto alignment = dpSimpleSemiGlobalEditDistance(readSequence, refSequence);
                alignment.readName = readName;
                alignment.refName = refSequenceName;
                alignment.refStart += indexer.getChunkNumber() * 320; // Adjust the reference start position based on the chunk number and chunk size (320 in this case)

                std::cout << "Read #" << readsProcessed << ": " << readName << std::endl;
                std::cout << "Number of seed hits: " << seedHits.size() << std::endl;
                std::cout << "Best alignment edit distance: " << alignment.editDistance << std::endl;
                std::cout << "Best alignment read start: " << alignment.readStart << std::endl;
                std::cout << "Best alignment reference start: " << alignment.refStart << std::endl;
                std::cout << "Best alignment length: " << alignment.alignmentLength << std::endl;
                std::cout << "Best alignment CIGAR string: " << alignment.cigarString << std::endl;

                readBestMatchesMap[readName].updateBestAlignmentAgainst(alignment);

                readsProcessed++;
            }
            
            // Expand Non-GAP filter. We try to expand the seed from right and count number of mismatches until we reach a gap or exceed a certain threshold of mismatches. This can help to filter out seed hits that are unlikely to lead to good alignments.
            // int MISMATCH_THRESHOLD = 10; // This is an example threshold, you can adjust it based on your needs.
            // int MAX_EXTENSION_LENGTH = 45; // This is an example maximum extension length, you can adjust it based on your needs.
            // std::vector<ReadSeedHit> filteredSeedHits;
            // for (const auto &hit : seedHits) {
            //     int readChunkId = std::get<0>(hit);
            //     int readPosition = std::get<1>(hit);
            //     int targetPosition = std::get<2>(hit);

            //     const std::string &readSequence = std::get<1>(reads[readChunkId][0]); // Get the read sequence from the first read in the chunk (you can modify this to get the correct read sequence based on your data structure)
            //     // Here you would implement the logic to expand the seed hit and count mismatches until you reach a gap or exceed the mismatch threshold.
            //     // If the seed hit passes the filter, you can add it to the filteredSeedHits vector for further processing in the DP alignment stage.
            //     int mismatches = 0;
            //     int extensionLength = 0;
            //     while (extensionLength < MAX_EXTENSION_LENGTH && mismatches <= MISMATCH_THRESHOLD) {
            //     }
                    
            // }
            // for now, copy all seed hits to filteredSeedHits without filtering to test the mapping.
            // for (const auto &hit : seedHits) {
            //     filteredSeedHits.push_back(hit);
            // }




            // FILTER STAGE: Filter candidate positions based on the number of shared seeds or other heuristics
            // DP ALIGN STAGE: Perform dynamic programming alignment (e.g., Smith-Waterman) for the remaining candidate positions to find the best alignment
            // WRITE STAGE: Write the alignment results to the SAM output file
            // std::cout << "Found " << seedHits.size() << " seed hits for REF chunk 0, #reads:" << reads.size() << std::endl;
            // // compute avarage hits per read and max hits per read
            // float avgHitsPerRead = (float)(seedHits.size()) / reads.size();
            // std::cout << "Average hits per read: " << avgHitsPerRead << std::endl;
            // // print amount of seed hits found for this chunk
            // for (int j = 0; j < seedHits.size(); j++) {
            //     ReadAlignment alignment = alignReadToReference(std::get<1>(reads[std::get<0>(seedHits[j])][0]), refSequence);
            // }
        }

        // std::cout << "Total processed seeds: " << processedSeeds << ", Total mismatches: " << mismatchCount << std::endl;
    }

    std::cout << "\n\n================================================================" << std::endl;
    std::cout << "               FINAL BEST ALIGNMENTS PER READ                  " << std::endl;
    std::cout << "================================================================" << std::endl;

    for (auto readsChunk : readChunks) {
        for (auto read : readsChunk) {
            const std::string &readName = std::get<0>(read);
            const std::string &readSequence = std::get<1>(read);
            if (readBestMatchesMap.find(readName) == readBestMatchesMap.end()) {
                std::cout << "No alignments found for read: " << readName << std::endl;
                continue;
            }
            const ReadBestMatches &bestMatches = readBestMatchesMap[readName];
            std::cout << "Read: " << readName << std::endl;
            std::cout << "-1st Edit Distance: " << bestMatches.bestAlignment.editDistance << std::endl;
            std::cout << "-1st Reference Name: " << bestMatches.bestAlignment.refName << std::endl;
            std::cout << "-1st Reference Start: " << bestMatches.bestAlignment.refStart << std::endl;
            std::cout << "-1st Length: " << bestMatches.bestAlignment.alignmentLength << std::endl;
            std::cout << "-1st CIGAR String: " << bestMatches.bestAlignment.cigarString << std::endl;

            if (bestMatches.secondBestAlignment.editDistance != -1) {
                std::cout << "-2nd Edit Distance: " << bestMatches.secondBestAlignment.editDistance << std::endl;
                std::cout << "-2nd Reference Name: " << bestMatches.secondBestAlignment.refName << std::endl;
                std::cout << "-2nd Reference Start: " << bestMatches.secondBestAlignment.refStart << std::endl;
                std::cout << "-2nd Length: " << bestMatches.secondBestAlignment.alignmentLength << std::endl;
                std::cout << "-2nd CIGAR String: " << bestMatches.secondBestAlignment.cigarString << std::endl;
            } else {
                std::cout << "No second-best alignment found for read: " << readName << std::endl;
            }

            // Here you would implement the logic to write the best alignment for this read to the SAM output file.
        }
    }

    writeAlignmentsToSAM(readChunks, readBestMatchesMap, samOutputPath);


    
    // FOR indexFilePath in indexFiles:
    //     LOAD index file into memory (e.g., using a hash table)
    //     FOR each readChunk in reads:
    //         SEED STAGE: Generate seeds for each read in the chunk and find candidate positions in the reference using the index
    //         FILTER STAGE: Filter candidate positions based on the number of shared seeds or other heuristics
    //         DP ALIGN STAGE: Perform dynamic programming alignment (e.g., Smith-Waterman) for the remaining candidate positions to find the best alignment
    //         WRITE STAGE: Write the alignment results to the SAM output file

    // ----------------- 1. SEED STAGE ----------------- //
    // ----------------- 2. FILTER STAGE ----------------- //
    // ----------------- 3. DP ALIGN STAGE ----------------- //
    // ----------------- 4. WRITE STAGE ----------------- //
}