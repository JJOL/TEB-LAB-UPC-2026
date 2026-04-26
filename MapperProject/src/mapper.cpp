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
#ifdef DEBUG
    std::cout << "Parent Directory: " << parentDirectory << std::endl;
#endif
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

#ifdef DEBUG
    std::cout << "Attempting to load reference sequence chunk: " << refSequenceName << ", chunk ID: " << chunkId << " from file: " << referencePath << std::endl;
#endif
    while (fastaReader.hasNextSequence()) {
        std::string sequenceName = fastaReader.readSequenceName();
        if (sequenceName.compare(refSequenceName) == 0) {
#ifdef DEBUG
            std::cout << "Found matching reference sequence: " << sequenceName << std::endl;
#endif
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
    int maxAllowedErrors = std::stoi(argsMap.at("max_errors"));
    int refChunkSize = std::stoi(argsMap.at("chunk_size"));
    std::string readsPath = argsMap.at("reads");
    std::string samOutputPath = argsMap.at("sam_output_path");
    // Implement the mapping logic here using the provided arguments
#ifdef DEBUG
    std::cout << "Mapping reads with the following parameters:" << std::endl;
    std::cout << "Reference Path: " << referencePath << std::endl;
    std::cout << "Index Path: " << indexPath << std::endl;
    std::cout << "Max Allowed Errors: " << maxAllowedErrors << std::endl;
    std::cout << "K-mer Size: " << kmerSize << std::endl;
    std::cout << "Reads Path: " << readsPath << std::endl;
    std::cout << "SAM Output Path: " << samOutputPath << std::endl;
#endif

    std::vector<std::string> whitelistSequences;
    if (!sequencesWhitelistString.empty()) {
        whitelistSequences = splitString(sequencesWhitelistString, ';');
    }

    // ------------------ Collect Index Files ----------------- //
    std::vector<std::string> indexFiles = collectIndexFilePaths(indexPath);
    indexFiles = filterIndexFilesByWhitelist(indexFiles, whitelistSequences);
#ifdef DEBUG
    std::cout << "Collected " << indexFiles.size() << " index files." << std::endl;
    for (const auto &file : indexFiles) {
        std::cout << "- Index File: " << file << std::endl;
    }
#endif

    // ------------------ Read FASTQ File in Chunks ----------------- //
    std::cout << "Reading reads from FASTQ file: " << readsPath << std::endl;
    const int READS_CHUNK_SIZE = 10000; // This is because 10,000 x 101 (read length) = 10,000,000 bases. Similar to the chunk size used for the reference sequences.
#ifdef DEBUG
    std::cout << "Starting to read FASTQ file: " << readsPath << std::endl;
#endif
    std::vector<std::vector<ReadRecord>> readChunks;
    {
        FASTQReader qReader(readsPath);
        while (qReader.hasNext()) {
            readChunks.push_back(qReader.readChunk(READS_CHUNK_SIZE));
        }
#ifdef DEBUG
        std::cout << "Finished reading FASTQ file." << std::endl;
        std::cout << "Read " << qReader.getNumReads() << " reads divided into " << qReader.getNumChunks() << " chunks." << std::endl;
#endif
    }

    std::string refSequenceName;
    std::string refSequence;

    std::unordered_map<std::string, ReadBestMatches> readBestMatchesMap; // Map to store the best and second-best alignments for each read
    
    int readsProcessed = 0;

    for (int i = 0; i < indexFiles.size(); ++i) {
    #ifdef DEBUG
        std::cout << "-- IndexChunk " << i << "/" << indexFiles.size()-1 << std::endl;
    #endif
        SimpleKmerIndexer indexer(kmerSize, "", i, -1, -1, -1);
        indexer.readIndexFromFile(indexFiles[i]);
        refSequenceName = indexer.getSequenceName();
        refSequence = loadReferenceSequenceChunk(referencePath, refSequenceName, indexer.getChunkNumber(), refChunkSize);

        // if (refSequence.empty() || refSequenceName.compare(indexer.getSequenceName()) != 0) {
        //     refSequenceName = indexer.getSequenceName();
        //     std::cout << "Reference sequence name: " << refSequenceName << std::endl;
        //     refSequence = loadReferenceSequence(referencePath, refSequenceName, 320);
        //     std::cout << "Loaded reference sequence with length: " << refSequence.size() << std::endl;
        // }
        std::unordered_map<int, std::tuple<std::string,std::string>> readSequencesMap; // Map to store readId->readSequence for all reads in the current chunk, to be used later when we have the seed hits and need to perform the DP alignment.
        for (int j = 0; j < readChunks.size(); ++j) {
            // std::cout << "Processing READ chunk " << j << "/" << readChunks.size()-1 << std::endl;            
            readsProcessed = 0;
            std::vector<ReadSeedHit> seedHits;
            for (const auto &read : readChunks[j]) {
                const std::string &readName = std::get<0>(read);
                const std::string &readSequence = std::get<1>(read);
                readSequencesMap[readsProcessed] = std::make_tuple(readName, readSequence); // Initialize with empty string for additional info if needed

                int qGramErrorThreshold = readSequence.size() - indexer.getKmerSize() + 1 - maxAllowedErrors * indexer.getKmerSize(); // This is the minimum number of shared k-mers required to consider a candidate position for alignment. It is derived from the fact that each error can affect at most k k-mers, so we need at least (readLength - k + 1) - (maxErrors * k) shared k-mers to have a chance of being within the max allowed errors.
                std::unordered_map<int, int> refStartvotes;

                // Seeding
                for (int k = 0; k <= readSequence.size() - kmerSize; ++k) {
                    std::string kmer = readSequence.substr(k, kmerSize);
                    for (int pos : indexer.kmerIndex[kmer]) {
                        int diagonal = pos - k;
                        refStartvotes[diagonal]++;
                        // seedHits.emplace_back(j, k, pos); // (readChunkId, readPosition, targetPosition)
                    }
                }
                // fast iterate over refStartVotes
                for (const auto &[diagonal, votes] : refStartvotes) {
                    int count = 0;
                    for (int offset = diagonal - maxAllowedErrors; offset <= diagonal + maxAllowedErrors; ++offset) {
                        if (refStartvotes.find(offset) != refStartvotes.end()) { // Weird way of counting the number of votes.
                            count += refStartvotes[offset];
                        }
                        // better than count += refStartVotes[offset]?
                        // idk because its more ifs, but perhaps is the same.
                    }
                    if (count >= qGramErrorThreshold) {
                        // seed.
                        seedHits.emplace_back(readsProcessed, i, diagonal); // (readId, chunkId, referenceStartPosition) 
                    }
                }
                
                readsProcessed++;
                // #ifdef DEBUG
                // std::cout << "Aligning with DP read '" << readName << "' with sequence: " << refSequenceName << std::endl;
                // #endif
                // auto alignment = dpSimpleSemiGlobalEditDistance(readSequence, refSequence);
                // #ifdef DEBUG
                // std::cout << "Done." << std::endl;
                // #endif
                // alignment.readName = readName;
                // alignment.refName = refSequenceName;
                // alignment.refStart += indexer.getChunkNumber() * refChunkSize; // Adjust the reference start position based on the chunk number and chunk size (320 in this case)

                // // std::cout << "Read #" << readsProcessed << ": " << readName << std::endl;
                // // std::cout << "Number of seed hits: " << seedHits.size() << std::endl;
                // // std::cout << "Best alignment edit distance: " << alignment.editDistance << std::endl;
                // // std::cout << "Best alignment read start: " << alignment.readStart << std::endl;
                // // std::cout << "Best alignment reference start: " << alignment.refStart << std::endl;
                // // std::cout << "Best alignment length: " << alignment.alignmentLength << std::endl;
                // // std::cout << "Best alignment CIGAR string: " << alignment.cigarString << std::endl;

                // if (alignment.editDistance != -1 && alignment.editDistance <= maxAllowedErrors) {
                //     readBestMatchesMap[readName].updateBestAlignmentAgainst(alignment);
                // }

            }

            // readId->sequence can be obtained from readSequences[readId] map.
            // chunkId->reference sequence is loaded in refSequence variable.
            // referencePosition has to be adjusted by pos-k to pos+len(read)+k to account for indels errors.
            float avgDPWindowSize = 0;
            for (const auto &seedHit : seedHits) {
                int readId, refChunkId, refStart;
                std::tie(readId, refChunkId, refStart) = seedHit;
                const std::string &readName = std::get<0>(readSequencesMap[readId]);
                const std::string &readSequence = std::get<1>(readSequencesMap[readId]);
                
                // build window around refStart - maxAllowedErrors to refStart + readSequence.size() + maxAllowedErrors
                int windowStart = std::max(0, refStart - maxAllowedErrors);
                int windowEnd = std::min(static_cast<int>(refSequence.size()), refStart + static_cast<int>(readSequence.size()) + maxAllowedErrors);
                std::string refWindow = refSequence.substr(windowStart, windowEnd - windowStart);
                avgDPWindowSize += refWindow.size();
                // Perform DP on a smaller window.
                auto alignment = dpSimpleSemiGlobalEditDistance(readSequence, refWindow);
                alignment.readName = readName;
                alignment.refName = refSequenceName;
                alignment.refStart += indexer.getChunkNumber() * refChunkSize + windowStart; // Adjust the reference start position based on the chunk number, chunk size and window start position.
                if (alignment.editDistance != -1 && alignment.editDistance <= maxAllowedErrors) {
                    readBestMatchesMap[readName].updateBestAlignmentAgainst(alignment);
                }
            }
#ifdef DEBUG
            // print how many seeds (aka DP executions) are compared to number of processedReads. I want to see the ratio per chunk.
            std::cout << "Processed " << readsProcessed << " reads in chunk " << j << "/" << readChunks.size()-1 << " with " << seedHits.size() << " seed hits (DP executions)." << std::endl;
            // print avg size DP window
            if (!seedHits.empty()) {
                avgDPWindowSize /= seedHits.size();
                std::cout << "Average DP window size: " << avgDPWindowSize << std::endl;
            }
#endif
        }
    }

#ifdef DEBUG
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
#endif

    std::cout << "Writing alignments to SAM output file: " << samOutputPath << std::endl;
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