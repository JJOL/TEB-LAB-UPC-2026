#include "reader.hpp"
#include <fstream>

FASTAReader::FASTAReader(const std::string &filePath) {
    file.open(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }
    nSeqs = 0;
    basesReadInChunk = 0;
    totalBasesRead = 0;

    std::getline(file, currentLine);
    currentLineNumber = 0;
    currentLineOffset = 0;
}

bool FASTAReader::hasNextSequence() const {
    return !currentLine.empty() && currentLine[0] == '>';
}

std::string FASTAReader::readSequenceName() {
    currentChunk = 0; // Reset chunk counter for the new sequence
    std::string readSequenceName = currentLine.substr(1); // Remove the '>' character

    std::getline(file, currentLine); // Move to the first line of the sequence
    currentLineNumber++; // Increment line number for the header line
    currentLineOffset = 0; // Reset line offset for the new sequence
    return readSequenceName;
}

bool FASTAReader::hasNextChunk() const {
    return !file.eof() && !currentLine.empty() && currentLine[0] != '>';
}

std::string FASTAReader::readChunk(size_t chunkSize) {
    // we read chunkSize or until we hit the next sequence header. If we hit the next sequence header, we stop and return what we have read so far.
    // currentLine ends at the next line or remaining part of the last line read. If currentLine is empty, we read the next line from the file.
    std::string chunk;

    while (chunk.size() < chunkSize && hasNextChunk()) {
        if (chunk.size() + currentLine.size() <= chunkSize) {
            chunk += currentLine;
            std::getline(file, currentLine); // Move to the next line for the next read
            currentLineNumber++; // Increment line number for the new line
            currentLineOffset = 0; // Reset line offset for the new line
        } else {
            size_t readSize = chunk.size();
            chunk += currentLine.substr(0, chunkSize - readSize);
            currentLine = currentLine.substr(chunkSize - readSize);
            currentLineOffset += chunkSize - readSize; // Update line offset for the next read
        }
    }

    return chunk;
}

FASTAReader::~FASTAReader() {
    if (file.is_open()) {
        file.close();
    }
}

FASTQReader::FASTQReader(const std::string &filePath) {
    file.open(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }
    nReads = 0;
    nChunks = 0;

    std::getline(file, currentLine);
}

bool FASTQReader::hasNext() const {
    return !currentLine.empty();
}

std::vector<ReadRecord> FASTQReader::readChunk(size_t chunkSize) {
    std::vector<ReadRecord> chunk;
    while (chunk.size() < chunkSize && hasNext()) {
        std::string readName = currentLine;
        std::string readSequence;
        std::string plusLine;
        std::string qualityScores;

        std::getline(file, readSequence);
        std::getline(file, plusLine);
        std::getline(file, qualityScores);

        chunk.emplace_back(readName, readSequence, qualityScores);

        nReads++;
        std::getline(file, currentLine); // Move to the next read
    }
    nChunks++;
    return chunk;
}

int FASTQReader::getNumReads() const {
    return nReads;
}

int FASTQReader::getNumChunks() const {
    return nChunks;
}

FASTQReader::~FASTQReader() {
    if (file.is_open()) {
        file.close();
    }
}


std::vector<ChromosomeInfo> findChromosomesInFASTA(const std::string &filePath) {
    std::vector<ChromosomeInfo> chromosomes;
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    // Chromosome start at the start of the line with > and until end of the line is the chromosome name. Then a bunch of lines follow with fixed size lines until the next chromosome header. 
    // we will read the file by scanning lines. We will count how many lines are in each chromosome and that will be the chromosome length. We will also keep track of the offset of the first base of each chromosome, which is the offset of the line after the chromosome header. We will return a vector of tuples with the chromosome name, length and offset.
    std::string line;
    std::string currentChromosomeName;
    size_t currentChromosomeOffset = 0;
    const size_t FIXED_LINE_SIZE = 60; // Assuming fixed line size of 60 bases per line

    while (std::getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            if (!currentChromosomeName.empty()) {
                chromosomes[chromosomes.size() - 1].length = currentChromosomeOffset - chromosomes[chromosomes.size() - 1].offset;
            }
            currentChromosomeName = line.substr(1); // Remove the '>' character
            chromosomes.emplace_back(currentChromosomeName, 0, currentChromosomeOffset); // Length will be updated later when we know how many bases are in the chromosome
        }
        currentChromosomeOffset += 1;
    }
    // After the loop, we need to save the last chromosome if it exists
    if (!currentChromosomeName.empty()) {
        chromosomes[chromosomes.size() - 1].length = currentChromosomeOffset - chromosomes[chromosomes.size() - 1].offset;
    }

    file.close();
    return chromosomes;
}