#include "reader.hpp"
#include <fstream>

FASTAReader::FASTAReader(const std::string &filePath) {
    file.open(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }
    nSeqs = 1;
    basesReadInChunk = 0;
    totalBasesRead = 0;

    std::getline(file, currentLine);
}

bool FASTAReader::hasNextSequence() const {
    return !currentLine.empty() && currentLine[0] == '>';
}

std::string FASTAReader::readSequenceName() {
    currentChunk = 0; // Reset chunk counter for the new sequence
    std::string readSequenceName = currentLine.substr(1); // Remove the '>' character
    std::getline(file, currentLine); // Move to the first line of the sequence
    return readSequenceName;
}

bool FASTAReader::hasNextChunk() const {
    return !file.eof() && !currentLine.empty() && currentLine[0] != '>';
}

std::string FASTAReader::readChunk(size_t chunkSize) {
    // we read chunkSize or until we hit the next sequence header. If we hit the next sequence header, we stop and return what we have read so far.
    // currentLine ends at the next line or remaining part of the last line read. If currentLine is empty, we read the next line from the file.
    std::string chunk;
    while (chunk.size() < chunkSize && !file.eof()) {
        if (currentLine.empty()) {
            std::getline(file, currentLine);
            if (currentLine.empty()) {
                continue; // Skip empty lines
            }
        }
        if (currentLine[0] == '>') {
            break; // Stop at the next sequence header
        }
        size_t remaining = chunkSize - chunk.size();
        if (currentLine.size() <= remaining) {
            chunk += currentLine;
            currentLine.clear();
        } else {
            chunk += currentLine.substr(0, remaining);
            currentLine = currentLine.substr(remaining);
        }
    }
    basesReadInChunk = chunk.size();
    totalBasesRead += chunk.size();
    return chunk;
}
