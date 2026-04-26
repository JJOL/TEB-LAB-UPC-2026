#ifndef WRITER_HPP
#define WRITER_HPP

#include "types.hpp"
#include <string>
#include <vector>
#include <unordered_map>

void writeAlignmentsToSAM(
    const std::vector<std::vector<ReadRecord>> &readChunks,
    const std::unordered_map<std::string, ReadBestMatches> &readBestMatchesMap,
    const std::string &samOutputPath
);

#endif // WRITER_HPP