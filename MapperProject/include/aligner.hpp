#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <string>
#include "types.hpp"

ReadAlignment alignReadToReference(const std::string &readSequence, const std::string &refSequence);

ReadAlignment dpSimpleSemiGlobalEditDistance(const std::string &readSequence, const std::string &refSequence);

#endif // ALIGNER_HPP