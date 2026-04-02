#include <zlib.h>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

static inline void strip_newlines(std::string& s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <file.fastq.gz> <PATTERN>\n";
        return 1;
    }

    const char* fname = argv[1];
    const std::string pattern = argv[2];

    gzFile gz = gzopen(fname, "rb");
    if (!gz) {
        std::cerr << "Error: cannot open " << fname << " (" << std::strerror(errno) << ")\n";
        return 1;
    }

    const int BUF_SIZE = 1 << 20; // 1MB per line
    std::string buf(BUF_SIZE, '\0');

    std::vector<std::string> sequences;
    sequences.reserve(1'000'000);

    std::uint64_t records = 0;
    std::size_t expected_len = 0;

    while (true) {
        // 1) header
        char* h = gzgets(gz, buf.data(), BUF_SIZE);
        if (!h) break; // EOF

        // If you want, you can validate header[0] == '@'
        // But we ignore header content as requested.

        // 2) sequence (COPY IT NOW, before reading more lines)
        char* s = gzgets(gz, buf.data(), BUF_SIZE);
        if (!s) {
            std::cerr << "Error: truncated FASTQ (missing sequence) near record " << records << "\n";
            gzclose(gz);
            return 1;
        }
        std::string seq(s);
        strip_newlines(seq);

        // 3) plus
        char* p = gzgets(gz, buf.data(), BUF_SIZE);
        // 4) quality
        char* q = gzgets(gz, buf.data(), BUF_SIZE);
        if (!p || !q) {
            std::cerr << "Error: truncated FASTQ (missing +/quality) near record " << records << "\n";
            gzclose(gz);
            return 1;
        }

        if (records == 0) {
            expected_len = seq.size();
            if (pattern.size() != expected_len) {
                std::cerr << "Error: pattern length " << pattern.size()
                          << " != read length " << expected_len << "\n";
                gzclose(gz);
                return 1;
            }
        } else if (seq.size() != expected_len) {
            std::cerr << "Error: read length mismatch at record " << records
                      << " (expected " << expected_len << ", got " << seq.size() << ")\n";
            gzclose(gz);
            return 1;
        }

        sequences.push_back(std::move(seq));
        ++records;
    }

    gzclose(gz);

    // Example: print the 4th sequence (index 3)
    if (sequences.size() > 3) {
        std::cerr << "Sequence[3] = " << sequences[3] << "\n";
    }

    // Brute-force exact matches
    std::uint64_t matches = 0;
    for (std::uint64_t i = 0; i < sequences.size(); ++i) {
        if (sequences[i] == pattern) ++matches;
    }

    std::cout << "Records read: " << sequences.size() << "\n";
    std::cout << "Read length: " << expected_len << "\n";
    std::cout << "Exact matches: " << matches << "\n";
    return 0;
}
