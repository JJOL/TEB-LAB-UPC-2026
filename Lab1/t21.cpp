#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>


struct Counts {
    // Index by unsigned char value directly (0..255) for O(1) increments
    std::string seq_identifier; // mandatory
    std::array<std::uint64_t, 256> c{};

    void add_line(std::string_view s) {
        // FASTA lines typically have no spaces; we still ignore '\r' safely.
        const unsigned char* p = reinterpret_cast<const unsigned char*>(s.data());
        const unsigned char* e = p + s.size();
        for (; p != e; ++p) {
            unsigned char ch = *p;
            if (ch == '\r') continue; // handle Windows CRLF
            ++c[ch];
        }
    }

    std::uint64_t total_bases() const {
        // If you only care about ACGTN/acgtn, sum those; otherwise sum all non-newline.
        return c['A'] + c['C'] + c['G'] + c['T'] + c['N'] +
               c['a'] + c['c'] + c['g'] + c['t'] + c['n'];
    }
};

static void print_counts(const Counts& cnt) {
    std::cout << "Sequence: " << cnt.seq_identifier << "\n";

    auto p = [&](char ch, const char* label) {
        std::cout << std::setw(6) << label << ": " << cnt.c[(unsigned char)ch] << "\n";
    };

    std::cout << "Uppercase:\n";
    p('A', "A"); p('C', "C"); p('G', "G"); p('T', "T"); p('N', "N");

    // std::cout << "\nLowercase:\n";
    // p('a', "a"); p('c', "c"); p('g', "g"); p('t', "t"); p('n', "n");

    // Optional: show non-ACGTN letters you might encounter
    std::uint64_t others = 0;
    for (int i = 0; i < 256; ++i) {
        unsigned char ch = static_cast<unsigned char>(i);
        if (ch == '\n' || ch == '\r') continue;
        if (ch=='A'||ch=='C'||ch=='G'||ch=='T'||ch=='N'||
            ch=='a'||ch=='c'||ch=='g'||ch=='t'||ch=='n') continue;
        others += cnt.c[ch];
    }
    std::cout << "Other (non ACGTN/acgtn, excluding newlines): " << others << "\n";

    std::cout << "Total ACGTN/acgtn bases counted: " << cnt.total_bases() << "\n\n";
}


void fna_read_impl_0(std::vector<Counts>& cnts, std::ifstream& in, const char* fname) {
    std::string line;
    line.reserve(128); // reserve a bit more than 80 to avoid reallocations

    std::uint64_t total_lines = 0;
    std::uint64_t header_lines = 0;
    std::uint64_t seq_lines = 0;

    auto t0 = std::chrono::steady_clock::now();

    int lastCountIndex = -1;
    while (std::getline(in, line)) {
        ++total_lines;
        if (!line.empty() && line[0] == '>') {
            lastCountIndex++;
            cnts.emplace_back(); // start a new count for this sequence
            // identifier is first work after >, so until first space or end of line
            line.resize(line.find_first_of(" \t\r\n", 1)); // trim to identifier only
            cnts.back().seq_identifier = line.substr(1); // store identifier without '>'
            ++header_lines;
            continue; // skip header
        }
        ++seq_lines;
        cnts[lastCountIndex].add_line(line);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dt = t1 - t0;

    std::cout << "File: " << fname << "\n";
    std::cout << "Total lines: " << total_lines
                << " | headers: " << header_lines
                << " | seq lines: " << seq_lines << "\n";
    std::cout << "Elapsed: " << std::fixed << std::setprecision(6) << dt.count() << " s\n\n";
}

void fna_read_impl_1(std::vector<Counts>& cnts, std::ifstream& in, const char* fname) {
    // We will now make use of the fact that FASTA files have lines of at most 80 characters. So we can read fixed-size blocks of 80 characters and process them directly, which may be more efficient than std::getline.

    std::string line;
    line.reserve(80); // reserve exactly 80 to avoid reallocations

    std::uint64_t total_lines = 0;
    std::uint64_t header_lines = 0;
    std::uint64_t seq_lines = 0;
    auto t0 = std::chrono::steady_clock::now();
    int lastCountIndex = -1;
    // to read exactly 80 characters, we can use in.read() into a char buffer of size 81 (80 + null terminator)
    char buffer[82]; // 80 chars + possible newline + null terminator
    while (in.read(buffer, 81) || in.gcount() > 0) {
        std::streamsize bytesRead = in.gcount();
        buffer[81] = '\0'; // reset buffer for next read
        buffer[80] = '\0';
        total_lines++;
        // buffer may not be null-terminated if we read exactly 80 chars, so we will construct a string_view from it
        std::string_view line(buffer, bytesRead);
        // std::cout << "Read line: '" << line << "'\n"; // debug print
        // if (total_lines > 10) break; // limit to first 10 lines for debugging
        if (!line.empty() && line[0] == '>') {
            lastCountIndex++;
            cnts.emplace_back(); // start a new count for this sequence
            // identifier is first work after >, so until first space or end of line
            size_t id_end = line.find_first_of(" \t\r\n", 1);
            cnts.back().seq_identifier = std::string(line.substr(1, id_end - 1)); // store identifier without '>'
            header_lines++;
        } else {
            seq_lines++;
            cnts[lastCountIndex].add_line(line);
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dt = t1 - t0;

    std::cout << "File: " << fname << "\n";
    std::cout << "Total lines: " << total_lines
                << " | headers: " << header_lines
                << " | seq lines: " << seq_lines << "\n";
    std::cout << "Elapsed: " << std::fixed << std::setprecision(6) << dt.count() << " s\n\n";
}


// what if we mmap the file to a whole block of ram and then process it as a single string_view?
// lets give it a try
void fna_read_impl_2(std::vector<Counts>& cnts, std::ifstream& in, const char* fname) {
    // mmap()
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    int fd = open(fname, O_RDONLY);
    // we need to include <fcntl.h> for O_RDONLY and <unistd.h> for close()
    if (fd == -1) {
        std::cerr << "Error: cannot open " << fname << " for mmap\n";
        return;
    }

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dt_open = t1 - t0;
    std::cout << "Time to open file for mmap: " << std::fixed << std::setprecision(6) << dt_open.count() << " s\n";
    
    close(fd); // we will reopen with ifstream, so close the fd immediately
}



int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fasta_count <file.fna | file.fna.gz>\n";
        return 1;
    }

    const char* fname = argv[1];

    auto fname_sv = std::string_view(fname);
    if (!fname_sv.ends_with(".fna") && !fname_sv.ends_with(".fna.gz")) {
        std::cerr << "Warning: file does not have .fna or .fna.gz extension, but will attempt to process anyway.\n";
    }

    std::ifstream in(fname, std::ios::in);
    if (!in) {
        std::cerr << "Error: cannot open " << fname << "\n";
        return 1;
    }

    if (fname_sv.ends_with(".fna")) {
        std::vector<Counts> cnts;

        fna_read_impl_1(cnts, in, fname);
    
        for (const auto& cnt : cnts) {
            print_counts(cnt);
        }
    }

    return 0;
}
