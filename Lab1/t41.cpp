#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

static inline int sym_index(unsigned char ch) {
    // Case-sensitive alphabet of 10 symbols:
    // 0..4  = A C G T N
    // 5..9  = a c g t n
    switch (ch) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'N': return 4;
        case 'a': return 5;
        case 'c': return 6;
        case 'g': return 7;
        case 't': return 8;
        case 'n': return 9;
        default:  return -1; // separator / invalid symbol
    }
}

static inline char sym_char(int idx) {
    static constexpr char map[10] = {'A','C','G','T','N','a','c','g','t','n'};
    return map[idx];
}

static bool safe_pow_u64(std::uint64_t base, std::uint64_t exp, std::uint64_t& out) {
    out = 1;
    for (std::uint64_t i = 0; i < exp; ++i) {
        if (out > std::numeric_limits<std::uint64_t>::max() / base) return false;
        out *= base;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: kmer_hist <file.fna> <K>\n";
        return 1;
    }

    const char* fname = argv[1];
    int K = 0;
    try {
        K = std::stoi(argv[2]);
    } catch (...) {
        std::cerr << "Error: K must be an integer.\n";
        return 1;
    }
    if (K <= 0) {
        std::cerr << "Error: K must be >= 1.\n";
        return 1;
    }

    std::ifstream in(fname, std::ios::in);
    if (!in) {
        std::cerr << "Error: cannot open " << fname << "\n";
        return 1;
    }

    constexpr std::uint64_t B = 10; // alphabet size
    std::uint64_t hist_size = 0;
    if (!safe_pow_u64(B, (std::uint64_t)K, hist_size)) {
        std::cerr << "Error: B^K overflows uint64_t for K=" << K << "\n";
        return 1;
    }

    // Practical memory guard (adjust if you want):
    // hist_size * 8 bytes. 200 million -> ~1.6GB
    const std::uint64_t MAX_BINS = 200'000'000ULL;
    if (hist_size > MAX_BINS) {
        std::cerr << "Error: histogram too large (10^" << K << " bins). "
                  << "Choose smaller K.\n";
        return 1;
    }

    std::vector<std::uint64_t> hist((size_t)hist_size, 0);

    // Precompute B^(K-1) for rolling removal
    std::uint64_t BK_1 = 1;
    for (int i = 1; i < K; ++i) BK_1 *= B;

    // Rolling state
    std::vector<int> ring((size_t)K, 0);
    int ring_pos = 0;
    int filled = 0;
    std::uint64_t code = 0;

    std::string line;
    line.reserve(4096);

    auto t0 = std::chrono::steady_clock::now();

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue; // header

        // Process sequence line
        const unsigned char* p = (const unsigned char*)line.data();
        const unsigned char* e = p + line.size();

        for (; p != e; ++p) {
            unsigned char ch = *p;
            if (ch == '\r') continue; // in case of CRLF

            int v = sym_index(ch);
            if (v < 0) {
                // break window on unknown symbol
                filled = 0;
                code = 0;
                ring_pos = 0;
                continue;
            }

            if (filled < K) {
                // grow window
                code = code * B + (std::uint64_t)v;
                ring[ring_pos] = v;
                ring_pos = (ring_pos + 1) % K;
                ++filled;

                if (filled == K) {
                    ++hist[(size_t)code];
                }
            } else {
                // slide window: remove oldest, add newest
                int oldest = ring[ring_pos];
                ring[ring_pos] = v;
                ring_pos = (ring_pos + 1) % K;

                code -= (std::uint64_t)oldest * BK_1; // remove leftmost
                code = code * B + (std::uint64_t)v;   // shift + append
                ++hist[(size_t)code];
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> dt = t1 - t0;

    // Output CSV to stdout
    std::cout << "kmer,count\n";
    std::string kmer((size_t)K, 'A');

    for (std::uint64_t idx = 0; idx < hist_size; ++idx) {
        std::uint64_t c = hist[(size_t)idx];
        if (c == 0) continue;

        // Decode idx -> base-B digits (length K)
        std::uint64_t x = idx;
        for (int pos = K - 1; pos >= 0; --pos) {
            int digit = (int)(x % B);
            x /= B;
            kmer[(size_t)pos] = sym_char(digit);
        }
        std::cout << kmer << "," << c << "\n";
    }

    // Timing to stderr so CSV stays clean on stdout
    std::cerr << "Elapsed: " << dt.count() << " s\n";
    return 0;
}
