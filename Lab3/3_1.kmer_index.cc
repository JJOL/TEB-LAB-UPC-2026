#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

/**
 * Load a FASTA file and return the concatenated sequence.
 */
string load_fasta(const string& filename) {
    string sequence = "";
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: File '" << filename << "' not found." << endl;
        return "";
    }
    
    string line;
    while (getline(file, line)) {
        // Skip header lines (starting with '>')
        if (!line.empty() && line[0] != '>') {
            sequence += line;
        }
    }
    
    file.close();
    return sequence;
}

/**
 * MapIndex class: builds and manages a k-mer index using a hashmap
 */
class MapIndex {
public:
    int k;
    string alphabet;
    unordered_map<string, vector<int>> index;
    
    MapIndex(const string& text, int k, const string& alphabet = "ACGT") 
        : k(k), alphabet(alphabet) {
        build_index(text);
    }
    
    /**
     * Build the k-mer index from the text
     */
    void build_index(const string& text) {
        for (size_t i = 0; i <= text.length() - k; i++) {
            string kmer = text.substr(i, k);
            index[kmer].push_back(i);
        }
    }
    
    /**
     * Get all k-mers sorted alphabetically
     */
    vector<string> get_sorted_kmers() const {
        vector<string> kmers;
        for (const auto& pair : index) {
            kmers.push_back(pair.first);
        }
        sort(kmers.begin(), kmers.end());
        return kmers;
    }
    
    /**
     * Find all positions of a k-mer in the index
     */
    vector<int> find_kmer(const string& kmer) const {
        auto it = index.find(kmer);
        if (it != index.end()) {
            return it->second;
        }
        return vector<int>();
    }
    
    /**
     * Get statistics about the index
     */
    void print_stats(size_t text_length) const {
        cout << "Loaded sequence length: " << text_length << endl;
        cout << "Total number of " << k << "-mers: " << (text_length - k + 1) << endl;
        cout << "Number of unique " << k << "-mers: " << index.size() << endl;
        
        if (!index.empty()) {
            auto sorted = get_sorted_kmers();
            string most_common = sorted[0];
            cout << "Most common " << k << "-mer: " << most_common 
                 << " (occurs " << index.at(most_common).size() << " times)" << endl;
        }
    }
};

void searchPattern(const string& pattern, const MapIndex& index) {
    cout << "\nSearching for pattern: " << pattern << endl;
    
    // 1. Divide pattern into distinict k-mers
    vector<string> pattern_kmers;
    for (size_t i = 0; i <= pattern.length(); i += index.k) {
        int substring_end = min(i + index.k, pattern.length()) ;
        pattern_kmers.push_back(pattern.substr(i, substring_end - i));
    }
    cout << "Pattern divided into " << pattern_kmers.size() << " k-mers." << endl;
    
    // 2. List all segments of the pattern
    for (const auto& kmer : pattern_kmers) {
        cout << "Segment: " << kmer << endl;
    }

    // 3. Get the k-mer from the pattern with the fewest matches in the index
    string best_kmer;
    size_t min_matches = SIZE_MAX;
    for (const auto& kmer : pattern_kmers) {
        size_t matches = index.find_kmer(kmer).size();
        cout << "K-mer '" << kmer << "' has " << matches << " matches in the index." << endl;
        if (matches < min_matches) {
            min_matches = matches;
            best_kmer = kmer;
        }
    }
    cout << "Best k-mer: " << best_kmer << " with " << min_matches << " matches." << endl;
}

int main() {
    string fasta_file = "Data/chr3.fna";
    
    // Load FASTA file
    auto start = chrono::high_resolution_clock::now();
    string text = load_fasta(fasta_file);
    auto end = chrono::high_resolution_clock::now();
    
    if (text.empty()) {
        return 1;
    }
    
    auto load_time = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "FASTA load time: " << load_time.count() << " ms" << endl;
    
    // Build index
    int k = 10;
    start = chrono::high_resolution_clock::now();
    MapIndex index(text, k);
    end = chrono::high_resolution_clock::now();
    
    auto index_time = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Index build time: " << index_time.count() << " ms" << endl;
    
    // Print statistics
    cout << "\n--- Index Statistics ---" << endl;
    index.print_stats(text.length());
    
    // Example: Query some k-mers
    cout << "\n--- K-mer Queries ---" << endl;
    string test_kmer = text.substr(1000, k);
    vector<int> positions = index.find_kmer(test_kmer);
    cout << "K-mer '" << test_kmer << "' found at " << positions.size() << " positions" << endl;
    if (positions.size() <= 10) {
        cout << "Positions: ";
        for (int pos : positions) {
            cout << pos << " ";
        }
        cout << endl;
    } else {
        cout << "First 10 positions: ";
        for (int i = 0; i < 10; i++) {
            cout << positions[i] << " ";
        }
        cout << "..." << endl;
    }
    
    // // Example: Test a k-mer not in sequence
    // cout << "\nSearching for k-mer 'ACCTCTTACG':" << endl;
    // positions = index.find_kmer("ACCTCTTACG");
    // cout << "Found at " << positions.size() << " positions" << endl;
    // // List an example:
    // if (!positions.empty()) {
    //     cout << "Example position: " << positions[0] << endl;
    // }
    string pattern = text.substr(2476511, 30); // Extract a 50bp pattern from the text
    searchPattern(pattern, index);
    
    return 0;
}
