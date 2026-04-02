#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

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

// Compare suffixes T[i..] and T[j..] lexicographically
bool suffix_less(const string& T, int i, int j) {
  int n = (int)T.size();
  while (i < n && j < n) {
    if (T[i] < T[j]) return true;
    if (T[i] > T[j]) return false;
    ++i;
    ++j;
  }
  return i == n && j != n;
}

// Naive suffix array construction
vector<int> build_suffix_array(const string& T) {
  int n = (int)T.size();
  vector<int> sa(n);
  for (int i = 0; i < n; ++i) sa[i] = i;

  sort(sa.begin(), sa.end(), [&](int a, int b) {
    return suffix_less(T, a, b);
  });

  return sa;
}

// Compare pattern P with suffix T[i..] lexicographically.
// Return -1 if P < suffix, 0 if P is a prefix of suffix, +1 if P > suffix.
int compare_pattern_to_suffix(const string& T, int i, const string& P) {
  int n = (int)T.size();
  int m = (int)P.size();

  int k = 0;
  while (k < m && i + k < n && P[k] == T[i + k]) ++k;

  if (k == m) return 0;                 // P fully matched, P is a prefix
  if (i + k == n) return 1;             // suffix ended, suffix < P, so P > suffix
  if (P[k] < T[i + k]) return -1;
  return 1;
}

// Binary search over suffix array to test if P occurs in T.
// comparisons counts how many SA midpoints were evaluated.
bool contains_pattern_sa(const string& T, const vector<int>& sa,
                         const string& P, int& comparisons) {
  comparisons = 0;
  int lo = 0;
  int hi = (int)sa.size() - 1;

  while (lo <= hi) {
    int mid = lo + (hi - lo) / 2;
    ++comparisons;

    int c = compare_pattern_to_suffix(T, sa[mid], P);
    if (c == 0) return true;
    if (c < 0) hi = mid - 1;
    else lo = mid + 1;
  }
  return false;
}

int main(int argc, char* argv[]) {
  string T;
  vector<string> patterns;
  
  if (argc > 1) {
    // Command-line argument: load chromosome from file
    string fasta_file = argv[1];
    cout << "Loading chromosome from: " << fasta_file << endl;
    T = load_fasta(fasta_file);
    
    if (T.empty()) {
      cerr << "Failed to load file or file is empty." << endl;
      return 1;
    }
    
    cout << "Loaded sequence length: " << T.size() << endl;
    
    // For large sequences, limit suffix array construction
    if (T.size() > 100000) {
      cout << "Warning: Sequence is very large (" << T.size() 
           << " bp). Truncating to first 100000 characters for demonstration." << endl;
      // T = T.substr(0, 100000);
    }
    
    // Extract some patterns from the loaded text for testing
    if (T.size() >= 100) {
      patterns.push_back(T.substr(10000, 10));    // 10-mer from position 10000
      patterns.push_back(T.substr(100000, 150));   // 15-mer from position 100000
      patterns.push_back(T.substr(500000, 20));   // 20-mer from position 500000
      patterns.push_back("GATTACA");            // Pattern unlikely to exist
      patterns.push_back(T.substr(5000000, 8));      // 8-mer from position 50
    } else {
      patterns = {"ACG", "CGT", "TGA", "GATTACA"};
    }
  } else {
    // No argument: use default test string
    T = "ACGTACGTGACG";
    patterns = {"ACG", "CGT", "TGA", "GATTACA", "ACGTG", "GAC"};
    cout << "Using test string: " << T << "\n\n";
  }

  cout << "Building suffix array..." << endl;
  vector<int> sa = build_suffix_array(T);
  cout << "Suffix array built (size: " << sa.size() << ")\n\n";

  cout << "Testing patterns using binary search on suffix array:\n";

  for (int idx = 0; idx < (int)patterns.size(); ++idx) {
    const string& P = patterns[idx];
    int comps = 0;
    bool ok = contains_pattern_sa(T, sa, P, comps);
    
    // Truncate pattern display if too long
    string display_P = P;
    if (P.size() > 50) {
      display_P = P.substr(0, 47) + "...";
    }
    
    cout << "P = " << display_P << "\toccurs = " << (ok ? "yes" : "no")
         << "\tSA comparisons = " << comps << "\n";
  }

  return 0;
}