#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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

// Print: rank, SA[rank], suffix text (tab-separated)
void print_suffix_array_view(const string& T, const vector<int>& sa) {
  cout << "rank\tSA\tsuffix\n";
  for (int r = 0; r < (int)sa.size(); ++r) {
    int i = sa[r];
    cout << r << '\t' << i << '\t';
    for (int j = i; j < (int)T.size(); ++j) cout << T[j];
    cout << '\n';
  }
}

int main(int argc, char* argv[]) {
  // Default: use a small test string, or read from a chromosome file
  string T;
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
      T = T.substr(0, 100000);
    }
  } else {
    // No argument: use default test string
    T = "ACGTACGTACGTACGT";
    cout << "Using test string: " << T << endl;
  }

  vector<int> sa = build_suffix_array(T);

  // Print raw suffix array
  cout << "SA: ";
  // Print all SA entries if small, otherwise just first and last 20
  if (sa.size() <= 40) {
    for (int i = 0; i < (int)sa.size(); ++i) {
      if (i) cout << ' ';
      cout << sa[i];
    }
  } else {
    for (int i = 0; i < 20; ++i) {
      if (i) cout << ' ';
      cout << sa[i];
    }
    cout << " ... ";
    for (int i = (int)sa.size() - 20; i < (int)sa.size(); ++i) {
      if (i != (int)sa.size() - 20) cout << ' ';
      cout << sa[i];
    }
  }
  cout << endl;

  // Pretty print (limit output for large arrays)
  if (sa.size() <= 100) {
    print_suffix_array_view(T, sa);
  } else {
    cout << "\nSuffix array is large (" << sa.size() 
         << " entries). Showing first 20 and last 20 entries:\n";
    cout << "rank\tSA\tsuffix\n";
    
    // First 20
    for (int r = 0; r < 20 && r < (int)sa.size(); ++r) {
      int i = sa[r];
      cout << r << '\t' << i << '\t';
      // Show first 50 characters of suffix
      for (int j = i; j < (int)T.size() && j < i + 50; ++j) 
        cout << T[j];
      if (i + 50 < (int)T.size()) cout << "...";
      cout << '\n';
    }
    
    cout << "...\n";
    
    // Last 20
    for (int r = max(20, (int)sa.size() - 20); r < (int)sa.size(); ++r) {
      int i = sa[r];
      cout << r << '\t' << i << '\t';
      for (int j = i; j < (int)T.size() && j < i + 50; ++j) 
        cout << T[j];
      if (i + 50 < (int)T.size()) cout << "...";
      cout << '\n';
    }
  }

  return 0;
}