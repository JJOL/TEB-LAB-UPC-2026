// Task 2.1 (Naive exact string matching) — Reference solution (C++)
// Constraints: simple constructs, explicit loops, 2-space indentation.

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

struct Result {
  vector<int> matches;
  long long comparisons;
};

Result naive_exact_search(const string& text, const string& pattern) {
  Result res;
  res.comparisons = 0;

  int n = (int)text.size();
  int m = (int)pattern.size();

  if (m == 0) return res;      // convention: empty pattern -> no matches
  if (m > n) return res;

  for (int i = 0; i <= n - m; ++i) {
    bool ok = true;
    for (int j = 0; j < m; ++j) {
      res.comparisons++;
      if (text[i + j] != pattern[j]) {
        ok = false;
        break;
      }
    }
    if (ok) res.matches.push_back(i);
  }

  return res;
}

void print_vector(const vector<int>& v) {
  cout << "[";
  for (int i = 0; i < (int)v.size(); ++i) {
    if (i) cout << ", ";
    cout << v[i];
  }
  cout << "]";
}

bool equal_vectors(const vector<int>& a, const vector<int>& b) {
  if ((int)a.size() != (int)b.size()) return false;
  for (int i = 0; i < (int)a.size(); ++i) {
    if (a[i] != b[i]) return false;
  }
  return true;
}

void run_tests() {
  struct Test {
    string text;
    string pattern;
    vector<int> expected;
  };

  vector<Test> tests;
  tests.push_back({"GATTACAGATTACA", "GATTACA", {0, 7}});
  tests.push_back({"AAAAA", "AAA", {0, 1, 2}});
  tests.push_back({"ACGTACGTACGT", "TAC", {3, 7}});
  tests.push_back({"ACGT", "ACGTA", {}});
  tests.push_back({"", "A", {}});
  tests.push_back({"ACGT", "", {}}); // empty pattern convention

  for (int t = 0; t < (int)tests.size(); ++t) {
    const Test& test = tests[t];
    Result r = naive_exact_search(test.text, test.pattern);

    cout << "text    : " << test.text << "\n";
    cout << "pattern : " << test.pattern << "\n";
    cout << "matches : ";
    print_vector(r.matches);
    cout << " expected: ";
    print_vector(test.expected);
    cout << "\n";
    cout << "comparisons: " << r.comparisons << "\n";
    cout << "OK: " << (equal_vectors(r.matches, test.expected) ? "true" : "false") << "\n";
    cout << "\n";
  }
}

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

int main() {
  // Run basic tests
  run_tests();
  
  // Load and search in chr3.fna
  cout << "Loading chr3.fna...\n";
  string text = load_fasta("../chr3.fna");
  
  if (!text.empty()) {
    cout << "Loaded sequence length: " << text.size() << " bp\n";
    // Example search
    string pattern = "GATTACA";
    Result r = naive_exact_search(text, pattern);
    cout << "\nSearching for '" << pattern << "' in chr3.fna\n";
    cout << "Found " << r.matches.size() << " matches\n";
    cout << "Number of comparisons: " << r.comparisons << "\n";
  }
  
  return 0;
}