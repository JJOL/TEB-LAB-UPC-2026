// Boyer–Moore preprocessing and search (Gusfield-based tables) — C++ translation
// Constraints: explicit loops, simple STL, no templates beyond necessity, 2-space indentation.

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

static vector<int> z_array(const string& s) {
  // Z algorithm (Gusfield theorem 1.4.1)
  int n = (int)s.size();
  if (n <= 1) {
    return vector<int>();
  }

  vector<int> z(n, 0);
  z[0] = n;

  // Initial comparison of s[1:] with prefix
  int z1 = 0;
  for (int i = 1; i < n; ++i) {
    if (s[i] == s[i - 1]) z1++;
    else break;
  }
  z[1] = z1;

  int r = 0;
  int l = 0;
  if (z[1] > 0) {
    r = z[1];
    l = 1;
  }

  for (int k = 2; k < n; ++k) {
    if (k > r) {
      // Case 1
      int zk = 0;
      for (int i = k; i < n; ++i) {
        if (s[i] == s[i - k]) zk++;
        else break;
      }
      z[k] = zk;
      r = k + z[k] - 1;
      l = k;
    } else {
      // Case 2
      int nbeta = r - k + 1;
      int zkp = z[k - l];

      if (nbeta > zkp) {
        // Case 2a
        z[k] = zkp;
      } else {
        // Case 2b
        int nmatch = 0;
        for (int i = r + 1; i < n; ++i) {
          if (s[i] == s[i - k]) nmatch++;
          else break;
        }
        l = k;
        r = r + nmatch;
        z[k] = r - k + 1;
      }
    }
  }

  return z;
}

static vector<int> n_array(const string& s) {
  // N array (Gusfield theorem 2.2.2) from Z array of reversed string
  string rev = s;
  reverse(rev.begin(), rev.end());

  vector<int> z = z_array(rev);
  reverse(z.begin(), z.end());
  return z;
}

static vector<int> big_l_prime_array(const string& p, const vector<int>& n) {
  // L' array (Gusfield theorem 2.2.2)
  int m = (int)p.size();
  vector<int> lp(m, 0);

  for (int j = 0; j < m - 1; ++j) {
    int i = m - n[j];
    if (i < m) {
      lp[i] = j + 1;
    }
  }
  return lp;
}

static vector<int> big_l_array(const string& p, const vector<int>& lp) {
  // L array (Gusfield theorem 2.2.2)
  int m = (int)p.size();
  vector<int> l(m, 0);
  if (m > 1) {
    l[1] = lp[1];
    for (int i = 2; i < m; ++i) {
      l[i] = max(l[i - 1], lp[i]);
    }
  }
  return l;
}

static vector<int> small_l_prime_array(const vector<int>& n) {
  // l' array (Gusfield theorem 2.2.4)
  int m = (int)n.size();
  vector<int> small_lp(m, 0);

  for (int i = 0; i < m; ++i) {
    if (n[i] == i + 1) {
      small_lp[m - i - 1] = i + 1;
    }
  }

  for (int i = m - 2; i >= 0; --i) {
    if (small_lp[i] == 0) small_lp[i] = small_lp[i + 1];
  }

  return small_lp;
}

static void good_suffix_table(const string& p,
                              vector<int>& big_l,
                              vector<int>& small_l_prime) {
  vector<int> n = n_array(p);
  vector<int> lp = big_l_prime_array(p, n);
  big_l = big_l_array(p, lp);
  small_l_prime = small_l_prime_array(n);
}

static vector< vector<int> > dense_bad_char_tab(const string& p,
                                                const vector<int>& amap) {
  // Dense bad character table indexed by pattern offset then by alphabet index.
  // Stores last occurrence position + 1 (0 means "not seen").
  int m = (int)p.size();
  int sigma = (int)amap.size();

  vector< vector<int> > tab;
  tab.reserve(m);

  vector<int> nxt(sigma, 0);
  for (int i = 0; i < m; ++i) {
    tab.push_back(nxt);
    unsigned char c = (unsigned char)p[i];
    int ci = amap[c];
    nxt[ci] = i + 1;
  }
  return tab;
}

class BoyerMoore {
public:
  string p;
  string alphabet;

  // Map ASCII char -> alphabet index, or -1 if not in alphabet
  vector<int> amap;
  vector< vector<int> > bad_char;
  vector<int> big_l;
  vector<int> small_l_prime;

  BoyerMoore(const string& pattern, const string& alpha = "ACGT") {
    p = pattern;
    alphabet = alpha;

    amap.assign(256, -1);
    for (int i = 0; i < (int)alphabet.size(); ++i) {
      unsigned char c = (unsigned char)alphabet[i];
      amap[c] = i;
    }

    bad_char = dense_bad_char_tab(p, amap);
    good_suffix_table(p, big_l, small_l_prime);
  }

  int bad_character_rule(int i, char c) const {
    unsigned char uc = (unsigned char)c;
    int ci = amap[uc];
    // If c is not in alphabet, treat as "not seen"
    if (ci < 0) return i + 1;

    int last_plus1 = bad_char[i][ci];  // 0 if not seen
    int last = last_plus1 - 1;         // -1 if not seen
    return i - last;
  }

  int good_suffix_rule(int i) const {
    int m = (int)big_l.size();
    if (i == m - 1) return 0;

    int ip1 = i + 1;
    if (big_l[ip1] > 0) return m - big_l[ip1];
    return m - small_l_prime[ip1];
  }

  int match_skip() const {
    int m = (int)small_l_prime.size();
    if (m <= 1) return 1;
    return m - small_l_prime[1];
  }
};

struct BoyerMooreResult {
  vector<int> matches;
  long long bad_rule_skips;
  long long goodsuffix_rule_skips;
};

static BoyerMooreResult boyer_moore_with_counts(const string& p, const BoyerMoore& bm, const string& t) {
  int n = (int)p.size();
  int m = (int)t.size();

  BoyerMooreResult result;
  result.matches = vector<int>();
  result.bad_rule_skips = 0;
  result.goodsuffix_rule_skips = 0;

  if (n == 0) return result;
  if (n > m) return result;

  int i = 0;
  while (i <= m - n) {
    int shift = 1;
    bool mismatched = false;

    for (int j = n - 1; j >= 0; --j) {
      if (p[j] != t[i + j]) {
        int skip_bc = bm.bad_character_rule(j, t[i + j]);
        int skip_gs = bm.good_suffix_rule(j);
        shift = max(shift, max(skip_bc, skip_gs));
        // Count which rule is applied
        if (skip_bc > skip_gs) {
          result.bad_rule_skips++;
        } else if (skip_gs > skip_bc) {
          result.goodsuffix_rule_skips++;
        }
        mismatched = true;
        break;
      }
    }

    if (!mismatched) {
      result.matches.push_back(i);
      int skip_gs = bm.match_skip();
      shift = max(shift, skip_gs);
    }

    i += shift;
  }

  return result;
}

static vector<int> boyer_moore(const string& p, const BoyerMoore& bm, const string& t) {
  int n = (int)p.size();
  int m = (int)t.size();

  vector<int> occ;
  if (n == 0) return occ;
  if (n > m) return occ;

  int i = 0;
  while (i <= m - n) {
    int shift = 1;
    bool mismatched = false;

    for (int j = n - 1; j >= 0; --j) {
      if (p[j] != t[i + j]) {
        int skip_bc = bm.bad_character_rule(j, t[i + j]);
        int skip_gs = bm.good_suffix_rule(j);
        shift = max(shift, max(skip_bc, skip_gs));
        mismatched = true;
        break;
      }
    }

    if (!mismatched) {
      occ.push_back(i);
      int skip_gs = bm.match_skip();
      shift = max(shift, skip_gs);
    }

    i += shift;
  }

  return occ;
}

static string load_fasta(const string& filename) {
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

static void print_vec(const vector<int>& v) {
  cout << "[";
  for (int i = 0; i < (int)v.size(); ++i) {
    if (i) cout << ", ";
    cout << v[i];
  }
  cout << "]";
}

int main() {
  // Example from the Python snippet
  //      012345678901234567
  //            TCTA    TCTA
  string t = "GCTAGCTCTACGAGTCTA";
  string p = "TCTA";

  BoyerMoore bm(p, "ACGT");
  vector<int> occ = boyer_moore(p, bm, t);

  print_vec(occ);
  cout << "\n";

  // Load and search in chr3.fna
  cout << "\nLoading chr3.fna...\n";
  string text = load_fasta("../chr3.fna");
  
  if (!text.empty()) {
    cout << "Loaded sequence length: " << text.size() << " bp\n";
    // Example search
    string pattern = "GATTACA";
    cout << "First segment of text: " << text.substr(0, 50) << " ...\n";
    BoyerMoore bm2(pattern, "ACGTN");
    BoyerMooreResult result = boyer_moore_with_counts(pattern, bm2, text);
    cout << "Found " << result.matches.size() << " matches for pattern '" << pattern << "' in chr3.fna\n";
    cout << "Bad character rule skips: " << result.bad_rule_skips << "\n";
    cout << "Good suffix rule skips: " << result.goodsuffix_rule_skips << "\n";
  }

  return 0;
}