# Task 2.1 (Naive exact string matching) — Reference solution (Python)
# Constraints: simple constructs, explicit loops, 2-space indentation.

def naive_exact_search(text, pattern):
  n = len(text)
  m = len(pattern)
  matches = []
  comparisons = 0

  if m == 0:
    return matches, comparisons
  if m > n:
    return matches, comparisons

  i = 0
  while i <= n - m:
    j = 0
    ok = True
    while j < m:
      comparisons += 1
      if text[i + j] != pattern[j]:
        ok = False
        break
      j += 1
    if ok:
      matches.append(i)
    i += 1

  return matches, comparisons


def load_fasta(filename):
  """Load a FASTA file and return the concatenated sequence."""
  sequence = []
  try:
    with open(filename, 'r') as f:
      for line in f:
        line = line.strip()
        # Skip header lines (starting with '>')
        if line and not line.startswith('>'):
          sequence.append(line)
  except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    return ""
  return ''.join(sequence)


def run_tests():
  # Basic tests, includes overlap cases
  tests = [
    ("GATTACAGATTACA", "GATTACA", [0, 7]),
    ("AAAAA", "AAA", [0, 1, 2]),
    ("ACGTACGTACGT", "TAC", [3, 7]),
    ("ACGT", "ACGTA", []),
    ("", "A", []),
    ("ACGT", "", []),  # by convention here, empty pattern returns no matches
  ]

  for text, pattern, expected in tests:
    matches, comps = naive_exact_search(text, pattern)
    print("text    :", text)
    print("pattern :", pattern)
    print("matches :", matches, "expected:", expected)
    print("comparisons:", comps)
    print("OK:", matches == expected)
    print()


if __name__ == "__main__":
  # Load chr3.fna and run basic tests
  fasta_file = "Data/chr3.fna"
  
  # Run basic tests
  run_tests()
  
  # Load and search in chr3.fna
#   print("\nLoading chr3.fna...")
#   text = load_fasta(fasta_file)
  
#   if text:
#     print(f"Loaded sequence length: {len(text)} bp")
#     # Example search
#     pattern = """CAGTCTGCATTGAGGGGCCAGGCTGAAGTCTGTGAGTGAAACCCCAATGTGGGAGAGAAGCG
# GAGAGTGTCGTTCAATGATCCACCTTTCCCTTGGGAAACTCTCCAACCAAGCCACTGGTGAGCACCCTGTCCCTCCCAAG
# CCCTGGACCTAACCTGGGGAGAGGCTGGGAGACTGTAAGAAGGAATGTCACTGGGAAGTGCTTCAAGCACTGATCCAGAC
# TTGGGACCCAGTTGGAGGACGCCATTCTCAATTCTAGCTCATAACAAGCTGGGAGGGTCCCTGCAAGCTAGCAGCTACAA
# CAGCCAAGGGGATTACAGAGTTTCAAGATGGAATTTAGAGTGCTGGGTCTGCTGGGGAAGAGGGCCCACAGCCAGAACTG
# ACAGGCAAGTGTGGTGTGTTTTCTAGCCATAAGCGCTGGAATTGGACTCTCCTGATGTAG"""
#     # pattern = "GATTACA"  # IGNORE
#     matches, comps = naive_exact_search(text, pattern)
#     print(f"\nSearching for '{pattern}' in chr3.fna")
#     print(f"Found {len(matches)} matches")
#     print(f"Number of comparisons: {comps}")