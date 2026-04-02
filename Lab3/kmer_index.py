


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


class MapIndex:
  def __init__(self, text, k, alphabet="ACGT"):
    self.k = k
    self.alphabet = alphabet
    self.index = {}
    for i in range(len(text) - k + 1):
      kmer = text[i:i+k]
      if kmer not in self.index:
        self.index[kmer] = []
      self.index[kmer].append(i)

  def get_sorted_kmers(self):
    return sorted(self.index.keys())

if __name__ == "__main__":
    fasta_file = "Data/chr3.fna"
    text = load_fasta(fasta_file)
    index = MapIndex(text, k=10)
    print(f"Loaded sequence length: {len(text)}")
    print(f"Total number of {index.k}-mers: {len(text) - index.k + 1}")
    print(f"Number of unique {index.k}-mers: {len(index.index)}")
    print(f"Most common {index.k}-mer: {index.get_sorted_kmers()[0]} (occurs {len(index.index[index.get_sorted_kmers()[0]])} times)")