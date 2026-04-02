"""
This is a master script. The objective is, given a text T, perform a mapping query on pattern P with allowance of k-edits.
The output is the list of positions in T where P matches with at most k edits, and the corresponding edit distance.

Implementation based on the Pigeonhole Principle:
- Partition pattern P into k+1 disjoint factors
- At least one factor must match exactly (pigeonhole principle)
- Use seed-and-verify strategy: search factors exactly, then verify candidates
"""
import random


# ============================================================================
# Hamming Distance
# ============================================================================

def hamming_distance(s1, s2):
    """
    Compute the Hamming distance between two strings of equal length.
    
    Args:
        s1: First string
        s2: Second string
    
    Returns:
        Number of positions where characters differ
    """
    if len(s1) != len(s2):
        raise ValueError("Strings must be of the same length")
    
    distance = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            distance += 1
    return distance


# ============================================================================
# Task 3.1: Partitioning the pattern
# ============================================================================

def partition_pattern(P, k):
    """
    Partition pattern P into k+1 consecutive factors of approximately equal length.
    
    Based on the pigeonhole principle: if P aligns with at most k mismatches,
    and P is divided into k+1 factors, then at least one factor must match exactly.
    
    Args:
        P: Pattern string
        k: Maximum number of mismatches allowed
    
    Returns:
        List of tuples (factor, start_position) where:
        - factor: substring of P
        - start_position: starting position of factor in P
    
    Example:
        >>> partition_pattern("ACGTACGT", 1)
        [('ACGT', 0), ('ACGT', 4)]  # k+1 = 2 factors
    """
    n = len(P)
    num_factors = k + 1
    
    if num_factors > n:
        # If we need more factors than characters, make single-character factors
        return [(P[i], i) for i in range(n)]
    
    # Calculate base length and number of factors that need an extra character
    base_length = n // num_factors
    extra_chars = n % num_factors
    
    factors = []
    pos = 0
    
    for i in range(num_factors):
        # First 'extra_chars' factors get one extra character
        length = base_length + (1 if i < extra_chars else 0)
        factor = P[pos:pos + length]
        factors.append((factor, pos))
        pos += length
    
    return factors


# ============================================================================
# Task 3.2: Exact factor search
# ============================================================================

def build_substring_index(T):
    """
    Build a hash table index mapping all substrings to their positions in text T.
    
    Args:
        T: Text string
    
    Returns:
        Dictionary mapping substrings to list of positions
    """
    index = {}
    n = len(T)
    
    # Index all substrings of various lengths
    # We'll index substrings up to a reasonable length to avoid explosion
    max_length = min(n, 50)  # Limit to avoid too many substrings
    
    for length in range(1, max_length + 1):
        for i in range(n - length + 1):
            substring = T[i:i + length]
            if substring not in index:
                index[substring] = []
            index[substring].append(i)
    
    return index


def search_factors_exact(factors, T):
    """
    Search for exact occurrences of each factor in text T.
    For each occurrence, compute candidate alignment positions for the full pattern.
    
    Args:
        factors: List of tuples (factor, start_pos_in_pattern)
        T: Text string
    
    Returns:
        Set of candidate positions where the full pattern might align
    """
    candidates = set()
    
    for factor, factor_start_in_pattern in factors:
        # Find all occurrences of this factor in T
        factor_len = len(factor)
        
        for i in range(len(T) - factor_len + 1):
            if T[i:i + factor_len] == factor:
                # If factor matches at position i in T,
                # and factor starts at position factor_start_in_pattern in P,
                # then P might start at position i - factor_start_in_pattern in T
                pattern_start = i - factor_start_in_pattern
                candidates.add(pattern_start)
    
    return candidates


# ============================================================================
# Task 3.3: Verification of candidate alignments
# ============================================================================

def verify_candidates(candidates, P, T, k):
    """
    Verify each candidate alignment using Hamming distance.
    
    Args:
        candidates: Set of candidate starting positions
        P: Pattern string
        T: Text string
        k: Maximum number of mismatches allowed
    
    Returns:
        List of tuples (position, distance) where distance <= k, sorted by position
    """
    matches = []
    n = len(P)
    
    for pos in sorted(candidates):
        # Check if position is valid
        if 0 <= pos <= len(T) - n:
            substring = T[pos:pos + n]
            distance = hamming_distance(P, substring)
            
            if distance <= k:
                matches.append((pos, distance))
    
    return matches


# ============================================================================
# Complete Pigeonhole-based Approximate Matching
# ============================================================================

def pigeonhole_approximate_match(P, T, k):
    """
    Find all approximate matches of pattern P in text T with at most k mismatches.
    Uses the pigeonhole principle for filtering.
    
    Strategy:
    1. Partition P into k+1 factors
    2. Search each factor exactly in T (seed phase)
    3. Verify candidates using Hamming distance (verify phase)
    
    Args:
        P: Pattern string
        T: Text string
        k: Maximum number of mismatches allowed
    
    Returns:
        List of tuples (position, distance) where pattern matches with distance <= k
    """
    # Task 3.1: Partition the pattern
    factors = partition_pattern(P, k)
    
    # Task 3.2: Search for exact occurrences of factors
    candidates = search_factors_exact(factors, T)
    
    # Task 3.3: Verify candidate alignments
    matches = verify_candidates(candidates, P, T, k)
    
    return matches


# ============================================================================
# Naive search for comparison
# ============================================================================

def naive_approximate_match(P, T, k):
    """
    Naive approximate search that checks every position.
    Used for verification and comparison.
    
    Args:
        P: Pattern string
        T: Text string
        k: Maximum number of mismatches allowed
    
    Returns:
        List of tuples (position, distance) where distance <= k
    """
    matches = []
    n = len(P)
    
    for i in range(len(T) - n + 1):
        substring = T[i:i + n]
        distance = hamming_distance(P, substring)
        if distance <= k:
            matches.append((i, distance))
    
    return matches


# ============================================================================
# Main demonstration
# ============================================================================

def main():
    print("=" * 80)
    print("SEQUENCE MAPPER - PIGEONHOLE PRINCIPLE APPROACH")
    print("=" * 80)
    
    # Configuration
    M = 1000  # Length of T
    N = 20   # Length of P
    k = 5    # Maximum mismatches allowed
    num_edits = 3  # Number of edits to make to pattern
    
    # Generate random text
    random.seed(42)  # For reproducibility
    T = ''.join(random.choices('ACGT', k=M))
    
    # Extract pattern from a random position in the text
    pattern_start = random.randint(0, M - N)
    P_original = T[pattern_start:pattern_start + N]
    
    # Make num_edits random substitutions to create the query pattern
    P = list(P_original)
    edit_positions = random.sample(range(N), num_edits)
    nucleotides = ['A', 'C', 'G', 'T']
    
    for pos in edit_positions:
        # Choose a different nucleotide
        original_char = P[pos]
        possible_chars = [c for c in nucleotides if c != original_char]
        P[pos] = random.choice(possible_chars)
    
    P = ''.join(P)
    
    print(f"\nConfiguration:")
    print(f"  Text length (M): {M}")
    print(f"  Pattern length (N): {N}")
    print(f"  Max mismatches (k): {k}")
    print(f"  Number of edits made: {num_edits}")
    print(f"\nText T (first 60 chars): {T[:60]}...")
    print(f"\nPattern extracted from position {pattern_start} with {num_edits} edits:")
    print(f"  Original: {P_original}")
    print(f"  Modified: {P}")
    print(f"  Edit positions: {sorted(edit_positions)}")
    
    # Show the differences
    diff_line = ""
    for i in range(N):
        if P_original[i] == P[i]:
            diff_line += "|"
        else:
            diff_line += "X"
    print(f"            {diff_line}")
    
    # ========================================================================
    # Task 3.1: Partition the pattern
    # ========================================================================
    print("\n" + "=" * 80)
    print("TASK 3.1: PARTITIONING THE PATTERN")
    print("=" * 80)
    
    factors = partition_pattern(P, k)
    print(f"\nPattern P divided into {k+1} factors (pigeonhole principle):")
    print(f"At least one factor must match exactly if P occurs with ≤{k} mismatches\n")
    
    for i, (factor, start_pos) in enumerate(factors, 1):
        print(f"  Factor {i}: '{factor}' (position {start_pos}, length {len(factor)})")
    
    # Verify partition
    reconstructed = ''.join([f[0] for f in factors])
    print(f"\nVerification: Reconstructed pattern = '{reconstructed}'")
    print(f"              Original pattern      = '{P}'")
    print(f"              Match: {reconstructed == P}")
    
    # ========================================================================
    # Task 3.2 & 3.3: Search and verify
    # ========================================================================
    print("\n" + "=" * 80)
    print("TASK 3.2 & 3.3: FACTOR SEARCH AND VERIFICATION")
    print("=" * 80)
    
    # Search using pigeonhole approach
    print(f"\nSearching for pattern in text...")
    print(f"Expected to find original position {pattern_start} with distance {num_edits}...")
    matches = pigeonhole_approximate_match(P, T, k)
    
    print(f"\nFound {len(matches)} approximate matches (≤{k} mismatches):")
    
    # Check if we found the expected position
    found_original = any(pos == pattern_start for pos, _ in matches)
    if found_original:
        print(f"✓ EXPECTED position {pattern_start} was found!")
    else:
        print(f"✗ WARNING: Expected position {pattern_start} NOT found!")
    
    if matches:
        # Show first 10 matches
        for pos, dist in matches[:10]:
            substring = T[pos:pos + N]
            marker = " ← ORIGINAL" if pos == pattern_start else ""
            print(f"\n  Position {pos}: distance = {dist}{marker}")
            print(f"    P: {P}")
            print(f"    T: {substring}")
            
            # Show alignment with markers
            alignment = ""
            for i in range(len(P)):
                if P[i] == substring[i]:
                    alignment += "|"
                else:
                    alignment += "X"
            print(f"       {alignment}")
        
        if len(matches) > 10:
            print(f"\n  ... and {len(matches) - 10} more matches")
    else:
        print("  No matches found.")
    
    # ========================================================================
    # Verification against naive search
    # ========================================================================
    print("\n" + "=" * 80)
    print("VERIFICATION: COMPARISON WITH NAIVE SEARCH")
    print("=" * 80)
    
    naive_matches = naive_approximate_match(P, T, k)
    
    print(f"\nResults:")
    print(f"  Pigeonhole approach: {len(matches)} matches")
    print(f"  Naive approach:      {len(naive_matches)} matches")
    print(f"  Results match:       {matches == naive_matches}")
    
    if matches != naive_matches:
        print("\n  WARNING: Results differ!")
        print(f"  Pigeonhole only: {set(matches) - set(naive_matches)}")
        print(f"  Naive only:      {set(naive_matches) - set(matches)}")
    
    # ========================================================================
    # Example with known pattern
    # ========================================================================
    print("\n" + "=" * 80)
    print("EXAMPLE WITH KNOWN PATTERN")
    print("=" * 80)
    
    # Create a text with known occurrences
    T_example = "ACGTACGTACGTACGTACGTACGT" + "GGGGGGGGGG" + "ACGTACGTACGTACGT"
    P_example = "ACGTACGT"
    k_example = 1
    
    print(f"\nText:    {T_example}")
    print(f"Pattern: {P_example}")
    print(f"k:       {k_example}")
    
    # Partition
    factors_ex = partition_pattern(P_example, k_example)
    print(f"\nFactors (k+1 = {k_example+1}):")
    for i, (factor, pos) in enumerate(factors_ex, 1):
        print(f"  Factor {i}: '{factor}' at position {pos}")
    
    # Search
    matches_ex = pigeonhole_approximate_match(P_example, T_example, k_example)
    
    print(f"\nMatches found: {len(matches_ex)}")
    for pos, dist in matches_ex:
        substring = T_example[pos:pos + len(P_example)]
        marker = "EXACT" if dist == 0 else f"{dist} mismatch(es)"
        print(f"  Position {pos:2d}: {substring} ({marker})")
    
    # ========================================================================
    # Performance insight
    # ========================================================================
    print("\n" + "=" * 80)
    print("ALGORITHM EFFICIENCY")
    print("=" * 80)
    
    total_positions = len(T) - len(P) + 1
    candidates = search_factors_exact(factors, T)
    
    print(f"\nNaive approach:")
    print(f"  Must check all {total_positions} positions")
    print(f"  Complexity: O(M * N) where M={M}, N={N}")
    
    print(f"\nPigeonhole approach:")
    print(f"  Generated {len(candidates)} candidate positions to verify")
    print(f"  Reduction: {100 * (1 - len(candidates)/total_positions):.1f}% fewer positions checked")
    print(f"  Uses {k+1} factors for filtering")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()