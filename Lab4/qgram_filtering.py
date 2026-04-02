"""
Q-gram Filtering for Approximate String Matching
Lab 5 - Tasks 2.1, 2.2, and 2.3

This module implements q-gram based filtering for approximate string matching
using the seed-and-verify strategy.
"""


def hamming_distance(s1, s2):
    """Compute the Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        raise ValueError("Strings must be of the same length")
    
    distance = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            distance += 1
    return distance


# ============================================================================
# Task 2.1: Enumerating q-grams
# ============================================================================

def enumerate_qgrams(P, q):
    """
    Generate all q-grams of pattern P.
    
    Given a pattern P of length n, there are n-q+1 overlapping q-grams,
    namely P[0..q-1], P[1..q], P[2..q+1], ..., P[n-q..n-1].
    
    Args:
        P: Pattern string
        q: Length of q-grams
    
    Returns:
        List of q-grams (substrings of length q)
    
    Example:
        >>> enumerate_qgrams("ACGT", 2)
        ['AC', 'CG', 'GT']
    """
    if q > len(P):
        return []
    
    qgrams = []
    for i in range(len(P) - q + 1):
        qgrams.append(P[i:i+q])
    return qgrams


# ============================================================================
# Task 2.2: Counting shared q-grams
# ============================================================================

def count_shared_qgrams(P, T, q):
    """
    Count how many q-grams of P appear in T.
    
    Args:
        P: Pattern string
        T: Text string
        q: Length of q-grams
    
    Returns:
        Number of q-grams from P that appear in T
    """
    # Get all q-grams from P
    pattern_qgrams = enumerate_qgrams(P, q)
    
    # Get all q-grams from T and store in a set for fast lookup
    text_qgrams = set(enumerate_qgrams(T, q))
    
    # Count how many pattern q-grams appear in text
    shared_count = 0
    for qgram in pattern_qgrams:
        if qgram in text_qgrams:
            shared_count += 1
    
    return shared_count


def verify_qgram_lemma(P, T, q, k):
    """
    Verify if P could match T with at most k mismatches using the q-gram lemma.
    
    Q-gram Lemma: If P aligns to a substring of T with at most k mismatches
    under the Hamming distance model, then the alignment must contain at least
    n - q + 1 - k*q exact matches of q-grams.
    
    Explanation: A single mismatch can affect at most q q-grams (all those that
    overlap the mismatched position). Therefore, k mismatches can invalidate
    at most k*q q-grams, leaving at least (n-q+1) - k*q exact q-gram matches.
    
    Args:
        P: Pattern string
        T: Text string
        q: Length of q-grams
        k: Maximum number of mismatches allowed
    
    Returns:
        Tuple (shared_count, required_count, is_possible)
        - shared_count: Number of q-grams from P that appear in T
        - required_count: Minimum q-grams required by the lemma
        - is_possible: Whether the condition is satisfied
    """
    n = len(P)
    
    # Calculate minimum number of shared q-grams required by lemma
    required_count = n - q + 1 - k * q
    
    # Count actual shared q-grams
    shared_count = count_shared_qgrams(P, T, q)
    
    # Check if the lemma condition is satisfied
    is_possible = shared_count >= required_count
    
    return shared_count, required_count, is_possible


# ============================================================================
# Task 2.3: Seed-and-verify approximate search
# ============================================================================

def build_qgram_index(T, q):
    """
    Build a hash table index mapping q-grams to their positions in text T.
    
    This index allows quick lookup of all positions where a given q-gram
    occurs in the text, which is essential for the seed-and-verify strategy.
    
    Args:
        T: Text string
        q: Length of q-grams
    
    Returns:
        Dictionary mapping q-grams to list of positions where they occur
    """
    index = {}
    for i in range(len(T) - q + 1):
        qgram = T[i:i+q]
        if qgram not in index:
            index[qgram] = []
        index[qgram].append(i)
    return index


def seed_and_verify_search(P, T, k, q=None):
    """
    Approximate search using seed-and-verify strategy with q-grams.
    
    Strategy:
    1. SEED: Use q-grams to identify candidate positions where P might match
    2. VERIFY: Check each candidate using Hamming distance
    
    The key insight is that if P occurs in T with at most k mismatches,
    then P must contain an exact substring of length at least floor(n/(k+1)).
    This is because k+1 disjoint substrings of length floor(n/(k+1)) cannot
    all be affected by only k mismatches.
    
    Args:
        P: Pattern string
        T: Text string
        k: Maximum number of mismatches allowed
        q: Length of q-grams (if None, use optimal q = floor(n/(k+1)))
    
    Returns:
        List of tuples (position, distance) where distance <= k, sorted by position
    """
    n = len(P)
    
    # Choose q based on the seed guarantee if not provided
    if q is None:
        q = n // (k + 1)
        if q == 0:
            q = 1  # Ensure q is at least 1
    
    # Build q-gram index for the text
    qgram_index = build_qgram_index(T, q)
    
    # Get all q-grams from pattern
    pattern_qgrams = enumerate_qgrams(P, q)
    
    # SEED PHASE: Collect candidate positions
    # For each q-gram in the pattern, find where it occurs in the text
    # and infer potential alignment positions for the whole pattern
    candidates = set()
    
    for i, qgram in enumerate(pattern_qgrams):
        if qgram in qgram_index:
            # Each occurrence of this q-gram suggests a candidate alignment
            # If q-gram at position i in pattern matches position j in text,
            # then pattern might start at position j - i
            for text_pos in qgram_index[qgram]:
                pattern_start = text_pos - i
                # Ensure candidate is within valid bounds
                if 0 <= pattern_start <= len(T) - len(P):
                    candidates.add(pattern_start)
    
    # VERIFY PHASE: Check each candidate position using Hamming distance
    matches = []
    for pos in sorted(candidates):
        substring = T[pos:pos+len(P)]
        distance = hamming_distance(substring, P)
        if distance <= k:
            matches.append((pos, distance))
    
    return matches


def naive_approximate_search(P, T, k):
    """
    Naive approximate search that checks every position.
    
    This is used for comparison to verify correctness of the seed-and-verify
    approach. It's O(n*m) where n is text length and m is pattern length.
    
    Args:
        P: Pattern string
        T: Text string
        k: Maximum number of mismatches allowed
    
    Returns:
        List of tuples (position, distance) where distance <= k
    """
    matches = []
    for i in range(len(T) - len(P) + 1):
        substring = T[i:i+len(P)]
        distance = hamming_distance(substring, P)
        if distance <= k:
            matches.append((i, distance))
    return matches


# ============================================================================
# Main demonstration
# ============================================================================

def main():
    print("=" * 80)
    print("Q-GRAM FILTERING FOR APPROXIMATE STRING MATCHING")
    print("=" * 80)
    
    # Example sequences
    T = "CACATGCTAGCTAGCTAGCTCAATGCAAATACTTATCGCATCACTAATCGACGAATATCTACATCATCGATGCAATATTATCGGATCGACTAGCTAGCGAATCGGCGACTACGATGCATGCGCGATGCGATGCATGCTAGCTAGCGCGCGGCATGACTGAGTCATCGTCAGCAT"
    P = "TCGATGCAATACCATCG"
    
    print(f"\nText T:    {T}")
    print(f"Pattern P: {P}")
    print(f"Length of T: {len(T)}, Length of P: {len(P)}")
    
    # ========================================================================
    # Task 2.1: Enumerate q-grams
    # ========================================================================
    print("\n" + "=" * 80)
    print("TASK 2.1: ENUMERATING Q-GRAMS")
    print("=" * 80)
    
    q = 3
    pattern_qgrams = enumerate_qgrams(P, q)
    print(f"\nQ-grams of pattern P with q={q}:")
    print(f"Number of q-grams: {len(pattern_qgrams)} (expected: {len(P) - q + 1})")
    for i, qgram in enumerate(pattern_qgrams):
        print(f"  Position {i}: '{qgram}'")
    
    # Test with different q values
    print(f"\nQ-grams with different values of q:")
    for test_q in [2, 3, 4, 5]:
        qgrams = enumerate_qgrams(P, test_q)
        print(f"  q={test_q}: {qgrams} (count: {len(qgrams)})")
    
    # ========================================================================
    # Task 2.2: Counting shared q-grams
    # ========================================================================
    print("\n" + "=" * 80)
    print("TASK 2.2: COUNTING SHARED Q-GRAMS AND Q-GRAM LEMMA")
    print("=" * 80)
    
    k = 2  # Maximum mismatches
    q = 3
    
    # Compare P with a substring of T
    test_substring = T[2:2+len(P)]  # Extract substring from T
    print(f"\nComparing pattern P with substring at position 2:")
    print(f"  P:         {P}")
    print(f"  Substring: {test_substring}")
    print(f"  Hamming distance: {hamming_distance(P, test_substring)}")
    
    shared, required, is_possible = verify_qgram_lemma(P, test_substring, q, k)
    print(f"\nQ-gram lemma verification (q={q}, k={k}):")
    print(f"  Total q-grams in P: {len(P) - q + 1}")
    print(f"  Shared q-grams: {shared}")
    print(f"  Required by lemma: {required} (formula: n-q+1-kq = {len(P)}-{q}+1-{k}*{q})")
    print(f"  Match is possible: {is_possible}")
    
    # Test with different k values
    print(f"\nTesting q-gram lemma with different k values:")
    print(f"  {'k':<5} {'Shared':<8} {'Required':<10} {'Possible':<10}")
    print(f"  {'-'*5} {'-'*8} {'-'*10} {'-'*10}")
    for test_k in range(0, 5):
        shared, required, is_possible = verify_qgram_lemma(P, test_substring, q, test_k)
        print(f"  {test_k:<5} {shared:<8} {required:<10} {str(is_possible):<10}")
    
    # ========================================================================
    # Task 2.3: Seed-and-verify approximate search
    # ========================================================================
    print("\n" + "=" * 80)
    print("TASK 2.3: SEED-AND-VERIFY APPROXIMATE SEARCH")
    print("=" * 80)
    
    k = 2
    q_optimal = len(P) // (k + 1)
    print(f"\nSearching for pattern P in text T with at most k={k} mismatches")
    print(f"Optimal q-gram length: q = floor(n/(k+1)) = floor({len(P)}/{k+1}) = {q_optimal}")
    
    # Perform seed-and-verify search
    matches = seed_and_verify_search(P, T, k)
    
    print(f"\nMatches found: {len(matches)}")
    for pos, dist in matches:
        substring = T[pos:pos+len(P)]
        print(f"\n  Position {pos}: distance={dist}")
        print(f"    P:      {P}")
        print(f"    T[{pos}]: {substring}")
        # Show mismatches
        mismatches = [i for i in range(len(P)) if P[i] != substring[i]]
        if mismatches:
            print(f"    Mismatches at positions: {mismatches}")
            # Visual representation of differences
            diff_line = "".join(["|" if P[i] == substring[i] else "X" for i in range(len(P))])
            print(f"            {diff_line}")
    
    # Compare with naive search
    print("\n" + "-" * 80)
    print("Verification: Comparison with naive search")
    print("-" * 80)
    naive_matches = naive_approximate_search(P, T, k)
    
    print(f"Naive search found:        {len(naive_matches)} matches")
    print(f"Seed-and-verify found:     {len(matches)} matches")
    print(f"Results match:             {naive_matches == matches}")
    
    # Show efficiency gain
    total_positions = len(T) - len(P) + 1
    candidates_checked = len(matches)  # In practice, more candidates are checked
    print(f"\nEfficiency:")
    print(f"  Total positions to check: {total_positions}")
    print(f"  Candidates verified:      varies (seeds filter most positions)")
    
    # ========================================================================
    # Additional example: Genomic sequence
    # ========================================================================
    print("\n" + "=" * 80)
    print("GENOMIC EXAMPLE: SEARCHING IN A LONGER SEQUENCE")
    print("=" * 80)
    
    # Simulated genomic sequence (150 bp)
    genome = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + \
             "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA" + \
             "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA"
    
    # Pattern with some mutations
    query = "ATCGAACGATCG"
    
    print(f"\nGenome (length {len(genome)} bp):")
    print(f"  {genome[:60]}...")
    print(f"\nQuery: {query} (length {len(query)} bp)")
    
    k = 2
    print(f"\nSearching with k={k} mismatches...")
    matches = seed_and_verify_search(query, genome, k)
    
    print(f"Found {len(matches)} approximate matches:")
    for pos, dist in matches[:10]:  # Show first 10 matches
        substring = genome[pos:pos+len(query)]
        print(f"  Position {pos:3d}: distance={dist} -> {substring}")
    
    # ========================================================================
    # Performance comparison example
    # ========================================================================
    print("\n" + "=" * 80)
    print("PERFORMANCE COMPARISON")
    print("=" * 80)
    
    import time
    
    # Test on a larger sequence
    large_text = genome * 10  # 1500 bp
    
    print(f"\nTesting on a text of length {len(large_text)} bp")
    print(f"Pattern length: {len(query)} bp, k={k}")
    
    # Naive search
    start = time.time()
    naive_result = naive_approximate_search(query, large_text, k)
    naive_time = time.time() - start
    
    # Seed-and-verify search
    start = time.time()
    seed_result = seed_and_verify_search(query, large_text, k)
    seed_time = time.time() - start
    
    print(f"\nNaive search:")
    print(f"  Time: {naive_time*1000:.3f} ms")
    print(f"  Matches: {len(naive_result)}")
    
    print(f"\nSeed-and-verify search:")
    print(f"  Time: {seed_time*1000:.3f} ms")
    print(f"  Matches: {len(seed_result)}")
    
    print(f"\nSpeedup: {naive_time/seed_time:.2f}x")
    print(f"Results match: {naive_result == seed_result}")
    
    # ========================================================================
    # Visualization (optional)
    # ========================================================================
    print("\n" + "=" * 80)
    print("HAMMING DISTANCE VISUALIZATION")
    print("=" * 80)
    
    T_vis = "CACATGCTAGCTAGCTAGCTCAATGCAAATACTTATCGCATCACTAATCGACGAATATCTACATCATCGATGCAATATTATCGGATCGACTAGCTAGCGAATCGGCGACTACGATGCATGCGCGATGCGATGCATGCTAGCTAGCGCGCGGCATGACTGAGTCTCGATGCAATACCATCGATCGTCAGCAT"
    P_vis = "TCGATGCAATACCATCG"
    # T_vis = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
    # P_vis = "CTCGCTAGC"
    
    distances = []
    for i in range(len(T_vis) - len(P_vis) + 1):
        substring = T_vis[i:i+len(P_vis)]
        distance = hamming_distance(substring, P_vis)
        distances.append((i, distance))
    
    print(f"\nHamming distances at each position (text length {len(T_vis)}):")
    print(f"Position  Distance  Visualization")
    print("-" * 50)
    for i, distance in distances[:30]:  # Show first 30
        bar = "█" * (10 - distance) + "░" * distance
        print(f"   {i:2d}        {distance}       {bar}")
    
    # Optional: Create heatmap visualization if matplotlib is available
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6))
        
        # Plot 1: Line plot of Hamming distances
        positions = [i for i, _ in distances]
        dists = [d for _, d in distances]
        
        ax1.plot(positions, dists, 'b-', linewidth=2, marker='o', markersize=4)
        ax1.axhline(y=k, color='r', linestyle='--', label=f'Threshold (k={k})')
        ax1.fill_between(positions, 0, dists, alpha=0.3)
        ax1.set_xlabel('Position in Text')
        ax1.set_ylabel('Hamming Distance')
        ax1.set_title('Hamming Distance at Each Position')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Heatmap
        heatmap_data = np.array([dists])
        im = ax2.imshow(heatmap_data, cmap='RdYlGn_r', aspect='auto', 
                       extent=[0, len(positions), 0, 1])
        ax2.set_xlabel('Position in Text')
        ax2.set_yticks([])
        ax2.set_title('Hamming Distance Heatmap (Red = high, Green = low)')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax2, orientation='vertical')
        cbar.set_label('Hamming Distance')
        
        plt.tight_layout()
        plt.savefig('/Users/juanjo/UPC/TEB/Lab5/hamming_visualization.png', dpi=150, bbox_inches='tight')
        print(f"\n✓ Visualization saved to 'hamming_visualization.png'")
        plt.show()
        
    except ImportError:
        print("\n(Matplotlib not available for visualization)")


if __name__ == "__main__":
    main()
