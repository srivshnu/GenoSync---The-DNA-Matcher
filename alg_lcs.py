# =============================================================================
# alg_lcs.py — LCS Similarity Score, Space-Optimised Sliding Window
# Returns similarity as a percentage — does NOT reconstruct the LCS string.
# See alg_hirschberg.py for string reconstruction.
#
# Time: O(m * n)  |  Space: O(n)  — vs standard LCS which is O(m * n) space
# =============================================================================


def lcs_percentage_only(seq1, seq2):
    """
    Computes LCS-based similarity between two sequences.
    Returns a float in [0.0, 100.0] relative to the LONGER sequence.

    Denominator is max(len) not min(len) — dividing by min inflates scores.
    e.g. "A" vs "AAAAAAAAAA" should be 10%, not 100%.
    """
    if not seq1 or not seq2:
        return 0.0

    # Capture max length BEFORE the swap below changes which is seq1/seq2
    original_max_len = max(len(seq1), len(seq2))

    # Put the shorter sequence as columns to minimise row width in memory.
    # LCS(A, B) == LCS(B, A) so swapping doesn't affect the result.
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1

    m = len(seq1)   # shorter — rows
    n = len(seq2)   # longer  — columns

    prev_row = [0] * (n + 1)
    curr_row = [0] * (n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                curr_row[j] = prev_row[j - 1] + 1              # match: extend diagonally
            else:
                curr_row[j] = max(curr_row[j - 1], prev_row[j])  # no match: best of left or above

        # [:] makes a value copy — without it prev_row and curr_row point to
        # the same list and overwrite each other on the next iteration
        prev_row = curr_row[:]

    return (prev_row[n] / original_max_len) * 100


if __name__ == "__main__":
    dna_A = "ATCGTACG"
    dna_B  = "ATGCTAC"
    print(f"Similarity    : {lcs_percentage_only(dna_A, dna_B):.2f}%")
    print(f"Identical     : {lcs_percentage_only('ATCG', 'ATCG'):.1f}%   (expect 100.0)")
    print(f"No overlap    : {lcs_percentage_only('AAAA', 'CCCC'):.1f}%   (expect 0.0)")
    print(f"Short vs long : {lcs_percentage_only('A', 'AAAAAAAAAA'):.1f}%  (expect 10.0)")