# =============================================================================
# alg_smith_waterman.py — Smith-Waterman Local Alignment Algorithm
#
# Finds the single best-matching REGION between two sequences, unlike
# Hirschberg/LCS which align end-to-end. More accurate for distantly
# related species where only a conserved region is shared.
#
# Scores never go below 0 — negative scores reset to 0, ending the current
# alignment region and allowing a new one to start elsewhere in the matrix.
# This floor is what makes it local rather than global.
#
# Scoring: +2 match  |  -1 mismatch  |  -1 gap
# Time: O(m * n)  |  Space: O(m * n) — full matrix kept for traceback
# =============================================================================


def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-1):
    """
    Returns (aligned_str, raw_score, percentage) for the best local region.
    Percentage is normalised against a perfect match across the shorter sequence.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Row 0 and col 0 stay zero — no gap penalties at the boundary.
    # This lets the best alignment start anywhere (local, not global).
    dp        = [[0] * cols for _ in range(rows)]
    best_score = 0
    best_pos   = (0, 0)

    # Fill the matrix
    for i in range(1, rows):
        for j in range(1, cols):
            diagonal   = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            score_up   = dp[i-1][j] + gap    # gap in seq2 (skip seq1 char)
            score_left = dp[i][j-1] + gap    # gap in seq1 (skip seq2 char)

            dp[i][j] = max(0, diagonal, score_up, score_left)

            if dp[i][j] > best_score:
                best_score = dp[i][j]
                best_pos   = (i, j)

    # Traceback from the peak back to the first 0.
    # Only exact matches contribute a character to the output string.
    # Mismatches are diagonal moves that don't record a character —
    # without this explicit branch they fall into the left-gap case (wrong).
    aligned = []
    i, j    = best_pos

    while i > 0 and j > 0 and dp[i][j] > 0:
        current  = dp[i][j]
        diag_val = dp[i-1][j-1]
        up_val   = dp[i-1][j]

        if seq1[i-1] == seq2[j-1] and current == diag_val + match:
            aligned.append(seq1[i-1])   # match — record and move diagonally
            i -= 1; j -= 1
        elif current == diag_val + mismatch:
            i -= 1; j -= 1              # mismatch — diagonal move, no character
        elif current == up_val + gap:
            i -= 1                      # gap in seq2
        else:
            j -= 1                      # gap in seq1

    aligned_str  = "".join(reversed(aligned))
    max_possible = match * min(len(seq1), len(seq2))
    percentage   = (best_score / max_possible * 100) if max_possible > 0 else 0.0

    return aligned_str, best_score, percentage


if __name__ == "__main__":
    s1 = "GGTTGACTA"
    s2 = "TGTTACGG"
    aligned, score, pct = smith_waterman(s1, s2)
    print(f"Classic test → aligned: {aligned}, score: {score}  (expect region GTTAC)")

    dna_A = "ATCGTACG"
    dna_B  = "ATGCTAC"
    aligned, score, pct = smith_waterman(dna_A, dna_B)
    print(f"DNA test     → aligned: {aligned}, score: {score}, similarity: {pct:.2f}%")