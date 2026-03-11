# 🧬 DNA Sequence Matching - An Application of LCS (Hirschberg's Algorithm and Smith-Waterman Algorithm)

> **How similar are two species, genetically?** This app fetches real DNA sequences from NCBI's public database and runs two classic bioinformatics algorithms to answer that question — with scores, alignments, and visual comparisons.

---
VISIT :
https://genosync---the-dna-matcher-uszpk6jhtdr2x3idzdiqgn.streamlit.app/
## 📖 Table of Contents

- [What This App Does](#what-this-app-does)
- [The Biology Behind It](#the-biology-behind-it)
  - [What is DNA?](#what-is-dna)
  - [Why Compare Sequences?](#why-compare-sequences)
  - [Gene Markers — What Are They?](#gene-markers--what-are-they)
- [Where the Data Comes From](#where-the-data-comes-from)
  - [NCBI & GenBank](#ncbi--genbank)
  - [How Sequences Are Fetched](#how-sequences-are-fetched)
  - [Marker Priority Chains](#marker-priority-chains)
- [The Three Algorithms](#the-three-algorithms)
  - [1. LCS — Longest Common Subsequence](#1-lcs--longest-common-subsequence)
  - [2. Hirschberg's Algorithm](#2-hirschbergs-algorithm)
  - [3. Smith-Waterman](#3-smith-waterman)
  - [Global vs Local Alignment — What's the Difference?](#global-vs-local-alignment--whats-the-difference)
- [How Similarity Score Is Calculated](#how-similarity-score-is-calculated)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Deployment](#deployment)
- [Limitations & Caveats](#limitations--caveats)

---

## What This App Does

1. You type in **two species** (common names like "human" and "chimpanzee")
2. The app resolves them to scientific names via the NCBI Taxonomy database
3. It fetches a **real DNA sequence** for each species — from the same gene, so the comparison is meaningful
4. It runs **two alignment algorithms** on those sequences
5. You get a **similarity percentage**, aligned sequences, and a score

---

## The Biology Behind It

### What is DNA?

DNA (Deoxyribonucleic Acid) is the molecule that carries the genetic instructions for all living organisms. It's made of four chemical **bases**:

| Base | Letter | Pairs With |
|------|--------|------------|
| Adenine | `A` | Thymine (T) |
| Thymine | `T` | Adenine (A) |
| Cytosine | `C` | Guanine (G) |
| Guanine | `G` | Cytosine (C) |

A DNA sequence is simply a long string of these letters — e.g. `ATCGGCTAGCTA...` — and the entire human genome contains about **3 billion** of them.

### Why Compare Sequences?

When two species share a common ancestor, their DNA sequences are similar. The more similar the sequences, the more recently they diverged. This is the foundation of **molecular phylogenetics** — building evolutionary trees from genetic data.

For example:
- Human vs Chimpanzee → ~98% similarity on coding genes
- Human vs Mouse → ~85% on conserved regions
- Human vs Banana → ~60% on fundamental cellular genes

By comparing specific gene regions, scientists can determine how closely related two species are — without needing fossils.

### Gene Markers — What Are They?

Not all parts of DNA mutate at the same rate. **Gene markers** are specific, well-studied regions chosen for comparison because they are:

- **Conserved enough** to exist in both species
- **Variable enough** to show differences between species
- **Short enough** to be practical (~500–1600 base pairs)

This app uses the following standard markers:

| Marker | Full Name | Used For | Approx. Length |
|--------|-----------|----------|----------------|
| **COI** | Cytochrome c Oxidase I | Animals (standard DNA barcode) | ~650 bp |
| **16S** | 16S Ribosomal RNA | Animals, bacteria | ~1500 bp |
| **12S** | 12S Ribosomal RNA | Animals (especially vertebrates) | ~950 bp |
| **cytb** | Cytochrome b | Animals (mammals especially) | ~1140 bp |
| **rbcL** | RuBisCO large subunit | Plants | ~550 bp |
| **matK** | Maturase K | Plants | ~850 bp |
| **ITS** | Internal Transcribed Spacer | Fungi, some plants | ~600 bp |
| **LSU** | Large Subunit rRNA | Fungi | ~2500 bp |
| **SSU** | Small Subunit rRNA | Fungi | ~1800 bp |

> ⚠️ It's critical that both species are compared on the **same gene**. Comparing COI from one species against cytb from another would be like comparing someone's height to another person's weight — the number is meaningless.

---

## Where the Data Comes From

### NCBI & GenBank

**NCBI** (National Center for Biotechnology Information) is a US government agency that maintains the world's largest public biological database. **GenBank** is their nucleotide sequence repository — containing over 200 million sequences from hundreds of thousands of organisms, all freely accessible.

This app uses NCBI's **Entrez API** (via Biopython) to query GenBank in real time.

### How Sequences Are Fetched

```
User types "dog" and "wolf"
        ↓
NCBI Taxonomy lookup → resolves to scientific names
  "dog"  → Canis lupus familiaris
  "wolf" → Canis lupus
        ↓
Try marker chain: COI → 16S → 12S → cytb ...
        ↓
Search: "Canis lupus familiaris[Organism] AND COI[Gene] AND 500[SLEN]:1600[SLEN]"
        ↓
Both species have COI? → fetch both sequences → proceed
One missing?          → try next marker in chain
        ↓
Return sequences to alignment algorithms
```

### Marker Priority Chains

The app tries markers **in order** until it finds one where **both** species have a result:

```python
"animal":  ["COI", "16S", "12S", "cytb"]
"plant":   ["rbcL", "matK", "ITS"]
"fungus":  ["ITS", "LSU", "SSU"]
"default": ["COI", "16S", "ITS", "rbcL", "12S"]
```

Results are **cached for 24 hours** — so fetching "Human vs Chimp" and then "Human vs Gorilla" reuses the Human sequence without a second API call.

---

## The Three Algorithms

### 1. LCS — Longest Common Subsequence

**What it is:**
LCS finds the longest sequence of characters that appear in the same order in both strings — but not necessarily contiguously. It's the mathematical foundation that both Hirschberg and (partially) Smith-Waterman are built on.

**DNA example:**
```
Seq A:  A T C G G A T
Seq B:  A T G G A C T
LCS  :  A T G G A T   (length 6)
```

**Why it matters:**
The length of the LCS divided by the sequence length gives a direct similarity score. A longer LCS = more shared evolutionary history.

**Pseudocode:**
```
function LCS(A, B):
    m = length of A
    n = length of B

    create table dp[0..m][0..n], all zeros

    for i from 1 to m:
        for j from 1 to n:
            if A[i] == B[j]:
                dp[i][j] = dp[i-1][j-1] + 1
            else:
                dp[i][j] = max(dp[i-1][j], dp[i][j-1])

    return dp[m][n]   ← this is the LCS length
```

**Time complexity:** O(m × n)
**Space complexity:** O(m × n) — the full table must be stored

> The space cost is the problem. For sequences of 1000bp each, this is a 1,000,000-cell table. Hirschberg solves this.

---

### 2. Hirschberg's Algorithm

**What it is:**
Hirschberg's algorithm computes the **global alignment** of two sequences — aligning them from start to finish. It produces the same result as naive LCS-based alignment but uses only **O(n) space** instead of O(m × n), by using a divide-and-conquer strategy.

**Global alignment** means both sequences are aligned in their entirety. Every base in both sequences must be accounted for — either matched, mismatched, or gapped.

**The key insight — divide and conquer:**
Instead of storing the entire DP table, Hirschberg:
1. Finds the midpoint of sequence A
2. Uses the **NeedlemanWunsch score** (just the last two rows, forward and backward) to find where that midpoint aligns in sequence B
3. Recursively aligns the left half and the right half

```
      A (full)
   ┌──────┬──────┐
   left   mid   right
     ↓           ↓
   align        align
  A_left       A_right
  with         with
  B[0..k]    B[k+1..n]
```

**Pseudocode:**
```
function Hirschberg(A, B):
    if length(A) == 0:
        return gap * length(B)
    if length(B) == 0:
        return gap * length(A)
    if length(A) == 1:
        return NeedlemanWunsch(A, B)   ← base case

    mid = length(A) / 2

    scoreL = NWScore(A[0..mid], B)          ← forward pass, last row only
    scoreR = NWScore(reverse(A[mid..end]),  ← backward pass, last row only
                     reverse(B))

    k = argmax over j of (scoreL[j] + scoreR[n-j])  ← optimal split point in B

    leftAlign  = Hirschberg(A[0..mid],  B[0..k])
    rightAlign = Hirschberg(A[mid..end], B[k..end])

    return leftAlign + rightAlign
```

**Time complexity:** O(m × n)
**Space complexity:** O(m + n) ← the big win over plain LCS

**What this app uses it for:**
Hirschberg gives the **global similarity percentage** — how similar the two sequences are across their entire shared length. This is the headline number.

---

### 3. Smith-Waterman

**What it is:**
Smith-Waterman finds the **local alignment** — the single best-matching region between two sequences, ignoring the rest. It doesn't force the whole sequence to align; it finds where they match *most*.

**Local alignment** means finding the substring of A and the substring of B that are most similar to each other. Mismatching flanks are ignored.

**Why this matters for DNA:**
Two distantly related species might share a highly conserved functional region (e.g., an active site in a gene) even if the surrounding sequence has diverged heavily. Smith-Waterman finds that conserved region.

**The key difference from LCS/Hirschberg:**
The DP cell is never allowed to go below zero. When the score would go negative, it resets to 0 — meaning "start a new alignment here rather than drag in a bad match."

**Pseudocode:**
```
function SmithWaterman(A, B, match=+2, mismatch=-1, gap=-2):
    m = length(A)
    n = length(B)

    create table dp[0..m][0..n], all zeros
    max_score = 0
    max_pos   = (0, 0)

    for i from 1 to m:
        for j from 1 to n:
            if A[i] == B[j]:
                diag = dp[i-1][j-1] + match
            else:
                diag = dp[i-1][j-1] + mismatch

            dp[i][j] = max(
                0,                    ← reset to zero (key difference!)
                diag,                 ← match or mismatch
                dp[i-1][j] + gap,     ← gap in B
                dp[i][j-1] + gap      ← gap in A
            )

            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_pos   = (i, j)

    ← Traceback from max_pos, stop when cell = 0
    aligned_A, aligned_B = traceback(dp, A, B, max_pos)

    similarity = matches_in_alignment / length_of_alignment * 100
    return aligned_A, aligned_B, max_score, similarity
```

**Time complexity:** O(m × n)
**Space complexity:** O(m × n)

**What this app uses it for:**
Smith-Waterman gives the **local similarity score and percentage** — showing the best-matching region and its score. A high SW score on a short region can indicate a conserved functional domain even in distantly related species.

---

### Global vs Local Alignment — What's the Difference?

| | Hirschberg (Global) | Smith-Waterman (Local) |
|--|--|--|
| **Aligns** | Entire sequence, start to finish | Best matching sub-region only |
| **Good for** | Closely related species | Distantly related, or finding conserved domains |
| **Penalises** | Every mismatch and gap | Ignores bad flanking regions |
| **Score of 0** | Impossible | Means "no local match at all" |
| **Output** | Full alignment | Substring alignment + score |

> **Rule of thumb:** Use Hirschberg % to answer "how similar are these species overall?" Use Smith-Waterman score to answer "do these species share any conserved region, and how strong is it?"

---

## How Similarity Score Is Calculated

Both sequences are trimmed to the same length before comparison:

```
compare_len = min(len(seq_A), len(seq_B), 1000)
```

This ensures neither species is penalised for having a longer sequence. The denominator is always `compare_len`.

**Hirschberg %:**
```
hirschberg_pct = (length of LCS returned by Hirschberg / compare_len) × 100
```

**Smith-Waterman %:**
```
sw_pct = (matching bases in local alignment / length of local alignment) × 100
```

These are two different questions — Hirschberg asks "how much do they share overall?" while SW asks "how good is their best matching region?"

---

## Project Structure

```
DNA-SEQUENCE-MATCHING/
│
├── app.py                  # Streamlit UI — all user interaction lives here
├── core.py                 # NCBI fetching, caching, and algorithm routing
├── alg_lcs.py              # Pure LCS implementation (used by Hirschberg)
├── alg_hirschberg.py       # Hirschberg global alignment
├── alg_smith_waterman.py   # Smith-Waterman local alignment
│
├── requirements.txt        # Python dependencies
├── .env                    # 🔒 Local secrets — NEVER commit this
├── .gitignore              # Ensures .env is excluded from git
└── README.md               # You are here
```

---

## Getting Started

### Prerequisites
- Python 3.9+
- A free [NCBI API key](https://www.ncbi.nlm.nih.gov/account/) (optional but recommended — raises rate limit from 3 to 10 requests/sec)

### Installation

```bash
# 1. Clone the repo
git clone https://github.com/srivshnu/GenoSync---The-DNA-Matcher.git
cd dna-sequence-matching

# 2. Install dependencies
pip install -r requirements.txt

# 3. Set up credentials
cp .env.example .env
# Edit .env and fill in your NCBI_EMAIL and NCBI_API_KEY

# 4. Run the app
streamlit run app.py
```

### `.env` file format
```
NCBI_EMAIL=you@example.com
NCBI_API_KEY=your_ncbi_api_key_here
```

## Limitations & Caveats

| Limitation | Explanation |
|---|---|
| **Sequence availability** | Not every species has every marker on NCBI. If no shared marker is found, no comparison is possible. |
| **Sequence length cap** | Sequences are capped at 1000bp for performance. Real barcodes are ~650bp so this rarely matters. |
| **Common name ambiguity** | "Dolphin" could match 30 species. The app picks the top NCBI taxonomy match — verify the scientific name shown. |
| **Not whole-genome** | This compares one gene marker, not the full genome. Results reflect that gene's evolutionary rate, not overall genomic similarity. |
| **NCBI rate limits** | Without an API key: 3 requests/sec. With API key: 10 requests/sec. The retry logic handles occasional failures. |
| **Cache is per-session** | Streamlit's `@st.cache_data` resets when the server restarts. This is a performance feature, not a data store. |

---

## References

- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
- [Biopython Entrez documentation](https://biopython.org/docs/latest/api/Bio.Entrez.html)
- Hirschberg, D.S. (1975). *A linear space algorithm for computing maximal common subsequences.* CACM.
- Smith, T.F. & Waterman, M.S. (1981). *Identification of common molecular subsequences.* Journal of Molecular Biology.
- [BOLD Systems — DNA Barcoding](https://www.boldsystems.org/)

---

<div align="center">
  Built with 🧬 Biopython · ⚡ Streamlit · 🐍 Python
</div>
