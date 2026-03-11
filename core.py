# =============================================================================
# core.py — NCBI Fetching, Caching, and Algorithm Orchestration
#
# Sits between app.py and the algorithms. Three responsibilities:
#   1. Fetch  — retrieves real DNA sequences from the NCBI API
#   2. Cache  — stores results for 24h so the same species is never re-fetched
#   3. Route  — runs both algorithms and returns a results dict to app.py
# =============================================================================

import time
import streamlit as st
from Bio import Entrez, SeqIO

from alg_hirschberg import hirschberg_reconstruct
from alg_smith_waterman import smith_waterman
import os
from dotenv import load_dotenv
load_dotenv()

Entrez.email   = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", "")

# ---------------------------------------------------------------------------
# Credentials — loaded from Streamlit Secrets (st.secrets) when deployed, or
# from a local .env file during development via the fallback block below.
# NEVER hardcode email or API keys directly in source code.
# ---------------------------------------------------------------------------
try:
    # Streamlit Cloud: set these under App Settings → Secrets
    Entrez.email   = st.secrets["NCBI_EMAIL"]
    Entrez.api_key = st.secrets.get("NCBI_API_KEY", "")

except (KeyError, FileNotFoundError):
    import os
    from dotenv import load_dotenv
    load_dotenv()
    Entrez.email   = os.environ.get("NCBI_EMAIL", "")
    Entrez.api_key = os.environ.get("NCBI_API_KEY", "")


# Gene marker priority chains — tried in order until BOTH species match the same gene.
# Mixing genes (e.g. COI for one species, cytb for another) produces meaningless scores.
# Standard barcodes: COI for animals (BOLD Systems), rbcL/matK for plants, ITS for fungi.
MARKER_CHAINS = {
    "animal":  ["COI", "16S", "12S", "cytb"],
    "plant":   ["rbcL", "matK", "ITS"],
    "fungus":  ["ITS", "LSU", "SSU"],
    "default": ["COI", "16S", "ITS", "rbcL", "12S"],
}

# Expected sequence length range per gene (base pairs).
# NCBI often returns full mitochondrial genomes (~16,000bp) when searching for COI
# because COI sits inside them. The SLEN filter and this range prevent that —
# without it we'd compare the wrong region entirely (D-loop instead of COI).
MARKER_LENGTH_RANGE = {
    "COI":  (500,  1600),   # barcode ~650bp, full gene ~1542bp
    "16S":  (400,  1800),
    "12S":  (300,  1100),
    "cytb": (800,  1200),
    "rbcL": (400,  600),
    "matK": (800,  900),
    "ITS":  (400,  800),
    "LSU":  (800,  3500),
    "SSU":  (1600, 2000),
}


def _fetch_with_retry(fn, retries=3, base_delay=1.0):
    """Calls fn() up to `retries` times with exponential backoff (1s → 2s → 4s).
    Returns the first non-None result, or None if all attempts fail."""
    for attempt in range(retries):
        try:
            result = fn()
            if result is not None:
                return result
        except Exception as e:
            print(f"[Retry {attempt + 1}/{retries}] {e}")
        time.sleep(base_delay * (2 ** attempt))
    return None


@st.cache_data(ttl=86400)
def search_species_options(common_name: str, max_results: int = 6):
    """
    Resolves a common name to a list of (common_name_label, scientific_name) tuples.
    Needed because common names are ambiguous — "dog" matches domestic dog,
    gray wolf, dingo, etc. Returns [] if nothing is found.
    """
    def _fetch():
        handle  = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_results)
        record  = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return []

        handle  = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        results = []
        for r in records:
            sci    = r.get("ScientificName", "")
            # GenbankCommonName is NCBI's standardised label; fall back to sci name
            common = r.get("GenbankCommonName") or r.get("CommonName") or sci
            results.append((common, sci))
        return results

    return _fetch_with_retry(_fetch) or []


@st.cache_data(ttl=86400)
def get_scientific_name(common_name: str):
    """Returns the single top scientific name for a common name. Kept for compatibility."""
    options = search_species_options(common_name, max_results=1)
    return options[0][1] if options else None


def _fetch_single_marker(species_name: str, gene: str):
    """
    Fetches up to 5 NCBI candidates for (species, gene) and returns the first
    sequence within MARKER_LENGTH_RANGE. The SLEN filter pre-filters on NCBI's
    side; the length check here catches anything that slips through.

    NCBI field tags used:
      [Organism] — restricts to species
      [Gene]     — restricts to gene name
      [SLEN]     — restricts sequence length (min:max)
    """
    min_len, max_len = MARKER_LENGTH_RANGE.get(gene, (300, 5000))
    search_term = (
        f"{species_name}[Organism] AND {gene}[Gene] "
        f"AND {min_len}[SLEN]:{max_len}[SLEN]"
    )

    def _fetch():
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            return None

        for seq_id in record["IdList"]:
            handle     = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()
            seq_str = str(seq_record.seq)

            if min_len <= len(seq_str) <= max_len:
                print(f"[NCBI] {species_name} {gene}: accepted {len(seq_str)}bp")
                return seq_str
            print(f"[NCBI] {species_name} {gene}: skipped {len(seq_str)}bp")

        return None

    return _fetch_with_retry(_fetch)


@st.cache_data(ttl=86400)
def get_marker_sequence(species_name: str, gene: str):
    """
    Cached wrapper around _fetch_single_marker.
    Caching at this level means individual (species, gene) pairs are reused
    across different comparisons — e.g. Human+COI is fetched once and shared
    between "Human vs Chimp" and "Human vs Gorilla".
    """
    return _fetch_single_marker(species_name, gene)


@st.cache_data(ttl=86400)
def fetch_both_sequences(species_A: str, species_B: str,
                         type_A="default", type_B="default"):
    """
    Iterates the marker chain and returns sequences only when BOTH species
    have a result for the SAME gene. Stops at the first shared marker found.

    Returns ((seq_A, marker), (seq_B, marker)) or ((None, None), (None, None)).

    Note: calls are sequential, not threaded — @st.cache_data requires the
    main Streamlit thread context and does not work inside ThreadPoolExecutor.
    """
    chain_key = type_A if type_A != "default" else type_B
    markers   = MARKER_CHAINS.get(chain_key, MARKER_CHAINS["default"])

    for gene in markers:
        seq_A = get_marker_sequence(species_A, gene)
        seq_B = get_marker_sequence(species_B, gene)

        if seq_A and seq_B:
            print(f"[NCBI] Both species matched on marker: {gene}")
            return (seq_A, gene), (seq_B, gene)
        if not seq_A:
            print(f"[NCBI] {species_A} has no {gene}, trying next...")
        if not seq_B:
            print(f"[NCBI] {species_B} has no {gene}, trying next...")

    return (None, None), (None, None)


def calculate_similarity_and_alignment(seq1: str, seq2: str, max_len: int = 1000):
    """
    Trims both sequences to min(len_A, len_B, max_len) so the comparison
    denominator is fair — neither sequence is penalised for extra length
    the other doesn't have. Then runs both algorithms.

    Returns a dict: hirschberg_pct, hirschberg_match, sw_pct, sw_match,
                    sw_score, compare_len.
    """
    if not seq1 or not seq2:
        return {
            "hirschberg_pct": 0.0, "hirschberg_match": "",
            "sw_pct": 0.0, "sw_match": "", "sw_score": 0, "compare_len": 0,
        }

    compare_len = min(len(seq1), len(seq2), max_len)
    s1 = seq1[:compare_len]
    s2 = seq2[:compare_len]

    matched_global = hirschberg_reconstruct(s1, s2)
    hirschberg_pct = (len(matched_global) / compare_len) * 100

    matched_local, sw_score, sw_pct = smith_waterman(s1, s2)

    return {
        "hirschberg_pct":   hirschberg_pct,
        "hirschberg_match": matched_global,
        "sw_pct":           sw_pct,
        "sw_match":         matched_local,
        "sw_score":         sw_score,
        "compare_len":      compare_len,
    }
