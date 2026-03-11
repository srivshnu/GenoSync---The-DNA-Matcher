# =============================================================================
# app.py — Streamlit UI for the Evolutionary DNA Matcher
#
# Handles all user interaction — inputs, dropdowns, progress, results display.
# Contains no algorithm logic; delegates everything to core.py.
#
# Streamlit re-runs this entire file on every user interaction.
# Widget state is preserved between reruns via the key= parameter.
# Run with: streamlit run app.py
# =============================================================================

import streamlit as st
from core import (
    search_species_options,
    fetch_both_sequences,
    calculate_similarity_and_alignment,
)

# Must be the first Streamlit call in the file
st.set_page_config(page_title="DNA Matcher", page_icon="🧬", layout="wide")

st.title("🧬 Evolutionary DNA Matcher")
st.markdown(
    "Compares two species using real DNA sequences from NCBI. "
    "Runs **Hirschberg** (global) and **Smith-Waterman** (local) alignment."
)
st.markdown("---")


def get_species_flow(label, key):
    """
    Renders a species search box for one side of the comparison (A or B).
    Resolves a typed common name to a scientific name via NCBI Taxonomy.
    Shows a disambiguation dropdown when multiple species match.

    key must be unique per instance ("A" / "B") — without it, both panels
    share the same Streamlit widget state and overwrite each other.
    """
    common   = st.text_input(label, placeholder="e.g. Dog, Hen, Tiger", key=f"input_{key}")
    org_type = st.selectbox(
        "Organism type",
        ["default", "animal", "plant", "fungus"],
        key=f"type_{key}",
        help="Helps pick the right genetic marker. Use 'default' if unsure.",
    )

    sci_name = None

    if common:
        with st.spinner(f"Searching for '{common}'..."):
            options = search_species_options(common)  # [(common_name, sci_name), ...]

        if not options:
            st.error("No species found. Try a different name.")

        elif len(options) == 1:
            common_label, sci_name = options[0]
            st.info(f"Found: **{common_label}** (*{sci_name}*)")

        else:
            # Show "Common Name — Scientific name" labels; return only the sci_name
            labels    = [f"{c}  —  {s}" for c, s in options]
            sci_names = [s for _, s in options]
            idx = st.selectbox(
                f"Multiple matches for '{common}' — pick one:",
                range(len(labels)),
                format_func=lambda i: labels[i],
                key=f"select_{key}",
            )
            sci_name = sci_names[idx]
            st.info(f"Selected: **{options[idx][0]}** (*{sci_name}*)")

    return sci_name, org_type


# Two-column layout — Species A left, Species B right
col1, col2 = st.columns(2)
with col1:
    st.subheader("Species A")
    species_A, type_A = get_species_flow("Common Name", "A")
with col2:
    st.subheader("Species B")
    species_B, type_B = get_species_flow("Common Name", "B")

st.markdown("---")

# Collapsed by default to keep the UI clean
with st.expander("⚙️ Advanced Settings", expanded=False):
    max_len = st.slider(
        "Sequence length to compare (base pairs)",
        min_value=200, max_value=2000, value=1000, step=100,
        help="Higher = more accurate but slower. COI markers are ~650bp.",
    )

st.markdown("---")

if st.button("🧬 Compare DNA Sequences", use_container_width=True):

    if not species_A or not species_B:
        st.warning("Please identify both species before running.")

    else:
        progress = st.progress(0)
        status   = st.empty()

        # Fetch DNA — only accepts a gene when BOTH species have it,
        # preventing meaningless cross-gene comparisons (e.g. COI vs cytb)
        status.text("Fetching DNA markers from NCBI...")
        (dna1, marker_A), (dna2, marker_B) = fetch_both_sequences(
            species_A, species_B, type_A, type_B
        )
        progress.progress(50)

        if not dna1 or not dna2:
            if not dna1:
                st.error(
                    f"Could not retrieve a DNA marker for **{species_A}**. "
                    "Try changing the organism type or using the scientific name directly."
                )
            if not dna2:
                st.error(
                    f"Could not retrieve a DNA marker for **{species_B}**. "
                    "Try changing the organism type or using the scientific name directly."
                )

        else:
            status.text(f"Running Hirschberg + Smith-Waterman on up to {max_len}bp...")
            results = calculate_similarity_and_alignment(dna1, dna2, max_len=max_len)
            progress.progress(100)
            status.text("Done ✅")

            st.balloons()

            # Show which gene was used and how many base pairs were compared
            st.markdown("#### Genetic Markers Used")
            m1, m2, m3 = st.columns(3)
            m1.success(f"**{species_A}** → `{marker_A}`")
            m2.success(f"**{species_B}** → `{marker_B}`")
            m3.info(f"Compared **{results['compare_len']} bp**")

            st.markdown("---")

            # Side-by-side results: global score drops for distant species,
            # local score stays higher by ignoring divergent regions
            st.markdown("#### Similarity Results")
            r1, r2 = st.columns(2)

            with r1:
                st.markdown("##### 🌍 Global Alignment (Hirschberg)")
                st.metric(
                    "Global Similarity",
                    f"{results['hirschberg_pct']:.2f}%",
                    help="End-to-end similarity across the full compared region.",
                )
                with st.expander("View Matched Sequence"):
                    st.code(results["hirschberg_match"] or "(no match)", language="text")

            with r2:
                st.markdown("##### 🔍 Local Alignment (Smith-Waterman)")
                st.metric(
                    "Local Similarity",
                    f"{results['sw_pct']:.2f}%",
                    help="Similarity of the best-matching conserved region.",
                )
                st.caption(f"Raw alignment score: **{results['sw_score']}**")
                with st.expander("View Best Local Region"):
                    st.code(results["sw_match"] or "(no match)", language="text")

            st.markdown("---")
            st.markdown("#### 📖 How to Read These Numbers")
            st.markdown("""
| Score Range | Likely Relationship |
|---|---|
| 95 – 100% | Same or nearly identical subspecies |
| 80 – 94% | Sister species (e.g. human vs chimpanzee) |
| 65 – 79% | Same family (e.g. lion vs tiger, orca vs blue whale) |
| 45 – 64% | Same order (e.g. whale vs dolphin vs cow) |
| Below 45% | Distant relatives (different class or phylum) |

> **Global** score drops for distant species because mismatches accumulate across the whole region.  
> **Local** score stays higher — it finds the best-preserved conserved region and ignores the rest.  
> ⚠️ LCS scores run ~10–15% lower than BLAST % identity. A 72% here ≈ ~85% in standard bioinformatics tools.
            """)