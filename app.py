#!/usr/bin/env python3
"""
CRISPR gRNA Designer – Streamlit Dashboard
================================================
User enters gene name + organism → tool fetches sequence from NCBI or Ensembl
→ designs ranked SpCas9 gRNAs for Knockout or Knockdown.
"""

import streamlit as st
import pandas as pd
from io import BytesIO
import plotly.express as px
from Bio.Seq import Seq

from grna_designer import (
    design_from_gene,
    design_guides,
    validate_guide,
    fetch_sequence,
)

st.set_page_config(
    page_title="CRISPR gRNA Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
st.sidebar.title("🧬 CRISPR gRNA Designer")
st.sidebar.markdown(
    """
**How it works**
1. Enter gene symbol + organism  
2. Choose NCBI or Ensembl  
3. Select Knockout or Knockdown  
4. Get ranked SpCas9 (NGG) guides  
"""
)

source = st.sidebar.selectbox("Sequence source", ["ncbi", "ensembl"], index=0)
application = st.sidebar.radio(
    "Application",
    ["knockout", "knockdown"],
    format_func=lambda x: "Knockout (NHEJ / frameshift)" if x == "knockout" else "Knockdown (CRISPRi-style, prefer 5′)",
)
max_guides = st.sidebar.slider("Max guides to show", 5, 50, 15)
min_score = st.sidebar.slider("Minimum score filter", 0, 80, 35)

st.sidebar.markdown("---")
st.sidebar.markdown(
    """
**Design rules implemented**
- PAM: NGG (SpCas9)
- Spacer: 20 nt
- GC 40–70 % preferred
- Avoid poly-G / poly-T
- Prefer G at pos 20
- Early-region bias for KO
- 5′-proximal bias for KD
"""
)

st.sidebar.info(
    "Off-target scanning requires a full genome index "
    "(Bowtie / Cas-OFFinder). This tool ranks on-target features. "
    "Always validate top guides with CRISPOR or CHOPCHOP for off-targets."
)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
st.title("CRISPR gRNA Designer")
st.markdown(
    "Enter a **gene name** and **organism**. The tool fetches the sequence "
    "from public databases and designs SpCas9 guide RNAs."
)

col1, col2, col3 = st.columns([2, 2, 1])
with col1:
    gene_name = st.text_input("Gene name / symbol", value="TP53", help="e.g. TP53, BRCA1, GFP, ACTB")
with col2:
    organism = st.text_input("Organism", value="Homo sapiens", help="e.g. Homo sapiens, Mus musculus, Danio rerio")
with col3:
    st.write("")  # spacer
    st.write("")
    run_btn = st.button("Design gRNAs", type="primary", use_container_width=True)

# Manual sequence option
with st.expander("Or paste your own DNA sequence"):
    custom_seq = st.text_area("DNA sequence (ATGC only)", height=120, placeholder="Paste sequence here…")
    use_custom = st.checkbox("Use pasted sequence instead of database fetch")

if run_btn:
    with st.spinner("Fetching sequence and designing guides…"):
        try:
            if use_custom and custom_seq.strip():
                seq = "".join(c for c in custom_seq.upper() if c in "ATGCU").replace("U", "T")
                if len(seq) < 50:
                    st.error("Sequence too short (need ≥ 50 bp).")
                    st.stop()
                if len(seq) > 20000:
                    st.error("Sequence too long. Please limit custom sequences to 20,000 bp to prevent timeouts.")
                    st.stop()
                acc, desc = "custom", "User-provided sequence"
                guides = design_guides(
                    seq,
                    application=application,
                    max_guides=max_guides,
                    min_score=min_score,
                )
            else:
                if not gene_name.strip() or not organism.strip():
                    st.error("Please provide both gene name and organism.")
                    st.stop()
                acc, desc, seq, guides = design_from_gene(
                    gene_name.strip(),
                    organism.strip(),
                    application=application,
                    source=source,
                    max_guides=max_guides,
                )
                guides = [g for g in guides if g.score >= min_score]

            st.session_state["run_results"] = {
                "guides": guides,
                "seq": seq,
                "acc": acc,
                "desc": desc,
                "gene_name": gene_name if not use_custom else "Custom",
            }

        except Exception as e:
            st.error(f"Error: {e}")
            st.info(
                "Tips:\n"
                "- Check gene symbol spelling\n"
                "- Try the other database\n"
                "- For non-model organisms paste the sequence manually\n"
                "- NCBI may rate-limit; wait a few seconds and retry"
            )

# Check if we have results stored in session state to render
if "run_results" in st.session_state:
    res = st.session_state["run_results"]
    guides = res["guides"]
    seq = res["seq"]
    acc = res["acc"]
    desc = res["desc"]
    g_name = res["gene_name"]

    # ---------- Results ----------
    st.success(f"Found **{len(guides)}** candidate guides")
    st.markdown(f"**Accession / ID:** `{acc}`  \n**Description:** {desc[:200]}…")
    st.markdown(f"**Sequence length:** {len(seq):,} bp")

    with st.expander("Sequence preview (first 300 bp)"):
        st.code(seq[:300] + ("…" if len(seq) > 300 else ""), language="text")

    if not guides:
        st.warning("No guides passed the score filter. Try lowering the minimum score.")
        st.stop()

    # Table
    df = pd.DataFrame([g.to_dict() for g in guides])
    st.subheader("Ranked Guide RNAs")
    st.dataframe(df, use_container_width=True, hide_index=True)

    fig = px.histogram(df, x="Score", nbins=15, title="On-target score distribution")
    st.plotly_chart(fig, use_container_width=True)

    # Validation of top 3
    st.subheader("Validation of top 3 guides")
    for i, g in enumerate(guides[:3]):
        checks = validate_guide(g)
        status = "✅ PASS" if checks["overall_pass"] else "⚠️ REVIEW"
        with st.expander(f"#{i+1} {g.sequence}  {status}"):
            st.json(checks)
            st.markdown(
                f"**Full target (spacer+PAM):** `{g.sequence}{g.pam}`  \n"
                f"**Strand:** {g.strand}  |  **Position:** {g.start+1}-{g.end}  \n"
                f"**GC%:** {g.gc_content:.1f}  |  **Score:** {g.score:.1f}"
            )
            if g.notes:
                st.warning("Notes: " + "; ".join(g.notes))

    # Download
    st.subheader("Download results")
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button("Download CSV", csv, file_name=f"{g_name}_gRNAs.csv", mime="text/csv")

    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="gRNAs")
        meta = pd.DataFrame(
            {
                "Field": ["Gene", "Organism", "Source", "Application", "Accession", "Seq length"],
                "Value": [g_name, organism, source, application, acc, len(seq)],
            }
        )
        meta.to_excel(writer, index=False, sheet_name="Metadata")
    
    st.download_button(
        "Download Excel",
        buffer.getvalue(),
        file_name=f"{g_name}_gRNAs.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )

# Footer
st.markdown("---")
st.caption(
    "This tool implements literature-derived on-target heuristics (Doench et al., Hsu et al., "
    "Moreno-Mateos et al.). It is intended for research / educational use. "
    "Always confirm off-target profiles with genome-wide tools (CRISPOR, Cas-OFFinder, GuideScan) "
    "before experimental use."
)