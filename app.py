#!/usr/bin/env python3
"""CRISPR Studio v3.1 polished Streamlit interface."""
from __future__ import annotations

import json
import os
import re
from statistics import median
from typing import List

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from Bio.Seq import Seq

from crispri import (
    CRISPRiGuide,
    CRISPRI_MAX,
    CRISPRI_MIN,
    CRISPRI_OPTIMAL_MAX,
    CRISPRI_OPTIMAL_MIN,
    design_crispri_from_sequence,
    design_crispri_guides,
    fetch_ensembl_tss_context,
)
from grna_designer import (
    GuideRNA,
    analyze_offtargets,
    clean_dna_sequence,
    design_guides,
    fetch_sequence,
    guidescan2_available,
    parse_reference,
    run_guidescan2,
    score_breakdown,
    validate_guide,
)

APP_VERSION = "3.1.0"
MAX_CUSTOM_BP = 100_000
MAX_LOCAL_REFERENCE_BP = 5_000_000

st.set_page_config(
    page_title="CRISPR Studio | gRNA Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer",
        "Report a bug": "https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer/issues",
        "About": "CRISPR Studio v3.1 — SpCas9 guide design with TSS-aware CRISPRi.",
    },
)

# ----------------------------- Theme -------------------------------------
dark_mode = st.sidebar.toggle("Dark mode", value=True, key="dark_mode")
P = (
    dict(app="#07110f", side="#0b1714", panel="#10211d", alt="#142a25", text="#eefbf6", muted="#9ebbb0", border="#24453b", accent="#34d399", accent2="#22d3ee", shadow="rgba(0,0,0,.28)")
    if dark_mode
    else dict(app="#f5faf8", side="#edf7f3", panel="#ffffff", alt="#f2f8f5", text="#12201b", muted="#5f756d", border="#d5e7df", accent="#059669", accent2="#0891b2", shadow="rgba(15,60,45,.10)")
)
st.markdown(
    f"""
<style>
:root{{--app:{P['app']};--side:{P['side']};--panel:{P['panel']};--alt:{P['alt']};--text:{P['text']};--muted:{P['muted']};--border:{P['border']};--accent:{P['accent']};--accent2:{P['accent2']};}}
.stApp{{background:var(--app);color:var(--text)}} [data-testid="stSidebar"]{{background:var(--side);border-right:1px solid var(--border)}}
.block-container{{max-width:1460px;padding-top:1.5rem;padding-bottom:4rem}} .stApp p,.stApp li,.stApp label{{color:var(--text)}}
.hero{{padding:2.2rem 2.35rem;border:1px solid var(--border);border-radius:24px;background:radial-gradient(circle at 92% 12%,rgba(34,211,238,.17),transparent 26%),radial-gradient(circle at 8% 100%,rgba(52,211,153,.17),transparent 30%),var(--panel);box-shadow:0 18px 55px {P['shadow']};margin-bottom:1.3rem}}
.hero .k{{color:var(--accent);font-size:.76rem;font-weight:800;letter-spacing:.14em;text-transform:uppercase}} .hero h1{{color:var(--text);font-size:clamp(2rem,4.6vw,4.1rem);letter-spacing:-.045em;line-height:.98;margin:.6rem 0 1rem}} .hero h1 span{{color:var(--accent)}} .hero .copy{{color:var(--muted);max-width:840px;font-size:1.05rem}}
.pills{{display:flex;flex-wrap:wrap;gap:.55rem;margin-top:1.25rem}} .pills span{{background:var(--alt);border:1px solid var(--border);border-radius:999px;padding:.42rem .72rem;font-size:.78rem;font-weight:650}}
div[data-testid="stMetric"]{{background:var(--panel);border:1px solid var(--border);padding:1rem;border-radius:16px;box-shadow:0 8px 24px {P['shadow']}}}
[data-testid="stForm"],[data-testid="stExpander"]{{background:var(--panel);border-color:var(--border)!important;border-radius:16px}} [data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:14px;overflow:hidden}}
.seq{{background:var(--alt);border:1px solid var(--border);padding:1rem 1.15rem;border-radius:14px;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;overflow-wrap:anywhere;letter-spacing:.035em}} .seq mark{{background:var(--accent);color:#04120d;padding:.1rem .25rem;border-radius:5px;font-weight:800}}
.method{{height:100%;background:var(--panel);border:1px solid var(--border);border-radius:15px;padding:1rem}} .method b{{color:var(--accent)}}
.stButton>button,.stDownloadButton>button{{border-radius:12px;min-height:2.8rem;font-weight:750;border-color:var(--border)}} .stButton>button[kind="primary"]{{background:linear-gradient(115deg,var(--accent),var(--accent2));color:#04120d;border:0}}
</style>
""",
    unsafe_allow_html=True,
)
plot_template = "plotly_dark" if dark_mode else "plotly_white"

# ----------------------------- Sidebar -----------------------------------
st.sidebar.markdown(f"### CRISPR Studio\n**SpCas9 workbench · v{APP_VERSION}**")
st.sidebar.subheader("Design settings")
max_guides = st.sidebar.slider("Maximum guides", 5, 50, 20, 5)
min_score = st.sidebar.slider("Minimum heuristic score", 0, 90, 35, 5)
max_mismatches = st.sidebar.slider("Off-target mismatches", 0, 4, 3)
with st.sidebar.expander("Models & backends"):
    st.markdown(
        "- **SpCas9:** 20 nt + NGG\n- **MIT/Hsu:** native\n- **Doench RS2 / CFD:** optional GuideMaker provider\n- **Whole genome:** GuideScan2 index\n- **CRISPRi:** Ensembl canonical TSS"
    )
st.sidebar.info("Research-use shortlist only. Validate genome build, target biology, off-targets and experimental controls before ordering guides.")

# ----------------------------- Helpers -----------------------------------
@st.cache_data(ttl=3600, show_spinner=False)
def cached_fetch(gene, organism, source):
    return fetch_sequence(gene, organism, source)

@st.cache_data(ttl=3600, show_spinner=False)
def cached_tss(gene, organism):
    return fetch_ensembl_tss_context(gene, organism)

def read_upload(f):
    if f is None:
        return ""
    b = f.getvalue()
    try:
        return b.decode("utf-8-sig")
    except UnicodeDecodeError:
        return b.decode("latin-1")

def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "crispr_guides"

def style_plot(fig, height=420):
    fig.update_layout(template=plot_template, height=height, margin=dict(l=20, r=20, t=55, b=20), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)")
    return fig

# ----------------------------- Hero + form -------------------------------
st.markdown(
    """
<section class="hero"><div class="k">Open-source CRISPR design workbench</div><h1>Design sharper <span>CRISPR guides.</span></h1><div class="copy">Rank SpCas9 guides, inspect Doench/MIT/CFD metrics, screen local references or indexed genomes, and design CRISPRi guides against a real annotated transcription start site.</div><div class="pills"><span>20 nt + NGG</span><span>Both strands</span><span>TSS-aware CRISPRi</span><span>MIT/Hsu</span><span>Doench RS2</span><span>GuideScan2</span><span>CSV · FASTA · JSON</span></div></section>
""",
    unsafe_allow_html=True,
)
st.subheader("Design workspace")
input_mode = st.radio("Target input", ["Gene lookup", "Paste sequence"], horizontal=True, label_visibility="collapsed")

with st.form("design_form", border=True):
    left, right = st.columns([1.25, 1])
    with left:
        if input_mode == "Gene lookup":
            gene_name = st.text_input("Gene symbol", value="TP53")
            organism = st.text_input("Organism", value="Homo sapiens")
            source = st.radio("Sequence database", ["Ensembl", "NCBI"], horizontal=True).lower()
            custom_sequence = ""
        else:
            custom_sequence = st.text_area("Target DNA / genomic FASTA", height=180, placeholder=">target\nATG...")
            gene_name, organism, source = "Custom target", "Not specified", "manual"
    with right:
        intent = st.radio("Design intent", ["Knockout", "CRISPRi repression"])
        application = "knockout" if intent == "Knockout" else "crispri"
        custom_tss = None
        if input_mode == "Paste sequence" and application == "crispri":
            custom_tss = st.number_input("TSS position in pasted sequence (1-based)", min_value=1, value=1, step=1, help="The pasted sequence must be genomic DNA in 5′→3′ transcriptional orientation. This avoids guessing a TSS.")
        screen_mode = st.selectbox("Specificity analysis", ["None", "Local reference", "Whole genome (GuideScan2)"])
        ref_upload = None; ref_text = ""; genome_index = ""
        if screen_mode == "Local reference":
            ref_upload = st.file_uploader("Reference FASTA / text", type=["fa", "fasta", "fna", "txt"])
            ref_text = st.text_area("Or paste reference sequence", height=80)
        elif screen_mode.startswith("Whole genome"):
            genome_index = st.text_input("GuideScan2 index path", value=os.getenv("GUIDESCAN_INDEX", ""))
            st.caption("GuideScan2 detected on this host: " + ("yes" if guidescan2_available() else "no"))
    submitted = st.form_submit_button("Design and analyze guides", type="primary", use_container_width=True)

if submitted:
    try:
        reference = None; raw_reference = ""
        if screen_mode == "Local reference":
            raw_reference = read_upload(ref_upload) or ref_text
            reference = parse_reference(raw_reference)
            total = sum(map(len, reference.values()))
            if not total: raise ValueError("Local reference mode is enabled but no reference was supplied.")
            if total > MAX_LOCAL_REFERENCE_BP: raise ValueError(f"Local reference exceeds {MAX_LOCAL_REFERENCE_BP:,} bp hosted limit.")

        crispri_rows: List[CRISPRiGuide] = []
        tss_context = None
        if application == "crispri":
            if input_mode == "Gene lookup":
                if not gene_name.strip() or not organism.strip(): raise ValueError("Enter both gene symbol and organism.")
                tss_context = cached_tss(gene_name.strip(), organism.strip())
                sequence = tss_context.sequence
                accession = tss_context.transcript_id
                description = tss_context.description
                crispri_rows = design_crispri_guides(tss_context, max_guides=max_guides, min_score=float(min_score), genome_context=reference, max_mismatches=max_mismatches)
            else:
                sequence = clean_dna_sequence(custom_sequence)
                if len(sequence) < 50: raise ValueError("The target sequence must contain at least 50 bp.")
                if len(sequence) > MAX_CUSTOM_BP: raise ValueError(f"Target exceeds {MAX_CUSTOM_BP:,} bp limit.")
                accession, description = "CUSTOM_TSS", "User-provided genomic sequence with declared TSS"
                crispri_rows = design_crispri_from_sequence(sequence, int(custom_tss), max_guides=max_guides, min_score=float(min_score), genome_context=reference, max_mismatches=max_mismatches)
            guides = [x.guide for x in crispri_rows]
        else:
            if input_mode == "Paste sequence":
                sequence = clean_dna_sequence(custom_sequence)
                if len(sequence) < 50: raise ValueError("The target sequence must contain at least 50 bp.")
                if len(sequence) > MAX_CUSTOM_BP: raise ValueError(f"Target exceeds {MAX_CUSTOM_BP:,} bp limit.")
                accession, description = "CUSTOM", "User-provided target sequence"
            else:
                if not gene_name.strip() or not organism.strip(): raise ValueError("Enter both gene symbol and organism.")
                accession, description, sequence = cached_fetch(gene_name.strip(), organism.strip(), source)
            guides = design_guides(sequence, application="knockout", max_guides=max_guides, min_score=float(min_score), genome_context=reference, max_mismatches=max_mismatches)

        genome_rows = []
        if screen_mode.startswith("Whole genome") and guides:
            genome_rows = run_guidescan2(guides, genome_index, max_mismatches=max_mismatches)

        st.session_state["analysis"] = dict(guides=guides, crispri_rows=crispri_rows, tss_context=tss_context, sequence=sequence, reference=reference, genome_rows=genome_rows, gene=gene_name.strip() or "Custom target", organism=organism.strip(), source=source, application=application, accession=accession, description=description, screen_mode=screen_mode, max_mismatches=max_mismatches)
        st.success(f"Analysis complete: {len(guides)} guide candidates passed the filter.")
    except Exception as exc:
        st.error(f"Could not complete the analysis: {exc}")

# ----------------------------- Results -----------------------------------
if "analysis" in st.session_state:
    r = st.session_state["analysis"]; guides: List[GuideRNA] = r["guides"]; sequence = r["sequence"]; reference = r["reference"]; crispri_rows: List[CRISPRiGuide] = r.get("crispri_rows", [])
    st.divider(); st.subheader(f"Analysis report · {r['gene']}")
    st.caption(f"{r['description'][:260]} | Accession: {r['accession']} | Analyzed DNA: {len(sequence):,} bp")
    if r["application"] == "crispri" and r.get("tss_context"):
        c = r["tss_context"]
        st.success(f"TSS-aware CRISPRi: {c.assembly} {c.chromosome}:{c.tss_coordinate:,} ({c.strand_label}) · canonical transcript {c.transcript_id} · accepted window {CRISPRI_MIN:+d} to {CRISPRI_MAX:+d} bp; preferred {CRISPRI_OPTIMAL_MIN:+d} to {CRISPRI_OPTIMAL_MAX:+d} bp.")
    elif r["source"] != "manual" and r["application"] == "knockout":
        st.warning("Knockout gene lookup uses representative transcript/cDNA for candidate discovery. Confirm genomic exon/assembly coordinates before experimental use.")

    if not guides:
        st.warning("No candidate guides passed the current filters. For CRISPRi, the annotated TSS window may simply contain no NGG site at the current threshold.")
    else:
        if crispri_rows:
            df = pd.DataFrame([x.to_dict() for x in crispri_rows])
        else:
            df = pd.DataFrame([g.to_dict() for g in guides])
        df.insert(0, "Rank", range(1, len(df) + 1))

        metrics = st.columns(5)
        metrics[0].metric("Guides retained", len(guides))
        metrics[1].metric("Best heuristic", f"{max(g.score for g in guides):.1f}")
        metrics[2].metric("Median GC", f"{median(g.gc_content for g in guides):.1f}%")
        ds = [g.doench_score for g in guides if g.doench_score is not None]
        metrics[3].metric("Best Doench RS2", f"{max(ds):.1f}" if ds else "Provider N/A")
        ms = [g.specificity_score for g in guides if g.specificity_score is not None]
        metrics[4].metric("Best MIT specificity", f"{max(ms):.1f}" if ms else ("GuideScan2 run" if r["genome_rows"] else "Not run"))

        tabs = st.tabs(["Ranked guides", "Design landscape", "Guide details", "Off-target screen", "Export", "Methods & limits"])
        with tabs[0]:
            st.markdown("#### Ranked candidate table")
            st.dataframe(df, use_container_width=True, hide_index=True, height=min(650, 100 + len(df) * 35))
        with tabs[1]:
            if crispri_rows:
                plot_df = df.copy()
                fig = px.scatter(plot_df, x="TSS distance (bp)", y="Heuristic", color="TSS band", size="GC%", hover_data=["Rank", "Spacer (20 nt)", "PAM", "Chromosome", "Genomic start"], title="CRISPRi candidates relative to the annotated TSS")
                fig.add_vrect(x0=50, x1=100, opacity=.10, line_width=0, annotation_text="preferred +50..+100")
                fig.add_vline(x=0, line_dash="dash", annotation_text="TSS")
                st.plotly_chart(style_plot(fig, 460), use_container_width=True)
            else:
                fig = px.scatter(df, x="Start", y="Heuristic", color="Strand", size="GC%", hover_data=["Rank", "Spacer (20 nt)", "PAM"], title="Candidate activity across the target")
                st.plotly_chart(style_plot(fig, 460), use_container_width=True)
            hist = px.histogram(df, x="GC%", nbins=12, title="GC-content distribution")
            st.plotly_chart(style_plot(hist, 330), use_container_width=True)
        with tabs[2]:
            idx = st.selectbox("Inspect a guide", range(len(guides)), format_func=lambda i: f"#{i+1} · {guides[i].sequence}")
            g = guides[idx]
            st.markdown(f'<div class="seq">5′—{g.sequence}<mark>{g.pam}</mark>—3′</div>', unsafe_allow_html=True)
            q = st.columns(5); q[0].metric("Heuristic", f"{g.score:.1f}"); q[1].metric("Doench RS2", f"{g.doench_score:.1f}" if g.doench_score is not None else "N/A"); q[2].metric("GC", f"{g.gc_content:.1f}%"); q[3].metric("MIT", f"{g.specificity_score:.1f}" if g.specificity_score is not None else "N/A"); q[4].metric("CFD spec.", f"{g.cfd_specificity:.1f}" if g.cfd_specificity is not None else "N/A")
            if crispri_rows:
                cg = crispri_rows[idx]
                st.info(f"TSS placement: {cg.tss_distance:+d} bp · {cg.tss_band}" + (f" · {cg.assembly} {cg.chromosome}:{cg.genomic_start}-{cg.genomic_end} ({cg.genomic_strand})" if cg.genomic_start else ""))
            checks = validate_guide(g); st.markdown("##### Validation checklist"); st.json(checks)
            with st.expander("Legacy heuristic breakdown"):
                st.json(score_breakdown(g.sequence, g.pam, application=None if g.application == "crispri" else g.application, start=g.start, sequence_length=len(sequence)))
            with st.expander("Example BbsI cloning oligos"):
                st.code(f"Forward: CACCG{g.sequence}\nReverse: AAAC{str(Seq(g.sequence).reverse_complement())}C")
        with tabs[3]:
            if reference:
                oi = st.selectbox("Inspect local off-targets", range(len(guides)), key="otguide", format_func=lambda i: f"#{i+1} · {guides[i].sequence}")
                rep = analyze_offtargets(guides[oi].sequence, reference, max_mismatches=r["max_mismatches"])
                c = st.columns(5); c[0].metric("MIT specificity", f"{rep.mit_specificity:.1f}"); c[1].metric("CFD specificity", f"{rep.cfd_specificity:.1f}" if rep.cfd_specificity is not None else "Provider N/A"); c[2].metric("Near-matches", len(rep.hits)); c[3].metric("PAM sites", rep.pam_sites_scanned); c[4].metric("Contigs", rep.reference_contigs)
                if rep.hits: st.dataframe(pd.DataFrame([h.to_dict() for h in rep.hits]), use_container_width=True, hide_index=True)
                else: st.success("No additional PAM-compatible near-matches within the selected mismatch radius.")
            elif r["genome_rows"]:
                st.markdown("#### GuideScan2 whole-genome results"); st.dataframe(pd.DataFrame(r["genome_rows"]), use_container_width=True, hide_index=True)
            else:
                st.info("No specificity analysis was run. Choose Local reference or Whole genome (GuideScan2).")
        with tabs[4]:
            base = safe_name(r["gene"]); csvb = df.to_csv(index=False).encode(); fasta = "\n".join(f">{base}_g{i}|pam={g.pam}|strand={g.strand}\n{g.sequence}" for i, g in enumerate(guides, 1)).encode(); js = json.dumps({"metadata": {"target": r["gene"], "organism": r["organism"], "application": r["application"], "accession": r["accession"], "version": APP_VERSION}, "guides": df.where(pd.notna(df), None).to_dict("records")}, indent=2).encode()
            d = st.columns(3); d[0].download_button("Download CSV", csvb, file_name=f"{base}_v31_guides.csv", use_container_width=True); d[1].download_button("Download FASTA", fasta, file_name=f"{base}_v31_spacers.fasta", use_container_width=True); d[2].download_button("Download JSON", js, file_name=f"{base}_v31_analysis.json", use_container_width=True)
        with tabs[5]:
            cards = st.columns(4)
            for col, (title, text) in zip(cards, [("Discovery", "Both-strand SpCas9 20 nt + NGG scanning."), ("On-target", "Doench Rule Set 2 when the optional provider is installed; transparent heuristic remains visible."), ("Specificity", "Native MIT/Hsu; optional CFD; local multi-contig or GuideScan2 whole-genome search."), ("CRISPRi", "Canonical Ensembl TSS, strand-aware genomic DNA, −50..+300 bp filter, +50..+100 preferred band.")]):
                col.markdown(f'<div class="method"><b>{title}</b><p>{text}</p></div>', unsafe_allow_html=True)
            st.markdown("#### Scientific boundaries")
            st.markdown("- **CRISPRi is now TSS-aware for gene lookup.** It is no longer ranked by the first 400 bp of cDNA.\n- Pasted-sequence CRISPRi requires an explicit TSS position; the app does not guess one.\n- Knockout gene lookup still begins from representative cDNA, so shortlisted guides should be mapped to the intended genomic coding exon/assembly.\n- Whole-genome specificity requires a GuideScan2 executable plus a matching prebuilt genome index.\n- Variant-aware filtering, chromatin state and DNA/RNA bulges are not inferred unless supplied by an external genome-aware workflow.\n- Computational scores prioritize candidates; they do not replace experimental validation.")

st.divider(); st.caption("CRISPR Studio v3.1 · research/educational use · MIT License")
