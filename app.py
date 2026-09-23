#!/usr/bin/env python3

# ============================================================================
# APP
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Builds the Streamlit user interface and connects the input, design, validation and reporting steps.
#
# HOW TO READ THIS FILE:
# 1. Start with the imports and constants.
# 2. Read each function separately; every function performs one part of the workflow.
# 3. Follow the function calls from the application/workflow rather than trying to
#    understand the entire file at once.
# 4. Scientific equations, thresholds, validation rules and public function names
#    are intentionally kept unchanged while the code is being humanized.
#
# MAIN TOP-LEVEL PARTS IN THIS FILE:
# - function: cached_fetch
# - function: cached_tss
# - function: cached_knockout
# - function: cached_accession
# - function: read_upload
# - function: safe_name
# - function: style_plot
# - function: render_validation
# ============================================================================

"""CRISPR Studio v3.3 professional Streamlit dashboard."""
from __future__ import annotations

import json
import os
import re
from html import escape
from statistics import median
from typing import List

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio.Seq import Seq

from knockout import fetch_ncbi_knockout_context, design_knockout_guides
from models import rs2_status
from reporting import make_metadata, export_json
from accession_lookup import AccessionRecord, fetch_accession
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
    screen_guides,
    ReferenceIndex,
)

APP_VERSION = "3.3.0"
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
        "About": "CRISPR Studio v3.3 — SpCas9 design, accession lookup and TSS-aware CRISPRi.",
    },
)

# ----------------------------- Theme -------------------------------------
dark_mode = st.sidebar.toggle("Dark mode", value=True, key="dark_mode")
P = (
    dict(
        app="#07111f", side="#0a1628", panel="#0f1d2e", panel2="#13243a",
        input="#162a42", text="#f3f7fb", muted="#9fb0c3", border="#29445f",
        accent="#2dd4bf", accent2="#60a5fa", accent3="#a78bfa",
        good="#34d399", warn="#fbbf24", soft="rgba(45,212,191,.11)",
        shadow="rgba(0,0,0,.30)", plot="plotly_dark",
    )
    if dark_mode
    else dict(
        app="#f6f8fc", side="#eef3f9", panel="#ffffff", panel2="#f7f9fc",
        input="#ffffff", text="#172033", muted="#64748b", border="#d6deea",
        accent="#0f766e", accent2="#2563eb", accent3="#7c3aed",
        good="#059669", warn="#b45309", soft="rgba(15,118,110,.08)",
        shadow="rgba(15,23,42,.08)", plot="plotly_white",
    )
)

st.markdown(
    f"""
<style>
:root{{--app:{P['app']};--side:{P['side']};--panel:{P['panel']};--panel2:{P['panel2']};--input:{P['input']};--text:{P['text']};--muted:{P['muted']};--border:{P['border']};--accent:{P['accent']};--accent2:{P['accent2']};--accent3:{P['accent3']};--good:{P['good']};--warn:{P['warn']};--soft:{P['soft']};}}
html,body,[data-testid="stAppViewContainer"],.stApp{{background:var(--app)!important;color:var(--text)!important}}
[data-testid="stHeader"],[data-testid="stToolbar"]{{background:transparent!important}}
[data-testid="stSidebar"]{{background:var(--side)!important;border-right:1px solid var(--border)!important}}
.block-container{{max-width:1480px;padding-top:1.35rem;padding-bottom:4rem}}
.stApp h1,.stApp h2,.stApp h3,.stApp h4,.stApp h5,.stApp p,.stApp li,.stApp label{{color:var(--text)!important}}
[data-testid="stCaptionContainer"] p,.stCaption{{color:var(--muted)!important}}

/* Core controls: force readable foreground/background in both modes. */
[data-baseweb="input"]>div,[data-baseweb="textarea"]>div,[data-baseweb="select"]>div{{background:var(--input)!important;border:1px solid var(--border)!important;color:var(--text)!important;box-shadow:none!important;border-radius:10px!important}}
[data-baseweb="input"] input,[data-baseweb="textarea"] textarea,.stApp input,.stApp textarea{{background:var(--input)!important;color:var(--text)!important;-webkit-text-fill-color:var(--text)!important;caret-color:var(--accent)!important}}
.stApp input::placeholder,.stApp textarea::placeholder{{color:var(--muted)!important;opacity:.85!important}}
[data-baseweb="select"] span,[data-baseweb="select"] div,[data-baseweb="select"] svg{{color:var(--text)!important;fill:var(--text)!important}}
[data-baseweb="popover"],[role="listbox"],ul[role="listbox"]{{background:var(--panel)!important;border-color:var(--border)!important}}
[role="option"]{{background:var(--panel)!important;color:var(--text)!important}}
[role="option"]:hover,[aria-selected="true"][role="option"]{{background:var(--panel2)!important;color:var(--text)!important}}
[data-testid="stNumberInput"] button{{background:var(--panel2)!important;color:var(--text)!important;border-color:var(--border)!important}}
[data-testid="stFileUploaderDropzone"]{{background:var(--panel2)!important;border-color:var(--border)!important;color:var(--text)!important}}
[data-testid="stFileUploaderDropzone"] *{{color:var(--text)!important}}
[data-testid="stRadio"] label,[data-testid="stCheckbox"] label,[data-testid="stToggle"] label{{color:var(--text)!important}}

.hero{{padding:2.25rem 2.4rem;border:1px solid var(--border);border-radius:22px;background:radial-gradient(circle at 90% 12%,rgba(96,165,250,.16),transparent 28%),radial-gradient(circle at 8% 100%,rgba(45,212,191,.14),transparent 32%),var(--panel);box-shadow:0 20px 55px {P['shadow']};margin-bottom:1.3rem}}
.hero .eyebrow{{color:var(--accent)!important;font-size:.74rem;font-weight:800;letter-spacing:.16em;text-transform:uppercase}}
.hero h1{{font-size:clamp(2.1rem,4.6vw,4.2rem);letter-spacing:-.048em;line-height:.98;margin:.62rem 0 .9rem}}
.hero h1 span{{background:linear-gradient(90deg,var(--accent),var(--accent2));-webkit-background-clip:text;color:transparent!important}}
.hero .copy{{color:var(--muted)!important;max-width:900px;font-size:1.04rem;line-height:1.65}}
.pills{{display:flex;flex-wrap:wrap;gap:.52rem;margin-top:1.2rem}} .pills span{{background:var(--panel2);border:1px solid var(--border);border-radius:999px;padding:.4rem .7rem;font-size:.77rem;font-weight:650;color:var(--text)}}

[data-testid="stForm"],[data-testid="stExpander"]{{background:var(--panel)!important;border:1px solid var(--border)!important;border-radius:16px!important}}
div[data-testid="stMetric"]{{background:linear-gradient(180deg,var(--panel),var(--panel2));border:1px solid var(--border);padding:1rem 1.05rem;border-radius:16px;box-shadow:0 7px 20px {P['shadow']}}}
div[data-testid="stMetricLabel"] p{{color:var(--muted)!important}} div[data-testid="stMetricValue"]{{color:var(--text)!important}}
button[data-baseweb="tab"]{{background:var(--panel2)!important;border:1px solid var(--border)!important;border-radius:9px!important;margin-right:5px;padding:.55rem .82rem!important}}
button[data-baseweb="tab"] p{{color:var(--text)!important}}
button[data-baseweb="tab"][aria-selected="true"]{{background:linear-gradient(90deg,var(--soft),rgba(96,165,250,.10))!important;border-color:var(--accent)!important}}
.stButton>button,.stDownloadButton>button{{border-radius:10px;min-height:2.8rem;font-weight:750;border:1px solid var(--border);background:var(--panel2);color:var(--text)}}
.stButton>button[kind="primary"]{{background:linear-gradient(115deg,var(--accent),var(--accent2));color:#07111f!important;border:0}}
[data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:13px;overflow:hidden;background:var(--panel)!important}}

.seq-card{{background:linear-gradient(90deg,var(--panel),var(--panel2));border:1px solid var(--border);padding:1.05rem 1.2rem;border-radius:14px;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:1.18rem;overflow-wrap:anywhere;letter-spacing:.035em;margin:.45rem 0 1rem}}
.seq-card .pam{{background:var(--accent);color:#06111f;padding:.12rem .34rem;border-radius:6px;font-weight:850}}
.summary-card{{padding:1rem 1.05rem;border-radius:14px;border:1px solid var(--border);background:linear-gradient(90deg,var(--soft),rgba(96,165,250,.06));margin:.8rem 0 1rem}}
.summary-card span{{color:var(--muted)}}
.check-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.72rem;margin:.55rem 0 1rem}}
.check-card{{border:1px solid var(--border);border-radius:13px;padding:.9rem 1rem;background:var(--panel)}}
.check-card.good{{border-color:rgba(52,211,153,.45);background:linear-gradient(180deg,rgba(52,211,153,.07),var(--panel))}}
.check-card.review{{border-color:rgba(251,191,36,.42);background:linear-gradient(180deg,rgba(251,191,36,.06),var(--panel))}}
.check-head{{display:flex;justify-content:space-between;gap:.8rem;font-weight:760;color:var(--text)}} .check-desc{{color:var(--muted);font-size:.88rem;margin-top:.2rem}}
.status-good{{color:var(--good);font-size:.85rem;white-space:nowrap}} .status-review{{color:var(--warn);font-size:.85rem;white-space:nowrap}}
.info-strip{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:.7rem;margin:.7rem 0 1rem}}
.info-chip{{background:var(--panel);border:1px solid var(--border);border-radius:13px;padding:.82rem .9rem}} .info-chip small{{display:block;color:var(--muted);margin-bottom:.18rem}}
.source-card{{background:var(--panel);border:1px solid var(--border);border-left:3px solid var(--accent2);border-radius:13px;padding:.9rem 1rem;margin:.7rem 0 1rem}}
.source-card small{{color:var(--muted)}}
.method{{height:100%;background:var(--panel);border:1px solid var(--border);border-radius:13px;padding:1rem}} .method b{{color:var(--accent)}}
.oligo{{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--panel2);border:1px solid var(--border);border-radius:10px;padding:.8rem 1rem;margin:.4rem 0;overflow-wrap:anywhere}}
hr{{border-color:var(--border)!important}}
@media(max-width:850px){{.check-grid,.info-strip{{grid-template-columns:1fr}}.hero{{padding:1.45rem}}}}
</style>
""",
    unsafe_allow_html=True,
)
plot_template = P["plot"]

# ----------------------------- Sidebar -----------------------------------
st.sidebar.markdown(f"### CRISPR Studio\n**SpCas9 workbench · v{APP_VERSION}**")
st.sidebar.subheader("Design settings")
max_guides = st.sidebar.slider("Maximum guides", 5, 50, 20, 5)
min_score = st.sidebar.slider("Minimum heuristic score", 0, 90, 35, 5)
max_mismatches = st.sidebar.slider("Off-target mismatches", 0, 4, 3)
with st.sidebar.expander("Models & backends"):
    st.markdown(
        "- **SpCas9:** 20 nt + NGG\n"
        "- **Online records:** NCBI Nucleotide + Ensembl\n"
        "- **MIT/Hsu:** native\n"
        "- **CFD:** bundled published weights\n"
        "- **Doench RS2:** optional GuideMaker provider\n"
        "- **Whole genome:** GuideScan2\n"
        "- **CRISPRi:** Ensembl canonical TSS"
    )
st.sidebar.info("Research-use shortlist only. Confirm genomic context, genome build, off-targets and experimental controls before ordering guides.")

# ----------------------------- Helpers -----------------------------------
@st.cache_data(ttl=3600, show_spinner=False)
def cached_fetch(gene, organism, source):
    return fetch_sequence(gene, organism, source)


@st.cache_data(ttl=3600, show_spinner=False)
def cached_tss(gene, organism, transcript_id=""):
    return fetch_ensembl_tss_context(gene, organism, transcript_id=transcript_id)


@st.cache_data(ttl=3600, show_spinner=False)
def cached_knockout(gene, organism, isoform=""):
    return fetch_ncbi_knockout_context(gene, organism, protein_id=isoform)


@st.cache_data(ttl=3600, show_spinner=False)
def cached_accession(accession, database):
    return fetch_accession(accession, database)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: read_upload
# ----------------------------------------------------------------------------
def read_upload(f):
    if f is None:
        return ""
    b = f.getvalue()
    try:
        return b.decode("utf-8-sig")
    except UnicodeDecodeError:
        return b.decode("latin-1")



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: safe_name
# ----------------------------------------------------------------------------
def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "crispr_guides"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: style_plot
# ----------------------------------------------------------------------------
def style_plot(fig, height=420):
    fig.update_layout(
        template=plot_template,
        height=height,
        margin=dict(l=20, r=20, t=55, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend_title_text="",
        font=dict(color=P["text"]),
    )
    return fig



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: render_validation
# ----------------------------------------------------------------------------
def render_validation(g: GuideRNA, min_score):
    checks = validate_guide(g, min_score=min_score)
    items = [
        ("Spacer length", "Exactly 20 nt for SpCas9.", checks.get("length_ok", False)),
        ("PAM compatibility", "Target uses a valid NGG PAM.", checks.get("pam_ok", False)),
        ("Preferred GC", "Preferred guide GC is 40–70%.", checks.get("gc_in_preferred_range", False)),
        ("Acceptable GC", "Broad acceptable GC is 30–80%.", checks.get("gc_in_acceptable_range", False)),
        ("Homopolymer check", "Avoids extreme same-base runs.", checks.get("no_extreme_homopolymer", False)),
        ("Heuristic threshold", "Passes the selected sequence heuristic threshold.", checks.get("score_above_threshold", False)),
        ("Poly-T check", "No TTTT motif that may affect U6 expression.", checks.get("no_poly_t", False)),
        ("Specificity screen", "A local off-target reference was screened; GuideScan2 output is reviewed separately.", checks.get("specificity_screened", False)),
        ("Screen scope", "Local search includes substitution sites through at least three mismatches.", checks.get("screen_scope_complete", False)),
        ("High-risk sites", "No non-intended Critical or High local hit was found.", checks.get("no_critical_or_high_hits", False)),
        ("Aggregate specificity", "MIT specificity is at least 50 within the stated local scope.", checks.get("aggregate_specificity_ok", False)),
        ("Locus evidence", "Intended locus is verified, unique and screened against an unambiguous reference.", checks.get("locus_evidence_ok", False)),
    ]
    quality = [x for x in items if x[0] != "Specificity screen"]
    passed = sum(bool(v) for _, _, v in quality)
    overall = checks.get("overall_pass", False)
    status = "LOCAL CHECKS MET" if overall else "REVIEW"
    status_class = "status-good" if overall else "status-review"
    st.markdown(
        f'<div class="summary-card"><b>Guide quality: <span class="{status_class}">{status}</span></b><br>'
        f'<span>{passed}/{len(quality)} quality checks passed. Specificity screening is reported separately.</span></div>',
        unsafe_allow_html=True,
    )
    cards = []
    for title, desc, ok in items:
        if title == "Specificity screen" and not ok:
            label, cls, scls = "Not run", "review", "status-review"
        else:
            label, cls, scls = ("Pass", "good", "status-good") if ok else ("Review", "review", "status-review")
        cards.append(
            f'<div class="check-card {cls}"><div class="check-head"><span>{title}</span>'
            f'<span class="{scls}">{label}</span></div><div class="check-desc">{desc}</div></div>'
        )
    st.markdown('<div class="check-grid">' + "".join(cards) + "</div>", unsafe_allow_html=True)


# ----------------------------- Workspace ---------------------------------
st.markdown(
    """<section class="hero"><div class="eyebrow">Open-source CRISPR design workbench</div>
    <h1>Design with <span>confidence and provenance.</span></h1>
    <div class="copy">Design and rank SpCas9 guides from a gene symbol, a verified online accession, or your own DNA. Inspect activity and specificity metrics, screen reference sequences, and use true TSS-aware CRISPRi for gene-based repression.</div>
    <div class="pills"><span>Gene lookup</span><span>NCBI / Ensembl accession</span><span>Paste DNA</span><span>TSS-aware CRISPRi</span><span>MIT/Hsu</span><span>GuideScan2</span></div></section>""",
    unsafe_allow_html=True,
)
st.subheader("Design workspace")
input_mode = st.radio(
    "Target input",
    ["Gene lookup", "Accession ID", "Paste sequence"],
    horizontal=True,
    label_visibility="collapsed",
)

with st.container(border=True):
    left, right = st.columns([1.25, 1])
    with left:
        gene_name = "Custom target"
        organism = "Not specified"
        source = "manual"
        custom_sequence = ""
        accession_query = ""
        accession_database = "Auto"
        if input_mode == "Gene lookup":
            gene_name = st.text_input("Gene symbol", value="TP53")
            organism = st.text_input("Organism", value="Homo sapiens")
            source = st.radio("Sequence database", ["NCBI", "Ensembl"], horizontal=True).lower()
            isoform = st.text_input("Isoform accession optional", help="For NCBI knockout: a versioned NP or NM accession from this gene. For CRISPRi: an Ensembl transcript ID. Leave blank for the stated representative selection.")
            st.caption("NCBI knockout uses genomic coding annotations. For annotated knockout design select NCBI; Ensembl gene lookup is reserved for CRISPRi. CRISPRi uses Ensembl human/mouse dCas9-KRAB placement rules.")
        elif input_mode == "Accession ID":
            accession_query = st.text_input(
                "Accession / stable ID",
                placeholder="Examples: ENSG00000141510, ENST00000269305, a small genomic accession",
            )
            accession_database = st.selectbox("Online database", ["Auto", "NCBI Nucleotide", "Ensembl"])
            st.caption("Auto routes ENS* stable IDs to Ensembl and other nucleotide accessions to NCBI Nucleotide. Only genomic DNA is designed: RNA/cDNA accessions are rejected and Ensembl IDs retrieve unspliced genomic intervals. Large chromosomes exceed the input limit.")
            gene_name = accession_query.strip() or "Accession target"
            source = "accession"
        else:
            custom_sequence = st.text_area("Target DNA / genomic FASTA", height=180, placeholder=">target\nATG...")

    with right:
        intent = st.radio("Design intent", ["Knockout", "CRISPRi repression"])
        application = "knockout" if intent == "Knockout" else "crispri"
        custom_tss = None
        if application == "crispri":
            st.caption("Mammalian dCas9-KRAB placement heuristic only. It does not model chromatin, the active cellular TSS, bacterial CRISPRi or plant repressors.")
        if input_mode == "Paste sequence" and application == "crispri":
            custom_tss = st.number_input(
                "TSS position in pasted sequence (1-based)",
                min_value=1,
                value=1,
                step=1,
                help="Sequence must be genomic DNA in 5′→3′ transcriptional orientation.",
            )
        if input_mode == "Accession ID" and application == "crispri":
            st.info("For scientifically TSS-aware CRISPRi, use Gene lookup. Accession mode currently designs knockout guides because many nucleotide accessions do not provide a reliable annotated TSS context.")

        screen_mode = st.selectbox("Specificity analysis", ["None", "Local reference", "Whole genome (GuideScan2)"])
        ref_upload = None
        ref_text = ""
        genome_index = ""
        target_reference_id = ""
        target_reference_start = 1
        target_reference_strand = "+"
        if screen_mode == "Local reference":
            ref_upload = st.file_uploader("Reference FASTA / text", type=["fa", "fasta", "fna", "txt"])
            ref_text = st.text_area("Or paste reference sequence", height=80)
            target_reference_id = st.text_input("Intended target record ID optional", help="Exact FASTA contig ID containing the intended target region. Leave blank to retain all exact hits.")
            target_reference_start = st.number_input("Target region start in reference (1-based left edge)", min_value=1, value=1, step=1)
            target_reference_strand = st.selectbox("Target region orientation in reference", ["+", "-"])
            st.caption("Local screening checks NGG sites in the supplied linear reference only. No alternative PAMs, variants, bulges or circular junctions are searched.")
        elif screen_mode.startswith("Whole genome"):
            genome_index = st.text_input("GuideScan2 index path", value=os.getenv("GUIDESCAN_INDEX", ""))
            st.caption("GuideScan2 detected on this host: " + ("yes" if guidescan2_available() else "no"))

    submitted = st.button("Design and analyze guides", type="primary", use_container_width=True)

# ----------------------------- Analysis ----------------------------------
if submitted:
    st.session_state.pop("analysis", None)
    try:
        reference = None
        if screen_mode == "Local reference":
            raw_reference = read_upload(ref_upload) or ref_text
            reference = parse_reference(raw_reference)
            total = sum(map(len, reference.values()))
            if not total:
                raise ValueError("Local reference mode is enabled but no reference was supplied.")
            if total > MAX_LOCAL_REFERENCE_BP:
                raise ValueError(f"Local reference exceeds {MAX_LOCAL_REFERENCE_BP:,} bp hosted limit.")

        intended_region = (target_reference_id.strip(), int(target_reference_start) - 1, target_reference_strand) if reference and target_reference_id.strip() else None
        crispri_rows: List[CRISPRiGuide] = []
        tss_context = None
        online_record = None
        knockout_context = None

        if input_mode == "Accession ID":
            if not accession_query.strip():
                raise ValueError("Enter an accession or stable ID.")
            if application == "crispri":
                raise ValueError("Accession-ID CRISPRi is not used because a nucleotide accession alone may not define a reliable TSS. Use Gene lookup for TSS-aware CRISPRi, or Paste sequence with an explicit TSS.")
            online_record: AccessionRecord = cached_accession(accession_query.strip(), accession_database)
            sequence = online_record.sequence
            accession = online_record.accession
            description = online_record.description
            gene_name = accession_query.strip()
            organism = online_record.database
            source = "accession"
            guides = design_guides(
                sequence,
                application="knockout",
                max_guides=max_guides,
                min_score=float(min_score),
                genome_context=reference,
                max_mismatches=max_mismatches,
                intended_region=intended_region,
            )
        elif application == "crispri":
            if input_mode == "Gene lookup":
                if not gene_name.strip() or not organism.strip():
                    raise ValueError("Enter both gene symbol and organism.")
                tss_context = cached_tss(gene_name.strip(), organism.strip(), isoform.strip())
                sequence = tss_context.sequence
                accession = tss_context.transcript_id
                description = tss_context.description
                crispri_rows = design_crispri_guides(
                    tss_context,
                    max_guides=max_guides,
                    min_score=float(min_score),
                    genome_context=reference,
                        max_mismatches=max_mismatches,
                    intended_region=intended_region,
                )
            else:
                sequence = clean_dna_sequence(custom_sequence)
                if len(sequence) < 50:
                    raise ValueError("The target sequence must contain at least 50 bp.")
                if len(sequence) > MAX_CUSTOM_BP:
                    raise ValueError(f"Target exceeds {MAX_CUSTOM_BP:,} bp limit.")
                accession = "CUSTOM_TSS"
                description = "User-provided genomic sequence with declared TSS"
                crispri_rows = design_crispri_from_sequence(
                    sequence,
                    int(custom_tss),
                    max_guides=max_guides,
                    min_score=float(min_score),
                    genome_context=reference,
                        max_mismatches=max_mismatches,
                    intended_region=intended_region,
                )
            guides = [x.guide for x in crispri_rows]
        else:
            if input_mode == "Gene lookup" and source == "ncbi":
                knockout_context = cached_knockout(gene_name.strip(), organism.strip(), isoform.strip())
                sequence = knockout_context.sequence
                accession = knockout_context.accession
                description = knockout_context.description
                guides = design_knockout_guides(knockout_context, max_guides=max_guides, min_score=float(min_score), genome_context=reference, max_mismatches=max_mismatches, intended_region=intended_region)
            else:
                if input_mode == "Paste sequence":
                    sequence = clean_dna_sequence(custom_sequence)
                    if len(sequence) < 50:
                        raise ValueError("The target sequence must contain at least 50 bp.")
                    if len(sequence) > MAX_CUSTOM_BP:
                        raise ValueError(f"Target exceeds {MAX_CUSTOM_BP:,} bp limit.")
                    accession = "CUSTOM"
                    description = "User-provided target sequence; coding annotation not supplied"
                else:
                    if not gene_name.strip() or not organism.strip():
                        raise ValueError("Enter both gene symbol and organism.")
                    raise ValueError("Annotated knockout gene design requires NCBI genomic CDS lookup. Select NCBI; Ensembl transcript discovery is disabled.")
                guides = design_guides(sequence, application="knockout", max_guides=max_guides,
                    min_score=float(min_score), genome_context=reference, max_mismatches=max_mismatches, intended_region=intended_region)

        genome_rows = (
            run_guidescan2(guides, genome_index, max_mismatches=max_mismatches)
            if screen_mode.startswith("Whole genome") and guides
            else []
        )
        st.session_state["analysis"] = dict(
            guides=guides,
            crispri_rows=crispri_rows,
            tss_context=tss_context,
            online_record=online_record,
            knockout_context=knockout_context,
            min_score=min_score,
            metadata=make_metadata(sequence, reference, version=APP_VERSION, target=gene_name, organism=organism, source=source, application=application, accession=accession, mismatch_limit=max_mismatches, min_score=min_score, max_guides=max_guides, screen_mode=screen_mode, intended_target_record=target_reference_id, intended_region=intended_region, model_status=rs2_status(), genomic_context=({"reference_accession": knockout_context.accession, "region_start": knockout_context.region_start, "selected_protein": knockout_context.selected_protein, "selection_rule": knockout_context.selection_rule} if knockout_context else None), tss_context=({k: v for k, v in vars(tss_context).items() if k != "sequence"} if tss_context else ({"custom_tss_1based": custom_tss} if application == "crispri" else None))),
            sequence=sequence,
            reference=reference,
            genome_rows=genome_rows,
            gene=gene_name.strip() or "Custom target",
            organism=organism.strip(),
            source=source,
            application=application,
            accession=accession,
            description=description,
            screen_mode=screen_mode,
            max_mismatches=max_mismatches,
            intended_region=intended_region,
        )
        st.success(f"Analysis complete: {len(guides)} guide candidates passed the filter.")
    except Exception as exc:
        st.error(f"Could not complete the analysis: {exc}")

# ----------------------------- Report ------------------------------------
if "analysis" in st.session_state:
    r = st.session_state["analysis"]
    guides: List[GuideRNA] = r["guides"]
    sequence = r["sequence"]
    reference = r["reference"]
    crispri_rows: List[CRISPRiGuide] = r.get("crispri_rows", [])

    st.divider()
    st.subheader(f"Analysis report · {r['gene']}")
    st.caption(f"{r['description'][:260]} | Accession: {r['accession']} | Analyzed DNA: {len(sequence):,} bp")

    if r.get("online_record"):
        rec: AccessionRecord = r["online_record"]
        st.markdown(
            f'<div class="source-card"><b>Verified online record</b><br>'
            f'<small>{escape(rec.database)} · {escape(rec.object_type)} · resolved accession {escape(rec.accession)} · {rec.length:,} bp</small></div>',
            unsafe_allow_html=True,
        )
    if r["application"] == "crispri" and r.get("tss_context"):
        c = r["tss_context"]
        st.success(
            f"TSS-aware CRISPRi: {c.assembly} {c.chromosome}:{c.tss_coordinate:,} ({c.strand_label}) · "
            f"annotated transcript {c.transcript_id} · accepted window {CRISPRI_MIN:+d} to {CRISPRI_MAX:+d} bp; "
            f"preferred {CRISPRI_OPTIMAL_MIN:+d} to {CRISPRI_OPTIMAL_MAX:+d} bp."
        )
    elif r.get("knockout_context"):
        context = r["knockout_context"]
        st.success(f"Genomic CDS filter applied: {context.accession}, selected protein {context.selected_protein}.")
        st.caption(context.selection_rule)
    elif r["application"] == "knockout":
        st.warning("Genomic sequence exploration: no CDS filter was applied. Confirm the coding or regulatory context; pasted DNA must be contiguous genomic sequence.")
    st.caption("Scores prioritize candidates; they are not probabilities of editing or experimental validation. " + r["metadata"]["model_status"])

    if not guides:
        st.warning("No candidate guides passed the current filters. Reduce the heuristic threshold or inspect whether the target region contains NGG PAM sites.")
    else:
        df = pd.DataFrame([x.to_dict() for x in crispri_rows]) if crispri_rows else pd.DataFrame([g.to_dict() for g in guides])
        df.insert(0, "Rank", range(1, len(df) + 1))

        metrics = st.columns(5)
        metrics[0].metric("Guides retained", len(guides))
        metrics[1].metric("Best heuristic", f"{max(g.score for g in guides):.1f}")
        metrics[2].metric("Median GC", f"{median(g.gc_content for g in guides):.1f}%")
        ds = [g.doench_score for g in guides if g.doench_score is not None]
        metrics[3].metric("Best Doench RS2", f"{max(ds):.1f}" if ds else "Not available")
        ms = [g.specificity_score for g in guides if g.specificity_score is not None]
        metrics[4].metric("Best MIT specificity", f"{max(ms):.1f}" if ms else ("GuideScan2 run" if r["genome_rows"] else "Not screened"))

        tabs = st.tabs(["Ranked guides", "Design landscape", "Guide details", "Off-target screen", "Export", "Methods & limits"])

        with tabs[0]:
            st.markdown("#### Ranked candidate table")
            st.dataframe(df, use_container_width=True, hide_index=True, height=min(650, 100 + len(df) * 35))

        with tabs[1]:
            if crispri_rows:
                fig = px.scatter(
                    df,
                    x="TSS distance (bp)",
                    y="Heuristic",
                    color="TSS band",
                    size="GC%",
                    hover_data=["Rank", "Spacer (20 nt)", "PAM", "Chromosome", "Genomic start"],
                    title="CRISPRi candidates relative to the annotated TSS",
                )
                fig.add_vrect(x0=50, x1=100, opacity=.10, line_width=0, annotation_text="preferred +50..+100")
                fig.add_vline(x=0, line_dash="dash", annotation_text="TSS")
            else:
                fig = px.scatter(
                    df,
                    x="Start",
                    y="Heuristic",
                    color="Strand",
                    size="GC%",
                    hover_data=["Rank", "Spacer (20 nt)", "PAM"],
                    title="Sequence heuristic across the target",
                )
            st.plotly_chart(style_plot(fig, 460), use_container_width=True)
            hist = px.histogram(df, x="GC%", nbins=12, title="GC-content distribution")
            st.plotly_chart(style_plot(hist, 330), use_container_width=True)

        with tabs[2]:
            idx = st.selectbox("Inspect a guide", range(len(guides)), format_func=lambda i: f"#{i+1} · {guides[i].sequence}")
            g = guides[idx]
            st.markdown(f'<div class="seq-card">5′—{g.sequence}<span class="pam">{g.pam}</span>—3′</div>', unsafe_allow_html=True)
            q = st.columns(5)
            q[0].metric("Heuristic", f"{g.score:.1f}")
            q[1].metric("Doench RS2", f"{g.doench_score:.1f}" if g.doench_score is not None else "N/A")
            q[2].metric("GC", f"{g.gc_content:.1f}%")
            q[3].metric("MIT", f"{g.specificity_score:.1f}" if g.specificity_score is not None else "N/A")
            q[4].metric("CFD spec.", f"{g.cfd_specificity:.1f}" if g.cfd_specificity is not None else "N/A")

            if crispri_rows:
                cg = crispri_rows[idx]
                chips = [
                    ("TSS distance", f"{cg.tss_distance:+g} bp"),
                    ("Placement", cg.tss_band),
                    ("Assembly", cg.assembly or "Custom"),
                    ("Genomic locus", f"{cg.chromosome}:{cg.genomic_start}-{cg.genomic_end}" if cg.genomic_start else "Custom sequence"),
                ]
                st.markdown(
                    '<div class="info-strip">' + "".join(f'<div class="info-chip"><small>{a}</small><strong>{b}</strong></div>' for a, b in chips) + "</div>",
                    unsafe_allow_html=True,
                )

            st.markdown("#### Validation & quality")
            render_validation(g, r["min_score"])
            st.markdown("#### Heuristic score breakdown")
            breakdown = score_breakdown(
                g.sequence,
                g.pam,
                application=None,
                start=g.start,
                sequence_length=len(sequence),
            )
            rows = [{"Component": k, "Contribution": v} for k, v in breakdown.items() if k != "Final score"]
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            st.markdown(
                f'<div class="summary-card"><b>Final heuristic score: {breakdown.get("Final score", "N/A")}</b><br>'
                '<span>This transparent heuristic is reported separately from Doench Rule Set 2.</span></div>',
                unsafe_allow_html=True,
            )
            st.caption("Exported sequences are 20 nt DNA spacers. Vector overhangs, an added U6 G and sgRNA scaffold depend on your expression system; these are not complete cloning oligos.")

        with tabs[3]:
            if reference:
                oi = st.selectbox("Inspect local off-targets", range(len(guides)), key="otguide", format_func=lambda i: f"#{i+1} · {guides[i].sequence}")
                rep = analyze_offtargets(guides[oi].sequence, reference, max_mismatches=r["max_mismatches"], intended_target=guides[oi].intended_target)
                c = st.columns(5)
                c[0].metric("MIT specificity", f"{rep.mit_specificity:.1f}")
                c[1].metric("CFD specificity", f"{rep.cfd_specificity:.1f}" if rep.cfd_specificity is not None else "N/A")
                c[2].metric("All local hits", rep.total_hits)
                c[3].metric("PAM sites", rep.pam_sites_scanned)
                c[4].metric("Contigs", rep.reference_contigs)
                st.caption(
                    f"Screened through {rep.screened_mismatch_radius} mismatch(es). Intended target status: {rep.intended_target_status}. "
                    f"Exact reference matches: {rep.exact_matches}. Risk counts across all hits: "
                    f"Critical {rep.risk_counts['Critical']}, High {rep.risk_counts['High']}, "
                    f"Moderate {rep.risk_counts['Moderate']}, Low {rep.risk_counts['Low']}. "
                    f"Maximum per-site risk: MIT {rep.max_mit_risk*100:.2f}%, CFD {rep.max_cfd_risk*100:.2f}%. "
                    f"Displayed {len(rep.hits)} of {rep.total_hits} hits; all hits contribute to scores and risk counts."
                )
                if rep.screened_mismatch_radius < 3:
                    st.warning("This local search covers fewer than three mismatches and cannot receive LOCAL CHECKS MET.")
                if rep.risk_counts['Critical'] or rep.risk_counts['High']:
                    st.warning("At least one non-intended Critical or High local hit requires review, regardless of the aggregate specificity score.")
                if rep.intended_target_status != "verified_locus":
                    st.warning("Intended locus was not verified. Exact matches were retained; a high score alone cannot confirm a valid or unique genomic target.")
                if rep.ambiguous_bases:
                    st.warning(f"Reference contains {rep.ambiguous_bases:,} ambiguous bases. Sites overlapping them were skipped.")
                if rep.hits:
                    st.dataframe(pd.DataFrame([h.to_dict() for h in rep.hits]), use_container_width=True, hide_index=True)
                else:
                    st.info("No reported local hits within this NGG-only search scope. This does not establish genome-wide specificity.")
            elif r["genome_rows"]:
                st.markdown("#### GuideScan2 whole-genome results")
                st.dataframe(pd.DataFrame(r["genome_rows"]), use_container_width=True, hide_index=True)
            else:
                st.info("Specificity was not screened for this run. Choose Local reference or Whole genome (GuideScan2) to calculate off-target evidence.")

        with tabs[4]:
            base = safe_name(r["gene"])
            csvb = df.to_csv(index=False).encode()
            fasta = "\n".join(f">{base}_g{i}|pam={g.pam}|strand={g.strand}\n{g.sequence}" for i, g in enumerate(guides, 1)).encode()
            js = export_json(r["metadata"], guides, rows=json.loads(df.to_json(orient="records")), genome_rows=r["genome_rows"])
            d = st.columns(3)
            d[0].download_button("Download CSV", csvb, file_name=f"{base}_v330_guides.csv", use_container_width=True)
            d[1].download_button("Download FASTA", fasta, file_name=f"{base}_v330_spacers.fasta", use_container_width=True)
            d[2].download_button("Download JSON", js, file_name=f"{base}_v330_analysis.json", use_container_width=True)

        with tabs[5]:
            cards = st.columns(4)
            content = [
                ("Input provenance", "Gene symbols, NCBI/Ensembl accessions, or user DNA with source metadata."),
                ("On-target", "Doench Rule Set 2 when the optional provider is installed; transparent heuristic stays separate."),
                ("Specificity", "Native MIT/Hsu and bundled CFD; local multi-contig or GuideScan2 whole-genome search. Local scores cover NGG sites only."),
                ("CRISPRi", "Human/mouse dCas9-KRAB heuristic around an annotated Ensembl TSS. Active TSS and chromatin are not measured."),
            ]
            for col, (title, text) in zip(cards, content):
                col.markdown(f'<div class="method"><b>{title}</b><p>{text}</p></div>', unsafe_allow_html=True)
            st.markdown("#### Scientific boundaries")
            st.markdown(
                "- **Accession lookup** resolves nucleotide targets from NCBI Nucleotide or Ensembl and records provenance.\n"
                "- **CRISPRi is TSS-aware for gene lookup.** Accession-only CRISPRi is intentionally not inferred because many nucleotide records do not define a reliable TSS.\n"
                "- NCBI knockout gene lookup scans genomic DNA and filters cuts to annotated CDS. Accession/pasted genomic sequence modes do not apply CDS annotations; Ensembl gene-symbol knockout is disabled.\n"
                "- Whole-genome specificity requires GuideScan2 plus a matching prebuilt genome index.\n"
                "- Variant-aware filtering, chromatin state and DNA/RNA bulges require external genome-aware resources.\n"
                "- Computational scores prioritize candidates; they do not replace experimental validation."
            )

st.divider()
st.caption("CRISPR Studio v3.3 · research/educational use · MIT License")
