#!/usr/bin/env python3
"""CRISPR Studio v3.1.1 polished Streamlit dashboard."""
from __future__ import annotations

import json
import os
import re
from statistics import median
from typing import List

import pandas as pd
import plotly.express as px
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

APP_VERSION = "3.1.1"
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
        "About": "CRISPR Studio v3.1.1 — SpCas9 design with TSS-aware CRISPRi.",
    },
)

dark_mode = st.sidebar.toggle("Dark mode", value=True, key="dark_mode")
P = (
    dict(app="#06110e", side="#091713", panel="#0d201b", alt="#112923", text="#eefbf6", muted="#9db8ae", border="#21483c", accent="#34d399", accent2="#22d3ee", soft="rgba(52,211,153,.10)", shadow="rgba(0,0,0,.30)")
    if dark_mode
    else dict(app="#f4faf7", side="#edf7f2", panel="#ffffff", alt="#f1f8f5", text="#11211b", muted="#5e756d", border="#d2e6dd", accent="#059669", accent2="#0891b2", soft="rgba(5,150,105,.08)", shadow="rgba(15,60,45,.10)")
)
st.markdown(f"""
<style>
:root{{--app:{P['app']};--side:{P['side']};--panel:{P['panel']};--alt:{P['alt']};--text:{P['text']};--muted:{P['muted']};--border:{P['border']};--accent:{P['accent']};--accent2:{P['accent2']};--soft:{P['soft']};}}
html,body,[data-testid="stAppViewContainer"],.stApp{{background:var(--app)!important;color:var(--text)!important}}
[data-testid="stHeader"],[data-testid="stToolbar"]{{background:transparent!important}}
[data-testid="stSidebar"]{{background:var(--side)!important;border-right:1px solid var(--border)}}
.block-container{{max-width:1460px;padding-top:1.45rem;padding-bottom:4rem}}
.stApp p,.stApp li,.stApp label,.stApp h1,.stApp h2,.stApp h3,.stApp h4,.stApp h5{{color:var(--text)}}
[data-testid="stCaptionContainer"] p{{color:var(--muted)!important}}
.hero{{padding:2.15rem 2.3rem;border:1px solid var(--border);border-radius:24px;background:radial-gradient(circle at 90% 10%,rgba(34,211,238,.16),transparent 27%),radial-gradient(circle at 9% 100%,rgba(52,211,153,.16),transparent 31%),var(--panel);box-shadow:0 18px 52px {P['shadow']};margin-bottom:1.25rem}}
.hero .k{{color:var(--accent);font-size:.76rem;font-weight:800;letter-spacing:.14em;text-transform:uppercase}}
.hero h1{{font-size:clamp(2rem,4.6vw,4.05rem);letter-spacing:-.045em;line-height:.98;margin:.6rem 0 1rem}}
.hero h1 span{{color:var(--accent)}} .hero .copy{{color:var(--muted);max-width:880px;font-size:1.04rem}}
.pills{{display:flex;flex-wrap:wrap;gap:.55rem;margin-top:1.2rem}} .pills span{{background:var(--alt);border:1px solid var(--border);border-radius:999px;padding:.42rem .72rem;font-size:.78rem;font-weight:650}}
[data-testid="stForm"],[data-testid="stExpander"]{{background:var(--panel)!important;border:1px solid var(--border)!important;border-radius:16px!important}}
[data-baseweb="input"]>div,[data-baseweb="textarea"]>div,[data-baseweb="select"]>div{{background:var(--alt)!important;border-color:var(--border)!important;color:var(--text)!important;box-shadow:none!important}}
[data-baseweb="select"] *{{color:var(--text)!important}}
.stApp input,.stApp textarea{{background:transparent!important;color:var(--text)!important;caret-color:var(--accent)!important}}
[data-baseweb="popover"],[role="listbox"]{{background:var(--panel)!important}} [role="option"]{{background:var(--panel)!important;color:var(--text)!important}} [role="option"]:hover{{background:var(--alt)!important}}
div[data-testid="stMetric"]{{background:linear-gradient(180deg,var(--panel),var(--alt));border:1px solid var(--border);padding:1rem 1.05rem;border-radius:17px;box-shadow:0 8px 22px {P['shadow']}}}
div[data-testid="stMetricLabel"] p{{color:var(--muted)!important}} div[data-testid="stMetricValue"]{{color:var(--text)!important}}
button[data-baseweb="tab"]{{background:var(--alt)!important;border:1px solid var(--border)!important;border-radius:10px!important;margin-right:5px;padding:.55rem .8rem!important}}
button[data-baseweb="tab"][aria-selected="true"]{{background:linear-gradient(90deg,var(--soft),rgba(34,211,238,.10))!important;border-color:var(--accent)!important}}
.stButton>button,.stDownloadButton>button{{border-radius:12px;min-height:2.8rem;font-weight:750;border-color:var(--border)}} .stButton>button[kind="primary"]{{background:linear-gradient(115deg,var(--accent),var(--accent2));color:#04120d;border:0}}
[data-testid="stDataFrame"]{{border:1px solid var(--border);border-radius:14px;overflow:hidden;background:var(--panel)!important}}
.seq-card{{background:linear-gradient(90deg,var(--panel),var(--alt));border:1px solid var(--border);padding:1.05rem 1.2rem;border-radius:16px;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:1.18rem;overflow-wrap:anywhere;letter-spacing:.035em;margin:.45rem 0 1rem}}
.seq-card .pam{{background:var(--accent);color:#04120d;padding:.12rem .34rem;border-radius:7px;font-weight:850}}
.summary-card{{padding:1rem 1.05rem;border-radius:15px;border:1px solid var(--border);background:linear-gradient(90deg,var(--soft),rgba(34,211,238,.07));margin:.8rem 0 1rem}} .summary-card span{{color:var(--muted)}}
.check-grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.75rem;margin:.55rem 0 1rem}} .check-card{{border:1px solid var(--border);border-radius:15px;padding:.9rem 1rem;background:var(--panel)}}
.check-card.good{{border-color:rgba(52,211,153,.48);background:linear-gradient(180deg,rgba(52,211,153,.09),var(--panel))}} .check-card.review{{border-color:rgba(245,158,11,.45);background:linear-gradient(180deg,rgba(245,158,11,.08),var(--panel))}}
.check-head{{display:flex;justify-content:space-between;gap:.8rem;font-weight:760;color:var(--text)}} .check-desc{{color:var(--muted);font-size:.88rem;margin-top:.2rem}} .status-good{{color:#34d399;font-size:.86rem;white-space:nowrap}} .status-review{{color:#f59e0b;font-size:.86rem;white-space:nowrap}}
.info-strip{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:.7rem;margin:.7rem 0 1rem}} .info-chip{{background:var(--panel);border:1px solid var(--border);border-radius:14px;padding:.82rem .9rem}} .info-chip small{{display:block;color:var(--muted);margin-bottom:.18rem}}
.method{{height:100%;background:var(--panel);border:1px solid var(--border);border-radius:15px;padding:1rem}} .method b{{color:var(--accent)}} .oligo{{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:var(--alt);border:1px solid var(--border);border-radius:12px;padding:.8rem 1rem;margin:.4rem 0;overflow-wrap:anywhere}}
hr{{border-color:var(--border)!important}} @media(max-width:850px){{.check-grid,.info-strip{{grid-template-columns:1fr}}.hero{{padding:1.45rem}}}}
</style>
""", unsafe_allow_html=True)
plot_template = "plotly_dark" if dark_mode else "plotly_white"

st.sidebar.markdown(f"### CRISPR Studio\n**SpCas9 workbench · v{APP_VERSION}**")
st.sidebar.subheader("Design settings")
max_guides = st.sidebar.slider("Maximum guides", 5, 50, 20, 5)
min_score = st.sidebar.slider("Minimum heuristic score", 0, 90, 35, 5)
max_mismatches = st.sidebar.slider("Off-target mismatches", 0, 4, 3)
with st.sidebar.expander("Models & backends"):
    st.markdown("- **SpCas9:** 20 nt + NGG\n- **MIT/Hsu:** native\n- **Doench RS2 / CFD:** optional provider\n- **Whole genome:** GuideScan2\n- **CRISPRi:** Ensembl canonical TSS")
st.sidebar.info("Research-use shortlist only. Validate genome build, target biology, off-targets and experimental controls before ordering guides.")

@st.cache_data(ttl=3600, show_spinner=False)
def cached_fetch(gene, organism, source): return fetch_sequence(gene, organism, source)
@st.cache_data(ttl=3600, show_spinner=False)
def cached_tss(gene, organism): return fetch_ensembl_tss_context(gene, organism)

def read_upload(f):
    if f is None: return ""
    b = f.getvalue()
    try: return b.decode("utf-8-sig")
    except UnicodeDecodeError: return b.decode("latin-1")

def safe_name(s): return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "crispr_guides"

def style_plot(fig, height=420):
    fig.update_layout(template=plot_template, height=height, margin=dict(l=20, r=20, t=55, b=20), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", legend_title_text="")
    return fig

def render_validation(g: GuideRNA):
    checks = validate_guide(g)
    items = [
        ("Spacer length", "Exactly 20 nt for SpCas9.", checks.get("length_ok", False)),
        ("PAM compatibility", "Target uses a valid NGG PAM.", checks.get("pam_ok", False)),
        ("Preferred GC", "Preferred guide GC is 40–70%.", checks.get("gc_in_preferred_range", False)),
        ("Acceptable GC", "Broad acceptable GC is 30–80%.", checks.get("gc_in_acceptable_range", False)),
        ("Homopolymer check", "Avoids extreme same-base runs.", checks.get("no_extreme_homopolymer", False)),
        ("Activity threshold", "Passes the selected activity threshold.", checks.get("score_above_threshold", False)),
        ("Poly-T check", "No TTTT motif that may affect U6 expression.", checks.get("no_poly_t", False)),
        ("Specificity screen", "An off-target reference/index was actually screened.", checks.get("specificity_screened", False)),
        ("Specificity quality", "Specificity is acceptable when a screen is available.", checks.get("specificity_ok", False)),
    ]
    biological = [x for x in items if x[0] != "Specificity screen"]
    passed = sum(bool(v) for _, _, v in biological)
    overall = checks.get("overall_pass", False)
    status = "PASS" if overall else "REVIEW"
    status_class = "status-good" if overall else "status-review"
    st.markdown(f'<div class="summary-card"><b>Guide quality: <span class="{status_class}">{status}</span></b><br><span>{passed}/{len(biological)} quality checks passed. Off-target screening is shown separately so an unscreened guide is not presented as inherently unsafe.</span></div>', unsafe_allow_html=True)
    cards = []
    for title, desc, ok in items:
        if title == "Specificity screen" and not ok: label, cls, scls = "Not run", "review", "status-review"
        else: label, cls, scls = ("Pass", "good", "status-good") if ok else ("Review", "review", "status-review")
        cards.append(f'<div class="check-card {cls}"><div class="check-head"><span>{title}</span><span class="{scls}">{label}</span></div><div class="check-desc">{desc}</div></div>')
    st.markdown('<div class="check-grid">' + ''.join(cards) + '</div>', unsafe_allow_html=True)

st.markdown("""<section class="hero"><div class="k">Open-source CRISPR design workbench</div><h1>Design sharper <span>CRISPR guides.</span></h1><div class="copy">Rank SpCas9 guides, inspect activity and specificity metrics, screen local references or indexed genomes, and design CRISPRi guides against a real annotated transcription start site.</div><div class="pills"><span>20 nt + NGG</span><span>Both strands</span><span>TSS-aware CRISPRi</span><span>MIT/Hsu</span><span>Doench RS2</span><span>GuideScan2</span><span>CSV · FASTA · JSON</span></div></section>""", unsafe_allow_html=True)
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
            custom_tss = st.number_input("TSS position in pasted sequence (1-based)", min_value=1, value=1, step=1, help="Sequence must be genomic DNA in 5′→3′ transcriptional orientation.")
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
        crispri_rows: List[CRISPRiGuide] = []; tss_context = None
        if application == "crispri":
            if input_mode == "Gene lookup":
                if not gene_name.strip() or not organism.strip(): raise ValueError("Enter both gene symbol and organism.")
                tss_context = cached_tss(gene_name.strip(), organism.strip()); sequence = tss_context.sequence; accession = tss_context.transcript_id; description = tss_context.description
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
        genome_rows = run_guidescan2(guides, genome_index, max_mismatches=max_mismatches) if screen_mode.startswith("Whole genome") and guides else []
        st.session_state["analysis"] = dict(guides=guides, crispri_rows=crispri_rows, tss_context=tss_context, sequence=sequence, reference=reference, genome_rows=genome_rows, gene=gene_name.strip() or "Custom target", organism=organism.strip(), source=source, application=application, accession=accession, description=description, screen_mode=screen_mode, max_mismatches=max_mismatches)
        st.success(f"Analysis complete: {len(guides)} guide candidates passed the filter.")
    except Exception as exc:
        st.error(f"Could not complete the analysis: {exc}")

if "analysis" in st.session_state:
    r = st.session_state["analysis"]; guides: List[GuideRNA] = r["guides"]; sequence = r["sequence"]; reference = r["reference"]; crispri_rows: List[CRISPRiGuide] = r.get("crispri_rows", [])
    st.divider(); st.subheader(f"Analysis report · {r['gene']}"); st.caption(f"{r['description'][:260]} | Accession: {r['accession']} | Analyzed DNA: {len(sequence):,} bp")
    if r["application"] == "crispri" and r.get("tss_context"):
        c = r["tss_context"]
        st.success(f"TSS-aware CRISPRi: {c.assembly} {c.chromosome}:{c.tss_coordinate:,} ({c.strand_label}) · canonical transcript {c.transcript_id} · accepted window {CRISPRI_MIN:+d} to {CRISPRI_MAX:+d} bp; preferred {CRISPRI_OPTIMAL_MIN:+d} to {CRISPRI_OPTIMAL_MAX:+d} bp.")
    elif r["source"] != "manual" and r["application"] == "knockout": st.warning("Knockout gene lookup uses representative transcript/cDNA for candidate discovery. Confirm genomic exon/assembly coordinates before experimental use.")
    if not guides:
        st.warning("No candidate guides passed the current filters. For CRISPRi, the annotated TSS window may simply contain no NGG site at the current threshold.")
    else:
        df = pd.DataFrame([x.to_dict() for x in crispri_rows]) if crispri_rows else pd.DataFrame([g.to_dict() for g in guides]); df.insert(0, "Rank", range(1, len(df)+1))
        metrics = st.columns(5); metrics[0].metric("Guides retained", len(guides)); metrics[1].metric("Best heuristic", f"{max(g.score for g in guides):.1f}"); metrics[2].metric("Median GC", f"{median(g.gc_content for g in guides):.1f}%")
        ds = [g.doench_score for g in guides if g.doench_score is not None]; metrics[3].metric("Best Doench RS2", f"{max(ds):.1f}" if ds else "Not available")
        ms = [g.specificity_score for g in guides if g.specificity_score is not None]; metrics[4].metric("Best MIT specificity", f"{max(ms):.1f}" if ms else ("GuideScan2 run" if r["genome_rows"] else "Not screened"))
        tabs = st.tabs(["Ranked guides", "Design landscape", "Guide details", "Off-target screen", "Export", "Methods & limits"])
        with tabs[0]:
            st.markdown("#### Ranked candidate table"); st.dataframe(df, use_container_width=True, hide_index=True, height=min(650, 100 + len(df)*35))
        with tabs[1]:
            if crispri_rows:
                fig = px.scatter(df, x="TSS distance (bp)", y="Heuristic", color="TSS band", size="GC%", hover_data=["Rank", "Spacer (20 nt)", "PAM", "Chromosome", "Genomic start"], title="CRISPRi candidates relative to the annotated TSS"); fig.add_vrect(x0=50, x1=100, opacity=.10, line_width=0, annotation_text="preferred +50..+100"); fig.add_vline(x=0, line_dash="dash", annotation_text="TSS")
            else:
                fig = px.scatter(df, x="Start", y="Heuristic", color="Strand", size="GC%", hover_data=["Rank", "Spacer (20 nt)", "PAM"], title="Candidate activity across the target")
            st.plotly_chart(style_plot(fig, 460), use_container_width=True); hist = px.histogram(df, x="GC%", nbins=12, title="GC-content distribution"); st.plotly_chart(style_plot(hist, 330), use_container_width=True)
        with tabs[2]:
            idx = st.selectbox("Inspect a guide", range(len(guides)), format_func=lambda i:f"#{i+1} · {guides[i].sequence}"); g = guides[idx]
            st.markdown(f'<div class="seq-card">5′—{g.sequence}<span class="pam">{g.pam}</span>—3′</div>', unsafe_allow_html=True)
            q = st.columns(5); q[0].metric("Heuristic", f"{g.score:.1f}"); q[1].metric("Doench RS2", f"{g.doench_score:.1f}" if g.doench_score is not None else "N/A"); q[2].metric("GC", f"{g.gc_content:.1f}%"); q[3].metric("MIT", f"{g.specificity_score:.1f}" if g.specificity_score is not None else "N/A"); q[4].metric("CFD spec.", f"{g.cfd_specificity:.1f}" if g.cfd_specificity is not None else "N/A")
            if crispri_rows:
                cg = crispri_rows[idx]; chips = [("TSS distance", f"{cg.tss_distance:+d} bp"), ("Placement", cg.tss_band), ("Assembly", cg.assembly or "Custom"), ("Genomic locus", f"{cg.chromosome}:{cg.genomic_start}-{cg.genomic_end}" if cg.genomic_start else "Custom sequence")]
                st.markdown('<div class="info-strip">' + ''.join(f'<div class="info-chip"><small>{a}</small><strong>{b}</strong></div>' for a,b in chips) + '</div>', unsafe_allow_html=True)
            st.markdown("#### Validation & quality"); render_validation(g)
            st.markdown("#### Heuristic score breakdown"); breakdown = score_breakdown(g.sequence, g.pam, application=None if g.application == "crispri" else g.application, start=g.start, sequence_length=len(sequence)); rows = [{"Component":k,"Contribution":v} for k,v in breakdown.items() if k != "Final score"]
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True); st.markdown(f'<div class="summary-card"><b>Final heuristic score: {breakdown.get("Final score", "N/A")}</b><br><span>This transparent heuristic is shown separately from Doench Rule Set 2.</span></div>', unsafe_allow_html=True)
            with st.expander("Example BbsI cloning oligos"):
                rev = str(Seq(g.sequence).reverse_complement()); st.markdown(f'<div class="oligo"><b>Forward</b><br>CACCG{g.sequence}</div><div class="oligo"><b>Reverse</b><br>AAAC{rev}C</div>', unsafe_allow_html=True)
        with tabs[3]:
            if reference:
                oi = st.selectbox("Inspect local off-targets", range(len(guides)), key="otguide", format_func=lambda i:f"#{i+1} · {guides[i].sequence}"); rep = analyze_offtargets(guides[oi].sequence, reference, max_mismatches=r["max_mismatches"])
                c = st.columns(5); c[0].metric("MIT specificity", f"{rep.mit_specificity:.1f}"); c[1].metric("CFD specificity", f"{rep.cfd_specificity:.1f}" if rep.cfd_specificity is not None else "N/A"); c[2].metric("Near-matches", len(rep.hits)); c[3].metric("PAM sites", rep.pam_sites_scanned); c[4].metric("Contigs", rep.reference_contigs)
                st.dataframe(pd.DataFrame([h.to_dict() for h in rep.hits]), use_container_width=True, hide_index=True) if rep.hits else st.success("No additional PAM-compatible near-matches within the selected mismatch radius.")
            elif r["genome_rows"]:
                st.markdown("#### GuideScan2 whole-genome results"); st.dataframe(pd.DataFrame(r["genome_rows"]), use_container_width=True, hide_index=True)
            else: st.info("Specificity was not screened for this run. Choose Local reference or Whole genome (GuideScan2) to calculate off-target evidence.")
        with tabs[4]:
            base = safe_name(r["gene"]); csvb = df.to_csv(index=False).encode(); fasta = "\n".join(f">{base}_g{i}|pam={g.pam}|strand={g.strand}\n{g.sequence}" for i,g in enumerate(guides,1)).encode(); js = json.dumps({"metadata":{"target":r["gene"],"organism":r["organism"],"application":r["application"],"accession":r["accession"],"version":APP_VERSION},"guides":df.where(pd.notna(df),None).to_dict("records")},indent=2).encode()
            d = st.columns(3); d[0].download_button("Download CSV", csvb, file_name=f"{base}_v311_guides.csv", use_container_width=True); d[1].download_button("Download FASTA", fasta, file_name=f"{base}_v311_spacers.fasta", use_container_width=True); d[2].download_button("Download JSON", js, file_name=f"{base}_v311_analysis.json", use_container_width=True)
        with tabs[5]:
            cards = st.columns(4); content = [("Discovery","Both-strand SpCas9 20 nt + NGG scanning."),("On-target","Doench Rule Set 2 when the optional provider is installed; transparent heuristic remains separate."),("Specificity","Native MIT/Hsu; optional CFD; local multi-contig or GuideScan2 whole-genome search."),("CRISPRi","Canonical Ensembl TSS, strand-aware genomic DNA, −50..+300 bp filter, +50..+100 preferred band.")]
            for col,(title,text) in zip(cards,content): col.markdown(f'<div class="method"><b>{title}</b><p>{text}</p></div>', unsafe_allow_html=True)
            st.markdown("#### Scientific boundaries"); st.markdown("- **CRISPRi is TSS-aware for gene lookup.** It no longer ranks by the first 400 bp of cDNA.\n- Pasted-sequence CRISPRi requires an explicit TSS position.\n- Knockout gene lookup still begins from representative cDNA; shortlisted guides should be mapped to the intended genomic coding exon/assembly.\n- Whole-genome specificity requires GuideScan2 plus a matching prebuilt genome index.\n- Variant-aware filtering, chromatin state and DNA/RNA bulges require external genome-aware resources.\n- Computational scores prioritize candidates; they do not replace experimental validation.")

st.divider(); st.caption("CRISPR Studio v3.1.1 · research/educational use · MIT License")
