#!/usr/bin/env python3
"""Interactive Streamlit dashboard for CRISPR gRNA design."""

from __future__ import annotations

import json
import re
from io import BytesIO
from statistics import median
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from Bio.Seq import Seq

from grna_designer import (
    GuideRNA,
    analyze_offtargets,
    clean_dna_sequence,
    design_guides,
    fetch_sequence,
    score_breakdown,
    validate_guide,
)


APP_VERSION = "2.0"
MAX_CUSTOM_BP = 50_000
MAX_REFERENCE_BP = 250_000

st.set_page_config(
    page_title="CRISPR Studio | gRNA Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer",
        "Report a bug": "https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer/issues",
        "About": "CRISPR Studio is an open-source SpCas9 guide-design workbench.",
    },
)


# ---------------------------------------------------------------------------
# Theme and presentation
# ---------------------------------------------------------------------------
dark_mode = st.sidebar.toggle("Dark mode", value=True, key="dark_mode")


def inject_theme(is_dark: bool) -> None:
    palette = (
        {
            "app": "#07110f",
            "sidebar": "#0b1714",
            "panel": "#10211d",
            "panel_alt": "#142a25",
            "text": "#eefbf6",
            "muted": "#9ebbb0",
            "border": "#24453b",
            "accent": "#34d399",
            "accent_2": "#22d3ee",
            "shadow": "rgba(0, 0, 0, 0.28)",
        }
        if is_dark
        else {
            "app": "#f5faf8",
            "sidebar": "#edf7f3",
            "panel": "#ffffff",
            "panel_alt": "#f2f8f5",
            "text": "#12201b",
            "muted": "#5f756d",
            "border": "#d5e7df",
            "accent": "#059669",
            "accent_2": "#0891b2",
            "shadow": "rgba(15, 60, 45, 0.10)",
        }
    )
    st.markdown(
        f"""
        <style>
        :root {{
            --app-bg: {palette['app']};
            --sidebar-bg: {palette['sidebar']};
            --panel-bg: {palette['panel']};
            --panel-alt: {palette['panel_alt']};
            --text: {palette['text']};
            --muted: {palette['muted']};
            --border: {palette['border']};
            --accent: {palette['accent']};
            --accent-2: {palette['accent_2']};
        }}
        .stApp {{ background: var(--app-bg); color: var(--text); }}
        [data-testid="stSidebar"] {{
            background: var(--sidebar-bg);
            border-right: 1px solid var(--border);
        }}
        [data-testid="stHeader"] {{ background: rgba(0,0,0,0); }}
        [data-testid="stSidebar"] *:not(svg),
        .stApp label, .stApp p, .stApp li {{ color: var(--text); }}
        .stCaption, [data-testid="stCaptionContainer"] p {{ color: var(--muted) !important; }}
        .block-container {{ max-width: 1460px; padding-top: 1.6rem; padding-bottom: 4rem; }}
        .hero {{
            padding: 2.2rem 2.35rem;
            border: 1px solid var(--border);
            border-radius: 24px;
            background:
                radial-gradient(circle at 92% 12%, rgba(34,211,238,.17), transparent 26%),
                radial-gradient(circle at 8% 100%, rgba(52,211,153,.17), transparent 30%),
                var(--panel-bg);
            box-shadow: 0 18px 55px {palette['shadow']};
            margin-bottom: 1.3rem;
        }}
        .hero-kicker {{
            color: var(--accent) !important;
            font-size: .76rem;
            font-weight: 800;
            letter-spacing: .14em;
            text-transform: uppercase;
            margin-bottom: .7rem;
        }}
        .hero h1 {{
            color: var(--text);
            font-size: clamp(2rem, 4.6vw, 4.15rem);
            letter-spacing: -.045em;
            line-height: .98;
            margin: 0 0 1rem 0;
        }}
        .hero h1 span {{ color: var(--accent); }}
        .hero-copy {{ color: var(--muted) !important; max-width: 760px; font-size: 1.06rem; }}
        .hero-pills {{ display: flex; flex-wrap: wrap; gap: .55rem; margin-top: 1.3rem; }}
        .hero-pills span {{
            color: var(--text);
            background: var(--panel-alt);
            border: 1px solid var(--border);
            border-radius: 999px;
            padding: .42rem .72rem;
            font-size: .78rem;
            font-weight: 650;
        }}
        .sidebar-brand {{
            padding: .45rem .1rem 1rem .1rem;
            border-bottom: 1px solid var(--border);
            margin-bottom: .9rem;
        }}
        .sidebar-brand strong {{ color: var(--text); font-size: 1.17rem; }}
        .sidebar-brand small {{ color: var(--muted); display: block; margin-top: .24rem; }}
        div[data-testid="stMetric"] {{
            background: var(--panel-bg);
            border: 1px solid var(--border);
            padding: 1rem 1.05rem;
            border-radius: 16px;
            box-shadow: 0 8px 24px {palette['shadow']};
        }}
        div[data-testid="stMetricLabel"] p {{ color: var(--muted) !important; }}
        div[data-testid="stMetricValue"] {{ color: var(--text); }}
        [data-testid="stForm"], [data-testid="stExpander"] {{
            background: var(--panel-bg);
            border-color: var(--border) !important;
            border-radius: 16px;
        }}
        [data-baseweb="input"] > div,
        [data-baseweb="textarea"] > div,
        [data-baseweb="select"] > div {{
            background: var(--panel-alt) !important;
            border-color: var(--border) !important;
            color: var(--text) !important;
        }}
        .stApp input, .stApp textarea {{
            background-color: var(--panel-alt) !important;
            color: var(--text) !important;
            caret-color: var(--accent) !important;
        }}
        [data-testid="stDataFrame"] {{
            border: 1px solid var(--border);
            border-radius: 14px;
            overflow: hidden;
        }}
        .stButton > button, .stDownloadButton > button {{
            border-radius: 12px;
            min-height: 2.8rem;
            font-weight: 750;
            border-color: var(--border);
        }}
        .stButton > button[kind="primary"] {{
            background: linear-gradient(115deg, var(--accent), var(--accent-2));
            color: #04120d;
            border: 0;
        }}
        button[data-baseweb="tab"] p {{ font-weight: 700; }}
        .sequence-card {{
            background: var(--panel-alt);
            border: 1px solid var(--border);
            padding: 1rem 1.15rem;
            border-radius: 14px;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            overflow-wrap: anywhere;
            letter-spacing: .035em;
            color: var(--text);
        }}
        .sequence-card mark {{
            background: var(--accent);
            color: #04120d;
            padding: .1rem .25rem;
            border-radius: 5px;
            font-weight: 800;
        }}
        .method-card {{
            height: 100%;
            background: var(--panel-bg);
            border: 1px solid var(--border);
            border-radius: 15px;
            padding: 1rem 1.05rem;
        }}
        .method-card b {{ color: var(--accent); }}
        code {{ color: var(--accent) !important; }}
        hr {{ border-color: var(--border) !important; }}
        @media (max-width: 700px) {{
            .hero {{ padding: 1.45rem; border-radius: 18px; }}
            .block-container {{ padding-left: 1rem; padding-right: 1rem; }}
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )


inject_theme(dark_mode)
plot_template = "plotly_dark" if dark_mode else "plotly_white"
text_color = "#eefbf6" if dark_mode else "#12201b"
grid_color = "#24453b" if dark_mode else "#d5e7df"

st.sidebar.markdown(
    f"""
    <div class="sidebar-brand">
      <strong>CRISPR Studio</strong>
      <small>SpCas9 guide workbench · v{APP_VERSION}</small>
    </div>
    """,
    unsafe_allow_html=True,
)

st.sidebar.subheader("Design settings")
max_guides = st.sidebar.slider("Maximum guides", 5, 50, 20, 5)
min_score = st.sidebar.slider("Minimum activity score", 0, 90, 35, 5)
max_mismatches = st.sidebar.slider(
    "Reference-screen mismatches",
    0,
    4,
    3,
    help="Maximum spacer substitutions allowed at an NGG-compatible reference site.",
)

with st.sidebar.expander("Design model", expanded=False):
    st.markdown(
        """
        - Nuclease: **SpCas9**
        - Spacer: **20 nt**
        - PAM: **NGG**
        - Strands: **both**
        - Preferred GC: **40-70%**
        """
    )

st.sidebar.info(
    "Research-use shortlist only. Confirm the genomic target, exon context, "
    "off-target profile, and experimental controls before ordering guides."
)


# ---------------------------------------------------------------------------
# Cached operations and export helpers
# ---------------------------------------------------------------------------
@st.cache_data(ttl=3600, show_spinner=False)
def cached_fetch(gene_name: str, organism: str, source: str) -> Tuple[str, str, str]:
    return fetch_sequence(gene_name, organism, source)


@st.cache_data(show_spinner=False)
def cached_design(
    sequence: str,
    application: str,
    guide_limit: int,
    score_floor: float,
    reference: Optional[str],
    mismatch_limit: int,
) -> List[GuideRNA]:
    return design_guides(
        sequence,
        application=application,
        max_guides=guide_limit,
        min_score=score_floor,
        genome_context=reference,
        max_mismatches=mismatch_limit,
    )


def read_uploaded_text(uploaded_file: Any) -> str:
    if uploaded_file is None:
        return ""
    payload = uploaded_file.getvalue()
    try:
        return payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        return payload.decode("latin-1")


def safe_filename(label: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", label.strip())
    return cleaned.strip("_") or "crispr_guides"


def validation_frame(guides: List[GuideRNA]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for rank, guide in enumerate(guides, start=1):
        checks = validate_guide(guide)
        rows.append(
            {
                "Rank": rank,
                "Spacer": guide.sequence,
                "Overall": "Pass" if checks.pop("overall_pass") else "Review",
                **checks,
            }
        )
    return pd.DataFrame(rows)


def off_target_frame(guides: List[GuideRNA]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for rank, guide in enumerate(guides, start=1):
        for hit in guide.off_target_details:
            rows.append({"Guide rank": rank, "Guide spacer": guide.sequence, **hit})
    return pd.DataFrame(rows)


def build_excel(
    guide_df: pd.DataFrame,
    guides: List[GuideRNA],
    metadata: Dict[str, Any],
) -> bytes:
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        guide_df.to_excel(writer, index=False, sheet_name="Ranked guides")
        validation_frame(guides).to_excel(writer, index=False, sheet_name="Validation")
        pd.DataFrame(
            {"Field": list(metadata.keys()), "Value": list(metadata.values())}
        ).to_excel(writer, index=False, sheet_name="Metadata")
        hits = off_target_frame(guides)
        if not hits.empty:
            hits.to_excel(writer, index=False, sheet_name="Reference hits")
    return buffer.getvalue()


def style_plot(figure: go.Figure, height: int = 430) -> go.Figure:
    figure.update_layout(
        template=plot_template,
        height=height,
        margin=dict(l=20, r=20, t=55, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color=text_color),
        legend_title_text="",
    )
    figure.update_xaxes(gridcolor=grid_color, zerolinecolor=grid_color)
    figure.update_yaxes(gridcolor=grid_color, zerolinecolor=grid_color)
    return figure


# ---------------------------------------------------------------------------
# Hero and design form
# ---------------------------------------------------------------------------
st.markdown(
    """
    <section class="hero">
      <div class="hero-kicker">Open-source guide design workbench</div>
      <h1>Design sharper <span>CRISPR guides.</span></h1>
      <p class="hero-copy">
        Fetch a representative transcript or paste your own target, rank SpCas9
        candidates on both strands, inspect every score contribution, and screen
        against an optional local reference sequence.
      </p>
      <div class="hero-pills">
        <span>20 nt spacers</span><span>NGG PAM</span><span>Both strands</span>
        <span>Explainable ranking</span><span>CSV · XLSX · FASTA</span>
      </div>
    </section>
    """,
    unsafe_allow_html=True,
)

st.subheader("Design workspace")
input_mode = st.radio(
    "Target input",
    ["Gene lookup", "Paste sequence"],
    horizontal=True,
    label_visibility="collapsed",
)

with st.form("design_form", border=True):
    left, right = st.columns([1.2, 1])
    with left:
        if input_mode == "Gene lookup":
            gene_name = st.text_input(
                "Gene symbol",
                value="TP53",
                placeholder="TP53",
                help="Use an official gene symbol when possible.",
            )
            organism = st.text_input(
                "Organism",
                value="Homo sapiens",
                placeholder="Homo sapiens",
            )
            source = st.radio(
                "Sequence database",
                ["Ensembl", "NCBI"],
                horizontal=True,
                help="Both modes retrieve a representative transcript sequence.",
            ).lower()
            custom_sequence = ""
        else:
            custom_sequence = st.text_area(
                "Target DNA or FASTA",
                height=190,
                placeholder=">target\nATG...",
                help=f"Accepted: DNA, RNA, or FASTA up to {MAX_CUSTOM_BP:,} bp.",
            )
            gene_name = "Custom target"
            organism = "Not specified"
            source = "manual"

    with right:
        application_label = st.radio(
            "Design intent",
            ["Knockout", "Knockdown / CRISPRi"],
            help=(
                "Knockout mildly favors the early sequence region. Knockdown "
                "strongly favors the first 400 bp as a simple 5-prime proxy."
            ),
        )
        application = "knockout" if application_label == "Knockout" else "knockdown"
        screen_reference = st.checkbox(
            "Add local-reference off-target screen",
            help=(
                "Upload or paste a genomic region, contig, or small reference. "
                "This is not a whole-genome aligner."
            ),
        )
        reference_upload = None
        reference_text = ""
        if screen_reference:
            reference_upload = st.file_uploader(
                "Reference FASTA / text",
                type=["fa", "fasta", "fna", "txt"],
            )
            reference_text = st.text_area(
                "Or paste reference sequence",
                height=100,
                placeholder="Paste a reference region here...",
            )

    submitted = st.form_submit_button(
        "Design and analyze guides",
        type="primary",
        width="stretch",
    )


if submitted:
    try:
        uploaded_reference = read_uploaded_text(reference_upload)
        raw_reference = uploaded_reference or reference_text
        reference = clean_dna_sequence(raw_reference) if screen_reference else ""
        if screen_reference and not reference:
            raise ValueError(
                "Reference screening is enabled, but no reference sequence was supplied."
            )
        if len(reference) > MAX_REFERENCE_BP:
            raise ValueError(
                f"The local reference is {len(reference):,} bp. Limit it to "
                f"{MAX_REFERENCE_BP:,} bp so the hosted app remains responsive."
            )

        if input_mode == "Paste sequence":
            sequence = clean_dna_sequence(custom_sequence)
            if len(sequence) < 50:
                raise ValueError("The target sequence must contain at least 50 bp.")
            if len(sequence) > MAX_CUSTOM_BP:
                raise ValueError(
                    f"The target sequence is {len(sequence):,} bp; the app limit is "
                    f"{MAX_CUSTOM_BP:,} bp."
                )
            accession = "CUSTOM"
            description = "User-provided target sequence"
        else:
            if not gene_name.strip() or not organism.strip():
                raise ValueError("Enter both a gene symbol and an organism.")
            with st.spinner(f"Fetching {gene_name.strip()} from {source.upper()}..."):
                accession, description, sequence = cached_fetch(
                    gene_name.strip(),
                    organism.strip(),
                    source,
                )

        with st.spinner("Scanning PAM sites, scoring candidates, and building the report..."):
            guides = cached_design(
                sequence,
                application,
                max_guides,
                float(min_score),
                reference or None,
                max_mismatches,
            )

        st.session_state["analysis"] = {
            "guides": guides,
            "sequence": sequence,
            "reference": reference,
            "accession": accession,
            "description": description,
            "gene": gene_name.strip() or "Custom target",
            "organism": organism.strip(),
            "source": source,
            "application": application,
            "max_mismatches": max_mismatches,
            "min_score": min_score,
        }
        st.success(f"Analysis complete: {len(guides)} guide candidates passed the filter.")
    except Exception as exc:
        st.error(f"Could not complete the analysis: {exc}")
        st.info(
            "Check the gene spelling and organism, try the other database, or use "
            "Paste sequence. Public database services can also be temporarily unavailable."
        )


# ---------------------------------------------------------------------------
# Results dashboard
# ---------------------------------------------------------------------------
if "analysis" in st.session_state:
    result = st.session_state["analysis"]
    guides: List[GuideRNA] = result["guides"]
    sequence: str = result["sequence"]
    reference: str = result["reference"]

    st.divider()
    st.subheader(f"Analysis report · {result['gene']}")
    st.caption(
        f"{result['description'][:240]} | Accession: {result['accession']} | "
        f"Input length: {len(sequence):,} bp"
    )

    if not guides:
        st.warning(
            "No candidates passed the current score threshold. Lower the minimum "
            "activity score in the sidebar or inspect a longer target sequence."
        )
    else:
        guide_df = pd.DataFrame([guide.to_dict() for guide in guides])
        guide_df.insert(0, "Rank", range(1, len(guide_df) + 1))

        metric_columns = st.columns(5)
        metric_columns[0].metric("Guides retained", len(guides))
        metric_columns[1].metric("Best activity", f"{guides[0].score:.1f}")
        metric_columns[2].metric(
            "Median GC",
            f"{median(guide.gc_content for guide in guides):.1f}%",
        )
        metric_columns[3].metric(
            "Strand balance",
            f"{sum(g.strand == '+' for g in guides)}+ / {sum(g.strand == '-' for g in guides)}-",
        )
        if reference:
            screened_scores = [
                guide.specificity_score
                for guide in guides
                if guide.specificity_score is not None
            ]
            best_specificity = max(screened_scores) if screened_scores else 0.0
            metric_columns[4].metric("Best specificity", f"{best_specificity:.1f}")
        else:
            metric_columns[4].metric("Reference screen", "Not run")

        tabs = st.tabs(
            [
                "Ranked guides",
                "Design landscape",
                "Guide details",
                "Off-target screen",
                "Export",
                "Methods & limits",
            ]
        )

        with tabs[0]:
            st.markdown("#### Ranked candidate table")
            visible_df = guide_df.copy()
            if not reference:
                visible_df = visible_df.drop(columns=["Specificity", "Off-target hits"])
            st.dataframe(
                visible_df,
                width="stretch",
                hide_index=True,
                height=min(620, 82 + len(visible_df) * 35),
                column_config={
                    "Rank": st.column_config.NumberColumn(width="small"),
                    "Spacer (20 nt)": st.column_config.TextColumn(width="large"),
                    "Score": st.column_config.ProgressColumn(
                        min_value=0,
                        max_value=100,
                        format="%.1f",
                    ),
                    "Specificity": st.column_config.ProgressColumn(
                        min_value=0,
                        max_value=100,
                        format="%.1f",
                    ),
                    "GC%": st.column_config.NumberColumn(format="%.1f%%"),
                },
            )
            with st.expander("Target sequence preview"):
                preview = sequence[:600]
                st.code(preview + ("..." if len(sequence) > 600 else ""), language="text")
                st.caption(
                    f"Showing {min(len(sequence), 600):,} of {len(sequence):,} bp. "
                    f"Ambiguous bases: {sequence.count('N'):,}."
                )

        with tabs[1]:
            chart_left, chart_right = st.columns([1.35, 1])
            with chart_left:
                landscape = px.scatter(
                    guide_df,
                    x="Start",
                    y="Score",
                    color="Strand",
                    size="GC%",
                    size_max=22,
                    color_discrete_map={"+": "#34d399", "-": "#22d3ee"},
                    hover_data=["Rank", "Spacer (20 nt)", "PAM", "GC%"],
                    title="Candidate activity across the target",
                )
                landscape.add_hrect(
                    y0=70,
                    y1=100,
                    fillcolor="#34d399",
                    opacity=0.07,
                    line_width=0,
                    annotation_text="high-score zone",
                    annotation_position="top left",
                )
                st.plotly_chart(
                    style_plot(landscape, 450),
                    width="stretch",
                    config={"displaylogo": False},
                )
            with chart_right:
                gc_chart = px.histogram(
                    guide_df,
                    x="GC%",
                    nbins=12,
                    color_discrete_sequence=["#34d399"],
                    title="GC-content distribution",
                )
                gc_chart.add_vrect(x0=40, x1=70, fillcolor="#22d3ee", opacity=0.08, line_width=0)
                st.plotly_chart(
                    style_plot(gc_chart, 450),
                    width="stretch",
                    config={"displaylogo": False},
                )

            strand_counts = (
                guide_df.groupby("Strand", as_index=False).size().rename(columns={"size": "Guides"})
            )
            strand_chart = px.bar(
                strand_counts,
                x="Strand",
                y="Guides",
                color="Strand",
                color_discrete_map={"+": "#34d399", "-": "#22d3ee"},
                title="Candidates by strand",
            )
            st.plotly_chart(
                style_plot(strand_chart, 330),
                width="stretch",
                config={"displaylogo": False},
            )

        with tabs[2]:
            selected_index = st.selectbox(
                "Inspect a guide",
                list(range(len(guides))),
                format_func=lambda index: f"#{index + 1} · {guides[index].sequence}",
            )
            selected = guides[selected_index]
            selected_rank = selected_index + 1
            checks = validate_guide(selected)

            details_left, details_right = st.columns([1.12, 1])
            with details_left:
                st.markdown(f"#### Guide #{selected_rank}")
                st.markdown(
                    f"""
                    <div class="sequence-card">
                      5′—{selected.sequence}<mark>{selected.pam}</mark>—3′
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
                st.caption("The highlighted 3 nt are the NGG PAM; the 20 nt spacer is ordered without the PAM.")
                detail_metrics = st.columns(4)
                detail_metrics[0].metric("Activity", f"{selected.score:.1f}")
                detail_metrics[1].metric("GC", f"{selected.gc_content:.1f}%")
                detail_metrics[2].metric("Strand", selected.strand)
                detail_metrics[3].metric("Coordinates", f"{selected.start + 1}-{selected.end}")

                st.markdown("##### Validation checklist")
                label_map = {
                    "length_ok": "20 nt spacer",
                    "pam_ok": "NGG PAM",
                    "gc_in_preferred_range": "Preferred GC range",
                    "gc_in_acceptable_range": "Acceptable GC range",
                    "no_extreme_homopolymer": "No long homopolymer",
                    "score_above_threshold": "Activity score ≥ 40",
                    "no_poly_t": "No poly-T terminator",
                    "specificity_screened": "Reference screened",
                    "specificity_ok": "Specificity acceptable",
                }
                check_columns = st.columns(2)
                visible_checks = [item for item in checks.items() if item[0] != "overall_pass"]
                for index, (check_name, passed) in enumerate(visible_checks):
                    label = label_map.get(check_name, check_name.replace("_", " ").title())
                    if check_name == "specificity_screened" and not reference:
                        check_columns[index % 2].info(f"○ {label}: not run")
                    elif passed:
                        check_columns[index % 2].success(f"✓ {label}")
                    else:
                        check_columns[index % 2].warning(f"! {label}")

                context_start = max(0, selected.start - 20)
                context_end = min(len(sequence), selected.end + 20)
                st.markdown("##### Sequence context")
                st.code(sequence[context_start:context_end], language="text")
                st.caption(f"Target coordinates shown with 20 bp flanks where available: {context_start + 1}-{context_end}.")

                with st.expander("Example BbsI cloning oligos"):
                    reverse_spacer = str(Seq(selected.sequence).reverse_complement())
                    st.code(
                        f"Forward: CACCG{selected.sequence}\nReverse: AAAC{reverse_spacer}C",
                        language="text",
                    )
                    st.caption("Overhangs are vector-specific. Confirm your plasmid protocol before ordering.")

            with details_right:
                gauge = go.Figure(
                    go.Indicator(
                        mode="gauge+number",
                        value=selected.score,
                        number={"suffix": "/100"},
                        title={"text": "Heuristic activity score"},
                        gauge={
                            "axis": {"range": [0, 100]},
                            "bar": {"color": "#34d399"},
                            "bgcolor": "rgba(0,0,0,0)",
                            "steps": [
                                {"range": [0, 40], "color": "rgba(244,63,94,.18)"},
                                {"range": [40, 70], "color": "rgba(245,158,11,.16)"},
                                {"range": [70, 100], "color": "rgba(52,211,153,.12)"},
                            ],
                        },
                    )
                )
                st.plotly_chart(
                    style_plot(gauge, 315),
                    width="stretch",
                    config={"displaylogo": False},
                )

                components = score_breakdown(
                    selected.sequence,
                    selected.pam,
                    application=selected.application,
                    start=selected.start,
                    sequence_length=len(sequence),
                )
                component_df = pd.DataFrame(
                    [
                        {"Feature": feature, "Contribution": value}
                        for feature, value in components.items()
                        if feature not in {"Baseline", "Final score"} and value != 0
                    ]
                )
                if not component_df.empty:
                    contribution_chart = px.bar(
                        component_df,
                        x="Contribution",
                        y="Feature",
                        orientation="h",
                        color="Contribution",
                        color_continuous_scale=["#fb7185", "#f1f5f9", "#34d399"],
                        color_continuous_midpoint=0,
                        title="Why this guide received its score",
                    )
                    contribution_chart.update_layout(coloraxis_showscale=False)
                    st.plotly_chart(
                        style_plot(contribution_chart, 360),
                        width="stretch",
                        config={"displaylogo": False},
                    )
                if selected.notes:
                    st.caption("Flags and preferences: " + " · ".join(selected.notes))

        with tabs[3]:
            st.markdown("#### PAM-aware local-reference similarity screen")
            if not reference:
                st.info(
                    "No reference was supplied for this run. Enable the local-reference "
                    "screen in the design form, then upload or paste a FASTA sequence."
                )
                st.markdown(
                    "This screen examines NGG-compatible sites on both strands and compares "
                    "their 20 nt spacers. For a real experiment, use a genome-indexed tool "
                    "such as CRISPOR, CHOPCHOP, Cas-OFFinder, or GuideScan."
                )
            else:
                screen_df = guide_df[
                    ["Rank", "Spacer (20 nt)", "Specificity", "Off-target hits"]
                ].copy()
                st.dataframe(
                    screen_df,
                    width="stretch",
                    hide_index=True,
                    column_config={
                        "Specificity": st.column_config.ProgressColumn(
                            min_value=0,
                            max_value=100,
                            format="%.1f",
                        )
                    },
                )
                screen_index = st.selectbox(
                    "Inspect reference hits",
                    list(range(len(guides))),
                    key="offtarget_guide",
                    format_func=lambda index: f"#{index + 1} · {guides[index].sequence}",
                )
                screened_guide = guides[screen_index]
                report = analyze_offtargets(
                    screened_guide.sequence,
                    reference,
                    max_mismatches=result["max_mismatches"],
                )
                screen_metrics = st.columns(4)
                screen_metrics[0].metric("Specificity", f"{report.specificity_score:.1f}")
                screen_metrics[1].metric("Near-matches", len(report.hits))
                screen_metrics[2].metric("PAM sites scanned", f"{report.pam_sites_scanned:,}")
                screen_metrics[3].metric("Reference length", f"{len(reference):,} bp")

                if report.hits:
                    hits_df = pd.DataFrame([hit.to_dict() for hit in report.hits])
                    st.dataframe(hits_df, width="stretch", hide_index=True)
                    risk_counts = (
                        hits_df.groupby("Risk", as_index=False)
                        .size()
                        .rename(columns={"size": "Hits"})
                    )
                    risk_chart = px.bar(
                        risk_counts,
                        x="Risk",
                        y="Hits",
                        color="Risk",
                        color_discrete_map={
                            "Critical": "#e11d48",
                            "High": "#f97316",
                            "Moderate": "#f59e0b",
                            "Low": "#34d399",
                        },
                        title="Near-matches by risk tier",
                    )
                    st.plotly_chart(
                        style_plot(risk_chart, 330),
                        width="stretch",
                        config={"displaylogo": False},
                    )
                else:
                    st.success(
                        "No additional PAM-compatible sites were found within the selected "
                        "mismatch limit in this reference."
                    )
                if report.on_target_excluded:
                    st.caption("One exact match was treated as the intended target and excluded from the hit count.")
                st.warning(
                    "A clean local-reference result does not prove genome-wide specificity. "
                    "This lightweight screen does not model bulges, chromatin, alternate PAMs, or variants."
                )

        with tabs[4]:
            st.markdown("#### Download a reproducible analysis bundle")
            metadata = {
                "Target": result["gene"],
                "Organism": result["organism"],
                "Source": result["source"],
                "Application": result["application"],
                "Accession": result["accession"],
                "Target length (bp)": len(sequence),
                "Minimum activity score": result["min_score"],
                "Reference screened": bool(reference),
                "Reference length (bp)": len(reference),
                "Mismatch limit": result["max_mismatches"],
                "App version": APP_VERSION,
            }
            base_name = safe_filename(result["gene"])
            csv_bytes = guide_df.to_csv(index=False).encode("utf-8")
            excel_bytes = build_excel(guide_df, guides, metadata)
            fasta_text = "\n".join(
                f">{base_name}_guide_{rank}|score={guide.score:.1f}|pam={guide.pam}|strand={guide.strand}\n{guide.sequence}"
                for rank, guide in enumerate(guides, start=1)
            )
            json_records = (
                guide_df.astype(object)
                .where(pd.notna(guide_df), None)
                .to_dict(orient="records")
            )
            json_bundle = json.dumps(
                {"metadata": metadata, "guides": json_records},
                indent=2,
                allow_nan=False,
            ).encode("utf-8")

            download_columns = st.columns(4)
            download_columns[0].download_button(
                "Download CSV",
                csv_bytes,
                file_name=f"{base_name}_guides.csv",
                mime="text/csv",
                width="stretch",
            )
            download_columns[1].download_button(
                "Download Excel",
                excel_bytes,
                file_name=f"{base_name}_analysis.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                width="stretch",
            )
            download_columns[2].download_button(
                "Download FASTA",
                fasta_text.encode("utf-8"),
                file_name=f"{base_name}_spacers.fasta",
                mime="text/plain",
                width="stretch",
            )
            download_columns[3].download_button(
                "Download JSON",
                json_bundle,
                file_name=f"{base_name}_analysis.json",
                mime="application/json",
                width="stretch",
            )
            st.caption(
                "The Excel workbook includes ranked guides, validation checks, metadata, "
                "and reference hits when screening was enabled."
            )

        with tabs[5]:
            st.markdown("#### What the dashboard calculates")
            method_columns = st.columns(4)
            cards = [
                (
                    "1 · Discovery",
                    "Scans both strands for a 20 nt spacer adjacent to an SpCas9 <b>NGG PAM</b>.",
                ),
                (
                    "2 · Activity rank",
                    "Combines GC window, position-specific bases, PAM context, homopolymers, and self-complementarity into an <b>explainable heuristic</b>.",
                ),
                (
                    "3 · Intent bias",
                    "Knockout mildly favors the early sequence; knockdown favors the first 400 bp as a <b>simple 5-prime proxy</b>.",
                ),
                (
                    "4 · Specificity screen",
                    "Compares spacers at NGG-compatible sites on both strands of a <b>user-supplied local reference</b>.",
                ),
            ]
            for column, (title, copy) in zip(method_columns, cards):
                column.markdown(
                    f'<div class="method-card"><b>{title}</b><p>{copy}</p></div>',
                    unsafe_allow_html=True,
                )

            st.markdown("#### Important scientific limitations")
            st.warning(
                "Database lookup returns a representative transcript/cDNA. A candidate can "
                "cross an exon-exon junction or miss isoform-specific and genomic context. "
                "Map every shortlisted spacer back to the intended genome build and coding exon."
            )
            st.markdown(
                """
                - The activity value is a literature-inspired ranking heuristic, **not** the trained Rule Set 2/Azimuth prediction.
                - The local screen is not a whole-genome alignment and does not evaluate DNA/RNA bulges, chromatin accessibility, SNPs, or non-NGG PAMs.
                - CRISPRi design normally needs an experimentally relevant TSS window; a transcript 5-prime position is only an approximation.
                - Use at least two independent guides, appropriate negative controls, and orthogonal validation.

                Useful primary references: [Jinek et al. 2012](https://doi.org/10.1126/science.1225829),
                [Hsu et al. 2013](https://doi.org/10.1038/nbt.2647),
                [Doench et al. 2016](https://doi.org/10.1038/nbt.3437), and
                [Moreno-Mateos et al. 2015](https://doi.org/10.1038/nmeth.3543).
                """
            )

st.divider()
st.caption(
    "CRISPR Studio is an educational research tool. It does not provide clinical advice "
    "or replace genome-wide specificity analysis and experimental validation."
)

