# 🧬 CRISPR gRNA Designer

[![Tests](https://github.com/abdulbasitbehlim/crispr-grna-designer/actions/workflows/tests.yml/badge.svg)](https://github.com/abdulbasitbehlim/crispr-grna-designer/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![Streamlit](https://img.shields.io/badge/built%20with-Streamlit-FF4B4B.svg)](https://streamlit.io/)
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io)

A lightweight, open-source Python tool that turns a **gene name + organism** into a ranked list of SpCas9 guide RNAs for **Knockout** or **Knockdown** experiments — no manual sequence lookup required.

```
Gene name + Organism  →  NCBI / Ensembl fetch  →  gRNA design & scoring  →  ranked table + download
```

> **Works on every OS (Windows / macOS / Linux) and on mobile** — just open the live web link in any modern browser. No installation required for the hosted version.

## 🚀 One-click live demo (recommended)

**https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/**

Anyone (including you on phone/tablet) can open that single link and use the tool immediately. No Python, no install, no account.

---

## Table of contents

- [Features](#features)
- [One-click deploy to Streamlit Cloud](#one-click-deploy-to-streamlit-cloud)
- [Upload to your GitHub profile](#upload-to-your-github-profile)
- [Local run (any OS)](#local-run-any-os)
- [Usage](#usage)
- [Project structure](#project-structure)
- [How guides are scored](#how-guides-are-scored)
- [Validation & limitations](#validation--limitations)
- [Running the tests](#running-the-tests)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [Scientific references](#scientific-references)
- [License](#license)

## Features

- 🔎 Automatic sequence retrieval from **NCBI** (Entrez) or the **Ensembl** REST API
- ✂️ SpCas9 (NGG) guide design scanned on both DNA strands
- 🎯 Application-aware ranking:
  - **Knockout** — mild preference for guides in the early region of the sequence
  - **Knockdown** — strong preference for 5′-proximal guides (near the TSS)
- 📊 Transparent, explainable on-target score inspired by Doench / Hsu / Moreno-Mateos rules
- ✅ Built-in validation checklist for each candidate guide
- 📥 CSV / Excel export
- 🧪 Pure Python core with unit tests — easy to embed, extend, or audit
- 📄 Paste-your-own-sequence mode for organisms not in public databases
- 📱 Fully usable on mobile browsers (responsive Streamlit layout)


## Local run (any OS)

Works on Windows, macOS, Linux, and even in cloud IDEs.

```bash
git clone https://github.com/abdulbasitbehlim/crispr-grna-designer.git
cd crispr-grna-designer
python -m venv venv

# Activate
#   macOS / Linux:
source venv/bin/activate
#   Windows (Command Prompt):
venv\Scripts\activate
#   Windows (PowerShell):
venv\Scripts\Activate.ps1

pip install -r requirements.txt
streamlit run app.py
```

Then open the URL shown in the terminal (usually http://localhost:8501).  
The same URL works on phones/tablets on the same Wi-Fi if you use `--server.address 0.0.0.0`.

## Usage

### Web dashboard (local or cloud)

1. Enter a gene symbol (e.g. `TP53`) and organism (e.g. `Homo sapiens`).
2. Choose the sequence source (NCBI or Ensembl) and the application (Knockout or Knockdown).
3. Click **Design gRNAs**.
4. Download results as CSV or Excel.

### Python API

```python
from grna_designer import design_from_gene

acc, desc, seq, guides = design_from_gene(
    "BRCA1", "Homo sapiens",
    application="knockout",
    source="ensembl",
    max_guides=10,
)

for g in guides:
    print(f"{g.sequence}  PAM={g.pam}  score={g.score:.1f}  strand={g.strand}")
```

You can also design guides from a sequence you already have:

```python
from grna_designer import design_guides, validate_guide

guides = design_guides(my_sequence, application="knockdown", max_guides=15)
for g in guides[:3]:
    print(g.to_dict(), validate_guide(g))
```

## Project structure

```
crispr-grna-designer/
├── app.py                  # Streamlit dashboard (entry point for Streamlit Cloud)
├── grna_designer.py        # Core library: fetch, design, score, validate
├── test_grna_designer.py   # Unit tests (no network required)
├── requirements.txt        # Runtime dependencies (used by Streamlit Cloud)
├── requirements-dev.txt    # Runtime + test dependencies
├── USER_GUIDE.md           # Plain-language walkthrough
├── ARCHITECTURE.md         # Design notes & scientific references
├── CONTRIBUTING.md         # Dev setup & PR guide
├── LICENSE
├── README.md
├── .gitignore
├── .streamlit/
│   └── config.toml         # Theme & server settings (mobile-friendly)
└── .github/
    └── workflows/
        └── tests.yml       # CI: runs the test suite on every push/PR
```

## How guides are scored

Each candidate spacer gets a transparent 0–100 score built from literature-derived features: GC content window, base preferences at specific positions (notably position 20), PAM identity, homopolymer penalties, and a simple self-complementarity check — plus an application-specific positional bias (early region for Knockout, 5′-proximal for Knockdown). Full details and citations are in [`ARCHITECTURE.md`](ARCHITECTURE.md).

## Validation & limitations

1. Every guide is checked with an internal `validate_guide()` function (length, PAM, GC range, homopolymers, score threshold).
2. **Always** cross-check top candidates on [CRISPOR](https://crispor.tefor.net) or [CHOPCHOP](https://chopchop.cbu.uib.no) — this tool does **not** perform full genome-wide off-target search by default.
3. Confirm experimentally (ICE / TIDE / NGS, Western blot, or RT-qPCR) before drawing conclusions.

Other known limitations are documented in [`ARCHITECTURE.md`](ARCHITECTURE.md#6-known-limitations).

## Running the tests

```bash
pip install -r requirements-dev.txt
pytest -v
```

The test suite covers PAM scanning, scoring, validation and off-target logic using synthetic sequences (no network access needed).

## Roadmap

- [ ] Integrate a proper ML on-target model (Rule Set 2 / Azimuth-style)
- [ ] Optional local off-target module via Cas-OFFinder / Bowtie2
- [ ] Support for Cas12a (TTTV PAM), base editors, and prime-editing pegRNAs
- [ ] Batch mode: gene list → full Excel guide library
- [ ] Docker image for one-command deployment

See [`ARCHITECTURE.md`](ARCHITECTURE.md#7-possible-extensions) for more detail.

## Contributing

Contributions are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to set up a dev environment, run tests, and submit a pull request.

## Scientific references

- Jinek et al., *Science* 2012 — Cas9 mechanism and PAM requirement
- Doench et al., *Nature Biotechnology* 2014 & 2016 — on-target scoring rules
- Hsu et al., *Nature Biotechnology* 2013 — off-target / seed-region considerations
- Moreno-Mateos et al., *Nature Methods* 2015 — CRISPRscan on-target features
- Community best practices from [CRISPOR](https://crispor.tefor.net) and [CHOPCHOP](https://chopchop.cbu.uib.no)

## License

Released under the [MIT License](LICENSE). If you use this tool in published research, please cite the scoring papers above and the public databases (NCBI, Ensembl) you queried.
