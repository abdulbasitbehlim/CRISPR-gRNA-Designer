# User Guide

A plain-language walkthrough of what this tool does, why it works the way
it does, and how to trust the results.

---

## What problem does this tool solve?

You want to edit a gene with CRISPR. You know the **gene name** (say,
`TP53`) and the **organism** (say, human). You don't want to dig through
NCBI or Ensembl yourself to find the DNA sequence, manually scan it for
PAM sites, and score each candidate by hand.

This tool automates that whole workflow in a simple web dashboard:

1. You type the gene name and organism.
2. The tool downloads the gene/transcript sequence automatically.
3. You choose **Knockout** (destroy the gene) or **Knockdown** (turn the
   gene down).
4. The tool returns a ranked list of candidate guide RNAs.

---

## What is a guide RNA (gRNA), in simple terms?

CRISPR-Cas9 works like a pair of molecular scissors. The **guide RNA**
is the address label that tells the scissors exactly where to cut.

- The address is a 20-letter DNA sequence (the **spacer**).
- Right next to it must be a short 3-letter code called the **PAM** —
  for the common SpCas9 enzyme this is `NGG`, where `N` can be any base.

The tool searches your gene for every location that has a 20-letter
sequence followed by `NGG` (on either strand) and ranks the results.

---

## Knockout vs. Knockdown

| Goal | What happens | Where the guide should sit |
|------|--------------|---------------------------|
| **Knockout** | Cut the DNA → the cell repairs it → this usually creates a frameshift → the gene is permanently broken | Ideally in an early exon, so the whole protein is disrupted |
| **Knockdown** | Typically uses a "dead" Cas9 that sits on the DNA and blocks transcription (CRISPRi) | Close to the start of the gene (near the promoter / transcription start site) |

The tool applies this as a ranking bias, not a hard rule:
- **Knockout** → mild preference for guides in the first ~30 % of the sequence.
- **Knockdown** → strong preference for guides in the first ~400 base pairs.

---

## How the scoring works

Each candidate guide gets a score from 0–100. Higher means better
predicted on-target activity. The score is built from published findings
(mainly from the Doench lab and related work):

- GC content between 40–70 % is favored.
- The last base of the 20-nt spacer is favored to be a G.
- Long runs of the same base (especially `GGGG` or `TTTT`) are penalized.
- A few specific positions (roughly positions 16–20) carry small
  preferences based on the literature.
- Guides that look like they might fold back on themselves get a small
  penalty.

This is a transparent, rule-based score — not a black-box model — so
it's easy to inspect, question, and improve.

**Important:** the score only predicts *on-target* activity. It does
**not** search the whole genome for off-targets. Always take your top
3–5 guides and check them on [CRISPOR](https://crispor.tefor.net) or
[CHOPCHOP](https://chopchop.cbu.uib.no) before ordering anything.

---

## Running the tool

### Easiest – one public link (works on phone, tablet, any computer)

Deploy the app once to Streamlit Community Cloud (free). You get a permanent
URL such as `https://crispr-grna-designer.streamlit.app`. Anyone can open
that single link in any modern browser — no installation required.

See the **“One-click deploy to Streamlit Cloud”** section in `README.md`
for the exact 5-step process.

### Local run (Windows / macOS / Linux)

```bash
git clone https://github.com/abdulbasitbehlim/crispr-grna-designer.git
cd crispr-grna-designer
python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Your browser will open automatically. Try gene `TP53`, organism
`Homo sapiens`, and click **Design gRNAs**.

You can also paste your own DNA sequence if the gene isn't available in
the public databases, using the "Or paste your own DNA sequence" panel.

---

## How to judge whether a result is good

1. Check the **Score** column — aim for guides scoring above ~50–60 where possible.
2. Check the **Notes** column — "OK" is best; "Homopolymer" or "Ends with T" are warnings worth a second look.
3. Open the validation panel for the top 3 guides — they should show **PASS**.
4. Copy the top guide's sequence + PAM into CRISPOR or CHOPCHOP and confirm the ranking is broadly consistent.
5. For real experiments, always test 2–3 guides in parallel — biology is never 100 % predictable from sequence alone.

---

## Project files

| File | What it is |
|------|------------|
| `app.py` | The Streamlit web dashboard |
| `grna_designer.py` | Sequence fetching, guide design, and scoring logic |
| `requirements.txt` | Runtime Python dependencies |
| `requirements-dev.txt` | Runtime + test dependencies |
| `test_grna_designer.py` | Unit tests for the core logic |
| `ARCHITECTURE.md` | Technical design notes and scientific references |
| `USER_GUIDE.md` | This document |
| `README.md` | Project overview and quick start |

---

## About similar tools

There are already excellent CRISPR design tools — CRISPOR, CHOPCHOP,
Benchling, GuideScan, crisprVerse (R), and others. This project is
intentionally smaller in scope: fully open-source, pure Python, and
designed to go directly from a gene symbol to ranked guides on your own
machine, without an account or a heavyweight installation.

---

## Final advice

- Treat the top 3–5 guides as strong candidates, not absolute answers.
- Always run an off-target check on a genome-wide tool before ordering oligos.
- For publication, cite the original scoring papers (Doench 2016, Hsu
  2013, Moreno-Mateos 2015) and the databases you used (NCBI, Ensembl).
