# Contributing to CRISPR Studio

Contributions are welcome for bug fixes, scientific improvements, UI/UX, documentation, testing and optional integrations.

## Current architecture

The application is split into focused modules:

- `app.py` — Streamlit presentation, state, plots and exports
- `grna_designer.py` — SpCas9 discovery, heuristic, MIT specificity, local screening and GuideScan2 adapter
- `crispri.py` — Ensembl TSS retrieval and TSS-aware CRISPRi
- `accession_lookup.py` — NCBI / Ensembl accession resolution

Please keep scientific logic out of the Streamlit layer where practical.

## Getting started

```bash
git clone https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer.git
cd CRISPR-gRNA-Designer
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements-dev.txt
```

## Running tests

```bash
pytest -v
```

GitHub Actions runs the suite on Python 3.10, 3.11 and 3.12.

Network-dependent API calls should be isolated so core unit tests remain deterministic. Streamlit smoke tests should continue to use local example sequence input rather than public services.

## Making a change

1. Create a feature branch.
2. Keep changes scoped to the appropriate module.
3. Add or update tests.
4. Run `pytest -v` locally.
5. Update README/USER_GUIDE/ARCHITECTURE/CHANGELOG when public behavior changes.
6. Open a pull request explaining both the code change and its scientific/UX motivation.

## Scientific contribution rules

- Do not label a custom heuristic as an established published model.
- Keep on-target activity and off-target specificity metrics separate.
- Document model/provider/version provenance.
- Do not present local-reference screening as whole-genome analysis.
- Do not infer a TSS from arbitrary accession records without validated annotation.
- Preserve assembly/transcript provenance whenever genomic coordinates are reported.

## Optional external backends

Large or specialized tools such as GuideScan2 should remain optional integrations. Do not commit whole-genome indexes or other large reference assets to the repository.

## Reporting issues

Please include:

- expected behavior;
- actual behavior;
- app version;
- input mode used (Gene lookup / Accession ID / Paste sequence);
- design intent (Knockout / CRISPRi);
- specificity mode;
- browser/OS or Python version;
- a non-sensitive reproducible example when possible.

Do not post confidential, patient-identifiable, proprietary or otherwise sensitive sequence data in public issues.
