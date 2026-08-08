# Contributing

Contributions are welcome — bug fixes, new scoring features, support for
additional PAMs/enzymes, or documentation improvements.

## Getting started

```bash
git clone https://github.com/abdulbasitbehlim/crispr-grna-designer.git
cd crispr-grna-designer
python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate
pip install -r requirements-dev.txt
```

## Running tests

```bash
pytest -v
```

The test suite covers the pure sequence-processing logic (PAM scanning,
scoring, validation) and does not require network access.

## Making a change

1. Fork the repository and create a branch: `git checkout -b feature/my-change`
2. Make your change, and add or update tests in `test_grna_designer.py`.
3. Run `pytest -v` and confirm everything passes.
4. Open a pull request describing the change and its motivation.

## Reporting issues

Please open a GitHub issue with:
- What you expected to happen
- What actually happened
- The gene/organism (or sequence) that reproduces the problem, if applicable

## Scope note

Off-target genome-wide search is intentionally out of scope for the core
tool (see `ARCHITECTURE.md`). If you'd like to contribute an
optional off-target module (e.g. via Cas-OFFinder or Bowtie2), please
open an issue first to discuss the approach.
