# CFD data provenance

`data/cfd_data.json` is copied unchanged from USDA-ARS-GBRU/GuideMaker,
commit `f5aac7743aa94d0a9ee4ca2f0fe44d000a6d2e99`, under its CC0 dedication
(`GuideMaker_CC0.txt`). Source:
https://github.com/USDA-ARS-GBRU/GuideMaker/blob/f5aac7743aa94d0a9ee4ca2f0fe44d000a6d2e99/guidemaker/data/cfd_data.json

These are the mismatch/PAM weights for Doench et al. (2016),
https://doi.org/10.1038/nbt.3437. CRISPR Studio implements the product formula
independently in `models.py`. Local enumeration searches NGG sites only;
including the full PAM weight table does not expand that search scope.
