# Genomic verification examples

NCBI genomic contexts retrieved 19 September 2026. ACTB was retrieved again after
completing the fixes. GAPDH is the saved retrieval from the earlier review in the
same work session. These are public reference sequences, not patient data.

- ACTB: NC_000007.14:5527144-5530605; representative NP_001092.1.
- GAPDH: NC_000012.12:6534513-6538375; representative NP_001276674.1.

`test_saved_ncbi_genomic_examples` verifies genomic target sequence, coding cut
position and successful shortlisting without requiring network access. These
fixtures test implementation, not experimental editing or genome-wide specificity.
