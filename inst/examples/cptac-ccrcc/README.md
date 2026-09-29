# CPTAC ccRCC prepared proteomic example

The worked guide covers 80 matched tumour/NAT pairs, the complete prepared
limma result, protein-to-gene coverage and the ordinary lisaR run:

```r
browseURL("https://olmedalab.org/lisaR/reader/cptac-ccrcc-proteomics.html?lang=en")
```

Obtain the six-file derived-input bundle through the public entry:

```r
bundle <- lisaR::install_lisa_example_bundle(
  "cptac-ccrcc", destination = "bundles/cptac"
)
bundle$next_command
```

The returned POSIX command enables `LISAR_CPTAC_CCRCC_PREPARE_ONLY=1` and
prepares a new project without running lisaR. The vignette supplies the R
alternative and separates preparation from execution. The generated study
uses `global proteomic`, one DE, no contrast and the complete 6,482-row input;
positive `tumor_minus_nat` means higher tumour abundance. It does not refit
limma, filter by significance or retrieve an original abundance matrix.

The versioned archive is publicly available from `DBM-OlmedaLab/lisaR-example-data`.
A verified cached archive can be reused offline. The external TERM2GENE
memberships must already be registered. Source data: Clark DJ et al., Cell
2019;179:964-983, CPTAC ccRCC PDC000127. The bundle contains derived inputs;
access to original proteomics and clinical material remains separate.
