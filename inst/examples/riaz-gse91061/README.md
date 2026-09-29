# Riaz GSE91061: prepared RNA example

Use the online walkthrough for the cohort, study design, prepared-input
commands, reconstruction route and selected integrated FULL report:

```r
browseURL("https://olmedalab.org/lisaR/reader/riaz-gse91061-worked-example.html?lang=en")
```

The prepared tutorial runs two DE inputs, PRE and ON, and one descriptive
LISA contrast, ON minus PRE. STANDARD and FULL use this same selection.
The original bundle also retains the other upstream study results for
reference; the tutorial does not run them.

Obtain the reviewed archive with the public installer:

```r
bundle <- lisaR::install_lisa_example_bundle(
  "riaz-gse91061", destination = "bundles/riaz"
)
bundle$next_command
```

The returned command prepares a new project without fitting DESeq2 or running
lisaR. The guide separates this step from the later standard run. The archive
is available from the public `DBM-OlmedaLab/lisaR-example-data` release; a verified
cached archive can be reused offline. Scientific dictionaries and the hierarchy
are bundled with lisaR. MSigDB memberships are obtained separately through the
resource guide. Original GEO inputs and the full upstream reconstruction are
not part of package installation.
