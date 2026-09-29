library(lisaR)

# Use a short, writable path. Change C:/lisa if needed.
root <- if (.Platform$OS.type == "windows") "C:/lisa" else "~/lisa"
root <- path.expand(root)
dir.create(root, recursive = TRUE, showWarnings = FALSE)

# Prepare the human gene-set membership resource once.
prepare_lisa_msigdb_resource(accept_terms = TRUE)

# Download the prepared Riaz data to a new directory.
riaz <- install_lisa_example_bundle(
  "riaz-gse91061",
  destination = file.path(root, "riaz-data")
)

# Analyse ON only: CR/PR responders versus PD during treatment.
result <- run_lisa(
  de = list(responders_vs_pd_on = file.path(
    riaz$bundle, "outputs", "lisa_inputs", "responders_vs_pd_on.tsv"
  )),
  output_dir = file.path(root, "riaz-on"),
  species = "Homo sapiens",
  method = "DESeq2",
  positive_direction =
    "Higher in CR/PR responders than in PD during treatment",
  mode = "standard"
)

stopifnot(identical(result$gate, "PASS"))
browseURL(result$report_index)
