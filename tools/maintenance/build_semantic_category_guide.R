#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Build the verified semantic-category catalog and its static vignette.
# Classification dictionaries are bundled; this catalogue emits category
# descriptions, aggregate counts and selected audited member identifiers.

parse_args <- function(args) {
  allowed <- c("core", "expanded", "category-map", "package-root", "audit")
  if (length(args) != 2L * length(allowed) ||
      any(!startsWith(args[seq.int(1L, length(args), 2L)], "--"))) {
    stop(
      paste(
        "usage: build_semantic_category_guide.R",
        "--core FILE --expanded FILE --category-map FILE",
        "--package-root DIR --audit FILE"
      ),
      call. = FALSE
    )
  }
  keys <- sub("^--", "", args[seq.int(1L, length(args), 2L)])
  if (anyDuplicated(keys) || !setequal(keys, allowed)) {
    stop("all arguments must be supplied exactly once", call. = FALSE)
  }
  stats::setNames(args[seq.int(2L, length(args), 2L)], keys)
}

read_tsv <- function(path) {
  utils::read.delim(
    path, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character()
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "", fileEncoding = "UTF-8"
  )
}

sha256 <- function(path) {
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
for (name in c("core", "expanded", "category-map")) {
  cli[[name]] <- normalizePath(cli[[name]], mustWork = TRUE)
}
package_root <- normalizePath(cli[["package-root"]], mustWork = TRUE)

# Pinned to the score-free runtime resources now published as the first LISA
# dictionaries (lisa_dictionary_core@1.0.0, lisa_dictionary_expanded@1.0.0) and
# the unchanged active category map (lisa_category_map@1.0.0). The digests
# below are the same bytes C6 produced: the first-publication release moved the
# identity and the filename only. The dictionary digests had changed earlier,
# when the runtime `LISA_score` column was removed; no assignment, label, tier
# or ordering changed then either.
expected_hashes <- c(
  core = "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994",
  expanded = "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736",
  `category-map` = "d61fcb2e1d40d0448d459daaab952975477203bb00b76cb53145aea8f61333b1"
)

# The resource identities those digests belong to, kept next to the hashes so
# the emitted provenance cannot drift from the bytes actually read.
source_identities <- c(
  core = "lisa_dictionary_core@1.0.0",
  expanded = "lisa_dictionary_expanded@1.0.0",
  `category-map` = "lisa_category_map@1.0.0"
)
observed_hashes <- vapply(cli[names(expected_hashes)], sha256, character(1L))
assert(identical(unname(observed_hashes), unname(expected_hashes)),
       "semantic source SHA-256 mismatch")

core <- read_tsv(cli[["core"]])
expanded <- read_tsv(cli[["expanded"]])
category_map <- read_tsv(cli[["category-map"]])
collections <- c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
assert(setequal(unique(core$universe), collections),
       "core dictionary collection set mismatch")
assert(setequal(unique(expanded$universe), collections),
       "expanded dictionary collection set mismatch")

core_key <- paste(core$universe, core$category_id, core$gene_set_id, sep = "\r")
expanded_key <- paste(
  expanded$universe, expanded$category_id, expanded$gene_set_id, sep = "\r"
)
assert(all(core_key %in% expanded_key),
       "a core category/gene-set membership is absent from expanded")
assert(!anyDuplicated(category_map$category_id),
       "category map contains duplicate category IDs")

represented_categories <- sort(unique(core$category_id))
map_only <- setdiff(category_map$category_id, represented_categories)
assert(length(represented_categories) == 158L, "expected 158 represented categories")
assert(!length(map_only), "active category map contains an unassigned category")
assert(!length(setdiff(represented_categories, category_map$category_id)),
       "a represented category is absent from the active category map")

category_map <- category_map[category_map$category_id %in% represented_categories,
                             , drop = FALSE]
map_index <- match(core$category_id, category_map$category_id)
assert(!anyNA(map_index), "a represented category is absent from category map")
assert(all(core$category_display_name == category_map$display_name[map_index]),
       "dictionary and category-map display names disagree")

supercategory_descriptions <- c(
  CC_IMMUNE_COMPLEXES = "Covers cellular complexes that assemble immune and inflammatory signaling machinery.",
  CC_METABOLIC_ORGANELLES = "Covers organelles that host energy production and intermediary metabolism, especially mitochondria and peroxisomes.",
  CC_NUCLEAR_GENE_EXPRESSION = "Covers nuclear, chromatin, transcriptional, RNA-processing and ribosomal compartments.",
  CC_PROTEOSTASIS_GENERIC_COMPLEXES = "Covers macromolecular complexes involved in protein folding, degradation and other general protein assemblies.",
  CC_SECRETORY_ENDOMEMBRANE = "Covers endoplasmic-reticulum, Golgi, vesicular, endolysosomal and secretory compartments.",
  CC_SURFACE_ADHESION_ECM = "Covers plasma-membrane, cell-junction, cytoskeletal and extracellular-matrix compartments.",
  CHEMOKINE_GPCR_SIGNALING = "Groups chemokine-receptor and related GPCR signaling axes that guide cell positioning and communication.",
  CYTOKINE_JAK_STAT = "Groups cytokine pathways that transmit signals through JAK kinases and STAT transcription factors.",
  DEATH_RECEPTOR_SIGNALING = "Groups extrinsic cell-death pathways initiated by death receptors and their proximal signaling complexes.",
  DEVELOPMENT_NEURAL = "Groups developmental, morphogenetic, neurogenic and synaptic biological programs.",
  GENOME_GENE_EXPRESSION = "Groups genome maintenance and the transcriptional, post-transcriptional and translational control of gene expression.",
  GPCR_SECOND_MESSENGER = "Groups GPCR and phospholipase pathways that use cyclic nucleotides, calcium and other second messengers.",
  HIPPO_MECHANOSIGNALING = "Groups Hippo/YAP/TAZ and mechanotransduction pathways that connect tissue mechanics to cell-state regulation.",
  IMMUNE_INFLAMMATION = "Groups adaptive, innate, cytokine, complement and antigen-presentation programs in immune and inflammatory biology.",
  INFLAMMATORY_NFKB = "Groups canonical and non-canonical NF-kB pathways downstream of inflammatory receptors.",
  INNATE_IMMUNE_SIGNALING = "Groups pattern-recognition, nucleic-acid sensing and inflammasome pathways of innate immunity.",
  METABOLISM_BIOENERGETICS = "Groups amino-acid, carbohydrate, lipid, nucleotide, mitochondrial and detoxification programs that support metabolism and energy balance.",
  MF_CATALYTIC_METABOLIC = "Groups molecular functions carried out by metabolic, redox, detoxification and other catalytic enzymes.",
  MF_CELL_CYCLE_DEVELOPMENT = "Groups molecular regulators of cell-cycle progression, differentiation and developmental decisions.",
  MF_GENERIC_OTHER = "Groups broad binding or host-interaction molecular functions that do not fit a more specific functional class.",
  MF_GENE_EXPRESSION = "Groups chromatin, transcription, RNA-processing and translation-related molecular functions.",
  MF_IMMUNE_LIGAND_RECEPTOR = "Groups immune recognition, complement, receptor and extracellular ligand molecular functions.",
  MF_PROTEOSTASIS_STRESS_DEATH = "Groups molecular functions involved in proteostasis, stress responses, proteolysis and regulated cell death.",
  MF_SIGNALING_REGULATION = "Groups receptors, kinases, phosphatases, adaptors and enzyme regulators that control signaling.",
  MF_STRUCTURE_ADHESION_ECM = "Groups structural, motor, adhesion and extracellular-matrix molecular functions.",
  MF_TRANSPORT_SENSING = "Groups channels, pumps, transporters, vesicle-trafficking factors and physicochemical sensors.",
  MICROENVIRONMENT_VASCULAR_ECM = "Groups extracellular-matrix remodeling, vascular, angiogenic and lymphatic programs in the tissue microenvironment.",
  NUCLEAR_RECEPTOR_SIGNALING = "Groups signaling and transcriptional programs controlled by steroid and metabolic nuclear receptors.",
  PI3K_AKT_MTOR = "Groups PI3K, AKT, mTOR and AMPK pathways that coordinate growth, survival and nutrient sensing.",
  PROLIFERATION_CELL_CYCLE = "Groups programs governing DNA replication, mitosis and cell proliferation.",
  RHO_GTPASE_CYTOSKELETON = "Groups RHO-family GTPase, adhesion and actin-regulatory pathways that organize cell shape and movement.",
  RTK_MATRIX_GUIDANCE_SIGNALING = "Groups receptor and guidance pathways that connect extracellular cues, matrix interactions and directed cell behavior.",
  RTK_RAS_MAPK = "Groups receptor-tyrosine-kinase signaling through RAS and the RAF-MEK-ERK cascade.",
  SIGNALING_COMMUNICATION = "Groups broad growth-factor, hormonal and developmental pathways used for intercellular communication.",
  STRESS_DEATH_PROTEOSTASIS = "Groups stress adaptation, autophagy, proteostasis and regulated cell-death programs.",
  STRESS_MAPK = "Groups stress-responsive p38 and JNK mitogen-activated protein kinase pathways.",
  STRESS_RESPONSE_SIGNALING = "Groups NRF2 antioxidant and unfolded-protein-response branches that manage oxidative and endoplasmic-reticulum stress.",
  STRUCTURE_MOTILITY_PLASTICITY = "Groups adhesion, polarity, cytoskeletal, migratory and cell-state-plasticity programs.",
  TGF_BETA_SMAD = "Groups TGF-beta-family ligands and their SMAD-dependent developmental and homeostatic pathways.",
  TRANSPORT_TRAFFICKING_HOMEOSTASIS = "Groups membrane transport, localization, vesicle traffic and ionic homeostasis programs.",
  T_CELL_B_CELL_RECEPTOR_SIGNALING = "Groups antigen-receptor and costimulatory signaling in T and B lymphocytes.",
  WNT_NOTCH_HEDGEHOG = "Groups WNT, NOTCH and Hedgehog developmental signaling, including canonical and non-canonical branches."
)
represented_macrogroups <- sort(unique(category_map$macrogroup_id))
assert(length(represented_macrogroups) == 42L, "expected 42 represented supercategories")
assert(setequal(names(supercategory_descriptions), represented_macrogroups),
       "supercategory description coverage mismatch")

collection_descriptions <- c(
  `GOBP-C2` = "Groups biological-process and curated C2 gene sets describing %s.",
  GOMF = "Groups Gene Ontology molecular-function gene sets describing %s.",
  GOCC = "Groups Gene Ontology cellular-component gene sets describing %s.",
  PATHWAYS = "Groups curated signaling and pathway gene sets describing %s."
)

tokens <- function(x) {
  x <- toupper(gsub("[^A-Za-z0-9]+", " ", x))
  out <- unlist(strsplit(x, "[[:space:]]+"), use.names = FALSE)
  out <- setdiff(out[nzchar(out)], c(
    "AND", "OR", "OF", "THE", "GENERAL", "CORE", "ACTIVITY", "PROCESS",
    "PATHWAY", "PATHWAYS", "SIGNALING", "REGULATION", "RESPONSE"
  ))
  ifelse(nchar(out) > 6L, substr(out, 1L, 6L), out)
}

candidate_rows <- unique(core[c(
  "universe", "category_id", "category_display_name", "gene_set_id"
)])
candidate_rows$lexical_overlap <- vapply(seq_len(nrow(candidate_rows)), function(i) {
  length(intersect(
    tokens(candidate_rows$category_display_name[[i]]),
    tokens(candidate_rows$gene_set_id[[i]])
  ))
}, integer(1L))
candidate_rows$gene_token_count <- vapply(
  candidate_rows$gene_set_id, function(x) max(1L, length(tokens(x))), integer(1L)
)
candidate_rows$lexical_density <-
  candidate_rows$lexical_overlap / candidate_rows$gene_token_count
candidate_rows$specificity_penalty <- vapply(
  candidate_rows$gene_set_id,
  function(x) sum(vapply(
    c("CANCER", "CARCINOMA", "TUMOR", "DISEASE", "INFECTION", "VIRUS",
      "SARS", "COVID", "MODEL", "PATIENT", "SYNDROME"),
    grepl, logical(1L), x = x, fixed = TRUE
  )),
  integer(1L)
)

sentence_topic <- function(x) {
  out <- tolower(sub("[.]$", "", x))
  replacements <- c(
    "dna" = "DNA", "rna" = "RNA", "rnp" = "RNP", "ecm" = "ECM",
    "emt" = "EMT", "er" = "ER", "hif" = "HIF", "ros" = "ROS",
    "prr" = "PRR", "rtk" = "RTK", "tgf-beta" = "TGF-beta",
    "smad" = "SMAD", "wnt" = "Wnt", "notch" = "Notch",
    "hippo" = "Hippo", "oxphos" = "OXPHOS", "etc" = "ETC",
    "gpcr" = "GPCR", "jak" = "JAK", "stat" = "STAT",
    "nf-kb" = "NF-kB", "pi3k" = "PI3K", "akt" = "AKT",
    "mtor" = "mTOR", "ampk" = "AMPK", "mapk" = "MAPK",
    "jnk" = "JNK", "nrf2" = "NRF2", "yap" = "YAP", "taz" = "TAZ"
  )
  for (from in names(replacements)) {
    out <- gsub(
      paste0("\\b", gsub("([-.])", "\\\\\\1", from), "\\b"),
      replacements[[from]], out, perl = TRUE
    )
  }
  out
}

selection_rows <- list()
catalog_rows <- list()
group_keys <- unique(candidate_rows[c("universe", "category_id")])
for (i in seq_len(nrow(group_keys))) {
  universe <- group_keys$universe[[i]]
  category_id <- group_keys$category_id[[i]]
  candidates <- candidate_rows[
    candidate_rows$universe == universe &
      candidate_rows$category_id == category_id,
    , drop = FALSE
  ]
  # The `-numeric_score` tie-break that used to sit between lexical overlap and
  # gene-set-id length has been removed with the runtime score column. It
  # was provably inert: every candidate here comes from the core dictionary,
  # where the score was constant (4) across the whole pool and within every
  # (universe, category_id) group, so it could never separate two candidates.
  # The remaining criteria are unchanged and still fully determine the order.
  candidates <- candidates[order(
    candidates$specificity_penalty,
    -candidates$lexical_density,
    -candidates$lexical_overlap,
    nchar(candidates$gene_set_id),
    candidates$gene_set_id,
    method = "radix"
  ), , drop = FALSE]
  member_count <- nrow(candidates)
  example_count <- if (member_count <= 5L) 1L else if (member_count <= 20L) 2L else 3L
  selected <- candidates[seq_len(min(example_count, member_count)), , drop = FALSE]
  selected$example_rank <- seq_len(nrow(selected))
  selected$core_member <- TRUE
  selected$expanded_member <- paste(
    selected$universe, selected$category_id, selected$gene_set_id, sep = "\r"
  ) %in% expanded_key
  selection_rows[[length(selection_rows) + 1L]] <- selected

  map_row <- category_map[category_map$category_id == category_id, , drop = FALSE]
  category_label <- map_row$display_name[[1L]]
  category_description <- sprintf(
    collection_descriptions[[universe]],
    sentence_topic(category_label)
  )
  examples <- rep("", 3L)
  examples[seq_len(nrow(selected))] <- selected$gene_set_id
  catalog_rows[[length(catalog_rows) + 1L]] <- data.frame(
    universe = universe,
    macrogroup_id = map_row$macrogroup_id[[1L]],
    macrogroup_name = map_row$macrogroup_name[[1L]],
    macrogroup_order = as.integer(map_row$macrogroup_order[[1L]]),
    supercategory_description = unname(
      supercategory_descriptions[[map_row$macrogroup_id[[1L]]]]
    ),
    category_id = category_id,
    category_display_name = category_label,
    category_order_within_macrogroup = as.integer(
      map_row$category_order_within_macrogroup[[1L]]
    ),
    category_description = category_description,
    core_member_gene_sets = member_count,
    example_count = nrow(selected),
    example_gene_set_1 = examples[[1L]],
    example_gene_set_2 = examples[[2L]],
    example_gene_set_3 = examples[[3L]],
    source_dictionary = source_identities[["core"]],
    source_dictionary_sha256 = observed_hashes[["core"]],
    expanded_dictionary_sha256 = observed_hashes[["expanded"]],
    category_map_sha256 = observed_hashes[["category-map"]],
    verification_status = "verified_exact_membership",
    stringsAsFactors = FALSE
  )
}

catalog <- do.call(rbind, catalog_rows)
catalog <- catalog[order(
  match(catalog$universe, collections), catalog$macrogroup_order,
  catalog$category_order_within_macrogroup, catalog$category_id,
  method = "radix"
), , drop = FALSE]
row.names(catalog) <- NULL
selection <- do.call(rbind, selection_rows)
selection <- selection[order(
  match(selection$universe, collections), selection$category_id,
  selection$example_rank, method = "radix"
), , drop = FALSE]
row.names(selection) <- NULL

assert(nrow(catalog) == 158L, "catalog row count mismatch")
assert(!anyDuplicated(paste(catalog$universe, catalog$category_id)),
       "duplicate collection/category catalog row")
assert(all(catalog$example_count >= 1L & catalog$example_count <= 3L),
       "example count must be between one and three")
assert(all(selection$core_member) && all(selection$expanded_member),
       "a selected example lacks exact core/expanded membership")
assert(all(nzchar(catalog$category_description)) &&
       all(nzchar(catalog$supercategory_description)),
       "a description is empty")

catalog_dir <- file.path(package_root, "inst", "extdata", "semantic-catalog")
catalog_path <- file.path(catalog_dir, "semantic_category_catalog.tsv")
supercategory_path <- file.path(catalog_dir, "semantic_supercategory_counts.tsv")
provenance_path <- file.path(catalog_dir, "SEMANTIC_CATALOG_PROVENANCE.tsv")
write_tsv(catalog, catalog_path)

# A gene set assigned to two categories in one supercategory counts once in
# the union, but twice as a category membership. Never sum category counts
# and label that sum as distinct gene sets or semantic coverage.
memberships <- unique(core[c("universe", "category_id", "gene_set_id")])
memberships$macrogroup_id <- category_map$macrogroup_id[
  match(memberships$category_id, category_map$category_id)
]
supercategory_keys <- unique(catalog[c("universe", "macrogroup_id")])
supercategory_counts <- do.call(rbind, lapply(seq_len(nrow(supercategory_keys)), function(i) {
  key <- supercategory_keys[i, , drop = FALSE]
  part <- memberships[
    memberships$universe == key$universe &
      memberships$macrogroup_id == key$macrogroup_id, , drop = FALSE
  ]
  data.frame(
    universe = key$universe, macrogroup_id = key$macrogroup_id,
    category_count = length(unique(part$category_id)),
    core_distinct_gene_sets = length(unique(part$gene_set_id)),
    core_member_assignments = nrow(part),
    source_dictionary_sha256 = observed_hashes[["core"]],
    category_map_sha256 = observed_hashes[["category-map"]],
    stringsAsFactors = FALSE
  )
}))
assert(nrow(supercategory_counts) == 42L, "supercategory count table mismatch")
assert(all(supercategory_counts$core_distinct_gene_sets <=
           supercategory_counts$core_member_assignments),
       "distinct supercategory gene sets exceed assignments")
write_tsv(supercategory_counts, supercategory_path)
write_tsv(data.frame(
  source = unname(source_identities[c("core", "expanded", "category-map")]),
  sha256 = unname(observed_hashes),
  rows = c(nrow(core), nrow(expanded), nrow(read_tsv(cli[["category-map"]]))),
  bundled = TRUE,
  use = c(
    "member counts and representative-example selection",
    "exact membership cross-check for selected examples",
    "category labels, order and supercategory membership"
  ),
  stringsAsFactors = FALSE
), provenance_path)
write_tsv(selection, cli[["audit"]])

escape_cell <- function(x) {
  x <- gsub("|", "\\|", x, fixed = TRUE)
  gsub("\n", " ", x, fixed = TRUE)
}
example_cell <- function(row) {
  examples <- unlist(row[c(
    "example_gene_set_1", "example_gene_set_2", "example_gene_set_3"
  )], use.names = FALSE)
  examples <- examples[nzchar(examples)]
  paste0("`", escape_cell(examples), "`", collapse = "<br>")
}

summary_rows <- lapply(collections, function(universe) {
  part <- catalog[catalog$universe == universe, , drop = FALSE]
  sprintf(
    "| `%s` | %d | %d | %d |",
    universe, length(unique(part$macrogroup_id)), nrow(part),
    sum(part$core_member_gene_sets)
  )
})

lines <- c(
  "---",
  "title: \"Semantic categories and supercategories\"",
  "lang: en-GB",
  "output:",
  "  rmarkdown::html_vignette:",
  "    mathjax: null",
  "    includes:",
  "      in_header: vignette-layout.inc",
  "    toc: true",
  "    toc_depth: 2",
  "vignette: >",
  "  %\\VignetteIndexEntry{Semantic categories and supercategories}",
  "  %\\VignetteEngine{knitr::rmarkdown}",
  "  %\\VignetteEncoding{UTF-8}",
  "---",
  "",
  "This catalogue describes the 158 LISA categories, their supercategories and",
  "representative member gene sets. It records dictionary membership, not the",
  "enrichment results of a particular study.",
  "",
  "The two small lookups below read the installed catalogue without running",
  "enrichment. The complete catalogue can also be read without running R.",
  "",
  "## Find a category",
  "",
  "A category is identified by its collection and category ID. This example",
  "finds the label, supercategory, member count and example gene sets for",
  "DNA damage repair and checkpoints in GOBP-C2:",
  "",
  "```{r semantic-lookup, eval=TRUE}",
  "catalog_path <- system.file(",
  "  \"extdata\", \"semantic-catalog\", \"semantic_category_catalog.tsv\",",
  "  package = \"lisaR\", mustWork = TRUE",
  ")",
  "semantic_catalog <- read.delim(",
  "  catalog_path, sep = \"\\t\", quote = \"\", comment.char = \"\",",
  "  colClasses = \"character\", check.names = FALSE",
  ")",
  "chosen <- semantic_catalog[",
  "  semantic_catalog$universe == \"GOBP-C2\" &",
  "    semantic_catalog$category_id == \"DNA_DAMAGE_REPAIR_CHECKPOINTS\",",
  "  c(\"category_display_name\", \"macrogroup_name\", \"core_member_gene_sets\",",
  "    \"example_gene_set_1\", \"example_gene_set_2\", \"example_gene_set_3\")",
  "]",
  "stopifnot(nrow(chosen) == 1L)",
  "knitr::kable(chosen, row.names = FALSE)",
  "```",
  "",
  "The count is membership in the frozen core dictionary. The three examples",
  "illustrate the category's scope; they are not the three strongest results of",
  "your study and are not its full membership list. The live evidence sheet,",
  "not this catalogue, supplies NES, FDR, leading edges and gene-level DE.",
  "",
  "## Member sets and category assignments",
  "",
  "A gene set may be assigned to two categories in one supercategory. It then",
  "contributes **one distinct set but two assignments**. The bundled count table",
  "lets you check that distinction directly:",
  "",
  "```{r semantic-distinct-support, eval=TRUE}",
  "semantic_counts <- read.delim(system.file(",
  "  \"extdata\", \"semantic-catalog\", \"semantic_supercategory_counts.tsv\",",
  "  package = \"lisaR\", mustWork = TRUE",
  "), sep = \"\\t\", check.names = FALSE)",
  "pi3k <- semantic_counts[",
  "  semantic_counts$macrogroup_id == \"PI3K_AKT_MTOR\",",
  "  c(\"universe\", \"macrogroup_id\", \"core_distinct_gene_sets\",",
  "    \"core_member_assignments\")",
  "]",
  "stopifnot(nrow(pi3k) == 1L,",
  "          pi3k$core_distinct_gene_sets == 59L,",
  "          pi3k$core_member_assignments == 88L)",
  "knitr::kable(pi3k, row.names = FALSE)",
  "```",
  "",
  "The 59 distinct gene sets have 88 assignments across categories in this",
  "supercategory. Assignment counts and distinct gene-set counts therefore",
  "measure different properties of the dictionary.",
  "",
  "To open the corresponding study evidence, use the category navigator or",
  "gene search described in the [report guide](scientific-report.html).",
  "",
  "## Catalogue sources and scope",
  "",
  "The four semantic collections are `GOBP-C2`, `GOMF`, `GOCC` and `PATHWAYS`.",
  "`HALLMARKS` is a separate direct gene-set collection, not another LISA category",
  "catalogue. Category descriptions and one to three exact member identifiers",
  "were checked against `lisa_dictionary_core@1.0.0`. The existing provenance also checks",
  "the same category vocabulary in the expanded 2.0.0 dictionary; this guide does",
  "not compare tiers or study results.",
  "",
  "Descriptions delimit a biological theme; they do not assert activation,",
  "inhibition or causality. Example selection uses semantic representativeness",
  "and stable tie-breaking on specificity, lexical match and gene-set identifier.",
  "Dictionary membership supplies no numerical weight or gene-priority score.",
  "Study evidence uses the member gene sets' enrichment and gene-level results.",
  "",
  "Unassigned gene sets remain in the complete tested universe. The",
  "`OTHER_UNCLASSIFIED` residual is tabular, not a biological category or figure.",
  "For `GOBP-C2`, `GOMF` and `GOCC` it has dedicated residual tables; `PATHWAYS`",
  "retains the complete GSEA ledger without a separate residual table.",
  "",
  "The core and expanded classification dictionaries and category map are",
  "included in lisaR. Gene-set-to-gene memberships are a separate external",
  "resource. Source conditions are recorded in `THIRD_PARTY_NOTICES.md`. See",
  "`vignette(\"resources-and-provenance\", package = \"lisaR\")` for those resource",
  "roles; no additional resource installation is needed to read this catalogue.",
  "",
  "Source versions and checksums accompany the installed catalogue in",
  "`extdata/semantic-catalog/SEMANTIC_CATALOG_PROVENANCE.tsv`.",
  "",
  "The active category map contains only categories with at least one exact",
  "assignment in the core or expanded dictionary. Supercategories without a",
  "retained category are absent from the active map.",
  "",
  "| Collection | Supercategories | Categories | Core member assignments |",
  "|---|---:|---:|---:|",
  unlist(summary_rows, use.names = FALSE),
  "",
  "# Catalogue",
  ""
)

for (universe in collections) {
  part <- catalog[catalog$universe == universe, , drop = FALSE]
  lines <- c(lines, paste0("## `", universe, "`"), "")
  macro_ids <- unique(part$macrogroup_id)
  for (macro_id in macro_ids) {
    block <- part[part$macrogroup_id == macro_id, , drop = FALSE]
    counts <- supercategory_counts[
      supercategory_counts$universe == universe &
        supercategory_counts$macrogroup_id == macro_id, , drop = FALSE
    ]
    macro_name <- block$macrogroup_name[[1L]]
    heading <- if (identical(macro_name, macro_id)) {
      display_overrides <- c(
        CHEMOKINE_GPCR_SIGNALING = "Chemokine GPCR signaling",
        CYTOKINE_JAK_STAT = "Cytokine/JAK/STAT signaling",
        DEATH_RECEPTOR_SIGNALING = "Death-receptor signaling",
        GPCR_SECOND_MESSENGER = "GPCR and second-messenger signaling",
        HIPPO_MECHANOSIGNALING = "Hippo and mechanosignaling",
        INFLAMMATORY_NFKB = "Inflammatory NF-kB signaling",
        INNATE_IMMUNE_SIGNALING = "Innate immune signaling",
        NUCLEAR_RECEPTOR_SIGNALING = "Nuclear-receptor signaling",
        PI3K_AKT_MTOR = "PI3K/AKT/mTOR signaling",
        RHO_GTPASE_CYTOSKELETON = "RHO GTPase and cytoskeletal signaling",
        RTK_MATRIX_GUIDANCE_SIGNALING = "RTK, matrix and guidance signaling",
        RTK_RAS_MAPK = "RTK/RAS/MAPK signaling",
        STRESS_MAPK = "Stress MAPK signaling",
        STRESS_RESPONSE_SIGNALING = "Stress-response signaling",
        TGF_BETA_SMAD = "TGF-beta/SMAD signaling",
        T_CELL_B_CELL_RECEPTOR_SIGNALING = "T- and B-cell receptor signaling",
        WNT_NOTCH_HEDGEHOG = "Wnt, Notch and Hedgehog signaling"
      )
      paste0(display_overrides[[macro_id]], " (`", macro_id, "`)")
    } else {
      paste0(escape_cell(macro_name), " (`", macro_id, "`)")
    }
    lines <- c(
      lines,
      paste0("### ", heading),
      "",
      block$supercategory_description[[1L]],
      "",
      sprintf(
        "Verified membership: %d %s, %d distinct core gene sets and %d category assignments.",
        nrow(block), if (nrow(block) == 1L) "category" else "categories",
        counts$core_distinct_gene_sets[[1L]],
        counts$core_member_assignments[[1L]]
      ),
      "",
      "| Category | Description | Core members | Representative member gene sets |",
      "|---|---|---:|---|"
    )
    for (j in seq_len(nrow(block))) {
      row <- block[j, , drop = FALSE]
      lines <- c(lines, sprintf(
        "| %s (`%s`) | %s | %d | %s |",
        escape_cell(row$category_display_name[[1L]]), row$category_id[[1L]],
        escape_cell(row$category_description[[1L]]),
        row$core_member_gene_sets[[1L]], example_cell(row)
      ))
    }
    lines <- c(lines, "")
  }
}

lines <- c(
  lines,
  "# Study evidence and custom categories",
  "",
  "The catalogue describes the dictionary. Study-specific NES, adjusted",
  "p-values and leading edges are available in the report's category evidence.",
  "See the [report guide](scientific-report.html) for that navigation.",
  "For custom resources and schema requirements, see",
  "`vignette(\"custom-dictionaries-and-supercategories\", package = \"lisaR\")`."
)

vignette_path <- file.path(package_root, "vignettes", "semantic-category-guide.Rmd")
writeLines(lines, vignette_path, useBytes = TRUE)

cat("SEMANTIC_CATALOG=PASS\n")
cat("collections=", length(collections), "\n", sep = "")
cat("supercategories=", length(unique(catalog$macrogroup_id)), "\n", sep = "")
cat("categories=", nrow(catalog), "\n", sep = "")
cat("selected_examples=", nrow(selection), "\n", sep = "")
cat("catalog=", catalog_path, "\n", sep = "")
cat("supercategory_counts=", supercategory_path, "\n", sep = "")
cat("vignette=", vignette_path, "\n", sep = "")
