# lisaR 1.0.0 (2026-09-29)

- Promote the package to stable release 1.0.0.
- Update package and pipeline release metadata while preserving configuration,
  report, manifest and scientific resource version identities.
- Require lisaR >= 1.0.0 and < 2.0.0 for the installed Riaz entry points.
- Preserve historical release notes and reproducibility references below.

# lisaR 0.99.1 (development, 2026-09-23)

- Show four highest-only adjusted category-P star levels in current Hommel
  views (`*` <= 0.05, `**` <= 0.01, `***` <= 0.001, `****` <= 0.0001),
  including portable SVG and offline help. Preserve exact P and d/N in the
  evidence card, and restore NES-axis ticks for large values in the separate
  member-NES panel. The historical S1–S5 and circle recipes remain available,
  but they are not the current default significance view.
- Present minimum category support as S1–S5 editorial bands, with exact d/N
  assignment and conservative percentages. Keep a single 95% confidence level,
  zero versus non-evaluable states and observed NES direction distinct. Add
  neutral support columns to both lollipop variants without changing historical
  geometry; retain historical downloads and reproducible source contracts.
  `refresh_lisa_support_grades()` refreshes only presentation in a new verified
  copy, preserving existing P values, counts, denominators and hypothesis families.
- Add robust, simultaneous category-enrichment inference to STANDARD/FULL
  reports, with a supported-count/NES figure and a complete coverage/results
  table. Category significance is adjusted P <= 0.05; older source recipes
  retain their original asterisk rendering.
  HALLMARKS and descriptive DE-profile contrasts are excluded. Shared
  gene-set nulls are conservatively deduplicated across LISA collections;
  eligible missing P values stay in the family as 1. Saved runs can be
  augmented only in a new verified copy, without repeating fgsea.
  Legacy blank category annotations remain unclassified, while their
  eligible gene sets are retained in the joint testing family.
- Start from a compact overview with analysis and contrast counts, result
  cards and a visible Analysis map. Summary tables and output counts remain
  available under Report details. Soft section colours distinguish analyses,
  contrasts and navigation without changing the scientific figures.
- `run_lisa()` accepts prepared DE results directly, resolves the bundled
  dictionaries and enables PNG, SVG, source tables and figure recipes by
  default. Advanced configuration files remain supported.
- `run_lisa(mode = "full")` now writes one report in the same `output_dir`:
  the standard report plus every FULL product in its analysis, collection and
  category. FULL products are rendered inside the run (under `artifacts/`),
  sealed by the same run manifest and checked by `verify_lisa_run()`; a failed
  FULL product fails the run. No `<output_dir>-full-report` folder is created
  and the result no longer has `extension_output_dir`.
- Native KEGG maps (`run_kegg_maps`) are rendered into the same `output_dir`
  when the declared `cache_only` snapshot is usable; otherwise the report
  cover and the log say so. No `<output_dir>-kegg-maps` folder is created.
- One source for figure formats: the run's `report_output_policy.tsv`. With
  `svg = TRUE` recurrent genes, contrast gene cards, paired heatmaps,
  gene-category networks and painted KEGG maps now also carry SVG; product
  cards flag a missing SVG. `build_LISA_report.R` accepts `--report-mode`.
- `render_lisa_categories()` keeps its own on-demand gallery.

# lisaR 0.99.1

- Install directly from GitHub; the dictionaries and synthetic sample project
  are included, with no second source-file installation.
- Use short example input paths, transaction directories and contrast figure
  names. Windows checks identify destinations exceeding the path budget.
- Preserve scientific IDs, input bytes, old resource versions and saved runs.

## First-publication dictionary identities (2026-09-13)

- The dictionaries lisaR ships are the first LISA dictionaries ever published,
  so they now carry version `1.0.0` in their own first-publication logical
  namespace: `lisa_dictionary_core@1.0.0`, `lisa_dictionary_expanded@1.0.0`,
  `lisa_dictionary_quickstart@1.0.0` and
  `lisa_dictionary_custom_example@1.0.0`. `lisa_category_map@1.0.0` was already
  at its first-publication version and is unchanged.

- This is an identity change only. Every recorded SHA-256 is unchanged, and no
  gene-set assignment, membership, tier, category label, ordering or scientific
  setting moved. `LISA_score` does not return to the active runtime. Schema
  identities (`lisa_dictionary@2`, `category_map@1`), the package version and
  the upstream MSigDB release `2026.1` are separate version spaces and are not
  renumbered.

- A new namespace is used rather than reusing `lisa_core@1.0.0` and
  `lisa_expanded@1.0.0`, because those names already identify the older
  score-carrying bytes. A published identity is never re-minted for different
  content, so those pins continue to read exactly what they always read.

- Configurations pinning `lisa_core@2.0.0`, `lisa_expanded@2.0.0`,
  `lisa_quickstart_dictionary@4.0.0` or
  `lisa_example_custom_dictionary@2.0.0`, and version-less references to
  `lisa_core` or `lisa_expanded`, keep working and read the identical bytes
  through a closed compatibility alias that verifies the declared digest and
  schema before returning. The alias is consulted only after ordinary
  resolution and the historical rows have failed, so it can never redirect a
  historical pin. Resolution reports the requested logical ID, the requested
  resource ID and the resolved one.

- The active artifacts were renamed so their filenames stop encoding retired
  internal numbers (`lisa_dictionary_core_runtime_v1_0.tsv` and siblings).
  Every historical artifact is retained byte for byte under its original name.
  `lisa_core@2.0.0`, `lisa_expanded@2.0.0` and
  `lisa_quickstart_dictionary@4.0.0` remain reserved against re-minting by an
  external registry; `lisa_example_custom_dictionary@2.0.0` is not reserved,
  because it was never a built-in registry row and users register it themselves
  in the custom-resource workflow.

## Resource-version migration (2026-09-11)

- The LISA project's own scientific resources now carry their
  first-publication versions: `lisa_core@1.0.0`, `lisa_expanded@1.0.0` and
  `lisa_category_map@1.0.0`, replacing `0.1.0`, `0.1.0` and `0.1.1`. The
  migration renames identities only. No dictionary, category map, assignment,
  tier, LISA_score or MSigDB source changed, every artifact keeps its exact
  bytes and SHA-256, and no file was renamed.

- Configurations pinning the previous identities keep working and read the
  identical artifact. Superseded identities are recorded in
  `extdata/dictionaries/RESOURCE_MIGRATION_MANIFEST.tsv`; resolution redirects
  an old pin when no registered row for that identity is eligible for the
  requested species and profile (exact modality first, then the
  profile-independent `all`), verifies the resolved bytes against the recorded
  digest, and reports both the requested and the active resource ID. A registry
  row that *is* eligible still wins, and an unknown version, an unknown species,
  an ambiguous match or a corrupt artifact still fails loudly instead of being
  redirected. Nothing is deleted and no cache is invalidated: the fgsea content
  cache keys on memberships and software versions, never on a resource label.

- Superseded identities stay reserved on the same key their built-in rows used:
  logical ID, version, species **and** modality. An external registry row for
  `lisa_core@0.1.0`, `lisa_expanded@0.1.0` or `lisa_category_map@0.1.1` under
  `Homo sapiens` and the profile-independent modality `all` is refused with the
  new `LISA-RESOURCE-030`, preserving exactly the protection
  `LISA-RESOURCE-020` gave those rows while they were built in. The reservation
  is deliberately **not** universal: as before this release, a registry may
  declare one of those identities under a *specific* profile such as
  `transcriptomic/genomic`, and that row then wins for that profile. Widening
  the guard to every profile would change existing custom-resource behaviour
  and is out of scope for this migration. Under any other profile the same pin
  still falls back to the renamed identity through the compatibility ledger.

- The quick-start fixtures deliberately keep `lisa_quickstart_dictionary@3.0.0`
  and `lisa_quickstart_term2gene@2.0.0`. Version `1.0.0` already exists for
  those same logical identities with different bytes, and repointing a
  published identity is not an acceptable price for uniform numbering.

- Unchanged by design: the package version (`0.99.0`), the internal schema
  versions (`lisa_dictionary@1`, `category_map@1`, `term2gene@1`) and the
  upstream `msigdb_term2gene@2026.1` provenance, which names an MSigDB release
  rather than a version of our resource.

## Full report integration (2026-09-08)

- Contrasts now get the same category-evidence navigator as Analyses: an
  inline SVG with one clickable row per category, linking straight into its
  exact contrast-evidence card (`build_contrast_navigation.R`, wired next to
  `build_contrast_evidence.R`). The plotted value is the already-computed
  `delta_mean_NES` (A minus B); no new delta, FDR or category test is
  introduced. Fixes the Contrasts category graph not being clickable and its
  presentation differing from Analyses.

- Leading-edge heatmap generation now automatically emits a self-contained,
  two-argument native R recipe with safely embedded exact source hash, category,
  sample/gene orders, captions and dimensions; PNG, PDF or native Cairo SVG
  is selected by output extension. Generation and reproduction use
  shared drawing code; scientific selection and saved values are unchanged.
  Dense matrices now receive independent bounded row/column dimensions and
  larger labels. No external context file or manual report repair is required.

- Attach full products to exact category evidence cards with portable,
  content-addressed asset paths and original download names.
- Treat complete contrast heatmap cards as the category-scoped FULL product,
  rather than requiring a second nonexistent multicategory heatmap gallery;
  incomplete PNG/source/recipe attachments continue to fail closed.
- Preserve saved evidence recipes and resolve contrast identities from the
  recorded contrast index in standalone extensions; keep nested navigation valid.
- Discover KEGG gene-set products using the source run's GSEA FDR threshold,
  matching rendering. These products are not painted KEGG pathway maps.

## Licence (2026-09-07)

- lisaR is released under GNU General Public License version 3 or later (`License: GPL (>= 3)`, SPDX `GPL-3.0-or-later`). Copyright holder: Agencia Estatal Consejo Superior de Investigaciones Científicas (CSIC); author: David Olmeda Casadomé.
- R, JavaScript and CSS sources carry copyright and licence headers.
- `inst/THIRD_PARTY_NOTICES.md` records the CC BY-SA 4.0 to GPLv3 compatibility route for KEGG MEDICUS annotations and the transitive-only status of `cowplot` and `fastmatch`.

## Report navigation

- Keep a hash-restored report section stable through reload or a fresh-tab
  deep link: the shell now re-applies `scrollIntoView` on the pending target
  until layout settles (window `load`, a `ResizeObserver` on the document
  body, or a short deadline) and holds the scroll-tracking guard the whole
  time, so the debounced scroll listener can no longer overwrite the
  restored section and hash once content above it finishes loading.
- Restore the full LISA logo and a persistent context/collection/section navigator.
- Unify LISA categories and NES views, retaining supercategories and downloads.
- Add `percentages` beside-point annotations; all four NES variants default on, any nonempty subset allowed. Existing `direction` and `dispersion` meanings are unchanged.
- Retain contrast all/same-direction/opposite-direction subsets in the unified viewer.
- New NES figures declare an isolated graphics context shared by rendering and reproduction; archived figures retain their recorded recipes.
- Bound generated NES basenames while preserving exact IDs in their metadata, so view/subset downloads remain portable in deep report folders.

# lisaR 0.99.0 — 1.0 beta

## Configuration schema 1.0

- Use `pipeline.schema_version: "1.0.0"` in YAML, JSON and R-list studies.
  Older configuration files require an explicit, reviewed migration; resource
  and package versions remain unchanged.
- Make category/gene/contrast evidence, legacy-product opt-in and initial
  matrix limits explicit in the installed templates.
- Keep `clean`, `percentages`, `direction` and `dispersion` as the default NES selection.
  Accept one to four unique variants, never an empty selection;
  no particular variant is compulsory.
- Start the Quick Start in `standard` mode. Full output remains an explicit
  choice after reviewing the output plan.

## Presentation and paired evidence

- Share one compact offline header and scientific context across report views.
- Integrate category navigation after LISA category figures without an iframe
  or redundant visible category table; keep full downloadable data.
- Place member-set figures at the end of each category evidence sheet.
- Add clean, direction and dispersion NES variants (all rendered by default),
  preserving canonical means/deltas and separating descriptive A/B statistics.
- Add paired category evidence with aligned gene sets/genes, explicit missing
  states, leading-edge overlap and exact roundtrip links through gene search.
- Preserve the original scientific outputs and dictionary assignments.

## Individual-analysis evidence

- Preserve GSEA category overview figures alongside standard evidence.
- Link individual category overview SVGs to scoped category evidence.
- Use red for positive and blue for negative NES/DE effects and expose dynamic
  set/gene ordering, explicit recurrence denominators and current-view exports.
- Add an offline shared gene search using full TERM2GENE membership, DE and
  GSEA, including significant member sets whose leading edge excludes the gene.

## Installation and examples

- Add an installation guide covering dependencies, an offline sample run
  and separately acquired scientific resources.
- Use the human release label **1.0 beta** and the R package version `0.99.0`.
  Configuration schema `1.0.0` is versioned separately. Resource versions are
  unchanged.
- Update Riaz example entry-point version checks and installed documentation.

## Corrections

- Accept the exact canonical category map without changing its scientific bytes.
- Preserve numeric configuration precision and explicit expression-matrix scale.
- Handle incomplete enrichment statistics and zero-weight sensitivity components.
- Reproduce lollipop figures through the original renderer and check manifest sizes.
- Correct unique gene-set counts, example paths, visible YAML and interpretation
  wording; count external indexes and reject unsupported native KEGG painting.

This is a beta release. Package and third-party resource licences remain unchanged.
## Hommel support presentation

- Integrate saved robust Hommel closed-testing results into the clickable
  category navigator and category-specific evidence sheets using exact
  analysis/collection/category keys and the raw-P-evaluable denominator.
- Earlier support-band presentation used one to five filled black circles,
  d/N and conservatively displayed percentages. Current default views use
  adjusted-P stars; the circle recipes remain compatible historical products.
- Share accessible offline interpretation help and the Goeman et al. (2019)
  reference; move complete inference statistics to a final collapsed technical
  reference. Contrasts remain descriptive and HALLMARKS remains excluded.
- Refresh presentation from verified saved results without recomputing inference;
  repeated refreshes retain previous figures and provenance archives.
