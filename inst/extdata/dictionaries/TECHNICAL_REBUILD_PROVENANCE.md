# Technical reconstruction and package boundary

lisaR contains exact, hash-pinned scientific core,
expanded and active category-map runtime TSVs. Their installed-package
identities and SHA-256 values are recorded in `RESOURCE_BUNDLE_MANIFEST.tsv`.
Resource-specific terms apply; provenance and
third-party attribution are recorded in `RESOURCE_PROVENANCE.tsv` and
`THIRD_PARTY_NOTICES.md`. The package also contains explicitly synthetic Quick
Start and minimal fixtures in their separate `extdata` directories.

The installed script `scripts/build_compact_dictionary_resources.R` creates
runtime tables from six frozen, already-scored inputs supplied
locally by an authorised maintainer. It performs no download and does not
replay the upstream annotation workflow.

## Restored dictionary versions

The historical resources `lisa_core@1.0.0` and `lisa_expanded@1.0.0`
were previously identified as `lisa_core@0.1.0` and `lisa_expanded@0.1.0`.
That migration renamed identities without changing their contents. The current
builder produces score-free tables under the `lisa_dictionary_*` identities
listed below. It retains the reviewed BioCarta, KEGG LEGACY and KEGG
MEDICUS gene-set-to-LISA-category assignments. The builder emits
a source-family inventory and must preserve every accepted frozen row.

## First-publication identities of the current contents

The resources the package actually ships today are the first LISA dictionaries
ever published, so they carry version `1.0.0` in their own first-publication
logical namespace: `lisa_dictionary_core@1.0.0`,
`lisa_dictionary_expanded@1.0.0`, `lisa_dictionary_quickstart@1.0.0` and
`lisa_dictionary_custom_example@1.0.0`.

They are deliberately **not** published as `lisa_core@1.0.0` or
`lisa_expanded@1.0.0`. Those names were already published, and they name the
older score-carrying bytes described in the previous section; a published
identity is never re-minted for different content. The identities these current
contents previously carried — `lisa_core@2.0.0`, `lisa_expanded@2.0.0`,
`lisa_quickstart_dictionary@4.0.0`, `lisa_example_custom_dictionary@2.0.0`,
and the version-less `lisa_core` and `lisa_expanded` — keep resolving to these
exact bytes through a closed, code-authenticated alias that re-checks the
declared digest and schema before returning. `RESOURCE_MIGRATION_MANIFEST.tsv`
records the same mapping as an audit ledger under `republish_first_version`.

Nothing scientific moved: gene-set assignments, membership, tier, category
labels, ordering and the score-free `lisa_dictionary@2` contract are byte-for-
byte identical to the artifacts this section already describes. The schema
numbers (`lisa_dictionary@2`, `category_map@1`), the package version and the
upstream MSigDB release `2026.1` are separate version spaces and are unchanged.

The active artifact filenames were renamed so they stop encoding retired
internal numbers: `lisa_dictionary_core_runtime_v1_0.tsv`,
`lisa_dictionary_expanded_runtime_v1_0.tsv`,
`lisa_category_map_runtime_v1_0.tsv`,
`lisa_dictionary_quickstart_v1_0.tsv` and
`lisa_dictionary_custom_example_v1_0.tsv`. The historical artifacts
`lisa_core_runtime_v0_1.tsv`, `lisa_expanded_runtime_v0_1.tsv`,
`example_dictionary.tsv`, `example_dictionary_v2_0.tsv`,
`example_dictionary_v3_0.tsv` and `example_custom_dictionary.tsv` are untouched
and keep serving their own pins. Note that the external frozen construction
ledgers named `lisa_dictionary_core_v0_1.tsv` and
`lisa_dictionary_expanded_v0_1.tsv` in `TECHNICAL_REBUILD_CONTRACT.tsv` are
inputs held in the source archive; they are a different kind of file from the
bundled `*_runtime_*` artifacts and remain outside the package.

## Historical category map

The historical map `lisa_category_map@0.1.0` remains byte-preserved
and unbundled, and is the immutable source for the active
`lisa_category_map@1.0.0` (previously identified as `0.1.1`). That identity did
not move in the first-publication release; only its filename did. The installed
`scripts/build_active_category_map.R` removes every category with zero
assignments across the exact core and expanded dictionaries; a supercategory
with no retained category consequently has no row in the active map. It emits
`CATEGORY_MAP_PRUNING.tsv` so every retained or removed row is auditable.

The historical frozen inputs and their build evidence remain outside the
package source tree. A provenance record does not grant redistribution
rights in its contents.

## Upstream resource identification

The frozen assignments use gene-set identifiers from the current technical
reconstruction, MSigDB release 2026.1 accessed through `msigdbr` 26.1.0. The
historical retrieval release remains unrecorded. These are two separate
provenance facts and neither determines redistribution authority.

## Reconstruction limit

The accepted score ledgers and upstream provenance archive remain external.
The builder can reproduce the complete 0.1 resources from those ledgers, but
it does not reproduce model responses, consensus formation or manual
scientific review. TERM2GENE is not bundled and public
distribution remains subject to the recorded resource terms and project
release decisions.
