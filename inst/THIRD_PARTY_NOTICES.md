# Third-party resources

## Licence of the lisaR package

The lisaR package (code, synthetic fixtures and documentation) is free
software under the GNU General Public License version 3 or (at the recipient's
option) any later version (GPL-3.0-or-later). The copyright
holder is the Agencia Estatal Consejo Superior de Investigaciones Científicas
(CSIC); the author is David Olmeda Casadomé. This notice explains which
third-party material is *not* covered by that licence, and how the package's
licence relates to its dependencies and to the LISA annotations.

GPL-3.0-or-later permits private and commercial use, copying, modification and
redistribution of the covered package. Distribution of covered copies or
modified versions carries GPL-3.0-or-later notice and corresponding-source obligations.
It does not grant rights in third-party resources.

### Compatibility of the LISA annotations with GPL-3

The LISA category annotations classify gene-set identifiers from several
sources. MSigDB identifies KEGG MEDICUS content as CC BY-SA 4.0. Creative
Commons declares a one-way compatibility from CC BY-SA 4.0 to GPLv3: adapted
material may be conveyed inside a GPLv3 work. This is why the package licence
is `GPL-3.0-or-later` (`License: GPL (>= 3)`), not GPL-2 or GPL (>= 2).
Retain MSigDB and KEGG attribution and any applicable collection-specific
conditions. The licence of the code does not relicense source material.

### Dependency licences

lisaR calls `ggrepel` (GPL-3), `rmarkdown` (GPL-3) and `DESeq2` (LGPL >= 3)
directly, which is compatible with GPL-3 and incompatible with GPL-2 only.
`cowplot` (GPL-2 only) and `fastmatch` (GPL-2) reach an installation only as
transitive dependencies of `fgsea` (MIT). An MIT-licensed intermediary does
not automatically neutralize any transitive copyleft obligations; the relevant
licences and the actual combination or distribution must be considered.
lisaR does not call these packages, and the source tarball verification fails
if `cowplot` reappears in lisaR's own code.
lisaR has no compiled code and bundles no third-party JavaScript, CSS or fonts.

lisaR bundles exact, hash-pinned LISA core and expanded
classification dictionaries plus the active category map. They are listed in
`extdata/dictionaries/RESOURCE_BUNDLE_MANIFEST.tsv` and retain the provenance
and attribution boundaries in `RESOURCE_PROVENANCE.tsv`. It does not bundle
MSigDB TERM2GENE memberships or a native KEGG snapshot for scientific
execution. Those remain external resources under their own terms. The
documentation links to the official KEGG pathway page rather than bundling
the native diagram, a coloured copy or its map-specific source files.

The documentation also includes selected derived Riaz GSE91061 and
CPTAC ccRCC figures, their plotting sources, and one Riaz report screenshot.
Their study citations and saved-output provenance are recorded in the
vignettes and the historical editorial asset manifest retained with the
project records (not included in this package). The package contains no complete
cohort result bundle or canonical model benchmark. Scientific input bundles
are obtained separately through the example installer. Their redistribution
status remains separate from that of the package code.

The small files under `extdata/quick-start/`, `extdata/custom-resources/` and
the minimal example TERM2GENE are original synthetic sample data and
format-validation fixtures. They do not contain MSigDB, KEGG or BioCarta
scientific content. The manifest under `extdata/custom-resources/` records the
reviewed SHA-256 values of that synthetic bundle; it does not certify or grant
rights over user-supplied resources.

## LISA annotation resources

lisaR includes `lisa_dictionary_core@1.0.0`,
`lisa_dictionary_expanded@1.0.0` and `lisa_category_map@1.0.0`.
Their hashes are recorded in the installed package manifest and checked when
resolving resources. The annotations and semantic catalogue
retain source-specific attribution and conditions. The package licence does
not grant rights in third-party content.

The fixed dictionaries retain BioCarta, KEGG LEGACY and KEGG MEDICUS
identifier-to-category classifications. They contain no third-party gene
memberships; the separate TERM2GENE resource remains external and subject to
its provider terms. Retaining identifier classifications does not grant or
alter third-party rights.

## MSigDB

External LISA dictionaries and TERM2GENE resources may use gene-set
identifiers or memberships originating from the Molecular Signatures Database
(MSigDB). MSigDB attribution and collection-specific conditions remain
applicable.

- MSigDB licence terms:
  <https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp>
- The exact release used for the original historical retrieval was not
  recorded.
- The current technical reconstruction is identified as MSigDB 2026.1,
  accessed through `msigdbr` 26.1.0. This records technical lineage; it does
  not prove the historical release or permission to redistribute a LISA
  resource.

Institutional ownership was resolved on 2026-09-07: CSIC is the rights holder
of the package, and no registration or prior authorisation is required to
publish it. Preserve MSigDB and KEGG attribution and the applicable terms for any
third-party content that is redistributed.

## KEGG

KEGG pathway diagrams, images and KGML are external resources. Configured
`run_kegg_maps = TRUE` accepts only `kegg_access_mode: cache_only` with an
explicit recipient-authorized immutable local snapshot. It does not retrieve
maps or accept terms. The documentation links to KEGG rather than providing a
runtime map snapshot. Enrichment figures for KEGG
gene sets remain available; they are not copies of native pathway diagrams.
These content distinctions do not grant additional redistribution rights.

The software licence of an R dependency does not replace the KEGG terms that
govern retrieved content. Reports containing KEGG-derived figures must retain
the required attribution and may be shared only as permitted by those terms.

- KEGG copyright and licence information:
  <https://www.kegg.jp/kegg/legal.html>

## User responsibility

Users must ensure that installation, retrieval, use and distribution of every
external scientific resource is permitted for their institution and intended
use. The GPL-3.0-or-later licence of the lisaR package grants no third-party rights.
