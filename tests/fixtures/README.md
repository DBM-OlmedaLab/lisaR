# Test Fixtures

Only small synthetic, non-study fixtures belong here. Stable fixture paths are
reserved for unit, integration, security, report-contract, and network-cache
tests. Runtime material and study data must remain outside the package.

`unit/transcriptomic_de.tsv` and `unit/proteomic_de.tsv` are synthetic Input contract
preflight fixtures. They contain no study data and intentionally have too few
features for their named profiles, so warning and strong-warning behavior can
be tested deterministically.
