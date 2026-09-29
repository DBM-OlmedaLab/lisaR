# Contributing to lisaR

lisaR is free software under the GNU General Public License version 3
or any later version (`License: GPL (>= 3)` in `DESCRIPTION`). The copyright holder is the Agencia
Estatal Consejo Superior de Investigaciones Científicas (CSIC). Contributions
are accepted under the same licence. Do not upload unpublished data, private
scientific resources or run outputs to the repository.

## Report an issue

Include:

- lisaR version and commit;
- operating system and R version;
- `sessionInfo()`;
- the smallest reproducible configuration with private paths and data removed;
- exact command and error;
- whether the failure occurs during validation, execution, reporting or run
  verification.

Never include credentials, patient data, unpublished expression matrices or
KEGG cache contents.

## Proposed changes

Open an issue before changing the public API, configuration schema, scientific
defaults, dictionary classification or report interpretation. Keep code,
tests, documentation and `NEWS.md` in the same change. A contribution is not
ready for review until the source tests and built-package checks pass.

Copyright in contributions made by CSIC staff belongs to CSIC. People whose
work is included are credited in `DESCRIPTION` (`Authors@R`), because CSIC
asks that everyone who took part in the development be named. Contributing
does not change the licence, which remains GPL-3.0-or-later.

## Documentation

Write documentation in clear, natural British English for scientists working
in R or RStudio. Preserve statistical,
privacy, licensing, and provenance boundaries while making the shortest safe
workflow easy to find.

Check code examples against the current package and review the rendered
guides for clarity. Edit the `.Rmd` sources in `vignettes/`. Build the package
with `R CMD build .` to compile its HTML vignettes; do not commit generated
copies. The website provides HTML guides and PDF downloads for readers.

Some optional integration tests use the distributed Riaz and CPTAC input
archives. Set `LISAR_RIAZ_ARCHIVE` and `LISAR_CPTAC_ARCHIVE` to their local
filenames to run those tests. Without those inputs, the relevant tests report
an explicit skip. The archives are not downloaded by the tests. Preparation
tests also accept `LISAR_TEST_CPTAC_BUNDLE` (an extracted bundle) and
`LISAR_TEST_RESOURCE_CACHE` (a prepared resource cache).

## Serial-versus-parallel release check

The worker protocol is a release invariant rather than a user workflow. Each
task must return a typed record with its task ID, success flag, value and
error. The coordinating process rejects a missing, `NULL`, malformed or failed
record. Linux workers use forked processes; lisaR does not expose a PSOCK
backend.

For a candidate release:

1. freeze one set of inputs, resources, package build and R environment;
2. run the same configuration with `workers: 1` and the candidate worker
   count;
3. require both runs and `verify_lisa_run()` to return `PASS`;
4. compare every shared scientific table exactly;
5. classify path, runtime and external-acquisition metadata separately; and
6. fail on every unexplained scientific difference.

On Linux, exercise `workers: 1`, `2` and `4`. The effective count may be
capped by the number of tasks, visible cores, scheduler allocations or CRAN's
`_R_CHECK_LIMIT_CORES_` rule. Test worker-record cardinality and structure as
well as final scientific equivalence.
