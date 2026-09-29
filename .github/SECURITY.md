# Security and private-data reporting

Do not open a public issue for a vulnerability, credential exposure or
unpublished-data leak. Contact the maintainer privately at
`dolmeda@iib.uam.es` with a minimal description and no sensitive attachment.

LISA configurations must not contain secrets. External-service credentials
must be managed through the user's normal environment or credential manager.
Reports and run directories can contain unpublished scientific data and must
be shared only through authorised channels.

## Private run and portable report boundary

The complete promoted run is the canonical private record. It may contain
absolute paths, input receipts and unpublished source material, so the run
directory must not be treated as a shareable report bundle. lisaR defines a
smaller portable report surface rooted at `report_index.html`. Copy only the
relative components enumerated in `shareable_report_manifest.tsv`; each is
bound to a SHA-256 digest. The manifest's own row uses the documented `SELF`
exception because a file cannot embed its own stable digest.

HTML uses only run-relative locations. Downloadable TSV/CSV report copies
rebase paths inside the run and hide the parent of external paths while
retaining the basename and role. The canonical source tables are not edited.
Anything outside the shareable manifest remains private by default.

This boundary does not decide whether an external basename is itself
sensitive, and it does not inspect text embedded in image pixels. Before
sharing, a human must review the listed filenames and rendered figures for
participant identifiers, unpublished labels or other confidential content.

## Managed paths and permissions

lisaR rejects pre-existing symbolic links inside managed run trees and
revalidates those trees before writes and final promotion. On Unix-like
systems, both staging and promoted results are private: directories are
`0700` and files are `0600`. lisaR does not change permissions on existing
ancestor directories. In a shared HPC project, the owner should first require
`verify_lisa_run(...)$gate` to return `"PASS"` and only then change result
permissions deliberately under the institutional data-sharing policy.

Base R does not provide an interface equivalent to
`openat(..., O_NOFOLLOW)`. Another process running as the same user and able to
modify the managed directory concurrently therefore retains a small
time-of-check/time-of-use race window. Do not execute lisaR in a tree that can
be modified concurrently by an untrusted process. Atomic promotion is not
promised for UNC paths or other network filesystems; use a local or
institutionally validated same-filesystem destination for transactional runs.

## Trusted RDS boundary

RDS deserialization is disabled by default and cannot be authorised from a
YAML/JSON study configuration. A caller using `run_lisa_de()` may opt in only
by naming `trusted_rds = TRUE`, which asserts that the file is local, stable
and fully trusted. `rds_max_bytes` checks serialized size before reading and
approximate object size afterwards, but it is diagnostic: it is neither a
sandbox nor a hard memory limit, and it does not protect against hostile
serialization, extreme compression, symlink changes or file replacement
between checks and `readRDS()`.

An accepted file produces a portable `lisa_rds_input_receipt_v2` receipt.
`source_sha256` identifies the serialized file; `object_sha256` hashes the
validated object using R serialization version 3 in XDR form. The object hash
is an internal branch-integrity fingerprint, not a stable semantic identity
across R or package versions. Receipts omit absolute source paths.
