# Environment records

`renv.lock` was generated from the observed R 4.4.1 release-preparation
environment using explicit package dependencies from `DESCRIPTION`.
`requirements.txt` records the relevant installed Python 3.10.18 packages
rather than an unfiltered environment dump.

The installed edgeR and limma versions came from a newer Bioconductor release
line than the default line detected for R 4.4.1. Their observed versions are
retained in the lockfile. A clean restoration must confirm that these records
remain obtainable and must not silently substitute different versions.

`session_info.txt` distinguishes observed release targets from exact historical
run-time evidence. When a historical run-specific session file exists, it has
priority for interpreting that run.
