# Release-candidate quality report

Date: 2026-10-26

## Completed checks

- Repository release checks passed, including required artifacts, UTF-8/LF
  text, local-path exclusion, credential patterns, file-size limits, script
  headers, Python syntax, and R parsing.
- R and Python implementations produced identical sparse4 scores on the
  synthetic fixture to an absolute tolerance of `1e-12`.
- The missing-gene test failed closed as required for a FULL score.
- The model lock, coefficient table, selected training-scale rows, and
  Supplementary Table 2 model definition agree on all four genes and full
  precision coefficients.
- Candidate provenance contains 539 historical Yellow symbols, of which 492
  are eligible in the locked feature set.
- The map contains 72 unique figure panels; 70 have at least one released
  binary asset. Final author-edited Fig. 7A and Fig. 7D are mapped without
  claiming an executable pixel generator or distributing an unavailable final
  asset.
- Sixty-two non-identifying aggregate TSV outputs are included. Participant-
  level score, prediction, donor, spot, section, and external patient ledgers
  are excluded.
- The requested restricted-term scan passed across text and binary content.

## Execution notes

- Data-intensive analyses require local retrieval of the public datasets and
  reconstruction of the inputs documented in `data/README.md` and the script
  headers.
- The released figure dispatcher verifies the panel map and existing assets;
  individual renderers use their corresponding local source objects.
- The public repository URL is recorded in `config/release.yml`; a persistent
  archive identifier can be added after a fixed release is created.
- The observed R lock records a mixed Bioconductor installation state. A clean
  restore must verify exact package availability and must not silently replace
  versions.
