
# Data acquisition and placement

The analyses use only public-access datasets. Raw and participant-level derived
files are not distributed in this repository. Dataset accessions, cohort roles,
and source pages are listed in `accessions.csv`.

Place downloaded material under `data/raw/<dataset>/` and derived objects under
`data/processed/<dataset>/`. The expected minimum filenames and schemas are
checked by `python3 tests/validate_repository.py --mode full-inputs`. Download
dates were not recoverable from the original analysis record and are therefore
reported as `not recorded`; no dates have been inferred.

Minimum reconstructed inputs currently checked by the guarded dispatcher are:

- `data/processed/TCGA-HNSC/vst_top5000.rds`
- `data/processed/TCGA-HNSC/clinical_survival.csv`
- `data/processed/TCGA-HNSC/yellow_module_members.csv`
- `data/processed/GSE164690/seurat_annotated.rds`
- the five CellChat contract/registry files named in
  `config/cellchat_input_authority.tsv`
- `data/processed/GSE208253/spatial_qc.h5ad`
- `data/processed/GSE208253/sample_manifest.tsv`
- `data/processed/GSE208253/rctd/` with one registered directory per section
- the three external-cohort deposited expression files named by the full-input
  validator.

Some workflow modules also use locally reconstructed regulon summaries,
clinical inclusion ledgers, annotation cross-checks, single-cell objects, and
flat source tables. Their expected paths and schemas are documented in the
corresponding script headers. Individual figure renderers use the matching
local source objects, while the bundled figure dispatcher verifies the
released panel assets and their map.

The contract files use `record_after_reconstruction` until the public inputs
have been rebuilt locally. Replace each marker with the SHA-256 digest of the
reconstructed file before formal execution; code must not bypass that check.

After download, record file digests with:

```bash
find data/raw -type f -print0 | sort -z | xargs -0 sha256sum > data/checksums.sha256
```

Do not commit raw data, participant-level score tables, expression matrices,
spatial barcodes, private metadata, or controlled-access files.
