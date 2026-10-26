
# sparse4 scoring specification

## Locked TCGA score

For a compatible TCGA-scale DESeq2-VST expression vector, the linear predictor is:

```text
0.0123106633921438 * CAMK2N1
- 0.0175920351043370 * MASP1
- 0.0039450538526462 * MIR9-3HG
- 0.0218626023395592 * MS4A1
```

Higher values denote higher modelled hazard. The primary analysis uses the
continuous score. A median split is descriptive and is not an optimized or
transportable clinical cutoff.

## External-cohort transport

FULL external scoring requires all four locked genes. Within each eligible
external cohort, each gene is standardized using that cohort's observed mean
and standard deviation, then weighted by `coefficient * TCGA training SD`.
There is no outcome-guided refit, sign reversal, or gene deletion. A cohort
missing any locked gene is not a FULL evaluation; any partial score is labelled
exploratory.

## Consistency authority

The gene order and coefficients must agree in all of:

- `config/model_lock.json`
- `config/coefficients.tsv`
- `config/training_scale.tsv`
- `results/tables/stable2_sparse4_model_definition.tsv`
- `scripts/10_sparse4_score.R`
- `scripts/10_sparse4_score.py`

Run `python3 tests/validate_repository.py --mode release` to check this contract.
