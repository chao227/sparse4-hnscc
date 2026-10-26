#!/usr/bin/env python3

# ==========================================================
# Script   : 10_sparse4_score.py
# Purpose  : Score compatible expression matrices with the locked four-gene coefficients
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: scoring specification; Supplementary Table 2
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/scoring/score_sparse4.py
# ==========================================================
"""Score the locked four-gene sparse4 model without outcome-guided refitting."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

COEFFICIENTS = {
    "CAMK2N1": 0.0123106633921438,
    "MASP1": -0.0175920351043370,
    "MIR9-3HG": -0.0039450538526462,
    "MS4A1": -0.0218626023395592,
}


def score_vst(samples_by_genes: pd.DataFrame) -> pd.Series:
    missing = [g for g in COEFFICIENTS if g not in samples_by_genes.columns]
    if missing:
        raise ValueError("FULL sparse4 score requires all four genes; missing: " + ", ".join(missing))
    x = samples_by_genes.loc[:, list(COEFFICIENTS)].apply(pd.to_numeric, errors="raise")
    return x.mul(pd.Series(COEFFICIENTS), axis=1).sum(axis=1).rename("sparse4_linear_predictor")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path, help="TSV with samples in rows and four genes in columns")
    p.add_argument("output", type=Path, help="Output TSV")
    p.add_argument("--sample-column", default="sample_id")
    args = p.parse_args()
    dat = pd.read_csv(args.input, sep="\t")
    if args.sample_column not in dat:
        raise ValueError(f"missing sample column: {args.sample_column}")
    out = pd.DataFrame({args.sample_column: dat[args.sample_column], "sparse4_linear_predictor": score_vst(dat)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
