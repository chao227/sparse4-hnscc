#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "synthetic_expression.tsv"
COEFFICIENTS = ROOT / "config" / "coefficients.tsv"


def read_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="	"))


with tempfile.TemporaryDirectory() as temporary:
    temporary = Path(temporary)
    py_out = temporary / "python.tsv"
    r_out = temporary / "r.tsv"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "10_sparse4_score.py"), str(FIXTURE), str(py_out)], check=True, cwd=ROOT)
    subprocess.run(["Rscript", str(ROOT / "scripts" / "10_sparse4_score.R"), str(FIXTURE), str(r_out)], check=True, cwd=ROOT)
    py_rows = read_rows(py_out)
    r_rows = read_rows(r_out)
    if len(py_rows) != len(r_rows):
        raise SystemExit("R/Python scoring row counts differ")
    for left, right in zip(py_rows, r_rows):
        if left["sample_id"] != right["sample_id"]:
            raise SystemExit("R/Python scoring sample order differs")
        if not math.isclose(float(left["sparse4_linear_predictor"]), float(right["sparse4_linear_predictor"]), rel_tol=0, abs_tol=1e-12):
            raise SystemExit("R/Python scores differ")

    missing_gene = temporary / "missing_gene.tsv"
    fixture_rows = read_rows(FIXTURE)
    with missing_gene.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "CAMK2N1", "MASP1", "MIR9-3HG"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in fixture_rows:
            writer.writerow({key: row[key] for key in writer.fieldnames})
    failure = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "10_sparse4_score.py"), str(missing_gene), str(temporary / "invalid.tsv")],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if failure.returncode == 0 or "requires all four genes" not in failure.stderr:
        raise SystemExit("Missing-gene gate did not fail closed")
print("Scoring parity: PASS")
