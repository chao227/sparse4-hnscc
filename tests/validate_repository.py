#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise SystemExit(message)


def tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def release_checks() -> None:
    required = [
        "README.md", "LICENSE", "CITATION.cff", "CHANGELOG.md", "run_all.sh",
        "config/paths.yml", "config/params.yml", "config/model_lock.json",
        "data/accessions.csv", "environment/renv.lock", "environment/requirements.txt",
        "docs/figure_table_map.md", "docs/sparse4_scoring_spec.md", "docs/known_limitations.md",
        "results/manifest.csv", "provenance/source_provenance.csv",
    ]
    missing = [item for item in required if not (ROOT / item).is_file()]
    if missing:
        fail("Missing required files: " + ", ".join(missing))

    coefficients = tsv(ROOT / "config" / "coefficients.tsv")
    lock = json.loads((ROOT / "config" / "model_lock.json").read_text(encoding="utf-8"))
    lock_pairs = [(row["gene_symbol"], float(row["coefficient"])) for row in lock["genes_and_coefficients"]]
    file_pairs = [(row["gene_symbol"], float(row["coefficient"])) for row in coefficients]
    if lock_pairs != file_pairs:
        fail("Model lock and coefficient table differ")

    table = tsv(ROOT / "results" / "tables" / "stable2_sparse4_model_definition.tsv")
    table_pairs = [(row["gene_symbol"], float(row["coefficient"])) for row in table]
    if lock_pairs != table_pairs:
        fail("Supplementary model-definition table differs from the model lock")

    scale = tsv(ROOT / "config" / "training_scale.tsv")
    selected_scale = [(row["gene_symbol"], float(row["coefficient"])) for row in scale if row["selected"] == "TRUE"]
    if lock_pairs != selected_scale:
        fail("Selected training-scale rows differ from the model lock")

    candidates = tsv(ROOT / "config" / "candidate_provenance.tsv")
    if len(candidates) != 539 or sum(row["eligible_for_model"] == "TRUE" for row in candidates) != 492:
        fail("Candidate provenance counts differ from the locked contract")

    mapping = list(csv.DictReader((ROOT / "docs" / "figure_table_map.csv").open(encoding="utf-8", newline="")))
    figure_rows = [row for row in mapping if row["item_type"] == "figure_panel"]
    if len(figure_rows) != 72 or len({row["paper_item"] for row in figure_rows}) != 72:
        fail("The figure map does not contain 72 unique panel rows")
    permitted_without_asset = {"Fig7A", "Fig7D"}
    for row in figure_rows:
        outputs = [item for item in row["output_file"].split(";") if item]
        if not outputs and row["paper_item"] not in permitted_without_asset:
            fail("Mapped panel lacks a released asset: " + row["paper_item"])
        for output in outputs:
            if not (ROOT / output).is_file():
                fail("Mapped panel output is missing: " + output)

    manifest = list(csv.DictReader((ROOT / "results" / "manifest.csv").open(encoding="utf-8", newline="")))
    manifest_files = {row["output_file"] for row in manifest}
    actual_result_files = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "results").rglob("*")
        if path.is_file() and path.name != "manifest.csv"
    }
    if manifest_files != actual_result_files:
        fail("Results manifest does not exactly cover released result files")
    for row in manifest:
        path = ROOT / row["output_file"]
        if sha256(path) != row["sha256"]:
            fail("Results manifest hash mismatch: " + row["output_file"])

    source_provenance = list(csv.DictReader((ROOT / "provenance" / "source_provenance.csv").open(encoding="utf-8", newline="")))
    for row in source_provenance:
        path = ROOT / row["release_path"]
        if not path.is_file() or sha256(path) != row["release_sha256"]:
            fail("Source provenance release hash mismatch: " + row["release_path"])

    for path in ROOT.rglob("*"):
        if ".git" in path.parts or not path.is_file():
            continue
        relative = path.relative_to(ROOT).as_posix()
        if any(part.startswith(".") and part not in {".gitignore"} for part in path.relative_to(ROOT).parts):
            fail("Unexpected hidden file: " + relative)
        if path.stat().st_size > 50 * 1024 * 1024:
            fail("File exceeds 50 MB release cap: " + relative)
        if path.suffix.lower() in {".rds", ".h5ad", ".h5", ".rdata", ".key", ".pem"}:
            fail("Disallowed file type: " + relative)
        text_extensions = {".r", ".py", ".md", ".txt", ".tsv", ".csv", ".yml", ".yaml", ".json", ".cff", ".sh"}
        if path.suffix.lower() in text_extensions or path.name in {"DESCRIPTION", "VERSION", "LICENSE"}:
            text = path.read_text(encoding="utf-8")
            if "\r" in text:
                fail("Non-LF text file: " + relative)
            if re.search(r"/(?:home|Users)/", text):
                fail("Local absolute path found: " + relative)
            if re.search(r"(?:AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9_]{20,}|BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY)", text):
                fail("Credential-like material found: " + relative)

    for path in (ROOT / "scripts").rglob("*.py"):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for path in (ROOT / "archive_exploratory").rglob("*.py"):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    for path in (ROOT / "scripts").iterdir():
        if path.is_file() and path.suffix.lower() in {".r", ".py"}:
            head = path.read_text(encoding="utf-8")[:1600]
            for label in ("Script   :", "Purpose  :", "Inputs   :", "Outputs  :", "Paper    :", "Seed     :", "Runtime  :", "Last run :"):
                if label not in head:
                    fail(f"Script header is incomplete in {path.name}: missing {label}")

    r_files = [str(path) for path in (ROOT / "scripts").rglob("*") if path.is_file() and path.suffix.lower() == ".r"]
    r_files += [str(path) for path in (ROOT / "archive_exploratory").rglob("*") if path.is_file() and path.suffix.lower() == ".r"]
    expression = "for (p in commandArgs(TRUE)) parse(file=p); cat('R parse: PASS\\n')"
    subprocess.run(["Rscript", "-e", expression, *r_files], check=True, cwd=ROOT)
    print("Repository release checks: PASS")


def full_input_checks() -> None:
    expected = [
        "data/processed/TCGA-HNSC/vst_top5000.rds",
        "data/processed/TCGA-HNSC/clinical_survival.csv",
        "data/processed/GSE164690/seurat_annotated.rds",
        "data/processed/GSE208253/spatial_qc.h5ad",
        "data/processed/GSE208253/sample_manifest.tsv",
        "data/processed/TCGA-HNSC/yellow_module_members.csv",
        "data/raw/GSE41613/GSE41613_series_matrix.txt.gz",
        "data/raw/GSE65858/GSE65858_series_matrix.txt.gz",
        "data/raw/GSE159067/GSE159067_expression.txt.gz",
    ]
    missing = [path for path in expected if not (ROOT / path).is_file()]
    if missing:
        fail("Required local workflow inputs are missing:\n" + "\n".join(missing))

    sparse_contract = json.loads((ROOT / "config" / "sparse4_analysis_contract.json").read_text(encoding="utf-8"))
    unresolved_sparse_hashes = [
        role for role, item in sparse_contract["inputs"].items()
        if item.get("sha256") == "record_after_reconstruction"
    ]
    if unresolved_sparse_hashes:
        fail("Sparse4 input hashes must be recorded before full execution: " + ", ".join(unresolved_sparse_hashes))

    cellchat_contract = tsv(ROOT / "config" / "cellchat_input_authority.tsv")
    unresolved_cellchat = [row["key"] for row in cellchat_contract if row["value"] == "record_after_reconstruction"]
    if unresolved_cellchat:
        fail("CellChat input hashes must be recorded before full execution: " + ", ".join(unresolved_cellchat))
    print("Full input gate: PASS")


parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["release", "full-inputs"], default="release")
args = parser.parse_args()
release_checks()
if args.mode == "full-inputs":
    full_input_checks()
