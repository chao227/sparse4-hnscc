#!/usr/bin/env python3

# ==========================================================
# Script   : 16_spatial_yellow_tls.py
# Purpose  : Compute patient-aware spatial Yellow-TLS associations
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: spatial analysis; Fig. 3; Supplementary Fig. 4 and Table 7
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/04_spatial/run_yellow_module_spatial.py
# ==========================================================
"""Produce the prespecified Yellow-module spatial projection analysis.

This script implements docs/spatial_analysis_contract.md. It uses the complete
historical Yellow module, performs no outcome-dependent filtering, and treats
patients (not spots) as the inferential units.
"""

from __future__ import annotations

import hashlib
import itertools
import os
import json
import platform
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.stats import rankdata, spearmanr


ROOT = Path(os.environ.get("SPARSE4_REPO_ROOT", Path.cwd())).resolve()
BASE = ROOT / "results" / "spatial_yellow_tls"
PROD = BASE
OUT = PROD / "source_data"
REPORTS = PROD / "reports"
LOGS = PROD / "logs"
for directory in (OUT, REPORTS, LOGS):
    directory.mkdir(parents=True, exist_ok=True)

H5AD = ROOT / "data" / "processed" / "GSE208253" / "spatial_qc.h5ad"
MANIFEST = ROOT / "data" / "processed" / "GSE208253" / "sample_manifest.tsv"
WGCNA = ROOT / "data" / "processed" / "TCGA-HNSC" / "yellow_module_members.csv"
RCTD_ROOT = ROOT / "data" / "processed" / "GSE208253" / "rctd"
CONTRACT = ROOT / "docs" / "spatial_analysis_contract.md"

TLS10 = ["CXCL13", "CCL19", "CCL21", "LTB", "LTA", "CXCR5", "CCR7", "CD79A", "ICAM1", "VCAM1"]
IKZF1 = "IKZF1"
SPARSE4 = ["CAMK2N1", "MASP1", "MIR9-3HG", "MS4A1"]
ANCHORS = ["CD79A", "MS4A1", "CD19", "CD79B", "JCHAIN"]
OFFSETS = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dense_cols(matrix, indices: list[int]) -> np.ndarray:
    selected = matrix[:, indices]
    return selected.toarray() if sparse.issparse(selected) else np.asarray(selected)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) != len(y) or len(x) < 2 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(spearmanr(x, y).statistic)


def rank_z_population(x: np.ndarray) -> np.ndarray:
    ranks = rankdata(np.asarray(x, float), method="average")
    sd = ranks.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        raise RuntimeError("Bivariate Moran input has zero rank variance")
    return (ranks - ranks.mean()) / sd


def visium_graph(rows: np.ndarray, cols: np.ndarray) -> tuple[csr_matrix, int]:
    coordinates = np.column_stack([rows, cols]).astype(int)
    if len({tuple(v) for v in coordinates}) != len(coordinates):
        raise RuntimeError("Duplicate Visium array coordinates within a section")
    lookup = {tuple(v): i for i, v in enumerate(coordinates)}
    edges: list[tuple[int, int]] = []
    for i, (row, col) in enumerate(coordinates):
        for drow, dcol in OFFSETS:
            j = lookup.get((row + drow, col + dcol))
            if j is not None and i < j:
                edges.append((i, j))
    ii = np.asarray([v for edge in edges for v in (edge[0], edge[1])], dtype=int)
    jj = np.asarray([v for edge in edges for v in (edge[1], edge[0])], dtype=int)
    graph = csr_matrix((np.ones(len(ii)), (ii, jj)), shape=(len(coordinates), len(coordinates)))
    return graph, len(edges)


def bivariate_moran(x: np.ndarray, y: np.ndarray, graph: csr_matrix) -> float:
    return float(rank_z_population(x) @ (graph @ rank_z_population(y)) / graph.sum())


def exact_signflip(values: np.ndarray) -> tuple[float, int, int]:
    values = np.asarray(values, float)
    if not np.isfinite(values).all():
        raise RuntimeError("Primary patient effects contain non-finite values")
    observed = abs(values.mean())
    extreme = 0
    total = 2 ** len(values)
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        permuted = abs((values * np.asarray(signs)).mean())
        extreme += int(permuted >= observed - 1e-15)
    return extreme / total, extreme, total


def natural_section(value: str) -> int:
    return int(str(value).upper().removeprefix("S"))


def average_sections_to_patients(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict] = []
    for patient_id, group in frame.groupby("patient_id", sort=True):
        group = group.sort_values("section_id", key=lambda x: x.map(natural_section))
        row = {
            "patient_id": patient_id,
            "n_sections": len(group),
            "section_ids": ";".join(group["section_id"].astype(str)),
        }
        for column in columns:
            row[column] = float(group[column].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def load_rctd(manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    qc_rows: list[dict] = []
    for record in manifest.itertuples(index=False):
        matches = sorted(RCTD_ROOT.glob(f"{record.GSM}_*/rctd_multi_weights.tsv"))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one RCTD weight table for {record.GSM}, found {len(matches)}")
        weight_path = matches[0]
        metadata_path = weight_path.parent / "spot_metadata.tsv"
        weights = pd.read_csv(weight_path, sep="\t", dtype={"barcode": str})
        metadata = pd.read_csv(metadata_path, sep="\t", dtype={"barcode": str})
        if "B_cell" not in weights or weights["barcode"].duplicated().any() or metadata["barcode"].duplicated().any():
            raise RuntimeError(f"Invalid RCTD schema or duplicate barcode for {record.GSM}")
        joined = weights[["barcode", "B_cell"]].merge(metadata[["barcode"]], on="barcode", validate="one_to_one")
        if len(joined) != len(weights) or len(joined) != len(metadata):
            raise RuntimeError(f"RCTD weights/metadata barcode mismatch for {record.GSM}")
        if joined["B_cell"].isna().any() or (joined["B_cell"] < 0).any():
            raise RuntimeError(f"Invalid RCTD B-cell weights for {record.GSM}")
        joined.insert(0, "GSM", str(record.GSM))
        rows.append(joined)
        qc_rows.append({
            "GSM": record.GSM,
            "patient_id": record.patient_id,
            "section_id": record.section_id,
            "rctd_directory": str(weight_path.parent.relative_to(ROOT)),
            "n_weight_rows": len(weights),
            "n_metadata_rows": len(metadata),
            "n_joined_rows": len(joined),
            "b_cell_min": float(joined["B_cell"].min()),
            "b_cell_max": float(joined["B_cell"].max()),
            "b_cell_mean": float(joined["B_cell"].mean()),
        })
    return pd.concat(rows, ignore_index=True), pd.DataFrame(qc_rows)


def main() -> None:
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)
    manifest = manifest.loc[manifest["eligible"].eq("YES")].copy()
    if len(manifest) != 12 or manifest["GSM"].nunique() != 12 or manifest["patient_id"].nunique() != 10:
        raise RuntimeError("Expected 12 eligible sections from 10 independent patients")

    hubs = pd.read_csv(WGCNA)
    yellow = hubs.loc[hubs["module"].astype(str).str.lower().eq("yellow"), ["gene", "module", "kME", "is_top20"]].copy()
    yellow["gene"] = yellow["gene"].astype(str)
    yellow["kME"] = pd.to_numeric(yellow["kME"], errors="coerce")
    if len(yellow) != 539 or yellow["gene"].nunique() != 539 or not np.isfinite(yellow["kME"]).all():
        raise RuntimeError("Historical Yellow membership is not exactly 539 unique genes with finite kME")
    if not set(ANCHORS).issubset(set(yellow["gene"])):
        raise RuntimeError("One or more orientation-anchor genes are absent from Yellow")
    anchor_kme = yellow.set_index("gene").loc[ANCHORS, "kME"]
    if not (anchor_kme > 0).all():
        raise RuntimeError("Historical eigengene orientation is not positive for every frozen anchor")

    a = ad.read_h5ad(H5AD)
    if a.n_obs != 23882 or a.obs["GSM"].astype(str).nunique() != 12:
        raise RuntimeError("Frozen spatial object dimensions changed")
    if "qc_pass" not in a.obs or not a.obs["qc_pass"].astype(bool).all():
        raise RuntimeError("Spatial object contains a non-QC-retained spot")
    if set(a.obs["GSM"].astype(str)) != set(manifest["GSM"]):
        raise RuntimeError("Spatial object and manifest GSM sets differ")
    if "gene_symbol" not in a.var or "counts" not in a.layers:
        raise RuntimeError("Spatial object lacks gene_symbol or raw counts")

    symbols = a.var["gene_symbol"].astype(str).to_numpy()
    symbol_to_hits = {symbol: np.flatnonzero(symbols == symbol) for symbol in set(yellow["gene"]) | set(TLS10) | {IKZF1}}
    mapping_rows: list[dict] = []
    mapped_genes: list[str] = []
    mapped_indices: list[int] = []
    for row in yellow.itertuples(index=False):
        hits = symbol_to_hits[row.gene]
        status = "UNIQUE_MAPPED" if len(hits) == 1 else ("UNMAPPED" if len(hits) == 0 else "AMBIGUOUS_MULTIMAP")
        mapping_rows.append({
            "gene_symbol": row.gene,
            "module": row.module,
            "kME": float(row.kME),
            "is_top20": bool(row.is_top20),
            "n_spatial_features": len(hits),
            "mapping_status": status,
            "var_index": int(hits[0]) if len(hits) == 1 else np.nan,
            "gene_id": str(a.var.iloc[hits[0]]["gene_id"]) if len(hits) == 1 else "",
            "included_primary": len(hits) == 1,
        })
        if len(hits) == 1:
            mapped_genes.append(row.gene)
            mapped_indices.append(int(hits[0]))
    mapping = pd.DataFrame(mapping_rows)
    if len(mapped_genes) != 520 or mapping["mapping_status"].eq("UNMAPPED").sum() != 19 or mapping["mapping_status"].eq("AMBIGUOUS_MULTIMAP").any():
        raise RuntimeError("Expected exactly 520 unique mappings and 19 unmapped Yellow genes")
    mapping.to_csv(OUT / "YELLOW_SPATIAL_FEATURE_MAPPING.tsv", sep="\t", index=False, float_format="%.12g")

    outcome_indices: dict[str, int] = {}
    for gene in TLS10 + [IKZF1]:
        hits = symbol_to_hits[gene]
        if len(hits) != 1:
            raise RuntimeError(f"Expected one spatial feature for {gene}, found {len(hits)}")
        outcome_indices[gene] = int(hits[0])

    yellow_expr = dense_cols(a.X, mapped_indices).astype(float)
    yellow_counts = dense_cols(a.layers["counts"], mapped_indices).astype(float)
    outcome_genes = TLS10 + [IKZF1]
    outcome_expr = dense_cols(a.X, [outcome_indices[g] for g in outcome_genes]).astype(float)
    outcome_counts = dense_cols(a.layers["counts"], [outcome_indices[g] for g in outcome_genes]).astype(float)
    all_expr = np.column_stack([yellow_expr, outcome_expr])
    all_counts = np.column_stack([yellow_counts, outcome_counts])
    totals = a.obs["total_counts"].astype(float).to_numpy()
    reconstructed = np.log1p(all_counts / totals[:, None] * 10000.0)
    transform_delta = float(np.max(np.abs(all_expr - reconstructed)))
    if transform_delta > 2e-5:
        raise RuntimeError(f"X is not registered CPM10000 log1p; maximum delta={transform_delta}")

    obs = a.obs[["GSM", "section_id", "barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]].copy().reset_index(drop=True)
    obs = obs.rename(columns={"pxl_col_in_fullres": "pixel_col_fullres", "pxl_row_in_fullres": "pixel_row_fullres"})
    for column in ("GSM", "section_id", "barcode"):
        obs[column] = obs[column].astype(str)
    obs["section_id"] = obs["section_id"].str.upper()
    obs = obs.merge(manifest[["GSM", "patient_id"]], on="GSM", how="left", validate="many_to_one")
    if obs["patient_id"].isna().any() or obs[["GSM", "barcode"]].duplicated().any():
        raise RuntimeError("Incomplete or duplicated spatial spot identity")

    rctd, rctd_qc = load_rctd(manifest)
    spot_identity = obs[["GSM", "barcode"]].copy()
    rctd_join = spot_identity.merge(rctd, on=["GSM", "barcode"], how="left", validate="one_to_one")
    if rctd_join["B_cell"].isna().any() or len(rctd_join) != len(obs):
        raise RuntimeError("RCTD B-cell weights do not map exactly to every frozen spot")
    rctd_qc.to_csv(OUT / "YELLOW_SPATIAL_RCTD_JOIN_QC.tsv", sep="\t", index=False, float_format="%.12g")

    kme = yellow.set_index("gene").loc[mapped_genes, "kME"].to_numpy(float)
    denominator = float(np.abs(kme).sum())
    signs = np.sign(kme)
    if denominator <= 0 or (signs == 0).any():
        raise RuntimeError("Invalid mapped Yellow kME weights")
    sparse4_excluded = [g for g in SPARSE4 if g in mapped_genes]
    keep_excluded = np.asarray([g not in SPARSE4 for g in mapped_genes], dtype=bool)
    excluded_denominator = float(np.abs(kme[keep_excluded]).sum())
    if set(sparse4_excluded) != set(SPARSE4) or excluded_denominator <= 0:
        raise RuntimeError("All four frozen model genes must be uniquely mapped before exclusion sensitivity")
    if "CD79A" not in mapped_genes:
        raise RuntimeError("Post-result shared-gene diagnostic requires uniquely mapped CD79A")
    keep_cd79a_excluded = np.asarray([g != "CD79A" for g in mapped_genes], dtype=bool)
    cd79a_excluded_denominator = float(np.abs(kme[keep_cd79a_excluded]).sum())
    if cd79a_excluded_denominator <= 0:
        raise RuntimeError("Invalid Yellow-minus-CD79A denominator")

    n_spots, n_genes = yellow_expr.shape
    yellow_z = np.zeros((n_spots, n_genes), dtype=np.float64)
    primary = np.empty(n_spots, dtype=float)
    sign_only = np.empty(n_spots, dtype=float)
    excluded = np.empty(n_spots, dtype=float)
    cd79a_excluded = np.empty(n_spots, dtype=float)
    tls10 = outcome_expr[:, :len(TLS10)].mean(axis=1)
    tls9_indices = [j for j, gene in enumerate(TLS10) if gene != "CD79A"]
    tls9_cd79a_excluded = outcome_expr[:, tls9_indices].mean(axis=1)
    ikzf1 = outcome_expr[:, -1]
    zero_rows: list[dict] = []
    detection_rows: list[dict] = []
    standardization_rows: list[dict] = []
    section_index: dict[tuple[str, str], np.ndarray] = {}

    grouped = obs.groupby(["GSM", "section_id"], sort=True).groups
    for (gsm, section_id), raw_idx in grouped.items():
        idx = np.asarray(list(raw_idx), dtype=int)
        section_index[(gsm, section_id)] = idx
        matrix = yellow_expr[idx, :]
        means = matrix.mean(axis=0)
        sds = matrix.std(axis=0, ddof=1)
        zero_var = ~np.isfinite(sds) | (sds <= 0)
        usable_sd = sds.copy()
        usable_sd[zero_var] = 1.0
        z = (matrix - means) / usable_sd
        z[:, zero_var] = 0.0
        yellow_z[idx, :] = z
        primary[idx] = (z @ kme) / denominator
        sign_only[idx] = (z @ signs) / len(mapped_genes)
        excluded[idx] = (z[:, keep_excluded] @ kme[keep_excluded]) / excluded_denominator
        cd79a_excluded[idx] = (z[:, keep_cd79a_excluded] @ kme[keep_cd79a_excluded]) / cd79a_excluded_denominator
        patient_id = str(obs.loc[idx[0], "patient_id"])
        detected = yellow_counts[idx, :] > 0
        for j, gene in enumerate(mapped_genes):
            standardization_rows.append({
                "GSM": gsm, "patient_id": patient_id, "section_id": section_id,
                "gene_symbol": gene, "kME": kme[j], "n_spots": len(idx),
                "section_mean_log1p": means[j], "section_sample_sd_log1p": sds[j],
                "zero_variance": bool(zero_var[j]),
            })
            detection_rows.append({
                "GSM": gsm, "patient_id": patient_id, "section_id": section_id,
                "gene_symbol": gene, "n_spots": len(idx),
                "n_detected": int(detected[:, j].sum()),
                "detection_rate": float(detected[:, j].mean()),
                "zero_variance": bool(zero_var[j]),
            })
            if zero_var[j]:
                zero_rows.append({
                    "GSM": gsm, "patient_id": patient_id, "section_id": section_id,
                    "gene_symbol": gene, "n_spots": len(idx),
                    "constant_log1p_value": float(means[j]),
                    "n_detected": int(detected[:, j].sum()),
                })

    if not all(np.isfinite(x).all() for x in (primary, sign_only, excluded, cd79a_excluded, tls10, tls9_cd79a_excluded, ikzf1)):
        raise RuntimeError("Non-finite activity or outcome values produced")
    standardization = pd.DataFrame(standardization_rows)
    detection = pd.DataFrame(detection_rows)
    zero_variance = pd.DataFrame(zero_rows, columns=["GSM", "patient_id", "section_id", "gene_symbol", "n_spots", "constant_log1p_value", "n_detected"])
    standardization.to_csv(OUT / "YELLOW_SPATIAL_SECTION_STANDARDIZATION.tsv", sep="\t", index=False, float_format="%.12g")
    detection.to_csv(OUT / "YELLOW_SPATIAL_DETECTION_SECTION.tsv", sep="\t", index=False, float_format="%.12g")
    zero_variance.to_csv(OUT / "YELLOW_SPATIAL_ZERO_VARIANCE_EVENTS.tsv", sep="\t", index=False, float_format="%.12g")

    detection_summary = detection.groupby("gene_symbol", as_index=False).agg(
        total_spots=("n_spots", "sum"),
        total_detected=("n_detected", "sum"),
        sections_detected=("n_detected", lambda x: int((x > 0).sum())),
        zero_variance_sections=("zero_variance", "sum"),
        median_section_detection_rate=("detection_rate", "median"),
        min_section_detection_rate=("detection_rate", "min"),
        max_section_detection_rate=("detection_rate", "max"),
    )
    detection_summary["overall_detection_rate"] = detection_summary["total_detected"] / detection_summary["total_spots"]
    detection_summary = mapping.loc[mapping["included_primary"], ["gene_symbol", "kME"]].merge(detection_summary, on="gene_symbol", validate="one_to_one")
    detection_summary.to_csv(OUT / "YELLOW_SPATIAL_DETECTION_SUMMARY.tsv", sep="\t", index=False, float_format="%.12g")

    spot = obs.copy()
    spot["yellow_activity_kME"] = primary
    spot["yellow_activity_sign_only"] = sign_only
    spot["yellow_activity_sparse4_excluded"] = excluded
    spot["yellow_activity_CD79A_excluded"] = cd79a_excluded
    spot["TLS10_score"] = tls10
    spot["TLS9_CD79A_excluded_score"] = tls9_cd79a_excluded
    spot["IKZF1_log1p_norm"] = ikzf1
    spot["RCTD_B_cell_weight"] = rctd_join["B_cell"].to_numpy(float)
    spot.to_csv(
        OUT / "YELLOW_SPATIAL_SPOT_LEVEL.tsv.gz", sep="\t", index=False,
        compression={"method": "gzip", "mtime": 0}, float_format="%.17g",
    )

    association_rows: list[dict] = []
    moran_rows: list[dict] = []
    adjacency_rows: list[dict] = []
    contribution_rows: list[dict] = []
    for (gsm, section_id), idx in section_index.items():
        patient_id = str(obs.loc[idx[0], "patient_id"])
        p = primary[idx]
        so = sign_only[idx]
        ex = excluded[idx]
        cdex = cd79a_excluded[idx]
        t = tls10[idx]
        t9 = tls9_cd79a_excluded[idx]
        ik = ikzf1[idx]
        bc = rctd_join.loc[idx, "B_cell"].to_numpy(float)
        association_rows.append({
            "GSM": gsm, "patient_id": patient_id, "section_id": section_id, "n_spots": len(idx),
            "yellow_TLS10_rho": safe_spearman(p, t),
            "yellow_IKZF1_rho": safe_spearman(p, ik),
            "yellow_RCTD_B_cell_rho": safe_spearman(p, bc),
            "IKZF1_TLS10_rho": safe_spearman(ik, t),
            "sign_only_TLS10_rho": safe_spearman(so, t),
            "sparse4_excluded_TLS10_rho": safe_spearman(ex, t),
            "yellow_minus_CD79A_TLS9_minus_CD79A_rho": safe_spearman(cdex, t9),
            "primary_sign_only_activity_rho": safe_spearman(p, so),
            "primary_sparse4_excluded_activity_rho": safe_spearman(p, ex),
            "inference_boundary": "section effect; patients are independent units",
        })
        graph, n_edges = visium_graph(obs.loc[idx, "array_row"].to_numpy(), obs.loc[idx, "array_col"].to_numpy())
        degree = np.asarray(graph.sum(axis=1)).ravel()
        adjacency_rows.append({
            "GSM": gsm, "patient_id": patient_id, "section_id": section_id, "n_spots": len(idx),
            "undirected_edges": n_edges, "directed_weight_sum": float(graph.sum()),
            "mean_degree": float(degree.mean()), "minimum_degree": int(degree.min()),
            "maximum_degree": int(degree.max()), "isolated_spots": int((degree == 0).sum()),
            "neighbor_offsets": "(0,+/-2);(+/-1,+/-1)",
        })
        moran_rows.append({
            "GSM": gsm, "patient_id": patient_id, "section_id": section_id, "n_spots": len(idx),
            "undirected_edges": n_edges,
            "yellow_focal_TLS10_neighbor_bivariate_Moran_I": bivariate_moran(p, t, graph),
        })
        variance = float(np.var(p, ddof=1))
        if not np.isfinite(variance) or variance <= 0:
            raise RuntimeError(f"Primary Yellow activity has zero variance in {gsm}")
        terms = yellow_z[idx, :] * (kme / denominator)
        covariance = ((terms - terms.mean(axis=0)) * (p - p.mean())[:, None]).sum(axis=0) / (len(idx) - 1)
        additive = covariance / variance
        for j, gene in enumerate(mapped_genes):
            contribution_rows.append({
                "GSM": gsm, "patient_id": patient_id, "section_id": section_id, "gene_symbol": gene,
                "kME": kme[j], "normalized_weight": kme[j] / denominator,
                "n_spots": len(idx), "term_variance": float(np.var(terms[:, j], ddof=1)),
                "term_activity_rho": safe_spearman(terms[:, j], p),
                "additive_variance_fraction": float(additive[j]),
                "detection_rate": float((yellow_counts[idx, j] > 0).mean()),
                "zero_variance": bool(np.var(yellow_z[idx, j], ddof=1) <= 0),
            })

    associations = pd.DataFrame(association_rows).sort_values("section_id", key=lambda x: x.map(natural_section))
    association_columns = [
        "yellow_TLS10_rho", "yellow_IKZF1_rho", "yellow_RCTD_B_cell_rho", "IKZF1_TLS10_rho",
        "sign_only_TLS10_rho", "sparse4_excluded_TLS10_rho", "primary_sign_only_activity_rho",
        "primary_sparse4_excluded_activity_rho", "yellow_minus_CD79A_TLS9_minus_CD79A_rho",
    ]
    if associations[association_columns].isna().any().any():
        raise RuntimeError("One or more section-level correlations are not estimable")
    associations.to_csv(OUT / "YELLOW_SPATIAL_SECTION_ASSOCIATIONS.tsv", sep="\t", index=False, float_format="%.12g")
    patients = average_sections_to_patients(associations, association_columns)
    patients.to_csv(OUT / "YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv", sep="\t", index=False, float_format="%.12g")

    moran = pd.DataFrame(moran_rows).sort_values("section_id", key=lambda x: x.map(natural_section))
    moran.to_csv(OUT / "YELLOW_SPATIAL_SECTION_MORAN.tsv", sep="\t", index=False, float_format="%.12g")
    patient_moran = average_sections_to_patients(moran, ["yellow_focal_TLS10_neighbor_bivariate_Moran_I"])
    patient_moran.to_csv(OUT / "YELLOW_SPATIAL_PATIENT_MORAN.tsv", sep="\t", index=False, float_format="%.12g")
    pd.DataFrame(adjacency_rows).sort_values("section_id", key=lambda x: x.map(natural_section)).to_csv(
        OUT / "YELLOW_SPATIAL_ADJACENCY_QC.tsv", sep="\t", index=False, float_format="%.12g"
    )

    contributions = pd.DataFrame(contribution_rows)
    additive_error = float(contributions.groupby(["GSM", "section_id"])["additive_variance_fraction"].sum().sub(1).abs().max())
    if additive_error > 5e-10:
        raise RuntimeError(f"Gene contribution decomposition failed; maximum error={additive_error}")
    contributions.to_csv(OUT / "YELLOW_SPATIAL_GENE_CONTRIBUTION_SECTION.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0}, float_format="%.10g")
    contribution_patient_rows: list[dict] = []
    contribution_columns = ["term_variance", "term_activity_rho", "additive_variance_fraction", "detection_rate"]
    for (patient_id, gene), group in contributions.groupby(["patient_id", "gene_symbol"], sort=True):
        group = group.sort_values("section_id", key=lambda x: x.map(natural_section))
        row = {
            "patient_id": patient_id, "gene_symbol": gene, "n_sections": len(group),
            "section_ids": ";".join(group["section_id"]), "kME": float(group["kME"].iloc[0]),
            "normalized_weight": float(group["normalized_weight"].iloc[0]),
            "zero_variance_sections": int(group["zero_variance"].sum()),
        }
        row.update({column: float(group[column].mean()) for column in contribution_columns})
        contribution_patient_rows.append(row)
    contribution_patients = pd.DataFrame(contribution_patient_rows)
    contribution_patients.to_csv(OUT / "YELLOW_SPATIAL_GENE_CONTRIBUTION_PATIENT.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0}, float_format="%.10g")
    contribution_summary = contribution_patients.groupby("gene_symbol", as_index=False).agg(
        kME=("kME", "first"), normalized_weight=("normalized_weight", "first"),
        median_patient_term_activity_rho=("term_activity_rho", "median"),
        median_patient_additive_variance_fraction=("additive_variance_fraction", "median"),
        median_patient_detection_rate=("detection_rate", "median"),
        total_zero_variance_patient_sections=("zero_variance_sections", "sum"),
    )
    contribution_summary["absolute_median_additive_fraction"] = contribution_summary["median_patient_additive_variance_fraction"].abs()
    contribution_summary = contribution_summary.sort_values("absolute_median_additive_fraction", ascending=False)
    contribution_summary.to_csv(OUT / "YELLOW_SPATIAL_GENE_CONTRIBUTION_SUMMARY.tsv", sep="\t", index=False, float_format="%.12g")
    never_detected_genes = int((detection_summary["total_detected"] == 0).sum())
    top_contribution_genes = contribution_summary.head(10)["gene_symbol"].astype(str).tolist()

    primary_values = patients["yellow_TLS10_rho"].to_numpy(float)
    p_value, extreme, total_permutations = exact_signflip(primary_values)
    primary_inference = {
        "analysis": "primary kME-weighted Yellow activity versus TLS10",
        "n_independent_patients": len(patients),
        "n_sections": len(associations),
        "mean_patient_rho": float(primary_values.mean()),
        "median_patient_rho": float(np.median(primary_values)),
        "minimum_patient_rho": float(primary_values.min()),
        "maximum_patient_rho": float(primary_values.max()),
        "positive_patients": int((primary_values > 0).sum()),
        "negative_patients": int((primary_values < 0).sum()),
        "zero_patients": int((primary_values == 0).sum()),
        "exact_signflip_extreme_assignments": extreme,
        "exact_signflip_total_assignments": total_permutations,
        "exact_signflip_two_sided_p": p_value,
        "formal_test_family": "one prespecified primary test only",
    }
    pd.DataFrame([primary_inference]).to_csv(OUT / "YELLOW_SPATIAL_PRIMARY_INFERENCE.tsv", sep="\t", index=False, float_format="%.12g")

    sensitivity_rows = []
    for analysis, column in [
        ("sign-only all mapped Yellow genes", "sign_only_TLS10_rho"),
        ("primary kME activity excluding frozen sparse4 genes", "sparse4_excluded_TLS10_rho"),
    ]:
        values = patients[column].to_numpy(float)
        sensitivity_rows.append({
            "analysis": analysis, "n_independent_patients": len(values),
            "mean_patient_rho": float(values.mean()), "median_patient_rho": float(np.median(values)),
            "minimum_patient_rho": float(values.min()), "maximum_patient_rho": float(values.max()),
            "positive_patients": int((values > 0).sum()), "negative_patients": int((values < 0).sum()),
            "formal_p_value": "NOT_CALCULATED__DESCRIPTIVE_SENSITIVITY",
        })
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(OUT / "YELLOW_SPATIAL_SENSITIVITY_SUMMARY.tsv", sep="\t", index=False, float_format="%.12g")

    diagnostic_values = patients["yellow_minus_CD79A_TLS9_minus_CD79A_rho"].to_numpy(float)
    shared_gene_diagnostic = {
        "analysis": "post-result shared-gene diagnostic: Yellow minus CD79A versus TLS9 minus CD79A",
        "protocol_status": "POST_RESULT_AMENDMENT__NOT_PRESPECIFIED",
        "n_independent_patients": len(diagnostic_values),
        "mean_patient_rho": float(diagnostic_values.mean()),
        "median_patient_rho": float(np.median(diagnostic_values)),
        "minimum_patient_rho": float(diagnostic_values.min()),
        "maximum_patient_rho": float(diagnostic_values.max()),
        "positive_patients": int((diagnostic_values > 0).sum()),
        "negative_patients": int((diagnostic_values < 0).sum()),
        "formal_p_value": "NOT_CALCULATED__DESCRIPTIVE_POST_RESULT_DIAGNOSTIC",
        "interpretation": "Removes direct CD79A overlap only; shared cell-composition effects remain possible.",
    }
    pd.DataFrame([shared_gene_diagnostic]).to_csv(
        OUT / "YELLOW_SPATIAL_SHARED_GENE_DIAGNOSTIC.tsv", sep="\t", index=False, float_format="%.12g"
    )

    context_rows = []
    for analysis, column in [
        ("Yellow activity versus IKZF1", "yellow_IKZF1_rho"),
        ("Yellow activity versus RCTD B-cell weight", "yellow_RCTD_B_cell_rho"),
        ("IKZF1 versus TLS10", "IKZF1_TLS10_rho"),
    ]:
        values = patients[column].to_numpy(float)
        context_rows.append({
            "analysis": analysis, "n_independent_patients": len(values),
            "mean_patient_effect": float(values.mean()), "median_patient_effect": float(np.median(values)),
            "minimum_patient_effect": float(values.min()), "maximum_patient_effect": float(values.max()),
            "positive_patients": int((values > 0).sum()), "negative_patients": int((values < 0).sum()),
            "inferential_status": "descriptive/exploratory; no formal p value",
        })
    moran_values = patient_moran["yellow_focal_TLS10_neighbor_bivariate_Moran_I"].to_numpy(float)
    context_rows.append({
        "analysis": "Yellow focal activity versus neighboring TLS10 rank-based bivariate Moran I",
        "n_independent_patients": len(moran_values),
        "mean_patient_effect": float(moran_values.mean()), "median_patient_effect": float(np.median(moran_values)),
        "minimum_patient_effect": float(moran_values.min()), "maximum_patient_effect": float(moran_values.max()),
        "positive_patients": int((moran_values > 0).sum()), "negative_patients": int((moran_values < 0).sum()),
        "inferential_status": "descriptive/exploratory; no formal p value",
    })
    pd.DataFrame(context_rows).to_csv(OUT / "YELLOW_SPATIAL_CONTEXT_SUMMARY.tsv", sep="\t", index=False, float_format="%.12g")

    input_paths = [H5AD, MANIFEST, WGCNA, CONTRACT]
    for record in manifest.itertuples(index=False):
        weight_path = next(iter(sorted(RCTD_ROOT.glob(f"{record.GSM}_*/rctd_multi_weights.tsv"))))
        input_paths.extend([weight_path, weight_path.parent / "spot_metadata.tsv"])
    input_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in input_paths}
    pd.DataFrame([{"path": path, "sha256": digest} for path, digest in input_hashes.items()]).to_csv(
        REPORTS / "INPUT_PROVENANCE.tsv", sep="\t", index=False
    )

    verdict = {
        "status": "PASS",
        "contract": "Frozen complete Yellow-module spatial activity",
        "n_yellow_input": len(yellow), "n_yellow_unique_mapped": len(mapped_genes),
        "n_yellow_unmapped": int(mapping["mapping_status"].eq("UNMAPPED").sum()),
        "n_spots": len(spot), "n_sections": len(associations), "n_independent_patients": len(patients),
        "zero_variance_gene_section_pairs": len(zero_variance),
        "mapped_genes_never_detected_in_spatial_counts": never_detected_genes,
        "transform_reconstruction_max_abs_delta": transform_delta,
        "contribution_additivity_max_abs_error": additive_error,
        "sparse4_genes_excluded_in_sensitivity": sparse4_excluded,
        "primary_inference": primary_inference,
        "sensitivity": sensitivity.to_dict(orient="records"),
        "post_result_shared_gene_diagnostic": shared_gene_diagnostic,
        "context": context_rows,
        "top_genes_by_absolute_median_additive_fraction": top_contribution_genes,
        "no_filtering_refit_retuning_smoothing_imputation_bootstrap": True,
        "claim_ceiling": "Descriptive independent-cohort spatial support for the historical Yellow module's B-cell/TLS context only; not spatial validation of sparse4 survival prediction, causal TLS biology, or histological mature-TLS confirmation.",
        "input_sha256": input_hashes,
        "software": {
            "python": platform.python_version(), "anndata": ad.__version__, "numpy": np.__version__,
            "pandas": pd.__version__, "scipy": scipy.__version__,
        },
    }
    (REPORTS / "PRODUCTION_VERDICT.json").write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")

    report = [
        "# Yellow-module spatial production report", "",
        "- Status: **PASS**.",
        f"- Frozen inputs: 539 unique Yellow genes; {len(mapped_genes)} uniquely mapped and {539-len(mapped_genes)} unmapped.",
        f"- Spatial cohort: {len(spot):,} QC-retained spots, {len(associations)} sections, {len(patients)} independent patients.",
        f"- Zero-variance gene-section pairs fixed to z=0 as prespecified: {len(zero_variance)}.",
        f"- Uniquely mapped Yellow genes never detected in any retained spatial spot: {never_detected_genes}; they remain in the fixed complete-module denominator and contribute z=0.",
        "- No outcome-driven filtering, smoothing, imputation, refitting, retuning, bootstrap, or method switching was performed.", "",
        "## Sole formal primary analysis", "",
        f"- Mean patient Yellow/TLS10 rho: {primary_inference['mean_patient_rho']:.6f}.",
        f"- Median patient Yellow/TLS10 rho: {primary_inference['median_patient_rho']:.6f}.",
        f"- Direction: {primary_inference['positive_patients']}/10 positive and {primary_inference['negative_patients']}/10 negative.",
        f"- Exact two-sided sign-flip p: {p_value:.9g} ({extreme}/{total_permutations} assignments).", "",
        "## Prespecified descriptive sensitivities", "",
    ]
    for row in sensitivity_rows:
        report.append(f"- {row['analysis']}: median patient rho {row['median_patient_rho']:.6f}; {row['positive_patients']}/10 positive. No additional formal p value was calculated.")
    report.extend(["", "## Post-result shared-gene diagnostic", "",
        f"- Yellow-minus-CD79A versus TLS9-minus-CD79A: median patient rho {shared_gene_diagnostic['median_patient_rho']:.6f}; {shared_gene_diagnostic['positive_patients']}/10 positive.",
        "- This was fixed only after Level-3 review identified the direct CD79A overlap. It is explicitly not prespecified, has no added formal p value, and does not remove shared cell-composition effects.",
        "", "## Interpretation boundary", "",
        f"- The largest absolute median contribution fractions are enriched for immunoglobulin/plasma-cell-associated genes ({', '.join(top_contribution_genes[:8])}); shared immune-cell composition is therefore a plausible explanation for the spatial association.",
        "- RCTD B-cell, IKZF1 and neighborhood/Moran results are descriptive context, not additional confirmatory tests.",
        "- Spots are measurement units; repeat sections for P7 and P8 were averaged only after section-level effects were computed.",
        "- The result cannot validate sparse4 survival prediction in space, establish causality, or identify histological mature TLS.",
    ])
    (REPORTS / "PRODUCTION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        failure = {"status": "FAIL", "error_type": type(exc).__name__, "error": str(exc)}
        (REPORTS / "PRODUCTION_VERDICT.json").write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(failure, indent=2), file=sys.stderr)
        raise
