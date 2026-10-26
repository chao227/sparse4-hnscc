#!/usr/bin/env python3

# ==========================================================
# Script   : 24_make_tables.py
# Purpose  : Build the supplementary table release from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Supplementary Tables 1-20
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/12_tables/build_supplementary_release.py
# ==========================================================
"""Build the unique manuscript-table release from guarded reback results."""

from __future__ import annotations

import hashlib
import re
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path.cwd().resolve()
OUT = ROOT / "final_audit_reconciliation_20261022" / "supplementary_release"
EDITABLE = OUT / "editable"
SOURCE = OUT / "source_data"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


WORKBOOK_METADATA = {
    "Main_Table1_current_cohorts": ("Cohort and data-source characteristics used in the revision.", "cohort/dataset", "Descriptive only.", "No multiplicity testing.", "Does not establish cohort comparability or validation."),
    "Main_Table2_common_roster_Cox": ("Common-roster Cox model coefficients for the current sparse4 analyses.", "patient", "Cox proportional-hazards estimates with model-based confidence intervals.", "P values are descriptive after model development.", "Same-cohort apparent association; not independent validation."),
    "STable1_Yellow_candidate_provenance": ("Yellow-module membership and availability in the fixed modelling feature set.", "gene", "Descriptive provenance.", "Not applicable.", "Module membership is a candidate-range definition, not model validation."),
    "STable2_sparse4_model_definition": ("Complete post hoc sparse4 definition and fixed scoring parameters.", "gene/model", "Penalized Cox model definition.", "Not applicable.", "Post hoc parsimonious model selected after exploratory results were inspected."),
    "STable3_TCGA_sparse4_scores": ("Patient-level TCGA sparse4 scores and grouping.", "patient", "Descriptive score application.", "Not applicable.", "Development-cohort scores; median groups are descriptive."),
    "STable4_GSE159067_partial_scores": ("Patient-level two-gene partial score in GSE159067.", "patient", "Exploratory partial-score application.", "Not applicable.", "Not the complete four-gene score and not LymphHot validation."),
    "STable5_VKO_unsigned_perturbation": ("Unsigned virtual-perturbation displacement results.", "gene/computational contrast", "Multiple-testing-adjusted computational displacement.", "BH correction as reported in the data sheet.", "No expression direction, causal knockout or mechanism is inferred."),
    "STable6_LymphHot_cutpoint_sensitivity": ("Exploratory LymphHot cutpoint sensitivity.", "patient/model", "Cox/log-rank sensitivity at fixed quantile cutpoints.", "Correction fields are reported in the data sheet where applicable.", "Cohort-dependent exploratory phenotype; no optimized clinical cutoff."),
    "STable7_spatial_score_definitions": ("Complete spatial Yellow and TLS10/TLS9 definitions and patient-level results.", "spot/section/patient as labelled by sheet", "Patient is the inferential unit; spots are nested observations.", "Exact sign-flip family and post-result diagnostics are labelled by sheet.", "Computational TLS-associated transcriptional context, not histological TLS."),
    "STable8_objective64_coefficients": ("Coefficients from the previously specified lambda.min procedure whose full-data fit retained 64 genes.", "gene/model", "Penalized survival model coefficients.", "Procedure-matched permutation belongs only to this model procedure.", "The procedure was specified; the 64-gene identity is the full-data fit result and is separate from post hoc sparse4."),
    "STable9_internal_validation_estimands": ("Internal-evaluation estimands for objective64 and sparse4 procedures.", "patient/repeat/procedure", "Repeated cross-validation estimands and non-bootstrap uncertainty where supplied.", "Permutation P=0.002 applies only to objective64.", "Conditional validation with fixed upstream discovery, not end-to-end discovery validation."),
    "STable10_LymphHot_TCGA_model_comparison": ("TCGA continuous and LymphHot model comparisons.", "patient", "Cox, likelihood-ratio and information-criterion comparisons.", "Families/corrections are reported in the table.", "Exploratory same-cohort comparisons."),
    "STable10b_LymphHot_formal_tests": ("Formal continuous/additive/interaction LymphHot tests and diagnostics.", "patient", "Cox interaction, likelihood-ratio, AIC and PH diagnostics.", "Families/corrections are reported in the table.", "Null incremental evidence does not prove equivalence."),
    "STable11_LymphHot_external_evaluation": ("Post hoc external reapplication of the LymphHot rule.", "patient/cohort", "Cohort-specific Cox estimates and diagnostics.", "Cohort-level families are reported in the data sheet.", "Previously inspected exploratory cohorts; not independent validation."),
    "STable12_sparse4_external_evaluation": ("Post hoc external full/partial sparse4 evaluations.", "patient/cohort", "Cohort-specific Cox estimates and BH correction.", "BH family is reported in the data sheet.", "Previously inspected cohorts; FULL and PARTIAL scores are not interchangeable."),
    "STable13_hypoxia_sensitivity": ("Post hoc sparse4-hypoxia sensitivity analyses.", "patient/gene/signature as labelled", "Spearman BH10 and fixed-roster Cox/LR/PH diagnostics.", "BH across the registered ten correlations.", "Does not prove B-cell specificity or eliminate confounding."),
    "STable14_nonbootstrap_uncertainty": ("Non-bootstrap uncertainty for conditional sparse4 performance and spatial effects.", "patient", "Patient-delete and model-based intervals; exact methods are labelled by sheet.", "No bootstrap was launched.", "Intervals match the registered conditional estimands."),
    "STable15_M4_missingness_selection": ("M4 complete-case selection and missingness diagnostics.", "patient", "Follow-up, events, denominators and standardized mean differences.", "Descriptive diagnostics; no MAR/MNAR test is claimed.", "Does not remove complete-case selection bias."),
    "STable16_B_state_annotation_evidence": ("Donor-aware reference and marker evidence for provisional B-lineage states.", "donor/state/cell as labelled", "Donor-aware confusion and marker summaries.", "Declared comparison families are reported by sheet.", "Labels remain provisional, especially GC-like and Transitional."),
    "STable17_expression_signature_dictionaries": ("Project expression-signature dictionaries.", "gene/signature", "Definitions only.", "Not applicable.", "Expression enrichment is not a measured cell fraction or official ESTIMATE score."),
    "STable18_external_mapping_endpoints": ("External gene coverage, endpoint mapping and score ledgers.", "patient/cohort/gene", "Mapping and descriptive endpoint provenance.", "Not applicable unless stated in a result sheet.", "Coverage does not itself establish transportability or validation."),
    "STable19_method_robustness_diagnostics": ("WGCNA, clustering, missing-data and continuous-functional-form diagnostics.", "gene/patient/donor as labelled", "Signed-network sensitivity, ARI, multiple-imputation sensitivity and restricted-cubic-spline diagnostics.", "BH families and feasibility gates are reported in the relevant sheets.", "Sensitivity analyses do not replace the historical discovery authority or make provisional labels definitive."),
    "STable20_composition_platform_sensitivities": ("Spatial composition, expression-signature and cross-platform external sensitivities.", "patient/section/spot/cohort as labelled", "Patient-level spatial inference, BH15 signature correlations and descriptive fixed-effect cohort synthesis.", "Exact sign-flip and BH families are labelled by sheet.", "No cell-fraction substitution, independent validation or pure platform effect is claimed."),
}


def workbook(name: str, sheets: list[tuple[str, Path | pd.DataFrame]]) -> Path:
    target = EDITABLE / f"{name}.xlsx"
    with pd.ExcelWriter(target, engine="openpyxl") as writer:
        purpose, unit, statistics, multiplicity, boundary = WORKBOOK_METADATA[name]
        notes = pd.DataFrame([
            {
                "workbook": name,
                "sheet": sheet[:31],
                "purpose": purpose,
                "independent_unit": unit,
                "statistics_or_content": statistics,
                "multiplicity": multiplicity,
                "claim_boundary": boundary,
                "source_file": str(value.relative_to(ROOT)) if isinstance(value, Path) and value.is_relative_to(ROOT) else (str(value) if isinstance(value, Path) else "generated from registered definitions"),
            }
            for sheet, value in sheets
        ])
        notes.to_excel(writer, sheet_name="Definitions", index=False)
        ws_notes = writer.book["Definitions"]
        ws_notes.freeze_panes = "A2"
        ws_notes.auto_filter.ref = ws_notes.dimensions
        for sheet, value in sheets:
            frame = read_table(value) if isinstance(value, Path) else value
            frame.to_excel(writer, sheet_name=sheet[:31], index=False)
            ws = writer.book[sheet[:31]]
            ws.freeze_panes = "A2"
            ws.auto_filter.ref = ws.dimensions
    return target


def export_flat_sources(name: str, sheets: list[tuple[str, Path | pd.DataFrame]]) -> None:
    """Export every workbook sheet as a directly inspectable TSV."""
    for sheet, value in sheets:
        frame = read_table(value) if isinstance(value, Path) else value
        frame.to_csv(SOURCE / f"{name}__{sheet}.tsv", sep="\t", index=False)


def parse_r_vector(script: Path, variable: str) -> list[str]:
    text = script.read_text(encoding="utf-8")
    match = re.search(rf"{re.escape(variable)}\s*<-\s*c\((.*?)\)", text, flags=re.S)
    if not match:
        raise RuntimeError(f"R vector not found: {variable} in {script}")
    return re.findall(r'[\"\']([^\"\']+)[\"\']', match.group(1))


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    EDITABLE.mkdir(parents=True)
    SOURCE.mkdir(parents=True)

    scan = ROOT / "scan_1019_remediation_20261019" / "output" / "tables"
    priority = ROOT / "priority_table_release_20261020" / "source_data"
    spatial = ROOT / "yellow_module_spatial_20261019" / "production" / "source_data"
    audit = ROOT / "final_audit_reconciliation_20261022" / "production"

    direct = {
        "Main_Table1_current_cohorts": priority / "Main_Table1_current_cohorts.tsv",
        "Main_Table2_common_roster_Cox": priority / "Main_Table2_common_roster_Cox.tsv",
        "STable1_Yellow_candidate_provenance": scan / "STable1_Yellow_candidate_provenance.tsv",
        "STable2_sparse4_model_definition": priority / "STable2_sparse4_model_definition.tsv",
        "STable3_TCGA_sparse4_scores": priority / "STable3_TCGA_sparse4_scores.tsv",
        "STable4_GSE159067_partial_scores": priority / "STable4_GSE159067_partial_scores.tsv",
        "STable5_VKO_unsigned_perturbation": scan / "STable5_VKO_unsigned_perturbation.tsv",
        "STable6_LymphHot_cutpoint_sensitivity": scan / "STable6_LymphHot_cutpoint_sensitivity.tsv",
        "STable10_LymphHot_TCGA_model_comparison": priority / "STable10_LymphHot_TCGA_model_comparison.tsv",
        "STable10b_LymphHot_formal_tests": priority / "STable10b_LymphHot_formal_tests.tsv",
        "STable11_LymphHot_external_evaluation": priority / "STable11_LymphHot_external_evaluation.tsv",
        "STable12_sparse4_external_evaluation": priority / "STable12_sparse4_external_evaluation.tsv",
    }
    for name, source in direct.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        copied = SOURCE / f"{name}.tsv"
        shutil.copy2(source, copied)
        workbook(name, [("data", copied)])

    # The historical source table carried two convenience flags referring to
    # retired submitted-model gene lists.  They are not part of the active
    # lambda.min sensitivity and are deliberately excluded from the release.
    objective64 = read_table(priority / "STable8_objective64_coefficients.tsv")[[
        "gene", "coefficient", "nestedCV_selection_frequency"
    ]]
    objective64.insert(0, "procedure", "previously specified lambda.min; full-data fit retained 64 genes")
    objective64_path = SOURCE / "STable8_objective64_coefficients.tsv"
    objective64.to_csv(objective64_path, sep="\t", index=False)
    workbook("STable8_objective64_coefficients", [("data", objective64_path)])

    internal_estimands = read_table(priority / "STable9_internal_validation_estimands.tsv")
    internal_estimands["procedure"] = internal_estimands["procedure"].replace({
        "Prespecified 64-gene lambda.min": "Previously specified lambda.min procedure; full-data fit retained 64 genes"
    })
    internal_path = SOURCE / "STable9_internal_validation_estimands.tsv"
    internal_estimands.to_csv(internal_path, sep="\t", index=False)
    workbook("STable9_internal_validation_estimands", [("data", internal_path)])

    tls10 = ["CXCL13", "CCL19", "CCL21", "LTB", "LTA", "CXCR5", "CCR7", "CD79A", "ICAM1", "VCAM1"]
    tls = pd.DataFrame({
        "signature": ["TLS10"] * 10 + ["TLS9_CD79A_excluded"] * 9,
        "gene_symbol": tls10 + [x for x in tls10 if x != "CD79A"],
        "weight": 1.0,
        "definition": ["unweighted mean of log1p(CPM10000) expression within each section; no section-wise gene z-scaling"] * 19,
    })
    workbook("STable7_spatial_score_definitions", [
        ("Yellow_mapping_kME", spatial / "YELLOW_SPATIAL_FEATURE_MAPPING.tsv"),
        ("TLS10_TLS9", tls),
        ("section_standardization", spatial / "YELLOW_SPATIAL_SECTION_STANDARDIZATION.tsv"),
        ("patient_associations", spatial / "YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv"),
        ("primary_inference", spatial / "YELLOW_SPATIAL_PRIMARY_INFERENCE.tsv"),
    ])
    tls.to_csv(SOURCE / "STable7_TLS10_TLS9_definitions.tsv", sep="\t", index=False)
    for src in ["YELLOW_SPATIAL_FEATURE_MAPPING.tsv", "YELLOW_SPATIAL_SECTION_STANDARDIZATION.tsv",
                "YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv", "YELLOW_SPATIAL_PRIMARY_INFERENCE.tsv"]:
        shutil.copy2(spatial / src, SOURCE / f"STable7_{src}")

    hyp = audit / "hypoxia"
    s13 = [
        ("signature_coverage", hyp / "SIGNATURE_COVERAGE.tsv"),
        ("score_correlations_BH10", hyp / "HYPOXIA_SPEARMAN_BH10.tsv"),
        ("M4_coefficients", hyp / "M4_HYPOXIA_COEFFICIENTS.tsv"),
        ("M4_likelihood_ratio", hyp / "M4_HYPOXIA_LR.tsv"),
        ("M4_PH", hyp / "M4_HYPOXIA_PH_DIAGNOSTICS.tsv"),
    ]
    workbook("STable13_hypoxia_sensitivity", s13)
    export_flat_sources("STable13_hypoxia_sensitivity", s13)

    unc = audit / "uncertainty"
    s14 = [
        ("sparse4_OOF_interval", unc / "SPARSE4_OOF_JACKKNIFE_INTERVAL.tsv"),
        ("sparse4_delete_patient", unc / "SPARSE4_OOF_DELETE_ONE_PATIENT_ESTIMATES.tsv"),
        ("spatial_intervals", unc / "SPATIAL_YELLOW_TLS10_INTERVALS.tsv"),
        ("spatial_delete_patient", unc / "SPATIAL_YELLOW_TLS10_DELETE_ONE_PATIENT_ESTIMATES.tsv"),
    ]
    workbook("STable14_nonbootstrap_uncertainty", s14)
    export_flat_sources("STable14_nonbootstrap_uncertainty", s14)

    miss = audit / "missing_annotation" / "results"
    s15 = [
        ("group_followup", miss / "M4_GROUP_EVENT_FOLLOWUP_SUMMARY.tsv"),
        ("continuous_SMD", miss / "M4_CONTINUOUS_SMD.tsv"),
        ("categorical_SMD", miss / "M4_CATEGORICAL_LEVEL_SMD.tsv"),
        ("missing_denominators", miss / "M4_MISSING_DENOMINATORS.tsv"),
        ("membership", miss / "M4_LOCKED_MEMBERSHIP_LEDGER.tsv"),
    ]
    workbook("STable15_M4_missingness_selection", s15)
    export_flat_sources("STable15_M4_missingness_selection", s15)
    s16 = [
        ("reference_confusion", miss / "B_STATE_REFERENCE_COMPACT_CONFUSION.tsv"),
        ("donor_agreement", miss / "B_STATE_DONOR_REFERENCE_AGREEMENT.tsv"),
        ("marker_contract", miss / "B_STATE_PREDEFINED_MARKER_CONTRACT.tsv"),
        ("marker_donor_state", miss / "B_STATE_MARKER_DONOR_BY_STATE.tsv"),
        ("marker_summary", miss / "B_STATE_MARKER_ACROSS_DONOR_SUMMARY.tsv"),
        ("object_scope", miss / "SINGLECELL_CURRENT_OBJECT_SCOPE.tsv"),
    ]
    workbook("STable16_B_state_annotation_evidence", s16)
    export_flat_sources("STable16_B_state_annotation_evidence", s16)

    tme_script = ROOT / "agent_work" / "production_path_repair" / "sandboxes" / "A_consistency" / "support" / "evidence_copies" / "tme" / "04a_feature_matrix.R"
    signature_rows = []
    signatures = {
        "project_immune_ssGSEA": parse_r_vector(tme_script, "immune_genes"),
        "project_stromal_ssGSEA": parse_r_vector(tme_script, "stromal_genes"),
        "simple_B_cell_naive": ["CD79A", "CD79B", "MS4A1", "BANK1", "PAX5", "FCER2", "CR2"],
        "simple_B_cell_memory": ["CD27", "AIM2", "SSPB1", "IGLC2", "JCHAIN"],
        "simple_plasma_cell": ["JCHAIN", "MZB1", "IGHG1", "IGHA1", "XBP1", "SDC1", "PRDM1"],
        "bulk_TLS12": ["CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13"],
    }
    signatures["bulk_TLS21"] = [
        "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
        "MS4A1", "CD79B", "IGHM", "IGHG1", "IGHG3", "IGKC", "JCHAIN", "MZB1", "CD38"
    ]
    for signature, genes in signatures.items():
        for order, gene in enumerate(genes, start=1):
            signature_rows.append({"signature": signature, "gene_order": order, "gene_symbol": gene,
                                   "weight": 1.0, "role": "project-defined expression enrichment; not a measured fraction"})
    signature_table = pd.DataFrame(signature_rows)
    signature_table.to_csv(SOURCE / "STable17_expression_signature_dictionaries.tsv", sep="\t", index=False)
    workbook("STable17_expression_signature_dictionaries", [("gene_sets", signature_table)])

    # STable18 must follow the current four-gene sparse4 external route.  The
    # retired d13_external_v2 package contains a different 23-gene model and
    # must never be exposed as sparse4 evidence.
    ext = ROOT / "agent_work" / "yellow_sparse_external_independent_audit_20261016"
    lymph = ROOT / "agent_work" / "lymphhot_formal_close_20261020" / "results"
    s18 = [
        ("sparse4_gene_coverage", ext / "INDEPENDENT_GENE_COVERAGE.tsv"),
        ("sparse4_endpoint_inventory", ext / "INDEPENDENT_ENDPOINT_SUMMARY.tsv"),
        ("sparse4_patient_ledger", ext / "INDEPENDENT_SAMPLE_ENDPOINT_SCORE_LEDGER.tsv"),
        ("sparse4_FULL_OS_BH", ext / "INDEPENDENT_FULL_OS_MULTIPLICITY.tsv"),
        ("LymphHot_gene_coverage", lymph / "LYMPHHOT_EXTERNAL_GENE_COVERAGE.tsv"),
        ("LymphHot_results", priority / "STable11_LymphHot_external_evaluation.tsv"),
    ]
    workbook("STable18_external_mapping_endpoints", s18)
    export_flat_sources("STable18_external_mapping_endpoints", s18)

    cluster = ROOT / "agent_work" / "cluster_resolution_concordance_20261016" / "formal_run_001" / "results"
    clinical = ROOT / "agent_work" / "sparse4_tcga_clinical_20261016"
    wgcna = ROOT / "agent_work" / "analytical_gap_wave1_20261016" / "formal_run_008" / "k_wgcna"
    hpv = ROOT / "agent_work" / "analytical_gap_wave1_20261016" / "formal_run_009" / "results"
    revision_primary = ROOT / "scan_1019_remediation_20261019" / "results" / "protocol_closure" / "REVISION_PRIMARY_3_COMPLETE.tsv"
    s19 = [
        ("cluster_resolutions", cluster / "RESOLUTION_INVENTORY.tsv"),
        ("cluster_pairwise_ARI", cluster / "PAIRWISE_ARI.tsv"),
        ("cluster_donor_coverage", cluster / "SELECTED_LABEL_DONOR_COVERAGE.tsv"),
        ("cluster_label_identity", cluster / "SELECTED_LABEL_IDENTITY.tsv"),
        ("stage_MI_pooled", clinical / "STAGE_ONLY_MI_POOLED.tsv"),
        ("stage_MI_diagnostics", clinical / "STAGE_ONLY_MI_DIAGNOSTICS.tsv"),
        ("WGCNA_target_module", wgcna / "K_TARGET_MODULE_RESOLUTION.tsv"),
        ("WGCNA_hub_definition", wgcna / "K_WGCNA_HUB_DEFINITION.tsv"),
        ("WGCNA_signed_sensitivity", wgcna / "K_WGCNA_SIGNED_SENSITIVITY.tsv"),
        ("HPV_score_gene_associations", hpv / "G_HPV_SCORE_GENE_ASSOCIATIONS.tsv"),
        ("HPV_LymphHot_association", hpv / "G_HPV_LYMPHHOT_ASSOCIATION.tsv"),
        ("HPV_LymphHot_counts", hpv / "G_HPV_LYMPHHOT_COUNTS.tsv"),
        ("HPV_interaction_feasibility", hpv / "G_HPV_INTERACTION_FEASIBILITY.tsv"),
        ("revision_primary_3", revision_primary),
    ]
    reporting = ROOT / "final_audit_reconciliation_20261022" / "wave2_reporting_closure"
    s19.extend([
        ("GEO_search_queries", reporting / "UPDATED_COHORT_SEARCH_QUERIES.tsv"),
        ("GEO_candidate_screen", reporting / "UPDATED_COHORT_CANDIDATE_SCREEN.tsv"),
    ])
    # One reader-facing map of every active multiplicity family.  This map does
    # not create new adjusted values; it identifies the estimand, family and
    # evidentiary status of values already generated in the source tables.
    multiplicity_map = pd.DataFrame([
        ("F01", "donor-aware B-state differential expression", "gene-by-state tests", "BH", "registered revision analysis", "inferential within family", "Fig.2G; Methods"),
        ("F02", "Yellow enrichment across provisional B states", "5 state tests", "BH", "registered revision analysis", "inferential within family", "Fig.2G; Methods"),
        ("F03", "complete-score external overall-survival evaluation", "2 FULL cohorts", "BH", "post hoc; cohorts inspected before model lock", "external evaluation, not validation", "STable12; STable18"),
        ("F04", "sparse4-hypoxia expression correlations", "10 gene/score-by-signature tests", "BH", "registered post hoc sensitivity", "same-cohort sensitivity", "STable13"),
        ("F05", "continuous functional-form diagnostics", "4 planned tests; feasibility gate applied", "BH", "registered before formal diagnostic run", "negative tests do not prove linearity", "STable19"),
        ("F06", "HPV relationship of sparse4 and four components", "5 tests", "BH", "registered descriptive analysis", "same-cohort description", "STable19 HPV_score_gene_associations"),
        ("F07", "paired-donor IKZF1 gene contrasts", "6 genes", "BH", "registered exploratory analysis", "association, not mechanism", "Fig.7B source table"),
        ("F08", "IKZF1 regulon state contrasts", "15 pairwise contrasts after one global test", "Tukey", "registered patient-level analysis", "method-dependent patient analysis", "Fig.5A source table"),
        ("F09", "LymphHot-versus-other sparse4 contrasts", "3 prespecified contrasts", "Holm", "registered exploratory analysis", "same-cohort description", "Fig.5F source table"),
        ("F10", "incremental delta-C comparisons", "2 comparisons", "Holm", "registered same-cohort analysis", "apparent performance", "SFig.2G source table"),
        ("F11", "marker-model delta-C comparisons", "4 comparisons", "Holm", "registered same-cohort analysis", "apparent performance", "SFig.2H source table"),
        ("F12", "LymphHot immune-feature comparisons", "4 features", "Holm", "registered exploratory analysis", "same-cohort description", "Fig.5G source table"),
        ("F13", "checkpoint-expression group comparisons", "6 genes", "BH", "registered exploratory analysis", "same-cohort description", "Fig.6B source table"),
        ("F14", "LymphHot external overall-survival evaluation", "2 evaluable cohorts", "BH", "post hoc rule reapplication", "external evaluation, not validation", "STable11"),
        ("F15", "patient-level Yellow-TLS10 primary spatial direction", "1 exact primary test", "none; single exact test", "registered primary spatial test", "patient-level association", "STable7; STable14"),
        ("F16", "CellChat direction-by-pathway endpoints", "144 descriptive endpoints", "no inferential P values", "descriptive by design", "not a cross-patient test family", "Fig.1G; Methods"),
        ("F17", "revision primary survival family", "sparse4 M1, sparse4 M4 and LymphHot M4", "Holm", "registered family restored; all estimates remain post-selection", "descriptive apparent estimates, not confirmatory inference", "Main Table2; STable19 revision_primary_3"),
        ("F18", "formal LymphHot continuous and binary comparisons", "5 prespecified Wald/LR tests", "Holm", "registered exploratory analysis", "null incremental evidence does not prove equivalence", "STable10b"),
    ], columns=[
        "family_id", "analysis", "family_members", "correction", "timing_or_status",
        "claim_boundary", "primary_report_location",
    ])
    s19.append(("multiplicity_map", multiplicity_map))
    functional = ROOT / "final_audit_reconciliation_20261022" / "wave2_functional_form_closure" / "production" / "C05_MODEL_RESULTS.tsv"
    if functional.is_file():
        functional_public = read_table(functional)
        functional_public.loc[
            functional_public["status"] == "COMPUTED_PENDING_DIFFERENT_ROLE_AUDIT", "status"
        ] = "COMPUTED"
        s19.append(("joint_marker_nonlinearity", functional_public))
    workbook("STable19_method_robustness_diagnostics", s19)
    export_flat_sources("STable19_method_robustness_diagnostics", s19)

    wave2 = ROOT / "final_audit_reconciliation_20261022" / "wave2_scientific_closure" / "production"
    s20 = [
        ("spatial_adjusted_inference", wave2 / "SPATIAL_COMPOSITION_ADJUSTED_INFERENCE.tsv"),
        ("spatial_adjusted_patient", wave2 / "SPATIAL_COMPOSITION_ADJUSTED_PATIENT.tsv"),
        ("spatial_adjusted_section", wave2 / "SPATIAL_COMPOSITION_ADJUSTED_SECTION.tsv"),
        ("spatial_join_QA", wave2 / "SPATIAL_RCTD_JOIN_QA.tsv"),
        ("fraction_availability", wave2 / "R1_11_01_DISPOSITION.tsv"),
        ("signature_coverage", wave2 / "B_LINEAGE_SIGNATURE_COVERAGE.tsv"),
        ("signature_correlations_BH15", wave2 / "SPARSE4_B_LINEAGE_SPEARMAN_BH15.tsv"),
        ("signature_patient_scores", wave2 / "B_LINEAGE_SIGNATURE_SCORES_476.tsv"),
        ("cross_platform_cohorts", wave2 / "CROSS_PLATFORM_COHORT_COMPARISON.tsv"),
        ("external_direction", wave2 / "EXTERNAL_FULL_SCORE_OS_DIRECTION.tsv"),
        ("external_heterogeneity", wave2 / "EXTERNAL_FULL_SCORE_OS_HETEROGENEITY.tsv"),
        ("external_ineligibility", wave2 / "EXTERNAL_HETEROGENEITY_INELIGIBILITY.tsv"),
    ]
    workbook("STable20_composition_platform_sensitivities", s20)
    export_flat_sources("STable20_composition_platform_sensitivities", s20)

    (OUT / "README.md").write_text(
        "# Unique supplementary table release\n\n"
        "This is the sole numbered table bundle for the current revision. STable7 is the spatial score dictionary "
        "and replaces the prior duplicate candidate-provenance numbering. STable13-16 contain the independently "
        "audited 2026-10-22 hypoxia, uncertainty, missingness, B-state annotation, method-robustness and composition/platform closures. "
        "Every workbook contains a Definitions sheet. XLSX files are "
        "submission attachments; source_data contains exact flat source tables.\n",
        encoding="utf-8",
    )
    manifest_rows = []
    for path in sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv"):
        manifest_rows.append({
            "relative_path": str(path.relative_to(OUT)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    pd.DataFrame(manifest_rows).to_csv(OUT / "MANIFEST.tsv", sep="\t", index=False)
    print(f"SUPPLEMENTARY_RELEASE_PASS files={sum(1 for p in OUT.rglob('*') if p.is_file())}")


if __name__ == "__main__":
    main()
