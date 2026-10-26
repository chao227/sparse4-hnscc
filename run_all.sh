#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$repo_root"
export SPARSE4_REPO_ROOT="$repo_root"
mkdir -p logs results/tables results/figures

mode="${1:---check-only}"
case "$mode" in
  --check-only)
    python3 tests/validate_repository.py --mode release
    ;;
  --score-fixture)
    python3 tests/validate_repository.py --mode release
    python3 tests/test_scoring.py
    ;;
  --full)
    python3 tests/validate_repository.py --mode full-inputs
    Rscript scripts/00_setup_and_download.R
    Rscript scripts/01_scrna_qc.R
    Rscript scripts/02_scrna_integration.R
    Rscript scripts/03_scrna_annotation.R
    Rscript scripts/04_bcell_cd8_states.R
    Rscript scripts/05_bcell_cross_sectional_ordering.R
    Rscript scripts/06_cellchat_donor_descriptive.R --mode formal --contract config/cellchat_input_authority.tsv --output results/cellchat_donor_descriptive
    Rscript scripts/07_tcga_wgcna_yellow.R
    Rscript scripts/08_yellow_module_context.R
    Rscript scripts/09_sparse4_construction.R full config/sparse4_analysis_contract.json
    Rscript scripts/11_sparse4_internal_validation.R
    Rscript scripts/13_survival_models_m1_m4.R
    Rscript scripts/14_calibration_nonbootstrap.R
    Rscript scripts/15_external_cohorts.R evaluate sparse
    python3 scripts/16_spatial_yellow_tls.py
    Rscript scripts/17_lymphhot_ikzf1.R
    Rscript scripts/18_scenic_patient_state.R
    Rscript scripts/20_hypoxia_sensitivity.R
    Rscript scripts/21_missingness_annotation.R "$repo_root" results/missingness_annotation
    Rscript scripts/22_nonbootstrap_uncertainty.R
    Rscript scripts/23_make_figures.R
    python3 scripts/24_make_tables.py
    ;;
  *)
    echo "Usage: bash run_all.sh [--check-only|--score-fixture|--full]" >&2
    exit 2
    ;;
esac
