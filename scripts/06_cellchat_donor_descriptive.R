#!/usr/bin/env Rscript

# ==========================================================
# Script   : 06_cellchat_donor_descriptive.R
# Purpose  : Run donor-separated descriptive CellChat analyses
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: donor-level communication analysis; Fig. 1G
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:agent_work/cellchat_donor_descriptive_nobsb_20261016/scripts/run_cellchat_donor_descriptive.R
# ==========================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(CellChat)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1L]]
}
mode <- arg_value("--mode")
if (!mode %in% c("preflight", "smoke", "formal")) stop("--mode must be preflight, smoke or formal")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Unable to resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
module_rel <- "results/cellchat_donor_descriptive"
contract_rel <- arg_value("--contract", "config/cellchat_input_authority.tsv")
output_rel <- arg_value("--output")
if (is.null(output_rel)) stop("--output is required")
launch_rel <- arg_value("--launch-contract")

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA")
}
read_tsv <- function(path) read.delim(path, sep = "\t", header = TRUE,
                                      stringsAsFactors = FALSE, check.names = FALSE,
                                      quote = "", comment.char = "", na.strings = NULL)
sha256_file <- function(path) {
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path)
  strsplit(out[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
resolve_regular <- function(rel_path, expected_sha) {
  path <- normalizePath(file.path(repo, rel_path), mustWork = TRUE)
  info <- file.info(path)
  if (!isTRUE(info$isdir == FALSE) || Sys.readlink(path) != "") {
    stop("Authority is not a standalone regular file: ", rel_path)
  }
  if (!identical(sha256_file(path), expected_sha)) stop("Authority hash mismatch: ", rel_path)
  path
}
collapse_ids <- function(x) paste(sort(unique(x)), collapse = ";")
parse_ids <- function(x) if (!nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1L]]

contract_path <- normalizePath(file.path(repo, contract_rel), mustWork = TRUE)
contract_tab <- read_tsv(contract_path)
if (!identical(names(contract_tab), c("key", "value")) || anyDuplicated(contract_tab$key)) {
  stop("Invalid input authority contract")
}
contract <- setNames(contract_tab$value, contract_tab$key)
contract_sha <- sha256_file(contract_path)
required_contract <- c(
  "runtime_db_manifest_path", "runtime_db_manifest_sha256", "metadata_contract_path",
  "metadata_contract_sha256", "direction_registry_path", "direction_registry_sha256",
  "pathway_registry_path", "pathway_registry_sha256", "donor_registry_path",
  "donor_registry_sha256", "donor_metadata_role", "celltype_metadata_role",
  "smoke_selection_rule", "algorithm_internal_nboot", "algorithm_internal_seed",
  "algorithm_internal_screen_only", "formal_parallel_workers", "donor_bootstrap_authorized",
  "donor_interval_authorized", "cross_donor_inference_authorized",
  "figure_generation_authorized", "manuscript_or_response_edit_authorized"
)
if (!all(required_contract %in% names(contract))) stop("Input contract keys missing")
for (key in c("donor_bootstrap_authorized", "donor_interval_authorized",
              "cross_donor_inference_authorized", "figure_generation_authorized",
              "manuscript_or_response_edit_authorized")) {
  if (!identical(contract[[key]], "false")) stop("Prohibited authorization enabled: ", key)
}
if (!identical(contract[["algorithm_internal_screen_only"]], "true")) {
  stop("CellChat internal resampling must be classified as algorithmic screening only")
}

runtime_manifest_path <- resolve_regular(contract[["runtime_db_manifest_path"]], contract[["runtime_db_manifest_sha256"]])
metadata_contract_path <- resolve_regular(contract[["metadata_contract_path"]], contract[["metadata_contract_sha256"]])
direction_path <- resolve_regular(contract[["direction_registry_path"]], contract[["direction_registry_sha256"]])
pathway_path <- resolve_regular(contract[["pathway_registry_path"]], contract[["pathway_registry_sha256"]])
donor_registry_path <- resolve_regular(contract[["donor_registry_path"]], contract[["donor_registry_sha256"]])

runtime_tab <- read_tsv(runtime_manifest_path)
if (!identical(names(runtime_tab), c("field", "value")) || anyDuplicated(runtime_tab$field)) stop("Invalid runtime DB manifest")
runtime <- setNames(runtime_tab$value, runtime_tab$field)
runtime_required <- c("source_object", "source_object_sha256", "database_rds", "database_rds_sha256",
                      "runtime_CellChat_version", "runtime_CellChat_git_sha", "database_interactions",
                      "database_pathways", "database_annotations")
if (!all(runtime_required %in% names(runtime))) stop("Runtime DB manifest keys missing")
source_path <- resolve_regular(runtime[["source_object"]], runtime[["source_object_sha256"]])
database_path <- resolve_regular(runtime[["database_rds"]], runtime[["database_rds_sha256"]])
source_sha <- runtime[["source_object_sha256"]]
database_sha <- runtime[["database_rds_sha256"]]

metadata_contract <- read_tsv(metadata_contract_path)
donor_candidates <- metadata_contract$field[metadata_contract$role == contract[["donor_metadata_role"]]]
celltype_candidates <- metadata_contract$field[metadata_contract$role == contract[["celltype_metadata_role"]]]
if (length(donor_candidates) != 1L || length(celltype_candidates) != 1L) stop("Metadata roles are not unique")
donor_field <- donor_candidates[[1L]]
celltype_field <- celltype_candidates[[1L]]

directions <- read_tsv(direction_path)
pathways <- read_tsv(pathway_path)
donor_registry <- read_tsv(donor_registry_path)
direction_fields <- c("direction_id", "direction_order", "sender", "receiver",
                      "external_min_cells_each_group", "filterCommunication_min_cells",
                      "minimum_successful_donors", "n_eligible", "eligible_donor_ids")
pathway_fields <- c("pathway_id", "pathway", "database_interactions")
donor_fields <- c("donor_id", "run_order", "n_cells", "n_cell_types_present", "cell_types_present")
if (!all(direction_fields %in% names(directions)) || !all(pathway_fields %in% names(pathways)) ||
    !all(donor_fields %in% names(donor_registry))) stop("Authority table schema mismatch")
if (anyDuplicated(directions$direction_id) || anyDuplicated(directions$direction_order) ||
    anyDuplicated(pathways$pathway_id) || anyDuplicated(pathways$pathway) ||
    anyDuplicated(donor_registry$donor_id) || anyDuplicated(donor_registry$run_order)) {
  stop("Authority identity/order fields must be unique")
}
directions <- directions[order(as.integer(directions$direction_order)), , drop = FALSE]
pathways <- pathways[order(pathways$pathway_id), , drop = FALSE]
donor_registry <- donor_registry[order(as.integer(donor_registry$run_order)), , drop = FALSE]
if (!identical(as.integer(directions$direction_order), seq_len(nrow(directions))) ||
    !identical(as.integer(donor_registry$run_order), seq_len(nrow(donor_registry)))) {
  stop("Authority order fields are not consecutive")
}

object <- readRDS(source_path)
database <- readRDS(database_path)
if (!inherits(object, "CellChat")) stop("Source object is not CellChat")
if (!identical(rownames(object@meta), colnames(object@data))) stop("Source metadata/data order mismatch")
if (!all(c(donor_field, celltype_field) %in% names(object@meta))) stop("Dynamic metadata fields absent")
meta_donor <- as.character(object@meta[[donor_field]])
meta_celltype <- as.character(object@meta[[celltype_field]])
if (anyNA(meta_donor) || anyNA(meta_celltype) || any(!nzchar(meta_donor)) || any(!nzchar(meta_celltype))) {
  stop("Missing dynamic donor or cell-type labels")
}
donors <- as.character(donor_registry$donor_id)
if (!identical(donors, sort(unique(meta_donor)))) stop("Donor registry does not exactly match source object")
celltypes <- sort(unique(meta_celltype))
if (!all(directions$sender %in% celltypes) || !all(directions$receiver %in% celltypes)) {
  stop("Direction authority contains cell types absent from source object")
}

runtime_description <- packageDescription("CellChat")
runtime_version <- as.character(packageVersion("CellChat"))
runtime_git_sha <- if (is.null(runtime_description$RemoteSha)) "NOT_RECORDED" else runtime_description$RemoteSha
if (!identical(runtime_version, runtime[["runtime_CellChat_version"]]) ||
    !identical(runtime_git_sha, runtime[["runtime_CellChat_git_sha"]])) stop("CellChat runtime drift")
if (!is.list(database) || is.null(database$interaction) ||
    !all(c("pathway_name", "annotation") %in% names(database$interaction))) stop("Retained DB schema mismatch")
if (nrow(database$interaction) != as.integer(runtime[["database_interactions"]]) ||
    length(unique(database$interaction$pathway_name)) != as.integer(runtime[["database_pathways"]])) {
  stop("Retained DB cardinality mismatch")
}
pathway_counts <- table(as.character(database$interaction$pathway_name))
if (!all(pathways$pathway %in% names(pathway_counts)) ||
    !identical(as.integer(pathway_counts[pathways$pathway]), as.integer(pathways$database_interactions))) {
  stop("Pathway authority is not exactly supported by retained DB")
}

cell_count_matrix <- table(factor(meta_donor, levels = donors), factor(meta_celltype, levels = celltypes))
registry_check <- vapply(seq_len(nrow(donor_registry)), function(i) {
  donor <- donors[[i]]
  observed_types <- celltypes[cell_count_matrix[donor, ] > 0L]
  identical(sum(cell_count_matrix[donor, ]), as.integer(donor_registry$n_cells[[i]])) &&
    identical(length(observed_types), as.integer(donor_registry$n_cell_types_present[[i]])) &&
    identical(collapse_ids(observed_types), donor_registry$cell_types_present[[i]])
}, logical(1))
if (!all(registry_check)) stop("Donor registry content mismatch")

eligibility_check <- vapply(seq_len(nrow(directions)), function(i) {
  sender <- as.character(directions$sender[[i]])
  receiver <- as.character(directions$receiver[[i]])
  threshold <- as.integer(directions$external_min_cells_each_group[[i]])
  ids <- donors[cell_count_matrix[, sender] >= threshold & cell_count_matrix[, receiver] >= threshold]
  identical(length(ids), as.integer(directions$n_eligible[[i]])) &&
    identical(collapse_ids(ids), directions$eligible_donor_ids[[i]])
}, logical(1))
if (!all(eligibility_check)) stop("Direction eligibility does not recompute from source object")
filter_values <- unique(as.integer(directions$filterCommunication_min_cells))
if (length(filter_values) != 1L) stop("Network-level filter threshold is not unique")
filter_min_cells <- filter_values[[1L]]

endpoint_rows <- vector("list", nrow(directions) * nrow(pathways))
counter <- 0L
for (i in seq_len(nrow(directions))) {
  for (j in seq_len(nrow(pathways))) {
    counter <- counter + 1L
    endpoint_rows[[counter]] <- data.frame(
      endpoint_order = counter,
      endpoint_id = paste(directions$direction_id[[i]], pathways$pathway_id[[j]], sep = "__"),
      direction_id = directions$direction_id[[i]], direction_order = directions$direction_order[[i]],
      sender = directions$sender[[i]], receiver = directions$receiver[[i]],
      pathway_id = pathways$pathway_id[[j]], pathway = pathways$pathway[[j]],
      minimum_successful_donors = as.integer(directions$minimum_successful_donors[[i]]),
      n_eligible_authority = as.integer(directions$n_eligible[[i]]),
      eligible_donor_ids_authority = directions$eligible_donor_ids[[i]],
      stringsAsFactors = FALSE
    )
  }
}
endpoints <- do.call(rbind, endpoint_rows)

prob_body <- paste(deparse(CellChat::computeCommunProb, width.cutoff = 500), collapse = "\n")
filter_body <- paste(deparse(CellChat::filterCommunication, width.cutoff = 500), collapse = "\n")
internal_nboot <- as.integer(contract[["algorithm_internal_nboot"]])
internal_seed <- as.integer(contract[["algorithm_internal_seed"]])

checks <- list()
add_check <- function(check_id, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = check_id, status = if (isTRUE(pass)) "PASS" else "FAIL",
    observed = paste(observed, collapse = ";"), expected = paste(expected, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
add_check("AUTHORITY_HASHES", TRUE, "all exact", "all exact")
add_check("SOURCE_CLASS_AND_ALIGNMENT", inherits(object, "CellChat") && identical(rownames(object@meta), colnames(object@data)), TRUE, TRUE)
add_check("DYNAMIC_METADATA_ROLES", length(donor_candidates) == 1L && length(celltype_candidates) == 1L,
          c(donor_field, celltype_field), "one donor and one cell-type role")
add_check("DONOR_REGISTRY_RECOMPUTED", all(registry_check), sum(registry_check), nrow(donor_registry))
add_check("DIRECTION_ELIGIBILITY_RECOMPUTED", all(eligibility_check), sum(eligibility_check), nrow(directions))
add_check("PATHWAYS_SUPPORTED_BY_DB", all(pathways$pathway %in% names(pathway_counts)), nrow(pathways), nrow(pathways))
add_check("PATHWAY_DB_COUNTS_EXACT", identical(as.integer(pathway_counts[pathways$pathway]), as.integer(pathways$database_interactions)), TRUE, TRUE)
add_check("RUNTIME_EXACT", identical(runtime_version, runtime[["runtime_CellChat_version"]]) && identical(runtime_git_sha, runtime[["runtime_CellChat_git_sha"]]), paste(runtime_version, runtime_git_sha), paste(runtime[["runtime_CellChat_version"]], runtime[["runtime_CellChat_git_sha"]]))
add_check("INTERNAL_LABEL_PERMUTATION_SEMANTICS", grepl("sample.int(nC, size = nC)", prob_body, fixed = TRUE) && grepl("p = nReject/nboot", prob_body, fixed = TRUE), TRUE, TRUE)
add_check("FILTER_OPERATOR_SEMANTICS", grepl("<= min.cells", filter_body, fixed = TRUE), TRUE, TRUE)
add_check("DONOR_RESAMPLING_PROHIBITED", identical(contract[["donor_bootstrap_authorized"]], "false"), contract[["donor_bootstrap_authorized"]], "false")
add_check("DONOR_INTERVALS_PROHIBITED", identical(contract[["donor_interval_authorized"]], "false"), contract[["donor_interval_authorized"]], "false")
add_check("CROSS_DONOR_INFERENCE_PROHIBITED", identical(contract[["cross_donor_inference_authorized"]], "false"), contract[["cross_donor_inference_authorized"]], "false")
add_check("FIGURES_AND_WRITING_PROHIBITED", identical(contract[["figure_generation_authorized"]], "false") && identical(contract[["manuscript_or_response_edit_authorized"]], "false"), TRUE, TRUE)
check_table <- do.call(rbind, checks)
if (any(check_table$status == "FAIL")) stop("Static/preflight checks failed")

out_root <- file.path(repo, output_rel)
if (dir.exists(out_root) && length(list.files(out_root, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Output root must be absent or empty")
}
dir.create(file.path(out_root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "reports"), recursive = TRUE, showWarnings = FALSE)
write_tsv(check_table, file.path(out_root, "reports", "PREFLIGHT_CHECKS.tsv"))
authority_paths <- c(contract_path, runtime_manifest_path, metadata_contract_path, direction_path,
                     pathway_path, donor_registry_path, source_path, database_path)
authority_snapshot <- data.frame(
  relative_path = sub(paste0("^", repo, "/"), "", authority_paths),
  bytes = file.info(authority_paths)$size,
  sha256 = vapply(authority_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_tsv(authority_snapshot, file.path(out_root, "reports", "INPUT_AUTHORITY_SNAPSHOT.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "reports", "SESSIONINFO.txt"))

if (identical(mode, "preflight")) {
  receipt <- list(
    schema = "CELLCHAT_DONOR_DESCRIPTIVE_PREFLIGHT_V1",
    status = "PASS_PREFLIGHT_NO_COMMUNICATION_PROBABILITY_COMPUTED",
    checks_passed = sum(check_table$status == "PASS"), checks_total = nrow(check_table),
    donor_field_read_from_role_contract = donor_field,
    celltype_field_read_from_role_contract = celltype_field,
    donors_read_from_registry = length(donors),
    directions_read_from_registry = nrow(directions),
    pathways_read_from_registry = nrow(pathways),
    endpoints_derived = nrow(endpoints),
    algorithm_internal_label_permutations = internal_nboot,
    algorithm_internal_screen_only = TRUE,
    communication_probability_computed = FALSE,
    donor_bootstrap_run = FALSE, donor_interval_computed = FALSE,
    cross_donor_inferential_contrast_run = FALSE, figure_generated = FALSE
  )
  write_json(receipt, file.path(out_root, "reports", "PREFLIGHT_RECEIPT.json"), auto_unbox = TRUE, pretty = TRUE)
  quit(save = "no", status = 0L)
}

if (identical(mode, "formal")) {
  if (is.null(launch_rel)) stop("Formal mode requires --launch-contract")
  launch_path <- normalizePath(file.path(repo, launch_rel), mustWork = TRUE)
  launch_tab <- read_tsv(launch_path)
  if (!identical(names(launch_tab), c("key", "value")) || anyDuplicated(launch_tab$key)) stop("Invalid launch contract")
  launch <- setNames(launch_tab$value, launch_tab$key)
  if (!all(c("smoke_receipt_path", "smoke_receipt_sha256", "formal_authorized") %in% names(launch)) ||
      !identical(launch[["formal_authorized"]], "true")) stop("Formal launch is not authorized")
  smoke_receipt_path <- resolve_regular(launch[["smoke_receipt_path"]], launch[["smoke_receipt_sha256"]])
  smoke_receipt <- fromJSON(smoke_receipt_path, simplifyVector = TRUE)
  if (!identical(smoke_receipt$status, "PASS_ISOLATED_DONOR_SMOKE") ||
      !identical(smoke_receipt$authority_contract_sha256, contract_sha) ||
      !identical(smoke_receipt$donor_bootstrap_run, FALSE) ||
      !identical(smoke_receipt$cross_donor_inferential_contrast_run, FALSE)) {
    stop("Smoke receipt is not valid for formal launch")
  }
}

eligible_counts <- setNames(integer(length(donors)), donors)
for (ids in directions$eligible_donor_ids) eligible_counts[parse_ids(ids)] <- eligible_counts[parse_ids(ids)] + 1L
if (!identical(contract[["smoke_selection_rule"]], "MAX_ELIGIBLE_DIRECTION_COVERAGE_THEN_RUN_ORDER")) stop("Unknown smoke selection rule")
smoke_candidates <- donors[eligible_counts == max(eligible_counts)]
smoke_donor <- smoke_candidates[which.min(match(smoke_candidates, donors))]
run_donors <- if (identical(mode, "smoke")) smoke_donor else donors

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
dir.create(file.path(out_root, "donor_networks"), showWarnings = FALSE)
dir.create(file.path(out_root, "donor_logs"), showWarnings = FALSE)

run_one_donor <- function(donor) {
  if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")
  cells <- rownames(object@meta)[meta_donor == donor]
  donor_log <- file.path(out_root, "donor_logs", paste0(donor, ".log"))
  donor_network <- file.path(out_root, "donor_networks", paste0(donor, "_netP.rds"))
  stage <- "DONOR_INPUT_PREPARATION"
  outcome <- tryCatch({
    donor_meta <- object@meta[cells, , drop = FALSE]
    donor_meta[[celltype_field]] <- droplevels(factor(as.character(donor_meta[[celltype_field]]),
                                                      levels = sort(unique(as.character(donor_meta[[celltype_field]])))))
    donor_data <- object@data[, cells, drop = FALSE]
    if (!identical(rownames(donor_meta), colnames(donor_data))) stop("donor data/metadata alignment failure")
    stage <- "CREATE_CELLCHAT"
    donor_object <- CellChat::createCellChat(object = donor_data, meta = donor_meta,
                                             group.by = celltype_field, datatype = "RNA")
    donor_object@DB <- database
    stage <- "SUBSET_DATA"
    donor_object <- CellChat::subsetData(donor_object)
    stage <- "IDENTIFY_OVEREXPRESSED_GENES"
    donor_object <- CellChat::identifyOverExpressedGenes(
      donor_object, only.pos = TRUE, thresh.pc = 0, thresh.fc = 0,
      thresh.p = 0.05, do.DE = TRUE, do.fast = TRUE, min.cells = 10
    )
    stage <- "IDENTIFY_OVEREXPRESSED_INTERACTIONS"
    donor_object <- CellChat::identifyOverExpressedInteractions(donor_object, variable.both = TRUE)
    stage <- "COMPUTE_COMMUNICATION_PROBABILITY"
    donor_object <- CellChat::computeCommunProb(
      donor_object, type = "triMean", trim = 0.1, raw.use = TRUE,
      population.size = FALSE, distance.use = NULL, contact.dependent = FALSE,
      nboot = internal_nboot, seed.use = internal_seed, Kh = 0.5, n = 1
    )
    stage <- "FILTER_COMMUNICATION"
    donor_object <- CellChat::filterCommunication(donor_object, min.cells = filter_min_cells)
    stage <- "COMPUTE_PATHWAY_PROBABILITY"
    donor_object <- CellChat::computeCommunProbPathway(donor_object, thresh = 0.05)
    stage <- "AGGREGATE_NETWORK"
    donor_object <- CellChat::aggregateNet(donor_object, thresh = 0.05)
    parameters <- donor_object@options$parameter
    parameter_ok <- identical(as.character(parameters$type.mean), "triMean") &&
      identical(parameters$raw.use, TRUE) && identical(parameters$population.size, FALSE) &&
      identical(parameters$seed.use, internal_seed) && is.null(parameters$distance.use) &&
      identical(parameters$contact.dependent, FALSE) &&
      identical(as.integer(parameters$nboot), internal_nboot) &&
      isTRUE(all.equal(as.numeric(parameters$trim), 0.1)) &&
      isTRUE(all.equal(as.numeric(parameters$Kh), 0.5)) &&
      isTRUE(all.equal(as.numeric(parameters$n), 1))
    if (!parameter_ok) stop(structure(list(message = "stored CellChat parameters violate active contract", call = NULL),
                                      class = c("systemic_contract_error", "error", "condition")))
    probability <- donor_object@netP$prob
    if (!is.array(probability) || length(dim(probability)) != 3L || is.null(dimnames(probability)) ||
        any(!is.finite(probability))) stop("netP probability output invalid")
    payload <- list(
      donor_id = donor, probability = probability, parameters = parameters,
      cell_counts = table(donor_meta[[celltype_field]]),
      source_object_sha256 = source_sha, retained_db_sha256 = database_sha,
      runtime_CellChat_version = runtime_version, runtime_CellChat_git_sha = runtime_git_sha,
      authority_contract_sha256 = contract_sha,
      internal_resampling_role = "ALGORITHM_LEVEL_LABEL_PERMUTATION_SCREEN_ONLY"
    )
    saveRDS(payload, donor_network, compress = "xz", version = 3)
    roundtrip <- readRDS(donor_network)
    if (!identical(roundtrip$donor_id, donor) || !identical(dim(roundtrip$probability), dim(probability)) ||
        !identical(dimnames(roundtrip$probability), dimnames(probability))) {
      stop(structure(list(message = "serialized donor payload roundtrip mismatch", call = NULL),
                     class = c("systemic_contract_error", "error", "condition")))
    }
    list(run_status = "SUCCESS", stage = "COMPLETE", probability = probability,
         message = NA_character_, systemic = FALSE)
  }, systemic_contract_error = function(e) {
    list(run_status = "SYSTEMIC_CONTRACT_FAILURE", stage = stage, probability = NULL,
         message = conditionMessage(e), systemic = TRUE)
  }, error = function(e) {
    list(run_status = "ISOLATED_DONOR_RUN_FAILURE", stage = stage, probability = NULL,
         message = conditionMessage(e), systemic = FALSE)
  })
  if (!identical(outcome$run_status, "SUCCESS") && file.exists(donor_network)) unlink(donor_network)
  writeLines(c(
    paste0("donor_id=", donor), paste0("run_status=", outcome$run_status),
    paste0("stage=", outcome$stage), paste0("message=", ifelse(is.na(outcome$message), "", outcome$message)),
    paste0("source_object_sha256=", source_sha),
    paste0("retained_db_sha256=", database_sha),
    paste0("internal_resampling_role=ALGORITHM_LEVEL_LABEL_PERMUTATION_SCREEN_ONLY")
  ), donor_log, useBytes = TRUE)
  outcome$donor_id <- donor
  outcome$n_cells <- length(cells)
  outcome$n_present_cell_types <- length(unique(meta_celltype[meta_donor == donor]))
  outcome$network_path <- donor_network
  outcome$log_path <- donor_log
  outcome
}

workers <- if (identical(mode, "smoke")) 1L else min(as.integer(contract[["formal_parallel_workers"]]),
                                                       length(run_donors), parallel::detectCores())
if (workers > 1L) {
  outcomes <- parallel::mclapply(run_donors, run_one_donor, mc.cores = workers,
                                 mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  outcomes <- lapply(run_donors, run_one_donor)
}
if (any(vapply(outcomes, inherits, logical(1), what = "try-error"))) stop("Worker process returned try-error")
if (any(vapply(outcomes, function(x) isTRUE(x$systemic), logical(1)))) stop("Systemic donor execution failure")

run_ledger <- do.call(rbind, lapply(seq_along(outcomes), function(i) {
  x <- outcomes[[i]]
  data.frame(
    donor_id = x$donor_id, attempt_order = match(x$donor_id, donors), n_cells = x$n_cells,
    n_present_cell_types = x$n_present_cell_types, attempted = TRUE, run_status = x$run_status,
    failure_reason_ref = if (identical(x$run_status, "SUCCESS")) NA_character_ else paste0("FAIL__", x$donor_id),
    network_rds_sha256 = if (identical(x$run_status, "SUCCESS")) sha256_file(x$network_path) else NA_character_,
    log_sha256 = sha256_file(x$log_path), source_object_sha256 = source_sha,
    retained_db_sha256 = database_sha, runtime_CellChat_version = runtime_version,
    runtime_CellChat_git_sha = runtime_git_sha, authority_contract_sha256 = contract_sha,
    stringsAsFactors = FALSE
  )
}))
write_tsv(run_ledger, file.path(out_root, "results", "CELLCHAT_DONOR_RUN_LEDGER.tsv"))

failed_outcomes <- outcomes[vapply(outcomes, function(x) !identical(x$run_status, "SUCCESS"), logical(1))]
failure_ledger <- if (length(failed_outcomes)) do.call(rbind, lapply(failed_outcomes, function(x) {
  affected <- directions$direction_id[vapply(directions$eligible_donor_ids,
                                              function(ids) x$donor_id %in% parse_ids(ids), logical(1))]
  data.frame(
    failure_id = paste0("FAIL__", x$donor_id), failure_scope = "ISOLATED_DONOR",
    donor_id = x$donor_id, pipeline_stage = x$stage,
    message_sha256 = digest::digest(x$message, algo = "sha256", serialize = FALSE),
    affected_direction_ids = collapse_ids(affected), propagation_action = "PROPAGATE_NA_AND_CONTINUE",
    stringsAsFactors = FALSE
  )
})) else data.frame(
  failure_id = character(), failure_scope = character(), donor_id = character(),
  pipeline_stage = character(), message_sha256 = character(), affected_direction_ids = character(),
  propagation_action = character(), stringsAsFactors = FALSE
)
write_tsv(failure_ledger, file.path(out_root, "results", "CELLCHAT_DONOR_FAILURE_LEDGER.tsv"))

outcome_by_donor <- setNames(outcomes, vapply(outcomes, `[[`, character(1), "donor_id"))
long_rows <- vector("list", nrow(endpoints) * length(run_donors)); row_counter <- 0L
for (e in seq_len(nrow(endpoints))) {
  endpoint <- endpoints[e, ]
  eligible_ids <- parse_ids(endpoint$eligible_donor_ids_authority)
  for (donor in run_donors) {
    row_counter <- row_counter + 1L
    eligible <- donor %in% eligible_ids
    outcome <- outcome_by_donor[[donor]]
    success <- identical(outcome$run_status, "SUCCESS")
    probability_value <- NA_real_; detected <- NA; value_status <- "INELIGIBLE_NA"
    dims_aligned <- FALSE
    if (success) {
      dims <- dimnames(outcome$probability)
      dims_aligned <- endpoint$sender %in% dims[[1L]] && endpoint$receiver %in% dims[[2L]]
    }
    if (eligible && !success) {
      value_status <- "TECHNICAL_FAILURE_NA"
    } else if (eligible && !dims_aligned) {
      stop("Successful eligible donor lacks requested sender/receiver dimensions")
    } else if (eligible && !(endpoint$pathway %in% dims[[3L]])) {
      probability_value <- 0; detected <- FALSE; value_status <- "ALGORITHMIC_NONDETECTION_ZERO"
    } else if (eligible) {
      probability_value <- as.numeric(outcome$probability[endpoint$sender, endpoint$receiver, endpoint$pathway])
      if (!is.finite(probability_value) || probability_value < 0) stop("Invalid endpoint probability")
      detected <- probability_value > 0
      value_status <- if (detected) "DETECTED_VALUE" else "ALGORITHMIC_NONDETECTION_ZERO"
    }
    row <- match(donor, rownames(cell_count_matrix))
    long_rows[[row_counter]] <- data.frame(
      donor_id = donor, endpoint_order = endpoint$endpoint_order, endpoint_id = endpoint$endpoint_id,
      direction_id = endpoint$direction_id, sender = endpoint$sender, receiver = endpoint$receiver,
      pathway_id = endpoint$pathway_id, pathway = endpoint$pathway,
      sender_cells = as.integer(cell_count_matrix[row, endpoint$sender]),
      receiver_cells = as.integer(cell_count_matrix[row, endpoint$receiver]),
      eligible = eligible, endpoint_attempted = eligible, donor_run_status = outcome$run_status,
      failure_reason_ref = if (success) NA_character_ else paste0("FAIL__", donor),
      dimensions_aligned = dims_aligned, probability = probability_value, detected = detected,
      value_status = value_status, stringsAsFactors = FALSE
    )
  }
}
donor_long <- do.call(rbind, long_rows)
write_tsv(donor_long, file.path(out_root, "results", "CELLCHAT_DONOR_ENDPOINT_LONG.tsv"))

if (identical(mode, "formal")) {
  summary_rows <- vector("list", nrow(endpoints))
  for (e in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[e, ]
    x <- donor_long[donor_long$endpoint_id == endpoint$endpoint_id, , drop = FALSE]
    eligible_ids <- x$donor_id[x$eligible]
    ineligible_ids <- x$donor_id[!x$eligible]
    attempted_ids <- x$donor_id[x$eligible & x$endpoint_attempted]
    successful_ids <- x$donor_id[x$eligible & x$donor_run_status == "SUCCESS" & !is.na(x$probability)]
    failed_ids <- setdiff(attempted_ids, successful_ids)
    detected_ids <- x$donor_id[x$eligible & x$donor_run_status == "SUCCESS" & x$value_status == "DETECTED_VALUE"]
    nondetected_ids <- x$donor_id[x$eligible & x$donor_run_status == "SUCCESS" & x$value_status == "ALGORITHMIC_NONDETECTION_ZERO"]
    values <- x$probability[match(successful_ids, x$donor_id)]
    evaluable <- length(successful_ids) >= endpoint$minimum_successful_donors
    stats <- rep(NA_real_, 7L)
    if (evaluable) {
      q1 <- unname(quantile(values, 0.25, type = 7)); q3 <- unname(quantile(values, 0.75, type = 7))
      min_value <- min(values); max_value <- max(values)
      stats <- c(median(values), q1, q3, q3 - q1, min_value, max_value, max_value - min_value)
    }
    summary_rows[[e]] <- data.frame(
      endpoint_order = endpoint$endpoint_order, endpoint_id = endpoint$endpoint_id,
      direction_id = endpoint$direction_id, sender = endpoint$sender, receiver = endpoint$receiver,
      pathway_id = endpoint$pathway_id, pathway = endpoint$pathway,
      endpoint_status = if (evaluable) "EVALUABLE_DESCRIPTIVE" else "NOT_EVALUABLE_INSUFFICIENT_SUCCESSFUL_DONORS",
      minimum_successful_donors = endpoint$minimum_successful_donors,
      n_total_donors = length(donors), total_donor_ids = collapse_ids(donors),
      n_eligible = length(eligible_ids), eligible_donor_ids = collapse_ids(eligible_ids),
      n_ineligible = length(ineligible_ids), ineligible_donor_ids = collapse_ids(ineligible_ids),
      n_attempted_eligible = length(attempted_ids), attempted_eligible_donor_ids = collapse_ids(attempted_ids),
      n_successful_eligible = length(successful_ids), successful_eligible_donor_ids = collapse_ids(successful_ids),
      n_failed_eligible = length(failed_ids), failed_eligible_donor_ids = collapse_ids(failed_ids),
      n_detected = length(detected_ids), detected_donor_ids = collapse_ids(detected_ids),
      n_nondetected_zero = length(nondetected_ids), nondetected_zero_donor_ids = collapse_ids(nondetected_ids),
      detection_fraction = if (evaluable) length(detected_ids) / length(successful_ids) else NA_real_,
      detection_fraction_numerator = length(detected_ids),
      detection_fraction_denominator = length(successful_ids),
      median_probability = stats[[1L]], q1_probability = stats[[2L]], q3_probability = stats[[3L]],
      iqr_probability = stats[[4L]], min_probability = stats[[5L]], max_probability = stats[[6L]],
      range_probability = stats[[7L]],
      source_object_sha256 = source_sha, retained_db_sha256 = database_sha,
      runtime_CellChat_version = runtime_version, runtime_CellChat_git_sha = runtime_git_sha,
      authority_contract_sha256 = contract_sha,
      claim_class = "DONOR_LEVEL_DESCRIPTIVE_CONDITIONAL_PROPENSITY",
      stringsAsFactors = FALSE
    )
  }
  summary_table <- do.call(rbind, summary_rows)
  forbidden_summary_fields <- grepl("bootstrap|interval|(^|_)ci(_|$)|p.?value", names(summary_table), ignore.case = TRUE)
  if (any(forbidden_summary_fields)) stop("Forbidden donor uncertainty/inference field present in endpoint summary")
  if (any(summary_table$n_total_donors != summary_table$n_eligible + summary_table$n_ineligible) ||
      any(summary_table$n_eligible != summary_table$n_successful_eligible + summary_table$n_failed_eligible) ||
      any(summary_table$n_successful_eligible != summary_table$n_detected + summary_table$n_nondetected_zero)) {
    stop("Endpoint denominator partition failure")
  }
  write_tsv(summary_table, file.path(out_root, "results", "CELLCHAT_ENDPOINT_DESCRIPTIVE_SUMMARY.tsv"))
}

result_files <- list.files(file.path(out_root, "results"), full.names = TRUE)
network_files <- list.files(file.path(out_root, "donor_networks"), full.names = TRUE)
receipt <- list(
  schema = "CELLCHAT_DONOR_DESCRIPTIVE_NOBSB_EVIDENCE_BOUNDARY_V1",
  status = if (identical(mode, "smoke")) "PASS_ISOLATED_DONOR_SMOKE" else "COMPUTED_PENDING_DIFFERENT_ROLE_RECOMPUTATION",
  mode = mode, authority_contract_sha256 = contract_sha,
  authority_snapshot_sha256 = sha256_file(file.path(out_root, "reports", "INPUT_AUTHORITY_SNAPSHOT.tsv")),
  donor_field_read_from_role_contract = donor_field, celltype_field_read_from_role_contract = celltype_field,
  donors_in_authority = length(donors), donors_attempted_this_run = length(run_donors),
  donors_successful = sum(run_ledger$run_status == "SUCCESS"),
  donors_failed = sum(run_ledger$run_status != "SUCCESS"),
  smoke_selection_rule = if (identical(mode, "smoke")) contract[["smoke_selection_rule"]] else NULL,
  smoke_donor_selected_from_authority = if (identical(mode, "smoke")) smoke_donor else NULL,
  directions_read_from_authority = nrow(directions), pathways_read_from_authority = nrow(pathways),
  endpoints_derived = nrow(endpoints), parallel_workers = workers,
  algorithm_internal_label_permutations = internal_nboot,
  algorithm_internal_seed = internal_seed,
  algorithm_internal_role = "ALGORITHM_LEVEL_LABEL_PERMUTATION_SCREEN_ONLY_NOT_DONOR_UNCERTAINTY",
  donor_bootstrap_run = FALSE, donor_interval_computed = FALSE,
  cross_donor_inferential_contrast_run = FALSE, multiplicity = "NONE_DESCRIPTIVE",
  cells_used_as_independent_biological_replicates = FALSE, independent_biological_unit = "donor",
  figure_generated = FALSE, manuscript_or_response_edited = FALSE,
  result_hashes = setNames(lapply(result_files, sha256_file), sub(paste0("^", out_root, "/"), "", result_files)),
  donor_network_hashes = setNames(lapply(network_files, sha256_file), basename(network_files)),
  claim_ceiling = "Donor-level descriptive CellChat pathway propensity conditional on captured cells, frozen labels, enriched fractions, retained DB, runtime and algorithmic screen; no cross-donor inference, validation, mechanism, causality or natural-tissue abundance claim."
)
write_json(receipt, file.path(out_root, "reports", "EVIDENCE_BOUNDARY_RECEIPT.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(paste0("CELLCHAT_", toupper(mode), "_COMPLETE donors_successful=", receipt$donors_successful,
           " donors_failed=", receipt$donors_failed, " workers=", workers, "\n"))
