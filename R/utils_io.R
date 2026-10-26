
require_regular_file <- function(path, label = basename(path)) {
  if (!file.exists(path) || dir.exists(path) || nzchar(Sys.readlink(path))) {
    stop(label, " is missing, not a regular file, or is a link: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

read_tsv_strict <- function(path, required_columns = character()) {
  path <- require_regular_file(path)
  value <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(value))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  value
}

write_tsv_atomic <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  write.table(value, temporary, sep = "	", row.names = FALSE, quote = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop("Atomic write failed: ", path, call. = FALSE)
  invisible(path)
}
