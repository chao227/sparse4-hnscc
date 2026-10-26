
write_session_info <- function(script_name, root = getwd()) {
  log_dir <- file.path(root, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- tools::file_path_sans_ext(basename(script_name))
  capture.output(sessionInfo(), file = file.path(log_dir, paste0(stem, "_session_info.txt")))
}
