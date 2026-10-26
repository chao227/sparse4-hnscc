
repo_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "config", "params.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Repository root could not be located", call. = FALSE)
    current <- parent
  }
}

read_project_config <- function(root = repo_root()) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("The yaml package is required", call. = FALSE)
  list(
    paths = yaml::read_yaml(file.path(root, "config", "paths.yml")),
    params = yaml::read_yaml(file.path(root, "config", "params.yml"))
  )
}

set_reproducible_seed <- function(seed) {
  RNGkind(kind = "L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
}
