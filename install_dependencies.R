get_entry_dir <- function() {
  ofiles <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]

  if (length(ofiles) > 0) {
    return(dirname(normalizePath(utils::tail(ofiles, 1), mustWork = FALSE)))
  }

  getwd()
}

repo_dir <- get_entry_dir()
source(file.path(repo_dir, "R", "dependencies.R"), local = TRUE)

message("Installing mortality app runtime dependencies if needed...")
install_project_dependencies(load_packages = FALSE)
message("Dependency check complete.")
