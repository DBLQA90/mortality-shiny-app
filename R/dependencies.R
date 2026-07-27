# =========================================================
# Runtime dependencies
# =========================================================

required_packages <- c(
  "glue",
  "PHEindicatormethods",
  "tidyverse",
  "shiny",
  "plotly",
  "forecast",
  "ineptr2",
  "strucchange",
  "memoise",
  "cachem",
  "later"
)

get_cran_repos <- function() {
  repos <- getOption("repos")

  if (
    is.null(repos) ||
      length(repos) == 0 ||
      identical(unname(repos[["CRAN"]]), "@CRAN@")
  ) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }

  repos
}

should_install_missing_packages <- function() {
  value <- tolower(Sys.getenv("MORTALITY_INSTALL_MISSING_PACKAGES", "true"))
  !value %in% c("0", "false", "no", "off")
}

install_project_dependencies <- function(
  packages = required_packages,
  install_missing = should_install_missing_packages(),
  load_packages = TRUE
) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    if (!isTRUE(install_missing)) {
      stop(
        "Missing required R package(s): ",
        paste(missing, collapse = ", "),
        ". Run source(\"install_dependencies.R\") first, or set ",
        "MORTALITY_INSTALL_MISSING_PACKAGES=true.",
        call. = FALSE
      )
    }

    install.packages(missing, repos = get_cran_repos())
  }

  still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "Could not install required R package(s): ",
      paste(still_missing, collapse = ", "),
      ". Check your internet connection, R library permissions, and system dependencies.",
      call. = FALSE
    )
  }

  if (isTRUE(load_packages)) {
    invisible(lapply(
      packages,
      function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    ))
  }

  invisible(packages)
}
