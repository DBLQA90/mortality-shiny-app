#!/usr/bin/env Rscript
# Run the app's unit tests:  Rscript run_tests.R
#
# The app is not an R package, so this loads the relevant R/ modules into a
# dedicated environment and runs the tests in tests/testthat against them.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[[1]])) else getwd()
setwd(root)

suppressMessages({
  library(testthat)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(glue)
})

for (pkg in c("PHEindicatormethods", "forecast")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required to run the tests.", call. = FALSE)
  }
}

# Load the modules under test into one environment. Only the network-free,
# calculation modules are needed; INE/Shiny wiring is intentionally excluded.
app_env <- new.env(parent = globalenv())

# The modules resolve committed data files (the NUTS lookups) relative to
# `app_dir`, exactly as they do when the app runs. Setting it here lets the
# tests exercise that resolution instead of reading the files by hand.
assign("app_dir", root, envir = app_env)

for (f in c("R/config.R", "R/regions.R", "R/metrics.R", "R/standardisation.R", "R/infant.R", "R/forecast_helpers.R")) {
  sys.source(file.path(root, f), envir = app_env)
}

reporter <- if (identical(Sys.getenv("TESTTHAT_REPORTER"), "")) "summary" else Sys.getenv("TESTTHAT_REPORTER")
test_dir(
  file.path(root, "tests", "testthat"),
  env = app_env,
  reporter = reporter,
  stop_on_failure = TRUE
)
