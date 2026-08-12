#!/usr/bin/env Rscript
# End-to-end smoke test: drives the real Shiny server without a browser.
#
#   Rscript tests/smoke_app.R
#
# run_tests.R covers the calculation modules in isolation. This exercises the
# server reactives themselves - inputs set, event observers fired, outputs read
# - against the committed snapshots, which is where wiring bugs actually live.
# It found two: a pooled window reaching into a year with no population
# estimate, and rate metrics offered for years that have none.
#
# Requires the snapshots to be present; makes no network calls.

Sys.setenv(MORTALITY_INSTALL_MISSING_PACKAGES="false", MORTALITY_DEFAULT_DATA_SOURCE="snapshot")
suppressMessages(library(shiny))
setwd(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), ".."))
app <- shiny::shinyAppDir(".")

testServer(app, {

  run <- function(metric, pooling, year = 2022, area = "Barrancos", ref = "Portugal") {
    session$setInputs(
      annual_year = year, annual_cause = "Todas as causas de morte", annual_sex = "HM",
      annual_area = area, annual_area_label = "", annual_metric = metric,
      annual_pooling = pooling, annual_smr_reference = ref,
      annual_data_source = "snapshot", go_annual_metrics = 1
    )
    out <- tryCatch(annual_metrics_long(), error = function(e) e)
    if (inherits(out, "error")) { cat("   ERROR:", conditionMessage(out), "\n"); return(NULL) }
    out
  }

  cat("=== SMR, Barrancos 2022, single year ===\n")
  r <- run("smr", "1")
  if (!is.null(r)) print(as.data.frame(r[, c("location","period","n_years","value","lower","upper")]))

  cat("\n=== SMR, Barrancos 2022, 5-year pooling ===\n")
  r5 <- run("smr", "5")
  if (!is.null(r5)) print(as.data.frame(r5[, c("location","period","n_years","value","lower","upper")]))

  cat("\n=== Crude, 2024 (no population -> must not invent a rate) ===\n")
  r24 <- run("crude", "1", year = 2024)
  if (!is.null(r24)) print(as.data.frame(r24[, c("location","period","value","lower","upper")]))

  cat("\n=== Deaths, 2024 (count metric -> must work) ===\n")
  d24 <- run("deaths", "1", year = 2024)
  if (!is.null(d24)) print(as.data.frame(d24[, c("location","period","value")]))
}, session = MockShinySession$new())
