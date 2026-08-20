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

  run <- function(metric, pooling, year = 2022, area = "Barrancos", ref = "Portugal",
                  vintage = "2024") {
    session$setInputs(
      nuts_vintage = vintage,
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

  # Infant mortality is the sparse case the flag exists for: Barrancos records
  # single-digit births, Portugal tens of thousands.
  cat("\n=== Infant rate, Barrancos 2024 (must be flagged) ===\n")
  i24 <- run("infant", "1", year = 2024)
  if (!is.null(i24)) {
    print(as.data.frame(i24[, c("location","period","value","lower","upper","flag")]))
    print(build_annual_metrics_table(i24, "infant"))
    cat("   footnote:", annual_metrics_footnotes(i24, "infant"), "\n")
  }

  cat("\n=== Infant rate, Portugal 2024 (must NOT be flagged) ===\n")
  ipt <- run("infant", "1", year = 2024, area = "Lisboa")
  if (!is.null(ipt)) {
    print(as.data.frame(ipt[, c("location","value","lower","upper","flag")]))
    cat("   footnotes:", length(annual_metrics_footnotes(ipt, "infant")), "\n")
  }

  cat("\n=== Infant deaths count, Barrancos 2024 + 5-year pooling ===\n")
  c24 <- run("infant_deaths", "1", year = 2024)
  if (!is.null(c24)) print(as.data.frame(c24[, c("location","period","n_years","value","upper","flag")]))
  c5 <- run("infant_deaths", "5", year = 2022)
  # Pooled counts stay window totals, so this must not be divided by 5.
  if (!is.null(c5)) print(as.data.frame(c5[, c("location","period","n_years","value")]))

  cat("\n=== Infant count reaches 1992, where the rate cannot ===\n")
  c92 <- run("infant_deaths", "1", year = 1992)
  if (!is.null(c92)) print(as.data.frame(c92[, c("location","period","value")]))
  cat("   rate 1992 (must be refused):\n")
  invisible(run("infant", "1", year = 1992))

  cat("\n=== AVPP, Portugal 2024 (0-4 band split by under-1 counts) ===\n")
  y24 <- run("ypll", "1", year = 2024, area = "Lisboa")
  if (!is.null(y24)) print(as.data.frame(y24[, c("location","value","source_detail")]))

  # The NUTS vintage regroups the same municipalities. Under NUTS 2013 the
  # municipal sums must reproduce INE's own published regional rows exactly.
  cat("\n=== Alentejo 2021 deaths under each vintage ===\n")
  for (v in c("2024", "2013")) {
    r <- run("deaths", "1", year = 2021, area = "Alentejo", vintage = v)
    if (!is.null(r)) {
      cat(sprintf("   NUTS %s: Alentejo = %.0f   (INE's own 0008206 row: 11742)\n",
                  v, r$value[r$location == "Alentejo"]))
    }
  }

  cat("\n=== AML (2013) vs Grande Lisboa + Península de Setúbal (2024) ===\n")
  aml <- run("deaths", "1", year = 2021, area = "Área Metropolitana de Lisboa", vintage = "2013")
  if (!is.null(aml)) cat(sprintf("   AML 2013           = %.0f\n", aml$value[3]))
  split <- run("deaths", "1", year = 2021,
               area = c("Grande Lisboa", "Península de Setúbal"), vintage = "2024")
  if (!is.null(split)) cat(sprintf("   GL + PS 2024       = %.0f\n", split$value[3]))

  cat("\n=== A 2024-only region selected under NUTS 2013 (must be refused) ===\n")
  invisible(run("deaths", "1", year = 2021, area = "Oeste e Vale do Tejo", vintage = "2013"))

  # The NUTS I level. The three units partition the country, so they must close
  # exactly against INE's own national row.
  cat("\n=== Continente + Açores + Madeira vs Portugal ===\n")
  for (y in c(2021, 2022, 2023)) {
    co <- run("deaths", "1", year = y, area = "Continente")
    ac <- run("deaths", "1", year = y, area = "Região Autónoma dos Açores")
    ma <- run("deaths", "1", year = y, area = "Região Autónoma da Madeira")
    total <- co$value[3] + ac$value[3] + ma$value[3]
    cat(sprintf("   %d  %.0f + %.0f + %.0f = %.0f vs Portugal %.0f  %s\n",
                y, co$value[3], ac$value[3], ma$value[3], total, co$value[1],
                if (isTRUE(all.equal(total, co$value[1]))) "EXACT" else "MISMATCH"))
  }
}, session = MockShinySession$new())
