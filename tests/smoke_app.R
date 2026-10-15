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
# Requires the snapshots to be present; makes no network calls. Every check
# prints its result, and the script exits non-zero if any failed, so CI
# (.github/workflows/tests.yml) runs it on every push.

Sys.setenv(MORTALITY_INSTALL_MISSING_PACKAGES="false", MORTALITY_DEFAULT_DATA_SOURCE="snapshot")
suppressMessages(library(shiny))
setwd(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), ".."))
app <- shiny::shinyAppDir(".")

smoke <- new.env()
smoke$failures <- character(0)
fail <- function(what) {
  smoke$failures <- c(smoke$failures, what)
  invisible(TRUE)
}
# TRUE, and recorded as a failure, when `x` is an error.
failed <- function(x, what = deparse(substitute(x))) inherits(x, "error") && fail(what)
verdict <- function(ok, pass = "OK", what = "value check") if (isTRUE(ok)) pass else { fail(what); "MISMATCH" }
# The opposite: the request must be refused.
expect_refusal <- function(x, what) {
  ok <- inherits(x, "error")
  if (!ok) fail(what)
  ok
}

# Deaths of unknown residence: in Portugal's all-ages total, in no municipality.
unknown_residence <- sapply(c("2021", "2022", "2023"), function(y) {
  indicator <- if (as.integer(y) >= 2022) "0013166" else "0008206"
  totals <- readRDS(file.path("data/snapshots/death_totals", indicator, paste0("year_", y, ".rds")))
  totals <- totals[totals$sex == "HM" & totals$cause == "Todas as causas de morte", ]
  totals$deaths[totals$area == "Portugal"] - sum(totals$deaths[totals$area %in% readRDS("data/nuts_lookup_2024.rds")$municipality])
})

testServer(app, {

  run <- function(metric, pooling, year = 2022, area = "Barrancos", ref = "Portugal",
                  vintage = "2024", refuse = FALSE) {
    session$setInputs(
      nuts_vintage = vintage,
      annual_year = year, annual_cause = "Todas as causas de morte", annual_sex = "HM",
      annual_area = area, annual_area_label = "", annual_metric = metric,
      annual_pooling = pooling, annual_smr_reference = ref,
      annual_data_source = "snapshot", go_annual_metrics = 1
    )
    out <- tryCatch(annual_metrics_long(), error = function(e) e)
    what <- paste(metric, pooling, year, paste(area, collapse = "+"), vintage)
    if (inherits(out, "error")) {
      if (refuse) {
        cat("   refused, as it must be:", conditionMessage(out), "\n")
      } else {
        cat("   ERROR:", conditionMessage(out), "\n")
        fail(what)
      }
      return(NULL)
    }
    if (refuse) fail(paste(what, "was not refused"))
    out
  }

  cat("=== SMR, Barrancos 2022, single year ===\n")
  r <- run("smr", "1")
  if (!is.null(r)) print(as.data.frame(r[, c("location","period","n_years","value","lower","upper")]))

  cat("\n=== SMR, Barrancos 2022, 5-year pooling ===\n")
  r5 <- run("smr", "5")
  if (!is.null(r5)) print(as.data.frame(r5[, c("location","period","n_years","value","lower","upper")]))

  cat("\n=== Crude, 2024 ===\n")
  r24 <- run("crude", "1", year = 2024)
  if (!is.null(r24)) print(as.data.frame(r24[, c("location","period","value","lower","upper")]))

  # Population reaches 2025 but deaths by cause stop at 2024: a window centred
  # on 2024 must drop 2025 rather than refuse the whole request.
  cat("\n=== DSR, 2024, 3-year pooling (window must stop at 2024) ===\n")
  p24 <- run("dsr", "3", year = 2024)
  if (!is.null(p24)) print(as.data.frame(p24[, c("location","period","n_years","value")]))

  # The mean population of 1991 needs 1990, which the series does not have.
  cat("\n=== DSR, 1991, mean population (first year keeps its own estimate) ===\n")
  d91 <- run("dsr", "1", year = 1991)
  if (!is.null(d91)) print(as.data.frame(d91[, c("location","period","value")]))

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
  invisible(run("infant", "1", year = 1992, refuse = TRUE))

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
  invisible(run("deaths", "1", year = 2021, area = "Oeste e Vale do Tejo", vintage = "2013", refuse = TRUE))

  # The NUTS I level. The three units partition the country: they must close
  # against INE's national row up to the deaths of unknown residence, which only
  # that row holds (Portugal's total less the sum of its municipalities).
  cat("\n=== Continente + Açores + Madeira vs Portugal ===\n")
  for (y in c(2021, 2022, 2023)) {
    co <- run("deaths", "1", year = y, area = "Continente")
    ac <- run("deaths", "1", year = y, area = "Região Autónoma dos Açores")
    ma <- run("deaths", "1", year = y, area = "Região Autónoma da Madeira")
    total <- co$value[3] + ac$value[3] + ma$value[3]
    unknown <- unknown_residence[[as.character(y)]]
    cat(sprintf("   %d  %.0f + %.0f + %.0f = %.0f vs Portugal %.0f, unknown residence %.0f  %s\n",
                y, co$value[3], ac$value[3], ma$value[3], total, co$value[1], unknown,
                verdict(abs(co$value[1] - total - unknown) < 0.5, "EXACT")))
  }

  # The planning tab, built around one location. Its figures must agree with
  # INE's own rows: the workbook's Continente 2020-2022 total is 356,333, and
  # the fertility index reproduces INE's published national series.
  cat("\n=== Planning indicators ===\n")
  session$setInputs(
    nuts_vintage = "2024", planning_area = "Matosinhos", planning_indicator = "infant_rate",
    planning_years = c(2015, 2025)
  )
  comparators <- tryCatch(planning_available_comparators(), error = function(e) e)
  if (!failed(comparators)) {
    cat("   Matosinhos comparators:", paste(paste0(comparators$level, "=", comparators$area), collapse = "; "), "\n")
  }
  session$setInputs(planning_comparators = comparators$area)
  series <- tryCatch(planning_series(), error = function(e) e)
  if (failed(series)) {
    cat("   ERROR:", conditionMessage(series), "\n")
  } else {
    cat("   infant rate series rows:", nrow(series), "areas:", length(unique(series$area)), "\n")
  }

  session$setInputs(planning_indicator = "pop_total")
  cat("   absolute indicator areas:", nrow(planning_areas()), "(location only)\n")

  session$setInputs(planning_area = "Continente", planning_indicator = "death_rate", planning_years = c(2020, 2022))
  prop <- tryCatch(planning_proportional_view(), error = function(e) e)
  if (!failed(prop)) {
    total <- prop$table$deaths[prop$table$area == "Continente" & prop$table$code == "C00"]
    # 356,355 on the current files: 22 above the workbook (0.006%) - 2022 is
    # read from 0013166, not 0008206.
    cat(sprintf("   Continente 2020-2022 all-cause deaths = %.0f (workbook I45: 356333, diff %+.0f) %s\n",
                total, total - 356333, verdict(abs(total - 356333) <= 100, "OK")))
  }

  # Proportional mortality under 75 (I46): Continente 2020-2022 against the
  # workbook's 98,546 deaths.
  session$setInputs(planning_proportional_ages = "under75", planning_years = c(2020, 2022))
  u75 <- tryCatch(planning_proportional_view(), error = function(e) e)
  if (!failed(u75)) {
    total <- u75$table$deaths[u75$table$area == "Continente" & u75$table$code == "C00"]
    cat(sprintf("   Continente 2020-2022 deaths under 75 = %.0f (workbook I46: 98546) %s\n",
                total, verdict(isTRUE(all.equal(total, 98546)), "EXACT")))
  }
  session$setInputs(planning_proportional_ages = "all")

  session$setInputs(planning_area = "Portugal", planning_indicator = "fertility_index", planning_years = c(2023, 2023))
  isf <- tryCatch(planning_series(), error = function(e) e)
  if (!failed(isf)) {
    value <- isf$value[isf$area == "Portugal" & isf$year == 2023]
    cat(sprintf("   Portugal 2023 fertility index = %.2f (INE 0001293: 1.32) %s\n", value, verdict(round(value, 2) == 1.32, "EXACT")))
  }
  # Life expectancy reproduces Eurostat for Portugal (2019: 82.0); INE's
  # published tables run about 0.9 years lower by method.
  session$setInputs(planning_area = "Portugal", planning_indicator = "life_expectancy", planning_years = c(2019, 2019))
  le <- tryCatch(planning_series(), error = function(e) e)
  if (!failed(le)) {
    value <- le$value[le$area == "Portugal" & le$year == 2019]
    cat(sprintf("   Portugal 2017-2019 life expectancy = %.1f (Eurostat 2019: 82.0) %s\n", value, verdict(abs(value - 82.0) <= 0.2, "OK")))
  }
  session$setInputs(planning_area = "ULS Guarda", planning_indicator = "ageing_index", planning_years = c(2020, 2025))
  ranking <- tryCatch(output$planningRankingPlot, error = function(e) e)
  cat("   ULS ranking chart:", if (failed(ranking)) conditionMessage(ranking) else "OK", "\n")

  # Portugal as the municipal sum replaces the published row everywhere.
  session$setInputs(planning_area = "Matosinhos", planning_indicator = "infant_rate", planning_years = c(2015, 2025),
                    planning_portugal = "municipal")
  session$setInputs(planning_comparators = planning_available_comparators()$area)
  shown <- planning_areas()$area
  cat("   municipal Portugal in comparators:", verdict("Portugal (soma dos municípios)" %in% shown && !"Portugal" %in% shown, "TRUE", "municipal Portugal in comparators"), "\n")
  series <- tryCatch(planning_series(), error = function(e) e)
  cat("   significance marks:", if (failed(series)) conditionMessage(series) else paste(names(table(series$significance)), table(series$significance), sep = "=", collapse = " "), "\n")
  funnel <- tryCatch(output$planningFunnelPlot, error = function(e) e)
  cat("   funnel (ULS):", if (failed(funnel)) conditionMessage(funnel) else "OK", "\n")
  session$setInputs(planning_funnel_units = "Município")
  funnel <- tryCatch(planning_funnel(), error = function(e) e)
  cat("   funnel (municipalities):", if (failed(funnel)) conditionMessage(funnel) else paste(nrow(funnel$data), "points"), "\n")
  session$setInputs(planning_indicator = "ageing_index")
  funnel <- tryCatch(planning_funnel(), error = function(e) e)
  cat("   funnel for an index refuses:", expect_refusal(funnel, "funnel for an index"), "\n")
  session$setInputs(planning_indicator = "pct_education_none", planning_education_age = "25")
  header <- tryCatch(output$planningIndicatorHeader, error = function(e) e)
  cat("   education header:", if (failed(header)) conditionMessage(header) else verdict(grepl("25 e mais anos", as.character(header$html)), "TRUE", "education header"), "\n")
  summary <- tryCatch(output$planningSummaryTable, error = function(e) e)
  cat("   summary table:", if (failed(summary)) conditionMessage(summary) else "OK", "\n")
  session$setInputs(planning_area = "ULS Matosinhos", planning_indicator = "smr_all", planning_years = c(2015, 2025), planning_education_age = "0")
  smr <- tryCatch(planning_series(), error = function(e) e)
  cat("   SMR series:", if (failed(smr)) conditionMessage(smr) else {
    own <- smr[smr$area == "ULS Matosinhos" & is.finite(smr$value), ]
    last <- own[own$year == max(own$year), ]
    paste(sprintf("ULS Matosinhos %s = %.1f (%s)", last$year, last$value, last$significance),
          verdict(nrow(last) == 1, "OK", "SMR series"))
  }, "\n")
  funnel <- tryCatch(planning_funnel(), error = function(e) e)
  cat("   SMR funnel:", if (failed(funnel)) conditionMessage(funnel) else paste(nrow(funnel$data), "ULS"), "\n")
  cause <- tryCatch(output$planningCauseTable, error = function(e) e)
  cat("   cause table:", if (failed(cause)) conditionMessage(cause) else "OK", "\n")
  session$setInputs(planning_cause_sex = "M")
  cause <- tryCatch(planning_cause_view(), error = function(e) e)
  cat("   cause view (women):", if (failed(cause)) conditionMessage(cause) else paste(nrow(cause$table), "rows"), "\n")
  # A ULS that shares a municipality, read whole and by parish weights.
  session$setInputs(planning_area = "ULS São José", planning_indicator = "pop_total", planning_years = c(2020, 2025),
                    planning_split_mode = "whole")
  whole <- tryCatch(planning_series(), error = function(e) e)
  session$setInputs(planning_split_mode = "parish")
  parish <- tryCatch(planning_series(), error = function(e) e)
  if (!failed(whole) && !failed(parish)) {
    last <- function(x) x$value[x$area == "ULS São José" & x$year == max(x$year)]
    cat(sprintf("   ULS São José population: whole %s, parish %s (%.0f%% of the whole)\n",
                format(last(whole), big.mark = " "), format(round(last(parish)), big.mark = " "), last(parish) / last(whole) * 100))
  } else {
    cat("   split ULS:", conditionMessage(if (inherits(whole, "error")) whole else parish), "\n")
  }
  notes <- tryCatch(output$planningIndicatorNotes, error = function(e) e)
  session$setInputs(planning_split_mode = "whole")

  session$setInputs(planning_area = "Matosinhos", planning_sns_indicator = "sns_mammography")
  sns <- tryCatch(planning_sns_view(), error = function(e) e)
  cat("   SNS view:", if (failed(sns)) conditionMessage(sns) else paste(paste(names(sns$levels), collapse = ", "), "| periods", length(unique(sns$table$period))), "\n")
  sns_table <- tryCatch(output$planningSnsTable, error = function(e) e)
  cat("   SNS table:", if (failed(sns_table)) conditionMessage(sns_table) else "OK", "\n")
  sns_rank <- tryCatch(output$planningSnsRanking, error = function(e) e)
  cat("   SNS ranking:", if (failed(sns_rank)) conditionMessage(sns_rank) else "OK", "\n")
  weekly <- tryCatch(planning_weekly_view(), error = function(e) e)
  cat("   weekly view:", if (failed(weekly)) conditionMessage(weekly) else paste(weekly$region$region, "| years", paste(unique(weekly$excess$year), collapse = ",")), "\n")
  session$setInputs(planning_weekly_age = "85plus")
  weekly_table <- tryCatch(output$planningWeeklyTable, error = function(e) e)
  cat("   weekly table (85+):", if (failed(weekly_table)) conditionMessage(weekly_table) else "OK", "\n")
  session$setInputs(planning_area = "Área Metropolitana do Porto")
  sns_nuts <- tryCatch(planning_sns_view(), error = function(e) e)
  cat("   SNS for a NUTS III not made of ULS refuses:", expect_refusal(sns_nuts, "SNS for a NUTS III"), "\n")
  # Empirical Bayes smoothing of the funnel's municipalities.
  session$setInputs(planning_area = "Matosinhos", planning_indicator = "infant_rate", planning_years = c(2015, 2024),
                    planning_funnel_units = "Município", planning_funnel_values = "smoothed")
  smoothed <- tryCatch(output$planningFunnelPlot, error = function(e) e)
  cat("   smoothed funnel:", if (failed(smoothed)) conditionMessage(smoothed) else "OK", "\n")
  session$setInputs(planning_funnel_values = "observed")

  session$setInputs(planning_area = "Matosinhos", planning_indicator = "pct_education_none", planning_education_age = "25")
  profile <- tryCatch(output$downloadPlanningProfile, error = function(e) e)
  cat("   profile download:", if (failed(profile)) conditionMessage(profile) else if (file.exists(profile)) paste(file.size(profile), "bytes") else "missing", "\n")

  # The latest year: all-cause deaths by age rebuilt from weekly deaths, for
  # whole regions only. Norte reaches it; a municipality stops a year before.
  cat("\n=== Latest year in the observed tab ===\n")
  for (area in c("Norte", "Lisboa")) {
    session$setInputs(area = area, area_label = "", cause = "Todas as causas de morte", sex = "HM", population = "Total",
                      rate_type = "dsr", years_import = c(2023, 2025), data_source = "snapshot",
                      go_rates = match(area, c("Norte", "Lisboa")))
    h <- tryCatch(observed_history(), error = function(e) e)
    if (!failed(h)) {
      has <- 2025 %in% h$series$year
      cat(sprintf("   %s: years %s - 2025 %s %s\n", area, paste(h$series$year, collapse = ", "),
                  if (has) "included" else "left out", verdict(has == identical(area, "Norte"), "OK", paste("latest year,", area))))
    }
  }

  # Infant mortality in the guided forecast, five-year periods, both methods.
  cat("\n=== Guided forecast of five-year infant mortality, Norte ===\n")
  for (method in c("annual", "rolling")) {
    session$setInputs(beginner_area = "Norte", beginner_area_label = "", beginner_indicator = "infant",
                      beginner_infant_window = "5", beginner_infant_method = method, beginner_data_source = "snapshot",
                      beginner_cause = "Todas as causas de morte", beginner_sex = "HM", beginner_horizon = 4,
                      beginner_years_fit = c(2004, 2025), beginner_mode = "recommended", beginner_validation = "rolling",
                      beginner_test_pct = 25, go_beginner_forecast = match(method, c("annual", "rolling")))
    f <- tryCatch(beginner_forecast(), error = function(e) e)
    if (!failed(f)) {
      last <- utils::tail(f$fc[f$fc$model == f$recommended_model, ], 1)
      cat(sprintf("   %s: %s, %d = %.2f (%.2f; %.2f) %s\n", method, f$recommended_model, last$year, last$mean, last$lower, last$upper,
                  verdict(is.finite(last$mean) && last$lower < last$mean && last$mean < last$upper, "OK", paste("infant forecast", method))))
    }
  }
}, session = MockShinySession$new())

if (length(smoke$failures) > 0) {
  cat("\nFAILED smoke checks:\n", paste0("  - ", smoke$failures, "\n"), sep = "")
  quit(status = 1)
}
cat("\nAll smoke checks passed.\n")
