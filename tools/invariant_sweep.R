#!/usr/bin/env Rscript
# Invariant sweep: run the planning module and the annual-metrics tab across
# areas, years, pooling windows and both population denominators, and flag what
# should never happen.
#
#   Rscript tools/invariant_sweep.R               # both parts
#   Rscript tools/invariant_sweep.R --planning    # planning module only (~3 min)
#   Rscript tools/invariant_sweep.R --server      # annual tab through the real server only
#   Rscript tools/invariant_sweep.R --update      # accept the current flags
#   Rscript tools/invariant_sweep.R --report=FILE # also write the report there
#
# regression_check.R compares a fixed list of validated numbers with a baseline.
# This looks at everything else: it found four bugs that the baseline could not
# see, because each lived in a year, area or option combination no tracked number
# used (census change blank depending on the years asked; Norte waste blank in
# 1999; 1991 failing under the mean population; every pooled 2024 window refused
# once 2025 population arrived).
#
# Planning checks, on Portugal, the Continente, the NUTS II regions and a few
# small areas, 1995 to the last year:
#   - missing values where the indicator has data, in the large areas;
#   - values outside a plausible range;
#   - year-on-year jumps of more than 15% in rates that move slowly;
#   - Portugal's published row against the sum of its municipalities;
#   - the Continente's published row against the aggregate of its regions;
#   - the mean against the end-of-year denominator, which must stay close.
# Server checks, through testServer on the committed snapshots: every metric of
# the annual tab, 1991 to the last year, with and without pooling, a region, a
# municipality and a ULS, both denominators - errors, refusals, missing or
# implausible values.
#
# Some flags are true of the data and stay (the infant rate before 2004 from
# incomplete municipal sums; a region blank because INE suppressed one of its
# municipalities). They are listed, with the reason, in
# data/invariant_sweep_expected.csv. Only flags not in that list fail the run.
#
# Exit codes: 0 no new flag, 1 a new flag, 2 the server raised an error or a
# part could not run.

args <- commandArgs(trailingOnly = TRUE)
update_expected <- "--update" %in% args
run_planning <- !"--server" %in% args
run_server <- !"--planning" %in% args
report_path <- sub("^--report=", "", grep("^--report=", args, value = TRUE))
expected_path <- "data/invariant_sweep_expected.csv"

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
suppressMessages({library(tibble); library(dplyr); library(tidyr); library(purrr); library(glue)})

flags <- list()
flag <- function(part, check, area, year, indicator, detail = "") {
  flags[[length(flags) + 1L]] <<- tibble(part, check, area, year = as.integer(year), indicator, detail = as.character(detail))
}
failed_parts <- character(0)

# ---------------------------------------------------------------------------
# Planning module
# ---------------------------------------------------------------------------
sweep_planning <- function() {
  e <- new.env(parent = globalenv())
  assign("app_dir", getwd(), envir = e)
  for (f in c("R/config.R", "R/helpers.R", "R/regions.R", "R/metrics.R", "R/standardisation.R", "R/avoidable.R",
              "R/infant.R", "R/regional_rows.R", "R/planning_parish.R", "R/planning_indicators.R",
              "R/life_expectancy.R", "R/planning_under75.R", "R/planning_standardised.R", "R/planning_sns.R",
              "R/data_versions.R")) {
    sys.source(f, envir = e)
  }
  lk <- e$get_nuts_lookup("2024")
  avail <- lapply(setNames(e$PLANNING_INDICATORS$id, e$PLANNING_INDICATORS$id), e$planning_indicator_years)
  last_year <- max(unlist(avail))
  years <- 1995:last_year
  big <- c("Portugal", "Continente", "Norte", "Centro", "Oeste e Vale do Tejo", "Grande Lisboa",
           "Península de Setúbal", "Alentejo", "Algarve")
  areas <- c(big, e$PLANNING_PORTUGAL_MUNICIPAL, "Lisboa", "Barrancos", "ULS Matosinhos")
  rate_ids <- c("death_rate", "birth_rate", "dsr_all", "dsr_premature", "smr_all")

  t0 <- Sys.time()
  tab <- e$planning_indicator_table(areas, years, lookup = lk, basis = "mean")
  tab_end <- e$planning_indicator_table(areas, years, ids = rate_ids, lookup = lk, basis = "end")
  message("planning tables: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")

  for (id in unique(tab$indicator)) {
    x <- tab %>% filter(indicator == id, area %in% big, year %in% avail[[id]], is.na(value))
    for (i in seq_len(nrow(x))) flag("planning", "NA where data exist", x$area[i], x$year[i], id)
  }

  ranges <- list(life_expectancy = c(70, 90), life_expectancy_65 = c(14, 25), dsr_all = c(500, 2500),
                 dsr_premature = c(150, 900), death_rate = c(5, 20), birth_rate = c(4, 15), infant_rate = c(0, 15),
                 fertility_index = c(0.9, 2.2), smr_all = c(60, 160))
  for (id in names(ranges)) {
    x <- tab %>% filter(indicator == id, area %in% big, !is.na(value), value < ranges[[id]][1] | value > ranges[[id]][2])
    for (i in seq_len(nrow(x))) flag("planning", "out of range", x$area[i], x$year[i], id, sprintf("%.2f", x$value[i]))
  }

  smooth <- c("dsr_all", "dsr_premature", "death_rate", "life_expectancy", "life_expectancy_65", "pop_total",
              "pct_65_plus", "ageing_index", "birth_rate", "smr_all")
  jumps <- tab %>% filter(indicator %in% smooth, area %in% big) %>% arrange(area, indicator, year) %>%
    group_by(area, indicator) %>% mutate(prev = lag(value), ratio = value / prev) %>% ungroup() %>%
    filter(is.finite(ratio), abs(ratio - 1) > 0.15)
  for (i in seq_len(nrow(jumps))) {
    flag("planning", "jump > 15%", jumps$area[i], jumps$year[i], jumps$indicator[i], sprintf("%.2f -> %.2f", jumps$prev[i], jumps$value[i]))
  }

  cmp <- tab %>% filter(area %in% c("Portugal", e$PLANNING_PORTUGAL_MUNICIPAL), indicator %in% c(smooth, "infant_rate"),
                        indicator != "pop_total") %>%
    select(area, year, indicator, value) %>% pivot_wider(names_from = area, values_from = value)
  names(cmp)[names(cmp) == e$PLANNING_PORTUGAL_MUNICIPAL] <- "municipal"
  bad <- cmp %>% filter(is.finite(Portugal), is.finite(municipal), Portugal != 0) %>%
    mutate(d = municipal / Portugal - 1) %>% filter(abs(d) > 0.015)
  for (i in seq_len(nrow(bad))) {
    flag("planning", "Portugal row vs municipal sum", "Portugal", bad$year[i], bad$indicator[i],
         sprintf("%.3f vs %.3f (%+.1f%%)", bad$Portugal[i], bad$municipal[i], 100 * bad$d[i]))
  }

  # The Continente as an ordinary aggregate of its regions, through the same engine.
  lk_c <- lk
  lk_c$nuts1[lk_c$nuts1 == "Continente"] <- "Continente (soma)"
  attr(lk_c, "vintage") <- "2024"
  alt <- e$planning_indicator_table("Continente (soma)", years, lookup = lk_c,
                                    ids = c("death_rate", "birth_rate", "dsr_all", "dsr_premature", "life_expectancy",
                                            "life_expectancy_65", "infant_rate", "pop_total", "deaths"))
  cc <- tab %>% filter(area == "Continente") %>% select(year, indicator, published = value) %>%
    inner_join(alt %>% select(year, indicator, aggregate = value), by = c("year", "indicator")) %>%
    filter(is.finite(published), is.finite(aggregate), aggregate != 0) %>% mutate(d = published / aggregate - 1) %>%
    filter(abs(d) > 0.01)
  for (i in seq_len(nrow(cc))) {
    flag("planning", "Continente row vs aggregate", "Continente", cc$year[i], cc$indicator[i],
         sprintf("%.3f vs %.3f (%+.1f%%)", cc$published[i], cc$aggregate[i], 100 * cc$d[i]))
  }

  mb <- tab %>% filter(indicator %in% rate_ids, area %in% big) %>% select(area, year, indicator, mean = value) %>%
    inner_join(tab_end %>% select(area, year, indicator, end = value), by = c("area", "year", "indicator")) %>%
    filter(is.finite(mean), is.finite(end), end != 0) %>% mutate(d = mean / end - 1) %>% filter(abs(d) > 0.06)
  for (i in seq_len(nrow(mb))) flag("planning", "mean vs end > 6%", mb$area[i], mb$year[i], mb$indicator[i], sprintf("%+.1f%%", 100 * mb$d[i]))
}

# ---------------------------------------------------------------------------
# Annual-metrics tab, through the real server
# ---------------------------------------------------------------------------
sweep_server <- function() {
  Sys.setenv(MORTALITY_INSTALL_MISSING_PACKAGES = "false", MORTALITY_DEFAULT_DATA_SOURCE = "snapshot")
  suppressMessages(library(shiny))
  app <- shiny::shinyAppDir(".")
  ranges <- list(crude = c(100, 5000), dsr = c(100, 5000), smr = c(20, 500), isr = c(100, 5000),
                 proportional = c(0, 100), infant = c(0, 60), deaths = c(0, 6e5), infant_deaths = c(0, 6000),
                 ypll = c(0, 1e7))
  server_errors <- 0L
  # The app's own globals are not visible inside testServer; the last year is
  # the newest population file, which the app's year list also reaches.
  last <- max(as.integer(sub("^year_(\\d{4})\\.rds$", "\\1", list.files("data/snapshots/population", "^year_\\d{4}\\.rds$"))))
  testServer(app, {
    i <- 0L
    run <- function(metric, pooling, year, area, basis = "mean") {
      i <<- i + 1L
      session$setInputs(
        nuts_vintage = "2024", population_basis = basis,
        annual_year = year, annual_cause = "Todas as causas de morte", annual_sex = "HM",
        annual_area = area, annual_area_label = "", annual_metric = metric,
        annual_pooling = as.character(pooling), annual_smr_reference = "Portugal",
        annual_data_source = "snapshot", go_annual_metrics = i
      )
      key <- paste0(metric, " p", pooling, " ", basis)
      out <- tryCatch(annual_metrics_long(), error = function(e) e)
      if (inherits(out, "error")) {
        refused <- inherits(out, "shiny.silent.error") || inherits(out, "validation")
        if (!refused) server_errors <<- server_errors + 1L
        flag("server", if (refused) "refused" else "ERROR", area, year, key, substr(conditionMessage(out), 1, 160))
        return(invisible())
      }
      r <- ranges[[metric]]
      bad <- out[is.na(out$value) | out$value < r[1] | out$value > r[2], , drop = FALSE]
      for (j in seq_len(nrow(bad))) flag("server", "NA or out of range", bad$location[j], year, key, sprintf("%.2f", bad$value[j]))
    }
    for (metric in names(ranges)) for (year in c(1991, 2001, 2014, 2021, last - 1L, last)) for (pooling in c(1, 3)) {
      for (area in c("Lisboa", "Alentejo", "ULS Matosinhos")) run(metric, pooling, year, area)
    }
    for (metric in c("crude", "dsr", "smr")) for (year in c(2001, 2021, last - 1L)) for (area in c("Lisboa", "Alentejo")) {
      run(metric, 1, year, area, basis = "end")
    }
    message("server runs: ", i)
  })
  server_errors
}

server_errors <- 0L
if (run_planning) {
  ok <- tryCatch({ sweep_planning(); TRUE }, error = function(e) { message("planning sweep failed: ", conditionMessage(e)); FALSE })
  if (!ok) failed_parts <- c(failed_parts, "planning")
}
if (run_server) {
  server_errors <- tryCatch(sweep_server(), error = function(e) { message("server sweep failed: ", conditionMessage(e)); NA })
  if (is.na(server_errors)) failed_parts <- c(failed_parts, "server")
}

found <- bind_rows(flags)
if (nrow(found) == 0) found <- tibble(part = character(), check = character(), area = character(), year = integer(),
                                      indicator = character(), detail = character())
key_cols <- c("part", "check", "area", "year", "indicator")
parts_run <- c(if (run_planning) "planning", if (run_server) "server")

expected <- if (file.exists(expected_path)) {
  utils::read.csv(expected_path, encoding = "UTF-8", stringsAsFactors = FALSE, colClasses = c(year = "integer"))
} else {
  tibble(part = character(), check = character(), area = character(), year = integer(), indicator = character(), reason = character())
}

if (update_expected) {
  # Keep the reasons already written; new rows get an empty reason to be filled in.
  kept <- expected %>% filter(!part %in% parts_run)
  new <- found %>% distinct(across(all_of(key_cols))) %>%
    left_join(expected %>% select(all_of(key_cols), reason), by = key_cols) %>%
    mutate(reason = ifelse(is.na(reason), "", reason))
  out <- bind_rows(kept, new) %>% arrange(part, check, indicator, area, year)
  utils::write.csv(out, expected_path, row.names = FALSE, fileEncoding = "UTF-8")
  message(sprintf("Recorded %d accepted flags in %s (%d without a reason yet).", nrow(out), expected_path, sum(out$reason == "")))
  quit(status = 0)
}

new_flags <- found %>% anti_join(expected, by = key_cols)
gone <- expected %>% filter(part %in% parts_run) %>% anti_join(found, by = key_cols)

lines <- c(
  sprintf("Invariant sweep, %s - parts: %s", format(Sys.time(), "%Y-%m-%d %H:%M"), paste(parts_run, collapse = ", ")),
  sprintf("Flags: %d (%d accepted, %d new). Accepted flags no longer raised: %d. Server errors: %s.",
          nrow(found), nrow(found) - nrow(new_flags), nrow(new_flags), nrow(gone), server_errors),
  if (length(failed_parts)) paste("FAILED to run:", paste(failed_parts, collapse = ", "))
)
fmt <- function(x) utils::capture.output(print(as.data.frame(x), row.names = FALSE, right = FALSE))
if (nrow(new_flags)) lines <- c(lines, "", "NEW FLAGS:", fmt(new_flags))
if (nrow(gone)) lines <- c(lines, "", "Accepted flags no longer raised (run --update to drop them):", fmt(gone[, key_cols]))
writeLines(lines)
if (length(report_path)) writeLines(lines, report_path)

status <- if (length(failed_parts) || isTRUE(server_errors > 0)) 2L else if (nrow(new_flags)) 1L else 0L
quit(status = status)
