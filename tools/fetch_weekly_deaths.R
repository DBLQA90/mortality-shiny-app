#!/usr/bin/env Rscript
# Fetch INE's weekly deaths by region of residence and age group into
# data/snapshots/weekly_deaths.
#
#   Rscript tools/fetch_weekly_deaths.R
#
#   0012100  NUTS 2024, from week 1 of 2021, updated weekly (a few weeks behind)
#   0010112  NUTS 2013, 2018 to 2024 - only for the baseline years before 2021,
#            and only for regions whose borders did not change
#
# Both sexes only. One request per indicator returns every week (about 290,000
# rows), so the whole series is refreshed each run; data_versions.R keeps what
# changed.
#
# Output: data/snapshots/weekly_deaths/<indicator>.rds
#   columns: year, week, code, region, age_group, deaths, source_indicator

suppressMessages(library(dplyr))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "fetch_weekly_deaths.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)
failed <- character(0)
for (indicator in c("0012100", "0010112")) {
  raw <- NULL
  for (attempt in 1:4) {
    raw <- tryCatch(client$get_data(indicator, dim3 = "T"), error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) break
    Sys.sleep(60 * attempt)
  }
  if (is.null(raw) || nrow(raw) == 0) {
    message(indicator, ": FAILED")
    failed <- c(failed, indicator)
    next
  }
  period <- as.character(raw$dim_1)
  tidy <- tibble::tibble(
    year = as.integer(sub(".* de (\\d{4})$", "\\1", period)),
    week = as.integer(sub("^(\\d+)\\..*$", "\\1", period)),
    code = as.character(raw$geocod),
    region = trimws(as.character(raw$geodsg)),
    age_group = trimws(as.character(raw$dim_4_t)),
    # A blank is not a zero; weekly counts by region and age are rarely blank.
    deaths = suppressWarnings(as.numeric(raw$valor)),
    source_indicator = indicator
  ) %>%
    dplyr::filter(!is.na(.data$year), !is.na(.data$week)) %>%
    dplyr::arrange(.data$region, .data$age_group, .data$year, .data$week)
  last <- tidy[tidy$year == max(tidy$year), ]
  message(sprintf("%s: %d rows, %d-W%02d to %d-W%02d, %d regions", indicator, nrow(tidy),
                  min(tidy$year), min(tidy$week[tidy$year == min(tidy$year)]), max(tidy$year), max(last$week), dplyr::n_distinct(tidy$region)))
  save_rds_atomic(tidy, file.path("data/snapshots/weekly_deaths", paste0(indicator, ".rds")))
  Sys.sleep(5)
}
message("Done. Failed: ", length(failed))
if (length(failed) > 0) quit(status = 1)
