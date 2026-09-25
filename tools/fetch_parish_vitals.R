#!/usr/bin/env Rscript
# Fetch INE's annual births and deaths by parish, for the ULS that share a
# municipality (R/planning_parish.R).
#
#   Rscript tools/fetch_parish_vitals.R
#
#   births  0012450 (NUTS 2024, 2021-), 0008234 (NUTS 2013), 0006239 (older)
#   deaths  0012542 (NUTS 2024, 2021-), 0008235 (NUTS 2013), 0005604 (older)
#
# Counts only - INE publishes no age or cause below the municipality - which is
# exactly what is needed to correct the census-based split of a municipality's
# births and deaths between its ULS, year by year. Both sexes.
#
# Output: data/snapshots/parish_vitals/<kind>.rds
#   columns: year, code (the census six-character code), parish, dico (4),
#   value, source_indicator
#
# These indicators write a parish as the seven-character municipal code plus
# two digits (nine characters); the census writes DICO plus the same two
# (six). The join key is the census form.

suppressMessages(library(dplyr))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())

EDITIONS <- list(
  births = c("0012450", "0008234", "0006239"),
  deaths = c("0012542", "0008235", "0005604")
)

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)
failed <- character(0)
for (kind in names(EDITIONS)) {
  parts <- list()
  for (indicator in EDITIONS[[kind]]) {
    raw <- NULL
    for (attempt in 1:3) {
      raw <- tryCatch(client$get_data(indicator, dim3 = "T"), error = function(e) NULL)
      if (!is.null(raw) && nrow(raw) > 0) break
      Sys.sleep(45 * attempt)
    }
    if (is.null(raw) || nrow(raw) == 0) {
      message("  ", indicator, ": FAILED")
      failed <- c(failed, indicator)
      next
    }
    tidy <- tibble::tibble(
      year = suppressWarnings(as.integer(as.character(raw$dim_1))),
      raw_code = as.character(raw$geocod),
      parish = trimws(as.character(raw$geodsg)),
      value = suppressWarnings(as.numeric(raw$valor)),
      source_indicator = indicator
    ) %>%
      dplyr::filter(!is.na(.data$year), nchar(.data$raw_code) %in% c(6L, 9L)) %>%
      dplyr::mutate(code = ifelse(nchar(.data$raw_code) == 9L,
                                  paste0(substr(.data$raw_code, 4, 7), substr(.data$raw_code, 8, 9)),
                                  .data$raw_code)) %>%
      dplyr::select(-raw_code)
    message(sprintf("  %s %s: %d rows, %d-%d, %d parishes", kind, indicator, nrow(tidy),
                    as.integer(min(tidy$year)), as.integer(max(tidy$year)), dplyr::n_distinct(tidy$code)))
    parts[[indicator]] <- tidy
    Sys.sleep(5)
  }
  if (length(parts) == 0) next
  # Newest edition wins for a year both cover.
  out <- dplyr::bind_rows(parts) %>%
    dplyr::mutate(rank = match(.data$source_indicator, EDITIONS[[kind]])) %>%
    dplyr::group_by(.data$year, .data$code) %>%
    dplyr::slice_min(.data$rank, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dico = substr(.data$code, 1, 4)) %>%
    dplyr::select(year, code, parish, dico, value, source_indicator) %>%
    dplyr::arrange(.data$year, .data$code)
  message(sprintf("%s: %d rows, %d-%d, %d parishes", kind, nrow(out), as.integer(min(out$year)), as.integer(max(out$year)), dplyr::n_distinct(out$code)))
  versioned_save_rds(out, file.path("data/snapshots/parish_vitals", paste0(kind, ".rds")),
                     tool = "fetch_parish_vitals.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))
}
message("Done. Failed: ", length(failed))
if (length(failed) == length(unlist(EDITIONS))) quit(status = 1)
