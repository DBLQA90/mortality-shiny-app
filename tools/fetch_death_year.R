#!/usr/bin/env Rscript
# Fetch one whole year of a death indicator into chunked RDS files.
#
#   Rscript tools/fetch_death_year.R indicator=0013166 year=2024 [minutes=60]
#
# Why this exists alongside build_death_snapshot_chunks.R
# ------------------------------------------------------
# That builder expands its request into one INE call per (year, area, cause).
# For a full year that is 309 x 66 = 20,394 calls, several hours of waiting.
#
# Measured against the live API, omitting the area dimension costs nothing:
# one call for a single area and a single cause returns 63 rows in ~3.8 s,
# while one call for a single cause across *all* areas returns 22,491 rows in
# ~3.8 s. Dropping the area filter therefore turns a year into 66 calls, a few
# minutes rather than most of a day.
#
# Omitting the year as well makes the response large enough that the request
# times out, so (year, cause) is the coarsest slice that reliably works.
#
# A further benefit: the all-areas response includes INE's own NUTS II and
# NUTS III rows, so regional geographies arrive with the municipalities instead
# of needing a separate backfill.

suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  value <- if (length(hit) == 0) Sys.getenv(toupper(name), unset = "") else sub(paste0("^", name, "="), "", hit[[1]])
  if (!nzchar(value)) default else value
}

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
repo_root <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

indicator <- get_arg("indicator", "0013166")
year <- as.integer(get_arg("year", ""))
out_dir <- get_arg("out", "data/snapshots")
budget_minutes <- suppressWarnings(as.numeric(get_arg("minutes", "60")))
if (!is.finite(budget_minutes) || budget_minutes <= 0) budget_minutes <- 60
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")

if (is.na(year)) stop("year= is required", call. = FALSE)

deadline <- Sys.time() + budget_minutes * 60
minutes_left <- function() as.numeric(difftime(deadline, Sys.time(), units = "mins"))

app_env <- new.env(parent = globalenv())
assign("app_dir", repo_root, envir = app_env)
assign("get_app_dir", function() repo_root, envir = app_env)
for (f in c("R/config.R", "R/helpers.R", "R/cache.R", "R/snapshots.R", "R/ine_client.R")) {
  sys.source(file.path(repo_root, f), envir = app_env)
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Could not move temporary file into ", path, call. = FALSE)
}

# Reuse the app's own tokeniser rather than reimplementing it. A second copy
# has to agree with snapshot_file_token() exactly or the reader will not find
# the files this tool writes, and the two had already drifted in the order of
# tolower() and iconv() - harmless for the current cause list, but only by luck.
cause_token <- function(cause) app_env$snapshot_file_token(cause)

causes <- app_env$get_dim_values_cached(indicator) %>%
  dplyr::filter(as.integer(dim_num) == 5) %>%
  dplyr::pull(categ_dsg) %>%
  as.character() %>%
  unique()

message("Fetching ", indicator, " ", year, ": ", length(causes), " causes, all areas (1 request each)")

# INE returns non-places alongside geographies; they must not enter the archive.
pseudo_areas <- c("Total", "Ignorado", "Estrangeiro")

written <- 0L
skipped <- 0L
failed <- character(0)

for (cause in causes) {
  if (minutes_left() < 1) {
    message("Time budget exhausted; stopping cleanly after ", written, " causes.")
    break
  }

  path <- file.path(out_dir, "deaths", indicator, paste0("year_", year), paste0("cause_", cause_token(cause), ".rds"))
  if (file.exists(path) && !overwrite) {
    skipped <- skipped + 1L
    next
  }

  fetched <- tryCatch(
    app_env$download_data(indicator, dims = list(dim1 = as.character(year), dim5 = cause), has_cause = TRUE),
    error = function(e) {
      message("  ", cause, ": FAILED (", conditionMessage(e), ")")
      failed <<- c(failed, cause)
      NULL
    }
  )

  if (is.null(fetched) || nrow(fetched) == 0) next

  chunk <- fetched %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(age_band = dplyr::case_when(
      age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
      TRUE ~ age_band
    )) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total"), !area %in% pseudo_areas) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::arrange(area, sex, age_band) %>%
    dplyr::select(year, area, sex, cause, age_band, deaths)

  save_rds_atomic(chunk, path)
  written <- written + 1L
  if (written %% 10 == 0) message("  ", written, " causes written ...")
}

message("Done: ", written, " written, ", skipped, " already present, ", length(failed), " failed.")
if (length(failed) > 0) message("Re-run to retry: ", paste(utils::head(failed, 5), collapse = "; "))
