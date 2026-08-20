#!/usr/bin/env Rscript
# Fetch live births by municipality into chunked RDS files.
#
#   Rscript tools/fetch_births.R [years=ALL] [minutes=60] [overwrite=false]
#
# Live births are the denominator of the infant mortality rate. The population
# indicators cannot serve: neither carries an under-1 age band, so infant deaths
# per under-1 population is not computable from them. Births are a separate INE
# dataset, and are published a year ahead of the population estimates.
#
# ---------------------------------------------------------------------------
# Why per-indicator dimension pinning is required
# ---------------------------------------------------------------------------
# These indicators carry dimensions the app does not model - mother's age,
# gestation duration, filiation - and the shared client flattens what it does
# not recognise. Summing the response naively multiplies the answer by the size
# of every unrequested dimension: for Portugal 2024 the raw sum is 1,354,272
# against a true 84,642, a factor of 16 built from sex (x2), mother's age
# groupings (x4) and gestation duration (x2), each contributing its own "Total"
# category alongside its components.
#
# Every dimension other than period and area is therefore pinned explicitly to
# its total category, per indicator, and the fetched total is checked against a
# plausible range before being written.
#
# ---------------------------------------------------------------------------
# Vintages
# ---------------------------------------------------------------------------
#   0000003  1995-2013  NUTS-2002  dims: period, area, mother age, sex, filiation
#   0008084  2011-2023  NUTS-2013  dims: period, area, sex, mother age, gestation
#   0012434  2021-2025  NUTS-2024  dims: period, area, sex, mother age, gestation
#
# Each is used only for the years the others do not cover, so the vintages never
# compete. They agree where they overlap: 0008084 and 0012434 both report 83,671
# live births for Portugal in 2022, which is the cross-check that the pinning is
# right in both.
#
# Deriving the middle years from the crude birth rate (0008264) times population
# was considered and rejected: validated against true counts in the 2021-2023
# overlap it was 1.8-4.4% out nationally and up to 8.9% out at the 90th
# percentile of municipalities. Real counts made that unnecessary.

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

years_arg <- get_arg("years", "ALL")
out_dir <- get_arg("out", "data/snapshots")
budget_minutes <- suppressWarnings(as.numeric(get_arg("minutes", "60")))
if (!is.finite(budget_minutes) || budget_minutes <= 0) budget_minutes <- 60
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")

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

# `pins` must cover EVERY dimension except period and area, so the response
# collapses to one value per area and can simply be summed.
#
# Do not filter on the client's own `sex` / `age_band` columns instead: it
# assumes dim3 is sex and dim4 is age, which holds for 0012434 but not for
# 0000003, whose dimensions run (mother's age, sex, filiation). There the
# client's `sex` column actually carries mother's age groups, and filtering it
# for "HM" silently selected the wrong cells - Portugal 2000 came out at
# 218,596 against a true figure near 120,000.
birth_sources <- list(
  list(
    indicator = "0000003",
    years = 1995:2013,
    # dim3 mother's age, dim4 sex, dim5 filiation
    pins = list(dim3 = "Total", dim4 = "HM", dim5 = "Total")
  ),
  list(
    indicator = "0008084",
    years = 2014:2020,
    # dim3 sex, dim4 mother's age, dim5 gestation duration
    pins = list(dim3 = "HM", dim4 = "Total", dim5 = "Total")
  ),
  list(
    indicator = "0012434",
    years = 2021:2100,
    # dim3 sex, dim4 mother's age, dim5 gestation duration
    pins = list(dim3 = "HM", dim4 = "Total", dim5 = "Total")
  )
)

pseudo_areas <- c("Total", "Ignorado", "Estrangeiro")

indicator_years <- function(indicator) {
  tryCatch(
    {
      values <- app_env$get_dim_values_cached(indicator)
      years <- values %>%
        dplyr::filter(as.integer(dim_num) == 1) %>%
        dplyr::pull(categ_dsg)
      sort(unique(suppressWarnings(as.integer(as.character(years)))))
    },
    error = function(e) {
      message("  ! metadata unreachable for ", indicator, ": ", conditionMessage(e))
      NULL
    }
  )
}

parse_years <- function(value, default_years) {
  if (identical(toupper(value), "ALL")) return(default_years)
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- suppressWarnings(as.integer(strsplit(value, ":", fixed = TRUE)[[1]]))
    bounds <- bounds[!is.na(bounds)]
    if (length(bounds) == 2) return(seq.int(min(bounds), max(bounds)))
  }
  years <- suppressWarnings(as.integer(unlist(strsplit(value, "[,|]"), use.names = FALSE)))
  years[!is.na(years)]
}

births_path <- function(year) file.path(out_dir, "births", paste0("year_", year, ".rds"))

written <- 0L
skipped <- 0L
failed <- character(0)

for (source in birth_sources) {
  available <- indicator_years(source$indicator)
  if (is.null(available)) {
    failed <- c(failed, source$indicator)
    next
  }

  years <- Reduce(intersect, list(available, source$years, parse_years(years_arg, available)))
  if (length(years) == 0) next

  message("\n== ", source$indicator, ": ", length(years), " year(s) ==")

  for (year in years) {
    if (minutes_left() < 1) {
      message("Time budget exhausted; stopping cleanly.")
      break
    }

    path <- births_path(year)
    if (file.exists(path) && !overwrite) {
      skipped <- skipped + 1L
      next
    }

    dims <- c(list(dim1 = as.character(year)), source$pins)
    fetched <- tryCatch(
      app_env$download_data(source$indicator, dims = dims, has_cause = FALSE),
      error = function(e) {
        message("  ", year, ": FAILED (", conditionMessage(e), ")")
        failed <<- c(failed, paste0(source$indicator, "/", year))
        NULL
      }
    )

    if (is.null(fetched) || nrow(fetched) == 0) next

    chunk <- fetched %>%
      dplyr::filter(!area %in% pseudo_areas) %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(births = sum(value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(year = as.integer(year), source_indicator = source$indicator) %>%
      dplyr::select(year, area, births, source_indicator) %>%
      dplyr::arrange(area)

    # A pinning mistake shows up as an implausible national total rather than as
    # an error, so check it before writing: Portugal has recorded between about
    # 70,000 and 130,000 live births a year across this period.
    national <- chunk$births[chunk$area == "Portugal"]
    if (length(national) == 1 && (national < 50000 || national > 200000)) {
      message("  ", year, ": REJECTED - Portugal total ", format(national),
              " is outside the plausible range; dimension pinning is probably wrong.")
      failed <- c(failed, paste0(source$indicator, "/", year))
      next
    }

    save_rds_atomic(chunk, path)
    written <- written + 1L
    message("  ", year, ": ok (", nrow(chunk), " areas, Portugal ", format(national), ")")
  }
}

message("\nDone: ", written, " written, ", skipped, " already present, ", length(failed), " failed.")
if (length(failed) > 0) message("Failures: ", paste(utils::head(failed, 8), collapse = "; "))
