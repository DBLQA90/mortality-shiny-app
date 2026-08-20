#!/usr/bin/env Rscript
# Fetch one whole year of a population indicator into data/snapshots/population.
#
#   Rscript tools/fetch_population_year.R indicator=0012918 year=2024 [overwrite=false]
#
# Why this exists alongside build_population_snapshot_chunks.R
# -----------------------------------------------------------
# That builder requests one (year, area) slice at a time and merges the result
# into whatever the year file already holds, which is right for topping up a
# partial archive. Replacing a year outright is a different job: the whole year
# has to come from one indicator, or the file ends up mixing two population
# series that disagree by several per cent.
#
# Omitting the area dimension costs nothing against the live API - the same
# measurement that made tools/fetch_death_year.R viable - so one request per
# year returns every area, sex and age band at once.
#
# The series this fetches
# -----------------------
# 0012918 is INE's NUTS-2024 population by sex and five-year age group,
# 2021-2025. It is not merely a continuation of 0008273 (NUTS-2013, 1991-2023):
# it carries a revised estimate that runs progressively higher, +1.71% for 2021,
# +3.93% for 2022 and +5.31% for 2023 on the national total. Both indicators
# were updated within the last two months and still disagree, so this is two
# published series, not stale data on one side.
#
# The archive therefore uses 0012918 for the whole of its range rather than only
# for the years 0008273 lacks. Splicing 2024 onto a 2023 from the older series
# would put a 5.3% step at the seam and read as a real fall in mortality; taking
# the revised series from 2021 moves the seam to 2020/2021, where the two
# differ by 1.7%. See METHODOLOGY.md.

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

indicator <- get_arg("indicator", "0012918")
years <- get_arg("years", get_arg("year", ""))
out_dir <- get_arg("out", "data/snapshots")
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")
dry_run <- tolower(get_arg("dry_run", "false")) %in% c("true", "1", "yes")

parse_years <- function(value) {
  if (!nzchar(value)) stop("year= or years= is required", call. = FALSE)
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(min(bounds), max(bounds)))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

years <- parse_years(years)

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

# INE returns non-places and roll-up categories alongside the real ones. A
# population file that keeps them double-counts every aggregate built from it.
pseudo_areas <- c("Total", "Ignorado", "Estrangeiro")
pseudo_bands <- c("Total", "Idade ignorada")

# Portugal's resident population has stayed between 9.5 and 12 million across
# every year any of these indicators covers. The guard is deliberately wide: it
# is not checking the estimate, it is checking that no dimension was left
# unpinned and silently summed. Fetching births taught this lesson expensively -
# an indicator publishing a Total beside its own components inflated the count
# sixteen-fold, and the figure still looked like a population.
plausible_range <- c(9.5e6, 12e6)

message("Fetching ", indicator, " for ", paste(years, collapse = ", "))

written <- 0L
for (year in years) {
  path <- file.path(out_dir, "population", paste0("year_", year, ".rds"))
  if (file.exists(path) && !overwrite) {
    message("  ", year, ": already present; pass overwrite=true to replace")
    next
  }

  fetched <- app_env$download_data(
    indicator,
    dims = list(dim1 = as.character(year)),
    has_cause = FALSE
  )

  chunk <- fetched %>%
    dplyr::rename(pop = value) %>%
    dplyr::filter(!area %in% pseudo_areas, !age_band %in% pseudo_bands) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(year = as.integer(year), source_indicator = indicator) %>%
    dplyr::arrange(area, sex, age_band) %>%
    dplyr::select(year, area, sex, age_band, pop, source_indicator)

  national <- sum(chunk$pop[chunk$area == "Portugal" & chunk$sex == "HM"], na.rm = TRUE)

  if (!is.finite(national) || national < plausible_range[[1]] || national > plausible_range[[2]]) {
    stop(
      "Implausible national total for ", year, ": ", format(round(national), big.mark = " "),
      ". Expected between ", format(plausible_range[[1]], big.mark = " "), " and ",
      format(plausible_range[[2]], big.mark = " "),
      " - a dimension is probably being summed instead of pinned.",
      call. = FALSE
    )
  }

  # Both-sexes must equal the sum of its parts, or HM is not what it claims.
  by_sex <- chunk %>%
    dplyr::filter(area == "Portugal") %>%
    dplyr::group_by(sex) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")
  hm <- by_sex$pop[by_sex$sex == "HM"]
  parts <- sum(by_sex$pop[by_sex$sex %in% c("H", "M")])
  if (length(hm) == 1 && abs(hm - parts) > 1) {
    stop("HM (", format(hm, big.mark = " "), ") does not equal H + M (",
         format(parts, big.mark = " "), ") for ", year, ".", call. = FALSE)
  }

  message(sprintf(
    "  %d: %s areas, Portugal = %s%s",
    year, format(dplyr::n_distinct(chunk$area)),
    format(round(national), big.mark = " "),
    if (dry_run) "  (dry run, not written)" else ""
  ))

  if (!dry_run) {
    save_rds_atomic(chunk, path)
    written <- written + 1L
  }
}

message("Done: ", written, " year files written.")
