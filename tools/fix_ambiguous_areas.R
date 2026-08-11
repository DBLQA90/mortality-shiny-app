#!/usr/bin/env Rscript
# Repair snapshot rows whose geography was resolved by an ambiguous INE label.
#
#   Rscript tools/fix_ambiguous_areas.R [dataset=all] [years=ALL] [minutes=120] [dry_run=false]
#
# ---------------------------------------------------------------------------
# The defect
# ---------------------------------------------------------------------------
# The snapshot builders resolve a geography through get_cat_code(), which
# matches on the INE category *label*. Several INE labels are not unique, and
# when a label matches more than one category the download returns all of them
# and they are summed into a single row. Three cases affect this archive:
#
#   0003182  "Lisboa"   code 17      NUTS-2002 region Lisboa
#                       code 1711106 municipio de Lisboa
#            -> the stored 1991-2013 population is region + municipio. Deaths
#               are always the municipio, so the crude rate for Lisboa reads
#               ~211/100k in 2013 against ~1231/100k in 2015: understated about
#               six-fold for 23 of the 32 years.
#
#   0008273  "Calheta"  codes 2004501 (R.A.A.) and 3003101 (R.A.M.)
#   0008206  "Lagoa"    codes 1500806 (Algarve) and 2004201 (R.A.A.)
#            -> two genuinely different municipalities, hundreds of kilometres
#               apart, added together. Population for "Lagoa" jumps 22,750 ->
#               37,027 across the 2013/2014 source seam for exactly this reason.
#
# 0008273 and 0013166 label the Lisbon region "Area Metropolitana de Lisboa",
# so "Lisboa" is unambiguous there and the 2014+ population is already correct.
#
# ---------------------------------------------------------------------------
# The repair
# ---------------------------------------------------------------------------
# Request each affected geography by its category *code*, which is unique, then
# write it back under an unambiguous label. Rows carrying the old ambiguous
# label are removed in the same write, so a chunk is never left holding both.
#
# Writes are atomic and per-chunk, so an interrupted run leaves a consistent
# archive and the next run redoes only what is outstanding.

suppressMessages({
  library(dplyr)
  library(tibble)
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

dataset_filter <- tolower(get_arg("dataset", "all"))
years_arg <- get_arg("years", "ALL")
budget_minutes <- suppressWarnings(as.numeric(get_arg("minutes", "120")))
if (!is.finite(budget_minutes) || budget_minutes <= 0) budget_minutes <- 120
dry_run <- tolower(get_arg("dry_run", "false")) %in% c("true", "1", "yes")
out_dir <- get_arg("out", "data/snapshots")

deadline <- Sys.time() + budget_minutes * 60
minutes_left <- function() as.numeric(difftime(deadline, Sys.time(), units = "mins"))

# ---------------------------------------------------------------------------
# What needs repairing.
#
# `stored_label` is the ambiguous label currently in the archive; `label` is the
# unambiguous name to write. Where a stored label conflated two municipalities,
# both replacements list the same `stored_label`, so the conflated row is
# removed once and both correct rows are written in its place.
# ---------------------------------------------------------------------------
# `year_min`/`year_max` pin each repair to the years that indicator is actually
# the archive's source for. 0003182 and 0008273 both publish 2011-2013, and the
# app resolves that overlap in favour of the legacy indicator, so the repair
# must follow the same split or the two rows would overwrite each other.
repairs <- tibble::tribble(
  ~dataset,     ~indicator, ~code,      ~stored_label, ~label,             ~year_min, ~year_max,
  "population", "0003182",  "1711106",  "Lisboa",      "Lisboa",                1991L,    2013L,
  "population", "0003182",  "2004501",  "Calheta",     "Calheta (R.A.A.)",      1991L,    2013L,
  "population", "0003182",  "3003101",  "Calheta",     "Calheta (R.A.M.)",      1991L,    2013L,
  "population", "0003182",  "2004201",  "Lagoa",       "Lagoa (R.A.A.)",        1991L,    2013L,
  "population", "0008273",  "2004501",  "Calheta",     "Calheta (R.A.A.)",      2014L,    2100L,
  "population", "0008273",  "3003101",  "Calheta",     "Calheta (R.A.M.)",      2014L,    2100L,
  "population", "0008273",  "1500806",  "Lagoa",       "Lagoa",                 2014L,    2100L,
  "population", "0008273",  "2004201",  "Lagoa",       "Lagoa (R.A.A.)",        2014L,    2100L,
  "deaths",     "0008206",  "2004501",  "Calheta",     "Calheta (R.A.A.)",      1991L,    2100L,
  "deaths",     "0008206",  "3003101",  "Calheta",     "Calheta (R.A.M.)",      1991L,    2100L,
  "deaths",     "0008206",  "1500806",  "Lagoa",       "Lagoa",                 1991L,    2100L,
  "deaths",     "0008206",  "2004201",  "Lagoa",       "Lagoa (R.A.A.)",        1991L,    2100L
)

if (!identical(dataset_filter, "all")) {
  repairs <- repairs %>% dplyr::filter(dataset == dataset_filter)
}

message("Ambiguous-area repair")
message("  Repairs queued: ", nrow(repairs))
message("  Dry run: ", dry_run)

# ---- App environment (INE client + config) --------------------------------
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
  if (!file.rename(tmp, path)) {
    stop("Could not move temporary file into ", path, call. = FALSE)
  }
}

population_path <- function(year) file.path(out_dir, "population", paste0("year_", year, ".rds"))
death_path <- function(indicator, year, cause_token) {
  file.path(out_dir, "deaths", indicator, paste0("year_", year), paste0("cause_", cause_token, ".rds"))
}

# Years the archive actually holds for a dataset, so the repair follows the
# files on disk rather than assuming a range.
archive_years <- function(dataset, indicator) {
  if (identical(dataset, "population")) {
    files <- list.files(file.path(out_dir, "population"), pattern = "^year_\\d+\\.rds$")
    years <- suppressWarnings(as.integer(sub("^year_(\\d+)\\.rds$", "\\1", files)))
  } else {
    dirs <- list.dirs(file.path(out_dir, "deaths", indicator), recursive = FALSE, full.names = FALSE)
    years <- suppressWarnings(as.integer(sub("^year_", "", dirs)))
  }
  sort(years[!is.na(years)])
}

# The two population indicators cover different periods; only touch a year the
# indicator in question is actually the source for.
indicator_years <- function(indicator) {
  tryCatch(
    {
      values <- app_env$get_dim_values_cached(indicator)
      years <- values %>%
        dplyr::filter(as.integer(dim_num) == 1) %>%
        dplyr::pull(categ_dsg)
      sort(unique(suppressWarnings(as.integer(as.character(years)))))
    },
    error = function(e) integer(0)
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

fetch_area <- function(indicator, year, code, has_cause = FALSE, cause = NULL) {
  dims <- list(dim1 = as.character(year), dim2 = code)
  if (has_cause) dims$dim5 <- cause
  app_env$download_data(indicator, dims = dims, has_cause = has_cause)
}

repaired <- 0L
skipped <- 0L
failed <- 0L

for (i in seq_len(nrow(repairs))) {
  if (minutes_left() < 2) {
    message("Time budget exhausted; stopping cleanly.")
    break
  }

  job <- repairs[i, ]
  present <- archive_years(job$dataset, job$indicator)
  years <- Reduce(intersect, list(
    present,
    indicator_years(job$indicator),
    parse_years(years_arg, present),
    seq.int(job$year_min, job$year_max)
  ))

  if (length(years) == 0) {
    message("  ", job$indicator, " ", job$label, ": no overlapping years in archive; skipped.")
    next
  }

  message("\n== ", job$dataset, " ", job$indicator, " code ", job$code,
          " -> '", job$label, "' (", length(years), " years) ==")

  for (year in years) {
    if (minutes_left() < 1) break

    if (identical(job$dataset, "population")) {
      path <- population_path(year)
      if (!file.exists(path)) next

      existing <- readRDS(path)
      # Already repaired: the corrected label is present and the ambiguous one
      # is gone.
      if (job$label %in% existing$area && !(job$stored_label %in% existing$area)) {
        skipped <- skipped + 1L
        next
      }

      if (dry_run) {
        message("   [dry-run] would repair ", year)
        next
      }

      fetched <- tryCatch(fetch_area(job$indicator, year, job$code), error = function(e) {
        message("   ", year, ": FAILED (", conditionMessage(e), ")")
        NULL
      })

      if (is.null(fetched) || nrow(fetched) == 0) {
        failed <- failed + 1L
        next
      }

      replacement <- fetched %>%
        dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
        dplyr::rename(pop = value) %>%
        dplyr::group_by(year, sex, age_band) %>%
        dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(area = job$label, year = as.integer(year)) %>%
        dplyr::select(dplyr::any_of(names(existing)))

      updated <- existing %>%
        dplyr::filter(!area %in% c(job$stored_label, job$label)) %>%
        dplyr::bind_rows(replacement) %>%
        dplyr::arrange(area, sex, age_band)

      save_rds_atomic(updated, path)
      repaired <- repaired + 1L
      message("   ", year, ": ok (", nrow(replacement), " rows)")
    } else {
      cause_files <- list.files(file.path(out_dir, "deaths", job$indicator, paste0("year_", year)),
                                pattern = "^cause_.*\\.rds$", full.names = TRUE)

      for (path in cause_files) {
        if (minutes_left() < 1) break

        existing <- readRDS(path)
        if (job$label %in% existing$area && !(job$stored_label %in% existing$area)) {
          skipped <- skipped + 1L
          next
        }

        cause_name <- unique(as.character(existing$cause))[[1]]

        if (dry_run) {
          message("   [dry-run] would repair ", year, " / ", cause_name)
          next
        }

        fetched <- tryCatch(
          fetch_area(job$indicator, year, job$code, has_cause = TRUE, cause = cause_name),
          error = function(e) {
            message("   ", year, " / ", cause_name, ": FAILED (", conditionMessage(e), ")")
            NULL
          }
        )

        if (is.null(fetched)) {
          failed <- failed + 1L
          next
        }

        replacement <- fetched %>%
          dplyr::rename(deaths = value) %>%
          dplyr::mutate(age_band = dplyr::case_when(
            age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
            TRUE ~ age_band
          )) %>%
          dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
          dplyr::group_by(year, sex, cause, age_band) %>%
          dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
          dplyr::mutate(area = job$label, year = as.integer(year)) %>%
          dplyr::select(dplyr::any_of(names(existing)))

        updated <- existing %>%
          dplyr::filter(!area %in% c(job$stored_label, job$label)) %>%
          dplyr::bind_rows(replacement) %>%
          dplyr::arrange(area, sex, age_band)

        save_rds_atomic(updated, path)
        repaired <- repaired + 1L
      }
      message("   ", year, ": done")
    }
  }
}

message("\nRepaired ", repaired, " chunk(s); ", skipped, " already correct; ", failed, " failed.")
if (failed > 0) message("Re-run to retry the failures; the repair is idempotent.")
