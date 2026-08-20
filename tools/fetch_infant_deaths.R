#!/usr/bin/env Rscript
# Fetch deaths under 1 year of age by municipality and cause.
#
#   Rscript tools/fetch_infant_deaths.R [years=ALL] [minutes=90] [overwrite=false]
#
# Both death indicators publish `Menos de 1 ano` as its own age band, but the
# main pipeline recodes it into `0 - 4 anos` on ingest, so the committed death
# archive cannot distinguish an infant death from a death at age four. For
# Portugal 2024 the archive holds 286 deaths in `0 - 4 anos`, which is 254
# infant deaths plus 32 at ages one to four.
#
# Rather than change the age bands everywhere - which would alter every existing
# rate - this writes a separate, parallel dataset holding only the under-1 band.
# The main pipeline is untouched; the infant mortality metric reads this instead.
#
# One request per year returns every area and every cause, so the whole series
# costs one call per year rather than one per area and cause.
#
# The by-cause indicators end at 2024. 0012540 publishes under-1 deaths for
# 2025, but only by municipality and month - no cause, no sex - so a year taken
# from it holds a single all-cause, both-sexes total per area. It is used only
# for years the by-cause indicators do not reach, and the rows carry its
# indicator id so the app can tell that the detail is absent rather than zero.
# It reconciles exactly with the by-cause archive where they overlap: 213 for
# 2023 and 254 for 2024, nationally.

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
budget_minutes <- suppressWarnings(as.numeric(get_arg("minutes", "90")))
if (!is.finite(budget_minutes) || budget_minutes <= 0) budget_minutes <- 90
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

# Same precedence the death pipeline uses: the current indicator wins where the
# two overlap, the historical one covers the years before it starts.
death_sources <- list(
  list(indicator = "0013166", priority = 2L),
  list(indicator = "0008206", priority = 1L)
)

# Under-1 deaths without a cause or sex breakdown. Reaches a year further than
# the by-cause indicators, and is used only where they cannot.
total_only_source <- list(indicator = "0012540", period_prefix = "S7A", total_code = "T")

pseudo_areas <- c("Total", "Ignorado", "Estrangeiro")
infant_band <- "Menos de 1 ano"
all_cause_label <- "Todas as causas de morte"

indicator_years <- function(indicator) {
  tryCatch(
    {
      app_env$get_dim_values_cached(indicator) %>%
        dplyr::filter(as.integer(dim_num) == 1) %>%
        dplyr::pull(categ_dsg) %>%
        as.character() %>%
        as.integer() %>%
        unique() %>%
        sort()
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

infant_path <- function(year) file.path(out_dir, "infant_deaths", paste0("year_", year, ".rds"))

# Which indicator should supply each year, highest priority first, so a year is
# fetched once from the best available source.
year_plan <- list()
for (source in death_sources) {
  available <- indicator_years(source$indicator)
  if (is.null(available)) next
  for (year in available) {
    key <- as.character(year)
    if (is.null(year_plan[[key]]) || year_plan[[key]]$priority < source$priority) {
      year_plan[[key]] <- list(indicator = source$indicator, priority = source$priority)
    }
  }
}

# Years only the total-only indicator reaches.
total_only_years <- setdiff(
  indicator_years(total_only_source$indicator),
  as.integer(names(year_plan))
)
if (length(total_only_years) > 0) {
  for (year in total_only_years) {
    year_plan[[as.character(year)]] <- list(
      indicator = total_only_source$indicator,
      priority = 0L,
      total_only = TRUE
    )
  }
}

requested <- parse_years(years_arg, sort(as.integer(names(year_plan))))
plan_years <- sort(intersect(as.integer(names(year_plan)), requested))

message("Fetching under-1 deaths for ", length(plan_years), " year(s), one request each")
if (length(intersect(plan_years, total_only_years)) > 0) {
  message("  ", paste(sort(intersect(plan_years, total_only_years)), collapse = ", "),
          " come from ", total_only_source$indicator, ": all causes, both sexes only.")
}

written <- 0L
skipped <- 0L
failed <- character(0)

for (year in plan_years) {
  if (minutes_left() < 1) {
    message("Time budget exhausted; stopping cleanly.")
    break
  }

  path <- infant_path(year)
  if (file.exists(path) && !overwrite) {
    skipped <- skipped + 1L
    next
  }

  plan <- year_plan[[as.character(year)]]
  indicator <- plan$indicator
  total_only <- isTRUE(plan$total_only)

  if (total_only) {
    # A different shape of response: no cause or sex dimension, and periods
    # keyed by code rather than label.
    fetched <- tryCatch(
      app_env$ine_client$get_data(
        indicator,
        dim1 = paste0(total_only_source$period_prefix, year),
        dim3 = total_only_source$total_code
      ),
      error = function(e) {
        message("  ", year, ": FAILED (", conditionMessage(e), ")")
        failed <<- c(failed, as.character(year))
        NULL
      }
    )

    if (is.null(fetched) || nrow(fetched) == 0) next

    chunk <- fetched %>%
      dplyr::transmute(
        year = as.integer(year),
        area = as.character(geodsg),
        sex = "HM",
        cause = all_cause_label,
        deaths = suppressWarnings(as.numeric(valor)),
        source_indicator = indicator
      ) %>%
      dplyr::filter(!area %in% pseudo_areas, !is.na(deaths)) %>%
      dplyr::group_by(year, area, sex, cause, source_indicator) %>%
      dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
      dplyr::select(year, area, sex, cause, deaths, source_indicator) %>%
      dplyr::arrange(area)

    national <- chunk$deaths[chunk$area == "Portugal"]
    save_rds_atomic(chunk, path)
    written <- written + 1L
    message("  ", year, " (", indicator, ", total only): ok - ", nrow(chunk),
            " rows, Portugal all-cause ", if (length(national) == 1) national else "?")
    next
  }

  fetched <- tryCatch(
    app_env$download_data(
      indicator,
      dims = list(dim1 = as.character(year), dim4 = infant_band),
      has_cause = TRUE
    ),
    error = function(e) {
      message("  ", year, ": FAILED (", conditionMessage(e), ")")
      failed <<- c(failed, as.character(year))
      NULL
    }
  )

  if (is.null(fetched) || nrow(fetched) == 0) next

  chunk <- fetched %>%
    dplyr::rename(deaths = value) %>%
    dplyr::filter(!area %in% pseudo_areas) %>%
    dplyr::group_by(area, sex, cause) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(year = as.integer(year), source_indicator = indicator) %>%
    dplyr::select(year, area, sex, cause, deaths, source_indicator) %>%
    dplyr::arrange(area, sex, cause)

  national <- chunk %>%
    dplyr::filter(area == "Portugal", sex == "HM", cause == "Todas as causas de morte") %>%
    dplyr::pull(deaths)

  save_rds_atomic(chunk, path)
  written <- written + 1L
  message("  ", year, " (", indicator, "): ok - ", nrow(chunk), " rows, Portugal all-cause ",
          if (length(national) == 1) national else "?")
}

message("Done: ", written, " written, ", skipped, " already present, ", length(failed), " failed.")
