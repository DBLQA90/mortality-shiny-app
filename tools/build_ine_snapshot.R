#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("build_ine_snapshot.R", mustWork = FALSE)
}
script_dir <- dirname(script_path)

trailing_args <- commandArgs(trailingOnly = TRUE)

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    out[[key]] <- value
  }
  out
}

cli <- parse_cli_args(trailing_args)

value_or_default <- function(x, default) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  x
}

parse_years <- function(value, default_years) {
  value <- value_or_default(value, "ALL")
  if (identical(toupper(value), "ALL")) {
    return(default_years)
  }
  if (grepl(":", value, fixed = TRUE)) {
    bounds <- suppressWarnings(as.integer(strsplit(value, ":", fixed = TRUE)[[1]]))
    bounds <- bounds[!is.na(bounds)]
    if (length(bounds) == 2) {
      return(seq.int(min(bounds), max(bounds)))
    }
  }
  years <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  years[!is.na(years)]
}

parse_values <- function(value, default_values) {
  value <- value_or_default(value, "ALL")
  if (identical(toupper(value), "ALL")) {
    return(default_values)
  }
  trimws(strsplit(value, "\\|")[[1]])
}

chunk_vector <- function(x, size) {
  size <- max(1L, as.integer(size))
  split(x, ceiling(seq_along(x) / size))
}

find_default_app_file <- function() {
  candidates <- c(
    file.path(script_dir, "mortality-shiny-app-v5.R"),
    file.path(script_dir, "mortality-shiny-app.R"),
    file.path(dirname(script_dir), "mortality-shiny-app.R")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Could not find mortality Shiny app file. Pass app=/path/to/mortality-shiny-app.R.", call. = FALSE)
  }
  found[[1]]
}

app_file <- value_or_default(cli$app, find_default_app_file())
out_dir <- value_or_default(cli$out, file.path(dirname(app_file), "data", "snapshots"))
year_batch_size <- as.integer(value_or_default(cli$year_batch_size, "5"))
area_batch_size <- as.integer(value_or_default(cli$area_batch_size, "25"))
cause_batch_size <- as.integer(value_or_default(cli$cause_batch_size, "5"))
population_year_batch_size <- as.integer(value_or_default(cli$population_year_batch_size, as.character(year_batch_size)))
population_area_batch_size <- as.integer(value_or_default(cli$population_area_batch_size, as.character(area_batch_size)))
death_year_batch_size <- as.integer(value_or_default(cli$death_year_batch_size, as.character(year_batch_size)))
death_area_batch_size <- as.integer(value_or_default(cli$death_area_batch_size, as.character(area_batch_size)))
death_cause_batch_size <- as.integer(value_or_default(cli$death_cause_batch_size, as.character(cause_batch_size)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

app <- new.env(parent = globalenv())
source(app_file, local = app)

years <- parse_years(cli$years, app$year_of_interest)
areas <- parse_values(cli$areas, app$local_area)
causes <- parse_values(cli$causes, app$diseases)

years <- sort(unique(as.integer(years)))
areas <- sort(unique(as.character(areas)))
causes <- sort(unique(as.character(causes)))

if (length(years) == 0 || length(areas) == 0 || length(causes) == 0) {
  stop("Snapshot requires at least one year, area, and cause.", call. = FALSE)
}

message("Building mortality snapshot")
message("App file: ", app_file)
message("Output directory: ", out_dir)
message("Years: ", min(years), " - ", max(years), " (", length(years), ")")
message("Areas: ", length(areas))
message("Causes: ", length(causes))

normalize_population <- function(indicator, years, areas, source_priority) {
  raw <- app$download_data(
    indicator,
    dims = list(dim1 = years, dim2 = areas),
    has_cause = FALSE
  )

  raw %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = source_priority)
}

normalize_deaths <- function(indicator, years, areas, causes, source_priority) {
  raw <- app$download_data(
    indicator,
    dims = list(dim1 = years, dim2 = areas, dim5 = causes),
    has_cause = TRUE
  )

  raw %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(
      age_band = dplyr::case_when(
        age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ age_band
      )
    ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = source_priority)
}

population_plan <- app$get_source_year_plan(
  indicators = c(app$population_indicator_current, app$population_indicator_legacy),
  priorities = c(1L, 2L),
  requested_years = years,
  year_order = "asc"
)

death_plan <- app$get_source_year_plan(
  indicators = c(app$death_indicator_legacy, app$death_indicator_current),
  priorities = c(1L, 2L),
  requested_years = years,
  year_order = "asc"
)

population_parts <- list()
for (plan_i in seq_len(nrow(population_plan))) {
  for (year_batch in chunk_vector(population_plan$years[[plan_i]], population_year_batch_size)) {
    for (area_batch in chunk_vector(areas, population_area_batch_size)) {
      message(
        "Population ",
        population_plan$indicator[[plan_i]],
        " years ", min(year_batch), "-", max(year_batch),
        " areas ", length(area_batch)
      )
      population_parts[[length(population_parts) + 1L]] <- normalize_population(
        indicator = population_plan$indicator[[plan_i]],
        years = year_batch,
        areas = area_batch,
        source_priority = population_plan$source_priority[[plan_i]]
      )
    }
  }
}

death_parts <- list()
for (plan_i in seq_len(nrow(death_plan))) {
  for (year_batch in chunk_vector(death_plan$years[[plan_i]], death_year_batch_size)) {
    for (area_batch in chunk_vector(areas, death_area_batch_size)) {
      for (cause_batch in chunk_vector(causes, death_cause_batch_size)) {
        message(
          "Deaths ",
          death_plan$indicator[[plan_i]],
          " years ", min(year_batch), "-", max(year_batch),
          " areas ", length(area_batch),
          " causes ", length(cause_batch)
        )
        death_parts[[length(death_parts) + 1L]] <- normalize_deaths(
          indicator = death_plan$indicator[[plan_i]],
          years = year_batch,
          areas = area_batch,
          causes = cause_batch,
          source_priority = death_plan$source_priority[[plan_i]]
        )
      }
    }
  }
}

population <- dplyr::bind_rows(population_parts) %>%
  dplyr::group_by(year, area, sex, age_band) %>%
  dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-source_priority)

deaths <- dplyr::bind_rows(death_parts) %>%
  dplyr::group_by(year, area, sex, cause, age_band) %>%
  dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-source_priority)

saveRDS(population, file.path(out_dir, "population.rds"), version = 2)
saveRDS(deaths, file.path(out_dir, "deaths.rds"), version = 2)
saveRDS(
  list(
    created_at = Sys.time(),
    app_file = normalizePath(app_file, mustWork = FALSE),
    years = years,
    areas = areas,
    causes = causes,
    population_rows = nrow(population),
    death_rows = nrow(deaths)
  ),
  file.path(out_dir, "snapshot_manifest.rds"),
  version = 2
)

message("Snapshot complete")
message("Population rows: ", nrow(population))
message("Death rows: ", nrow(deaths))
message("Files written to: ", out_dir)
