#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("build_0008206_snapshot_chunks.R", mustWork = FALSE)
}
script_dir <- dirname(script_path)

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
    file.path(script_dir, "mortality-shiny-app.R"),
    file.path(dirname(script_dir), "mortality-shiny-app.R")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Could not find mortality Shiny app file. Pass app=/path/to/mortality-shiny-app.R.", call. = FALSE)
  }
  found[[1]]
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
app_file <- value_or_default(cli$app, find_default_app_file())
out_dir <- value_or_default(cli$out, file.path(dirname(app_file), "data", "snapshots"))
max_chunks <- as.integer(value_or_default(cli$max_chunks, "20"))
population_area_batch_size <- as.integer(value_or_default(cli$population_area_batch_size, "50"))
death_area_batch_size <- as.integer(value_or_default(cli$death_area_batch_size, "50"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

app <- new.env(parent = globalenv())
source(app_file, local = app)

all_death_years <- app$get_source_year_plan(
  indicators = c(app$death_indicator_legacy, app$death_indicator_current),
  priorities = c(1L, 2L),
  requested_years = app$year_of_interest,
  year_order = "asc"
) %>%
  dplyr::filter(indicator == app$death_indicator_legacy) %>%
  dplyr::pull(years) %>%
  unlist(use.names = FALSE) %>%
  as.integer() %>%
  sort()

years <- parse_years(cli$years, all_death_years)
areas <- parse_values(cli$areas, app$local_area)
causes <- parse_values(cli$causes, app$diseases)

years <- sort(unique(as.integer(years)))
areas <- sort(unique(as.character(areas)))
causes <- sort(unique(as.character(causes)))

if (length(years) == 0 || length(areas) == 0 || length(causes) == 0) {
  stop("Chunked snapshot requires at least one year, area, and cause.", call. = FALSE)
}

population_path <- function(year) {
  file.path(out_dir, "population", paste0("year_", year, ".rds"))
}

death_path <- function(year, cause) {
  file.path(
    out_dir,
    "deaths",
    "0008206",
    paste0("year_", year),
    paste0("cause_", app$snapshot_file_token(cause), ".rds")
  )
}

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

normalize_deaths <- function(indicator, years, areas, cause, source_priority) {
  raw <- app$download_data(
    indicator,
    dims = list(dim1 = years, dim2 = areas, dim5 = cause),
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

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  invisible(file.rename(tmp, path))
}

build_population_year <- function(year) {
  path <- population_path(year)
  if (file.exists(path)) {
    message("Population year ", year, " already exists")
    return(invisible(FALSE))
  }

  population_plan <- app$get_source_year_plan(
    indicators = c(app$population_indicator_current, app$population_indicator_legacy),
    priorities = c(1L, 2L),
    requested_years = year,
    year_order = "asc"
  )

  parts <- list()
  for (plan_i in seq_len(nrow(population_plan))) {
    for (area_batch in chunk_vector(areas, population_area_batch_size)) {
      message(
        "Population ",
        population_plan$indicator[[plan_i]],
        " year ", year,
        " areas ", length(area_batch)
      )
      parts[[length(parts) + 1L]] <- normalize_population(
        indicator = population_plan$indicator[[plan_i]],
        years = year,
        areas = area_batch,
        source_priority = population_plan$source_priority[[plan_i]]
      )
    }
  }

  population <- dplyr::bind_rows(parts) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)

  save_rds_atomic(population, path)
  message("Saved ", path, " rows ", nrow(population))
  invisible(TRUE)
}

build_death_chunk <- function(year, cause) {
  path <- death_path(year, cause)
  if (file.exists(path)) {
    message("Deaths year ", year, " cause '", cause, "' already exists")
    return(invisible(FALSE))
  }

  parts <- list()
  for (area_batch in chunk_vector(areas, death_area_batch_size)) {
    message(
      "Deaths 0008206 year ", year,
      " cause '", cause,
      "' areas ", length(area_batch)
    )
    parts[[length(parts) + 1L]] <- normalize_deaths(
      indicator = app$death_indicator_legacy,
      years = year,
      areas = area_batch,
      cause = cause,
      source_priority = 1L
    )
  }

  deaths <- dplyr::bind_rows(parts) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)

  save_rds_atomic(deaths, path)
  message("Saved ", path, " rows ", nrow(deaths))
  invisible(TRUE)
}

grid <- tidyr::expand_grid(year = years, cause = causes) %>%
  dplyr::arrange(year, cause) %>%
  dplyr::mutate(path = vapply(seq_len(dplyr::n()), function(i) death_path(.data$year[[i]], .data$cause[[i]]), character(1))) %>%
  dplyr::filter(!file.exists(.data$path))

if (nrow(grid) == 0) {
  message("All requested 0008206 chunks already exist.")
  quit(save = "no", status = 0)
}

todo <- utils::head(grid, max_chunks)
message("Building ", nrow(todo), " missing 0008206 chunk(s) out of ", nrow(grid), " remaining.")
message("Output directory: ", out_dir)

for (year in sort(unique(todo$year))) {
  build_population_year(year)
}

built <- 0L
for (i in seq_len(nrow(todo))) {
  if (isTRUE(build_death_chunk(todo$year[[i]], todo$cause[[i]]))) {
    built <- built + 1L
  }
}

manifest <- list(
  created_at = Sys.time(),
  app_file = normalizePath(app_file, mustWork = FALSE),
  years_requested = years,
  areas_requested = areas,
  causes_requested = causes,
  chunks_built_this_run = built,
  chunks_remaining_after_selection = nrow(grid) - nrow(todo)
)
save_rds_atomic(manifest, file.path(out_dir, "chunked_0008206_manifest.rds"))

message("Built death chunks this run: ", built)
message("Remaining missing chunks after this selection: ", nrow(grid) - nrow(todo))
