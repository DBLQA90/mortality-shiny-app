#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("build_population_snapshot_chunks.R", mustWork = FALSE)
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

env_or_default <- function(name, default) {
  value_or_default(Sys.getenv(name, unset = ""), default)
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
  years <- suppressWarnings(as.integer(unlist(strsplit(value, "[,|]", perl = TRUE), use.names = FALSE)))
  years[!is.na(years)]
}

parse_values <- function(value, default_values) {
  value <- value_or_default(value, "ALL")
  if (identical(toupper(value), "ALL")) {
    return(default_values)
  }
  values <- unlist(strsplit(value, "\\|", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  values[nzchar(values)]
}

chunk_vector <- function(x, size) {
  size <- max(1L, as.integer(size))
  split(x, ceiling(seq_along(x) / size))
}

find_default_app_file <- function() {
  candidates <- c(
    file.path(dirname(script_dir), "mortality-shiny-app.R"),
    file.path(script_dir, "mortality-shiny-app.R")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Could not find mortality-shiny-app.R. Pass app=/path/to/mortality-shiny-app.R.", call. = FALSE)
  }
  found[[1]]
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) {
    stop(paste0("Could not move temporary file into ", path, "."), call. = FALSE)
  }
}

missing_areas_for_path <- function(path, required_areas) {
  required_areas <- unique(as.character(required_areas))
  if (!file.exists(path)) {
    return(required_areas)
  }

  existing <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(existing) || !"area" %in% names(existing)) {
    return(required_areas)
  }

  setdiff(required_areas, unique(as.character(existing$area)))
}

collapse_population <- function(parts) {
  dplyr::bind_rows(parts) %>%
    dplyr::group_by(.data$year, .data$area, .data$sex, .data$age_band) %>%
    dplyr::arrange(dplyr::desc(.data$source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
app_file <- value_or_default(cli$app, find_default_app_file())
out_dir <- value_or_default(cli$out, file.path(dirname(app_file), "data", "snapshots"))
area_batch_size <- as.integer(value_or_default(cli$area_batch_size, env_or_default("AREA_BATCH_SIZE", "50")))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

app <- new.env(parent = globalenv())
source(app_file, local = app)

available_years <- sort(unique(app$get_indicator_years(app$population_indicators)))
default_years <- sort(intersect(app$year_of_interest, available_years))
years <- parse_years(value_or_default(cli$years, env_or_default("YEARS", "ALL")), default_years)
areas <- parse_values(value_or_default(cli$areas, env_or_default("AREAS", "Portugal|Norte")), app$local_area)

years <- sort(unique(as.integer(years)), decreasing = TRUE)
years <- years[years %in% available_years]
areas <- unique(as.character(areas))

if (length(years) == 0 || length(areas) == 0) {
  stop("Population snapshot chunks require at least one valid year and area.", call. = FALSE)
}

population_path <- function(year) {
  file.path(out_dir, "population", paste0("year_", year, ".rds"))
}

normalize_population <- function(indicator, years, areas, source_priority) {
  app$download_data(
    indicator,
    dims = list(dim1 = years, dim2 = areas),
    has_cause = FALSE
  ) %>%
    dplyr::filter(!.data$age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value) %>%
    dplyr::group_by(.data$year, .data$area, .data$sex, .data$age_band) %>%
    dplyr::summarise(pop = sum(.data$pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = source_priority)
}

message("Population snapshot chunks")
message("  Years: ", paste(years, collapse = ", "))
message("  Areas: ", paste(areas, collapse = " | "))

built <- 0L
for (year in years) {
  path <- population_path(year)
  missing_areas <- missing_areas_for_path(path, areas)
  if (length(missing_areas) == 0) {
    message("Skipping covered year ", year)
    next
  }

  existing <- if (file.exists(path)) {
    readRDS(path) %>% dplyr::mutate(source_priority = 0L)
  } else {
    NULL
  }

  population_plan <- app$get_source_year_plan(
    indicators = c(app$population_indicator_current, app$population_indicator_legacy),
    priorities = c(1L, 2L),
    requested_years = year,
    year_order = "asc"
  )

  parts <- list()
  for (plan_i in seq_len(nrow(population_plan))) {
    for (area_batch in chunk_vector(missing_areas, area_batch_size)) {
      message(
        "Building population ",
        population_plan$indicator[[plan_i]],
        " / ", year,
        " / areas ", length(area_batch)
      )
      parts[[length(parts) + 1L]] <- normalize_population(
        indicator = population_plan$indicator[[plan_i]],
        years = year,
        areas = area_batch,
        source_priority = population_plan$source_priority[[plan_i]]
      )
    }
  }

  population <- collapse_population(c(list(existing), parts)[!vapply(c(list(existing), parts), is.null, logical(1))])
  save_rds_atomic(population, path)
  message("Saved ", path, " rows ", nrow(population))
  built <- built + 1L
}

manifest <- list(
  created_at = Sys.time(),
  app_file = app_file,
  years_requested = years,
  areas_requested = areas,
  years_built_this_run = built
)
save_rds_atomic(manifest, file.path(out_dir, "chunked_population_manifest.rds"))

message("Done. Population year chunks built: ", built, ".")
