#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("build_death_snapshot_chunks.R", mustWork = FALSE)
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

parse_limit <- function(value) {
  value <- value_or_default(value, "Inf")
  if (tolower(value) %in% c("inf", "infinite", "all")) {
    return(Inf)
  }
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 0) {
    return(Inf)
  }
  out
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

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
app_file <- value_or_default(cli$app, find_default_app_file())
out_dir <- value_or_default(cli$out, file.path(dirname(app_file), "data", "snapshots"))
indicator <- value_or_default(cli$indicator, env_or_default("INDICATOR", "0013166"))
area_batch_size <- as.integer(value_or_default(cli$area_batch_size, env_or_default("AREA_BATCH_SIZE", "12")))
cause_batch_size <- as.integer(value_or_default(cli$cause_batch_size, env_or_default("CAUSE_BATCH_SIZE", "66")))
year_batch_size <- as.integer(value_or_default(cli$year_batch_size, env_or_default("YEAR_BATCH_SIZE", "1")))
max_batches <- parse_limit(value_or_default(cli$max_batches, env_or_default("MAX_BATCHES", "Inf")))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

app <- new.env(parent = globalenv())
source(app_file, local = app)

available_years <- sort(unique(app$get_indicator_years(indicator)))
default_years <- sort(intersect(app$year_of_interest, available_years))
cause_values <- app$get_dimension_categories(indicator, target_dim_num = 5) %>%
  dplyr::arrange(.data$categ_ord_num, .data$categ_dsg) %>%
  dplyr::pull(.data$categ_dsg) %>%
  unique()

years <- parse_years(value_or_default(cli$years, env_or_default("YEARS", "ALL")), default_years)
areas <- parse_values(value_or_default(cli$areas, env_or_default("AREAS", "Portugal|Norte")), app$local_area)
causes <- parse_values(value_or_default(cli$causes, env_or_default("CAUSES", "ALL")), cause_values)

years <- sort(unique(as.integer(years)), decreasing = TRUE)
years <- years[years %in% available_years]
areas <- unique(as.character(areas))
causes <- unique(as.character(causes))

if (length(years) == 0 || length(areas) == 0 || length(causes) == 0) {
  stop("Death snapshot chunks require at least one valid year, area, and cause.", call. = FALSE)
}

source_priority <- dplyr::case_when(
  indicator == app$death_indicator_current ~ 2L,
  indicator == app$death_indicator_legacy ~ 1L,
  TRUE ~ 1L
)

death_path <- function(year, cause) {
  file.path(
    out_dir,
    "deaths",
    indicator,
    paste0("year_", year),
    paste0("cause_", app$snapshot_file_token(cause), ".rds")
  )
}

batch_needs_build <- function(years, areas, causes) {
  for (year in years) {
    for (cause in causes) {
      if (length(missing_areas_for_path(death_path(year, cause), areas)) > 0) {
        return(TRUE)
      }
    }
  }
  FALSE
}

normalize_deaths <- function(years, areas, causes) {
  app$download_data(
    indicator,
    dims = list(dim1 = years, dim2 = areas, dim5 = causes),
    has_cause = TRUE
  ) %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(
      age_band = dplyr::case_when(
        .data$age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ .data$age_band
      )
    ) %>%
    dplyr::filter(!.data$age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::group_by(.data$year, .data$area, .data$sex, .data$cause, .data$age_band) %>%
    dplyr::summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = source_priority)
}

save_death_chunks <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(0L)
  }

  age_order <- if (exists("age_levels", envir = app, inherits = FALSE)) app$age_levels else unique(data$age_band)
  chunks <- data %>%
    dplyr::group_by(.data$year, .data$cause) %>%
    dplyr::group_split()

  written <- 0L
  for (chunk in chunks) {
    year <- unique(chunk$year)
    cause <- unique(chunk$cause)
    chunk_areas <- unique(chunk$area)
    path <- death_path(year, cause)

    existing <- if (file.exists(path)) readRDS(path) else NULL
    if (!is.null(existing) && nrow(existing) > 0) {
      existing <- existing %>%
        dplyr::filter(!(.data$year %in% year & .data$cause %in% cause & .data$area %in% chunk_areas))
      chunk <- dplyr::bind_rows(existing, chunk)
    }

    chunk <- chunk %>%
      dplyr::mutate(age_band = factor(.data$age_band, levels = age_order, ordered = TRUE)) %>%
      dplyr::arrange(.data$area, .data$sex, .data$age_band) %>%
      dplyr::mutate(age_band = as.character(.data$age_band)) %>%
      dplyr::select(year, area, sex, cause, age_band, deaths)

    save_rds_atomic(chunk, path)
    written <- written + 1L
  }

  written
}

message("Death snapshot chunks")
message("  Indicator: ", indicator)
message("  Years: ", paste(years, collapse = ", "))
message("  Areas: ", paste(areas, collapse = " | "))
message("  Causes: ", if (identical(sort(causes), sort(cause_values))) "ALL" else paste(causes, collapse = " | "))

year_batches <- chunk_vector(years, year_batch_size)
area_batches <- chunk_vector(areas, area_batch_size)
cause_batches <- chunk_vector(causes, cause_batch_size)

batches_done <- 0L
chunks_written <- 0L

for (year_batch in year_batches) {
  for (area_batch in area_batches) {
    for (cause_batch in cause_batches) {
      if (!batch_needs_build(year_batch, area_batch, cause_batch)) {
        next
      }
      if (batches_done >= max_batches) {
        break
      }

      message(
        "Building ",
        indicator,
        " / years ", paste(year_batch, collapse = ", "),
        " / areas ", paste(area_batch, collapse = " | "),
        " / causes ", length(cause_batch)
      )

      data <- tryCatch(
        normalize_deaths(year_batch, area_batch, cause_batch),
        error = function(e) {
          warning(conditionMessage(e), call. = FALSE)
          NULL
        }
      )

      chunks_written <- chunks_written + save_death_chunks(data)
      batches_done <- batches_done + 1L
    }

    if (batches_done >= max_batches) {
      break
    }
  }

  if (batches_done >= max_batches) {
    break
  }
}

manifest <- list(
  created_at = Sys.time(),
  app_file = app_file,
  indicator = indicator,
  years_requested = years,
  areas_requested = areas,
  causes_requested = causes,
  batches_built_this_run = batches_done,
  chunks_written_this_run = chunks_written
)
save_rds_atomic(manifest, file.path(out_dir, paste0("chunked_", indicator, "_manifest.rds")))

message("Done. Batches built: ", batches_done, ". RDS chunks written: ", chunks_written, ".")
