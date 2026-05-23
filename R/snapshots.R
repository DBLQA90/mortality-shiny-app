# =========================================================
# Snapshot files
# =========================================================

get_snapshot_dir <- function() {
  Sys.getenv(
    "MORTALITY_SNAPSHOT_DIR",
    file.path(get_app_dir(), "data", "snapshots")
  )
}

get_combined_snapshot_file <- function() {
  Sys.getenv(
    "MORTALITY_SNAPSHOT_RDS",
    file.path(get_snapshot_dir(), "mortality_ine_snapshot.rds")
  )
}

get_snapshot_file <- function(kind) {
  env_name <- switch(
    kind,
    population = "MORTALITY_POPULATION_SNAPSHOT_RDS",
    deaths = "MORTALITY_DEATHS_SNAPSHOT_RDS",
    stop("Unknown snapshot kind: ", kind, call. = FALSE)
  )

  Sys.getenv(
    env_name,
    file.path(get_snapshot_dir(), paste0(kind, ".rds"))
  )
}

get_snapshot_inventory_file <- function() {
  Sys.getenv(
    "MORTALITY_SNAPSHOT_INVENTORY_RDS",
    file.path(get_snapshot_dir(), "snapshot_inventory.rds")
  )
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", as.character(path))
}

snapshot_inventory_full_path <- function(path) {
  path <- as.character(path)
  ifelse(is_absolute_path(path), path, file.path(get_snapshot_dir(), path))
}

normalize_data_source <- function(data_source = "ine") {
  if (is.null(data_source) || length(data_source) == 0 || is.na(data_source[[1]])) {
    return("ine")
  }

  data_source <- as.character(data_source[[1]])
  if (!data_source %in% unname(data_source_choices)) {
    return("ine")
  }
  data_source
}

get_data_source_label <- function(data_source = "ine") {
  data_source <- normalize_data_source(data_source)
  label <- names(data_source_choices)[match(data_source, unname(data_source_choices))]
  ifelse(is.na(label), data_source, label)
}

get_data_source_progress_message <- function(data_source = "ine", context = "data") {
  data_source <- normalize_data_source(data_source)
  if (identical(data_source, "snapshot")) {
    return(if (identical(context, "annual")) "A obter métricas anuais dos ficheiros RDS..." else "A obter dados dos ficheiros RDS...")
  }

  if (identical(context, "annual")) "A obter métricas anuais do INE..." else "A obter dados do INE..."
}

snapshot_file_token <- function(x) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

get_chunked_death_indicators <- function() {
  tibble::tibble(
    indicator = c(death_indicator_legacy, death_indicator_current),
    source_priority = c(1L, 2L)
  )
}

get_chunked_death_snapshot_path <- function(indicator, year, cause) {
  file.path(
    get_snapshot_dir(),
    "deaths",
    indicator,
    paste0("year_", year),
    paste0("cause_", snapshot_file_token(cause), ".rds")
  )
}

read_snapshot_object <- memoise::memoise(function(path) {
  readRDS(path)
})

read_snapshot_dataset <- memoise::memoise(function(kind, separate_path, combined_path) {
  if (file.exists(separate_path)) {
    return(read_snapshot_object(separate_path))
  }

  if (file.exists(combined_path)) {
    snapshot <- read_snapshot_object(combined_path)
    if (is.list(snapshot) && !is.null(snapshot[[kind]])) {
      return(snapshot[[kind]])
    }
  }

  stop(
    glue::glue(
      "Snapshot RDS file not found for '{kind}'. Expected '{separate_path}' or combined snapshot '{combined_path}'."
    ),
    call. = FALSE
  )
})

get_snapshot_dataset <- function(kind) {
  read_snapshot_dataset(
    kind = kind,
    separate_path = get_snapshot_file(kind),
    combined_path = get_combined_snapshot_file()
  )
}

validate_snapshot_columns <- function(data, required_cols, label) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      glue::glue(
        "{label} snapshot is missing required column(s): {paste(missing_cols, collapse = ', ')}."
      ),
      call. = FALSE
    )
  }

  data
}

read_snapshot_inventory_object <- memoise::memoise(function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  inventory <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(inventory)) {
    warning(
      glue::glue("Snapshot inventory '{path}' could not be read; falling back to direct chunk lookup."),
      call. = FALSE
    )
    return(NULL)
  }

  required_cols <- c("dataset", "indicator", "year", "cause", "path", "rows", "areas")
  missing_cols <- setdiff(required_cols, names(inventory))
  if (length(missing_cols) > 0) {
    warning(
      glue::glue("Snapshot inventory is missing required column(s): {paste(missing_cols, collapse = ', ')}; falling back to direct chunk lookup."),
      call. = FALSE
    )
    return(NULL)
  }

  inventory %>%
    dplyr::mutate(
      dataset = as.character(.data$dataset),
      indicator = as.character(.data$indicator),
      year = as.integer(.data$year),
      cause = as.character(.data$cause),
      path = as.character(.data$path),
      rows = as.integer(.data$rows)
    )
})

get_snapshot_inventory <- function() {
  read_snapshot_inventory_object(get_snapshot_inventory_file())
}

snapshot_inventory_has_area <- function(area_lists, requested_areas = NULL) {
  if (is.null(requested_areas)) {
    return(rep(TRUE, length(area_lists)))
  }

  requested_areas <- as.character(requested_areas)
  vapply(
    area_lists,
    function(areas) any(requested_areas %in% as.character(areas)),
    logical(1)
  )
}

snapshot_relative_path <- function(path, snapshot_dir = get_snapshot_dir()) {
  root <- normalizePath(snapshot_dir, mustWork = FALSE)
  full_path <- normalizePath(path, mustWork = FALSE)

  if (startsWith(full_path, paste0(root, .Platform$file.sep))) {
    return(substring(full_path, nchar(root) + 2L))
  }

  full_path
}

summarize_snapshot_chunk <- function(path, dataset, indicator = NA_character_) {
  chunk <- read_snapshot_object(path)
  cause <- if ("cause" %in% names(chunk)) sort(unique(as.character(chunk$cause))) else NA_character_

  tibble::tibble(
    dataset = dataset,
    indicator = indicator,
    year = if ("year" %in% names(chunk)) as.integer(sort(unique(chunk$year))[[1]]) else NA_integer_,
    cause = if (length(cause) > 0) cause[[1]] else NA_character_,
    path = snapshot_relative_path(path),
    rows = nrow(chunk),
    areas = list(if ("area" %in% names(chunk)) sort(unique(as.character(chunk$area))) else character(0)),
    sexes = list(if ("sex" %in% names(chunk)) sort(unique(as.character(chunk$sex))) else character(0)),
    age_bands = list(if ("age_band" %in% names(chunk)) sort(unique(as.character(chunk$age_band))) else character(0))
  )
}

build_snapshot_inventory <- function(snapshot_dir = get_snapshot_dir()) {
  population_paths <- list.files(
    file.path(snapshot_dir, "population"),
    pattern = "^year_[0-9]+\\.rds$",
    full.names = TRUE
  )
  death_paths <- list.files(
    file.path(snapshot_dir, "deaths"),
    pattern = "^cause_.*\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  population_inventory <- purrr::map_dfr(
    population_paths,
    summarize_snapshot_chunk,
    dataset = "population",
    indicator = "population"
  )

  death_inventory <- purrr::map_dfr(death_paths, function(path) {
    relative <- snapshot_relative_path(path, snapshot_dir)
    parts <- strsplit(relative, .Platform$file.sep, fixed = TRUE)[[1]]
    indicator <- if (length(parts) >= 3) parts[[2]] else NA_character_
    summarize_snapshot_chunk(path, dataset = "deaths", indicator = indicator)
  })

  dplyr::bind_rows(population_inventory, death_inventory) %>%
    dplyr::left_join(get_chunked_death_indicators(), by = "indicator") %>%
    dplyr::mutate(
      source_priority = dplyr::if_else(
        .data$dataset == "population",
        1L,
        dplyr::coalesce(.data$source_priority, 0L)
      )
    ) %>%
    dplyr::arrange(.data$dataset, .data$indicator, .data$year, .data$cause, .data$path)
}

write_snapshot_inventory <- function(path = get_snapshot_inventory_file(), snapshot_dir = get_snapshot_dir()) {
  inventory <- build_snapshot_inventory(snapshot_dir = snapshot_dir)
  ensure_cache_dir(dirname(path))
  saveRDS(inventory, path, version = 2)
  inventory
}

read_chunked_population_snapshot <- function(years, areas = NULL) {
  years <- sort(unique(as.integer(years)))
  inventory <- get_snapshot_inventory()

  if (!is.null(inventory)) {
    selected <- inventory %>%
      dplyr::filter(.data$dataset == "population", .data$year %in% .env$years)

    if (nrow(selected) > 0) {
      selected <- selected[snapshot_inventory_has_area(selected$areas, areas), , drop = FALSE]
      missing_years <- setdiff(years, unique(selected$year))

      if (length(missing_years) > 0) {
        stop(
          glue::glue("Missing population snapshot chunk(s) for year(s): {paste(missing_years, collapse = ', ')}."),
          call. = FALSE
        )
      }

      paths <- snapshot_inventory_full_path(selected$path)
      missing_paths <- paths[!file.exists(paths)]
      if (length(missing_paths) > 0) {
        stop(
          glue::glue("Missing population snapshot chunk file(s): {paste(missing_paths, collapse = ', ')}."),
          call. = FALSE
        )
      }

      return(dplyr::bind_rows(lapply(paths, read_snapshot_object)))
    }
  }

  paths <- file.path(get_snapshot_dir(), "population", paste0("year_", years, ".rds"))

  if (!any(file.exists(paths))) {
    return(NULL)
  }

  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      glue::glue("Missing population snapshot chunk(s): {paste(missing_paths, collapse = ', ')}."),
      call. = FALSE
    )
  }

  dplyr::bind_rows(lapply(paths, read_snapshot_object))
}

collapse_death_snapshot_chunks <- function(selected) {
  rows <- purrr::map2(
    selected$path,
    selected$source_priority,
    function(path, source_priority) {
      chunk <- read_snapshot_object(path)
      if (!"source_priority" %in% names(chunk)) {
        chunk$source_priority <- source_priority
      }
      chunk
    }
  )

  dplyr::bind_rows(rows) %>%
    dplyr::group_by(.data$year, .data$area, .data$sex, .data$cause, .data$age_band) %>%
    dplyr::arrange(dplyr::desc(.data$source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)
}

read_chunked_death_snapshot <- function(years, causes, areas = NULL) {
  years <- sort(unique(as.integer(years)))
  causes <- sort(unique(as.character(causes)))

  requested <- tidyr::expand_grid(year = years, cause = causes)
  inventory <- get_snapshot_inventory()

  if (!is.null(inventory)) {
    inventory_selected <- inventory %>%
      dplyr::filter(
        .data$dataset == "deaths",
        .data$year %in% .env$years,
        .data$cause %in% .env$causes
      )

    if (nrow(inventory_selected) > 0) {
      available <- inventory_selected %>%
        dplyr::distinct(.data$year, .data$cause)

      missing <- dplyr::anti_join(requested, available, by = c("year", "cause"))
      if (nrow(missing) > 0) {
        missing_labels <- paste0(missing$year, " / ", missing$cause)
        stop(
          glue::glue("Missing deaths snapshot chunk(s) for: {paste(missing_labels, collapse = ', ')}."),
          call. = FALSE
        )
      }

      selected <- inventory_selected[snapshot_inventory_has_area(inventory_selected$areas, areas), , drop = FALSE]
      if (nrow(selected) == 0) {
        return(NULL)
      }

      selected <- selected %>%
        dplyr::select(-dplyr::any_of("source_priority")) %>%
        dplyr::left_join(get_chunked_death_indicators(), by = "indicator") %>%
        dplyr::mutate(
          source_priority = dplyr::coalesce(as.integer(.data$source_priority), 0L),
          path = snapshot_inventory_full_path(.data$path)
        )

      return(collapse_death_snapshot_chunks(selected))
    }
  }

  indicator_priority <- get_chunked_death_indicators()
  grid <- tidyr::expand_grid(requested, indicator_priority) %>%
    dplyr::mutate(
      path = purrr::pmap_chr(
        list(.data$indicator, .data$year, .data$cause),
        get_chunked_death_snapshot_path
      ),
      exists = file.exists(.data$path)
    )

  if (!any(grid$exists)) {
    return(NULL)
  }

  available <- grid %>%
    dplyr::filter(.data$exists) %>%
    dplyr::distinct(.data$year, .data$cause)

  missing <- dplyr::anti_join(requested, available, by = c("year", "cause"))
  if (nrow(missing) > 0) {
    missing_labels <- paste0(missing$year, " / ", missing$cause)
    stop(
      glue::glue("Missing deaths snapshot chunk(s) for: {paste(missing_labels, collapse = ', ')}."),
      call. = FALSE
    )
  }

  # Keep all available indicators for overlapping year/cause pairs. The final
  # row-level collapse lets 0008206 fill areas missing from 0013166, while still
  # preferring 0013166 wherever it exists.
  selected <- grid %>%
    dplyr::filter(.data$exists)

  collapse_death_snapshot_chunks(selected)
}

validate_snapshot_combinations <- function(data, required, by, label) {
  present <- data %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(by)))
  missing <- dplyr::anti_join(required, present, by = by)

  if (nrow(missing) == 0) {
    return(invisible(TRUE))
  }

  missing_labels <- apply(
    missing,
    1,
    function(row) paste(paste(names(row), row, sep = "="), collapse = " / ")
  )
  if (length(missing_labels) > 10) {
    missing_labels <- c(utils::head(missing_labels, 10), glue::glue("... and {nrow(missing) - 10} more"))
  }

  stop(
    glue::glue("{label} snapshot is missing requested combination(s): {paste(missing_labels, collapse = '; ')}."),
    call. = FALSE
  )
}

get_flat_or_chunked_population_snapshot <- function(years, areas = NULL) {
  tryCatch(
    get_snapshot_dataset("population"),
    error = function(e) read_chunked_population_snapshot(years, areas = areas)
  )
}

get_flat_or_chunked_death_snapshot <- function(years, causes, areas = NULL) {
  tryCatch(
    get_snapshot_dataset("deaths"),
    error = function(e) read_chunked_death_snapshot(years, causes, areas = areas)
  )
}

get_snapshot_population_data <- function(years, area) {
  year_filter <- as.integer(years)
  area_filter <- as.character(area)
  pop_snapshot <- get_flat_or_chunked_population_snapshot(years, areas = area_filter)
  if (is.null(pop_snapshot)) {
    stop(
      glue::glue(
        "Snapshot RDS file not found for 'population'. Expected '{get_snapshot_file('population')}', combined snapshot '{get_combined_snapshot_file()}', or chunked files in '{file.path(get_snapshot_dir(), 'population')}'."
      ),
      call. = FALSE
    )
  }

  pop <- pop_snapshot %>%
    validate_snapshot_columns(
      required_cols = c("year", "area", "sex", "age_band", "pop"),
      label = "Population"
    ) %>%
    dplyr::mutate(
      year = as.integer(year),
      area = as.character(area),
      sex = as.character(sex),
      age_band = as.character(age_band),
      pop = as.numeric(pop)
    ) %>%
    dplyr::filter(
      .data$year %in% .env$year_filter,
      .data$area %in% .env$area_filter
    ) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop")

  validate_snapshot_combinations(
    data = pop,
    required = tidyr::expand_grid(year = year_filter, area = area_filter),
    by = c("year", "area"),
    label = "Population"
  )

  if (nrow(pop) == 0) {
    stop(
      "The population snapshot has no rows for the selected years/areas. Choose INE or rebuild the snapshot.",
      call. = FALSE
    )
  }

  pop
}

get_snapshot_death_data <- function(years, area, cause) {
  year_filter <- as.integer(years)
  area_filter <- as.character(area)
  cause_filter <- as.character(cause)
  death_snapshot <- get_flat_or_chunked_death_snapshot(years, cause, areas = area_filter)
  if (is.null(death_snapshot)) {
    stop(
      glue::glue(
        "Snapshot RDS file not found for 'deaths'. Expected '{get_snapshot_file('deaths')}', combined snapshot '{get_combined_snapshot_file()}', or chunked files in '{file.path(get_snapshot_dir(), 'deaths', '<indicator>')}'."
      ),
      call. = FALSE
    )
  }

  deaths <- death_snapshot %>%
    validate_snapshot_columns(
      required_cols = c("year", "area", "sex", "cause", "age_band", "deaths"),
      label = "Deaths"
    ) %>%
    dplyr::mutate(
      year = as.integer(year),
      area = as.character(area),
      sex = as.character(sex),
      cause = as.character(cause),
      age_band = as.character(age_band),
      deaths = as.numeric(deaths)
    ) %>%
    dplyr::filter(
      .data$year %in% .env$year_filter,
      .data$area %in% .env$area_filter,
      .data$cause %in% .env$cause_filter
    ) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

  validate_snapshot_combinations(
    data = deaths,
    required = tidyr::expand_grid(year = year_filter, area = area_filter, cause = cause_filter),
    by = c("year", "area", "cause"),
    label = "Deaths"
  )

  if (nrow(deaths) == 0) {
    stop(
      "The deaths snapshot has no rows for the selected years/areas/causes. Choose INE or rebuild the snapshot.",
      call. = FALSE
    )
  }

  deaths
}
