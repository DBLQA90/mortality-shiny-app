# =========================================================
# Snapshot files
# =========================================================

get_default_remote_snapshot_dir <- function() {
  Sys.getenv(
    "MORTALITY_REMOTE_SNAPSHOT_DIR",
    "https://raw.githubusercontent.com/DBLQA90/mortality-shiny-app/main/data/snapshots"
  )
}

is_url_path <- function(path) {
  grepl("^https?://", as.character(path), ignore.case = TRUE)
}

snapshot_path_join <- function(...) {
  args <- lapply(list(...), as.character)

  if (length(args) == 0) {
    return("")
  }

  max_length <- max(lengths(args))
  args <- lapply(args, rep_len, length.out = max_length)

  do.call(
    mapply,
    c(
      list(FUN = function(...) {
        parts <- as.character(c(...))
        parts <- parts[nzchar(parts)]

        if (length(parts) == 0) {
          return("")
        }

        if (is_url_path(parts[[1]])) {
          return(paste(c(sub("/+$", "", parts[[1]]), gsub("^/+|/+$", "", parts[-1])), collapse = "/"))
        }

        do.call(file.path, as.list(parts))
      }),
      args,
      list(USE.NAMES = FALSE)
    )
  )
}

snapshot_dir_has_local_data <- function(path) {
  if (!dir.exists(path)) {
    return(FALSE)
  }

  if (file.exists(file.path(path, "mortality_ine_snapshot.rds"))) {
    return(TRUE)
  }

  has_population <- file.exists(file.path(path, "population.rds"))
  has_deaths <- file.exists(file.path(path, "deaths.rds"))

  for (subdir in "population") {
    chunk_dir <- file.path(path, subdir)
    if (
      dir.exists(chunk_dir) &&
        length(list.files(chunk_dir, pattern = "\\.rds$", recursive = TRUE, full.names = FALSE)) > 0
    ) {
      has_population <- TRUE
    }
  }

  for (subdir in "deaths") {
    chunk_dir <- file.path(path, subdir)
    if (
      dir.exists(chunk_dir) &&
        length(list.files(chunk_dir, pattern = "\\.rds$", recursive = TRUE, full.names = FALSE)) > 0
    ) {
      has_deaths <- TRUE
    }
  }

  isTRUE(has_population && has_deaths)
}

get_snapshot_dir <- function() {
  configured <- Sys.getenv("MORTALITY_SNAPSHOT_DIR", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }

  local_snapshot_dir <- file.path(get_app_dir(), "data", "snapshots")
  use_local <- !tolower(Sys.getenv("MORTALITY_USE_LOCAL_SNAPSHOTS", unset = "true")) %in% c("0", "false", "no", "n")
  if (isTRUE(use_local) && snapshot_dir_has_local_data(local_snapshot_dir)) {
    return(local_snapshot_dir)
  }

  get_default_remote_snapshot_dir()
}

get_combined_snapshot_file <- function() {
  Sys.getenv(
    "MORTALITY_SNAPSHOT_RDS",
    snapshot_path_join(get_snapshot_dir(), "mortality_ine_snapshot.rds")
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
    snapshot_path_join(get_snapshot_dir(), paste0(kind, ".rds"))
  )
}

get_snapshot_inventory_file <- function() {
  Sys.getenv(
    "MORTALITY_SNAPSHOT_INVENTORY_RDS",
    snapshot_path_join(get_snapshot_dir(), "snapshot_inventory.rds")
  )
}

is_absolute_path <- function(path) {
  is_url_path(path) | grepl("^(/|[A-Za-z]:[/\\\\])", as.character(path))
}

snapshot_inventory_full_path <- function(path) {
  path <- as.character(path)
  ifelse(is_absolute_path(path), path, snapshot_path_join(get_snapshot_dir(), path))
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

get_default_data_source <- function() {
  configured <- Sys.getenv("MORTALITY_DEFAULT_DATA_SOURCE", unset = "")
  if (nzchar(configured)) {
    return(normalize_data_source(configured))
  }

  snapshot_dir <- get_snapshot_dir()
  if (is_url_path(snapshot_dir) || snapshot_dir_has_local_data(snapshot_dir)) {
    return("snapshot")
  }

  "ine"
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
  snapshot_path_join(
    get_snapshot_dir(),
    "deaths",
    indicator,
    paste0("year_", year),
    paste0("cause_", snapshot_file_token(cause), ".rds")
  )
}

snapshot_path_exists <- function(path) {
  path <- as.character(path)
  ifelse(is_url_path(path), TRUE, file.exists(path))
}

read_snapshot_object <- memoise::memoise(function(path) {
  if (is_url_path(path)) {
    con <- url(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(readRDS(con))
  }

  readRDS(path)
})

read_snapshot_dataset <- memoise::memoise(function(kind, separate_path, combined_path) {
  if (snapshot_path_exists(separate_path)) {
    return(read_snapshot_object(separate_path))
  }

  if (snapshot_path_exists(combined_path)) {
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

ensure_source_indicator_column <- function(data, source_indicator) {
  if (!"source_indicator" %in% names(data)) {
    data$source_indicator <- source_indicator
  }

  data
}

read_snapshot_inventory_object <- memoise::memoise(function(path) {
  if (!snapshot_path_exists(path)) {
    return(NULL)
  }

  inventory <- tryCatch(read_snapshot_object(path), error = function(e) NULL)
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

format_snapshot_values <- function(values, max_values = 6) {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]

  if (length(values) == 0) {
    return("Nenhum")
  }

  shown <- utils::head(values, max_values)
  suffix <- if (length(values) > max_values) {
    paste0(" +", length(values) - max_values)
  } else {
    ""
  }

  paste0(paste(shown, collapse = ", "), suffix)
}

format_snapshot_year_span <- function(years) {
  years <- sort(unique(as.integer(years)))
  years <- years[!is.na(years)]

  if (length(years) == 0) {
    return("Nenhum")
  }

  if (length(years) == 1) {
    return(as.character(years))
  }

  paste0(min(years), " - ", max(years), " (", length(years), ")")
}

build_snapshot_inventory_summary <- function(inventory = get_snapshot_inventory()) {
  inventory_path <- get_snapshot_inventory_file()

  if (is.null(inventory)) {
    return(tibble::tibble(
      Item = c("Inventário", "Caminho"),
      Valor = c("Não encontrado", inventory_path)
    ))
  }

  population <- inventory %>%
    dplyr::filter(.data$dataset == "population")
  deaths <- inventory %>%
    dplyr::filter(.data$dataset == "deaths")
  inventory_mtime <- if (!is_url_path(inventory_path) && file.exists(inventory_path)) {
    format(file.info(inventory_path)$mtime, "%Y-%m-%d %H:%M")
  } else if (is_url_path(inventory_path)) {
    "GitHub"
  } else {
    "N/D"
  }
  max_area_count <- if (nrow(inventory) > 0) {
    max(vapply(inventory$areas, length, integer(1)), na.rm = TRUE)
  } else {
    NA_integer_
  }

  tibble::tibble(
    Item = c(
      "Inventário",
      "Actualizado",
      "Chunks totais",
      "Chunks de população",
      "Anos de população",
      "Chunks de óbitos",
      "Indicadores de óbitos",
      "Anos de óbitos",
      "Causas de óbito",
      "Máximo de áreas por chunk"
    ),
    Valor = c(
      inventory_path,
      inventory_mtime,
      as.character(nrow(inventory)),
      as.character(nrow(population)),
      format_snapshot_year_span(population$year),
      as.character(nrow(deaths)),
      format_snapshot_values(deaths$indicator),
      format_snapshot_year_span(deaths$year),
      as.character(dplyr::n_distinct(deaths$cause)),
      as.character(max_area_count)
    )
  )
}

get_snapshot_inventory_areas <- function(inventory_rows) {
  sort(unique(unlist(inventory_rows$areas, use.names = FALSE)))
}

get_relevant_snapshot_indicators <- function(inventory_rows, requested_areas) {
  if (nrow(inventory_rows) == 0) {
    return(character(0))
  }

  has_requested_area <- snapshot_inventory_has_area(inventory_rows$areas, requested_areas)
  sort(unique(inventory_rows$indicator[has_requested_area]))
}

summarize_snapshot_availability_row <- function(inventory_rows, requested_areas) {
  available_areas <- get_snapshot_inventory_areas(inventory_rows)
  present_areas <- intersect(requested_areas, available_areas)
  missing_areas <- setdiff(requested_areas, available_areas)

  status <- dplyr::case_when(
    length(requested_areas) == 0 ~ "Sem locais seleccionados",
    length(missing_areas) == 0 ~ "Disponível",
    length(present_areas) > 0 ~ "Parcial",
    TRUE ~ "Indisponível"
  )

  tibble::tibble(
    Estado = status,
    `Áreas pedidas` = format_snapshot_values(requested_areas, max_values = 4),
    `Áreas cobertas` = length(present_areas),
    `Áreas em falta` = format_snapshot_values(missing_areas, max_values = 4),
    `Fonte(s)` = format_snapshot_values(get_relevant_snapshot_indicators(inventory_rows, requested_areas)),
    Ficheiros = nrow(inventory_rows),
    Linhas = sum(inventory_rows$rows, na.rm = TRUE)
  )
}

build_population_snapshot_availability <- function(inventory, years, areas) {
  purrr::map_dfr(years, function(current_year) {
    inventory_rows <- inventory %>%
      dplyr::filter(
        .data$dataset == "population",
        .data$year == .env$current_year
      )

    summarize_snapshot_availability_row(inventory_rows, areas) %>%
      dplyr::mutate(
        Ano = current_year,
        Conjunto = "População",
        .before = 1
      )
  })
}

build_death_snapshot_availability <- function(inventory, years, areas, causes) {
  requested <- tidyr::expand_grid(
    Ano = years,
    `Causa de morte` = causes
  )

  purrr::pmap_dfr(requested, function(Ano, `Causa de morte`) {
    inventory_rows <- inventory %>%
      dplyr::filter(
        .data$dataset == "deaths",
        .data$year == .env$Ano,
        .data$cause == .env$`Causa de morte`
      ) %>%
      dplyr::arrange(dplyr::desc(.data$source_priority), .data$indicator)

    summarize_snapshot_availability_row(inventory_rows, areas) %>%
      dplyr::mutate(
        Ano = Ano,
        `Causa de morte` = `Causa de morte`,
        .before = 1
      )
  })
}

build_snapshot_availability_table <- function(
  dataset = "deaths",
  years = year_of_interest,
  areas = get_default_area_selection(),
  causes = NULL,
  show_missing = TRUE,
  inventory = get_snapshot_inventory()
) {
  if (is.null(inventory)) {
    return(tibble::tibble(
      Estado = "Inventário RDS não encontrado",
      Detalhe = get_snapshot_inventory_file()
    ))
  }

  dataset <- if (identical(dataset, "population")) "population" else "deaths"
  years <- sort(unique(as.integer(years)))
  years <- years[!is.na(years)]
  areas <- sort(unique(as.character(areas)))
  areas <- areas[!is.na(areas) & nzchar(areas)]

  out <- if (identical(dataset, "population")) {
    build_population_snapshot_availability(inventory, years, areas)
  } else {
    causes <- sort(unique(as.character(causes)))
    causes <- causes[!is.na(causes) & nzchar(causes)]

    if (length(causes) == 0) {
      return(tibble::tibble(
        Estado = "Sem causas seleccionadas",
        Detalhe = "Seleccione pelo menos uma causa de morte."
      ))
    }

    build_death_snapshot_availability(inventory, years, areas, causes)
  }

  if (!isTRUE(show_missing)) {
    out <- out %>%
      dplyr::filter(.data$Estado != "Indisponível")
  }

  out
}

format_snapshot_availability_issue <- function(row) {
  pieces <- character(0)
  if ("Ano" %in% names(row)) {
    pieces <- c(pieces, as.character(row[["Ano"]][[1]]))
  }
  if ("Causa de morte" %in% names(row)) {
    cause <- as.character(row[["Causa de morte"]][[1]])
    if (!is.na(cause) && nzchar(cause)) {
      pieces <- c(pieces, cause)
    }
  }

  issue_label <- if (length(pieces) > 0) paste(pieces, collapse = " / ") else "Selecção"
  issue <- paste0(issue_label, ": ", as.character(row[["Estado"]][[1]]))

  missing_areas <- if ("Áreas em falta" %in% names(row)) {
    as.character(row[["Áreas em falta"]][[1]])
  } else {
    ""
  }
  if (nzchar(missing_areas) && !identical(missing_areas, "Nenhum")) {
    issue <- paste0(issue, " (áreas em falta: ", missing_areas, ")")
  }

  sources <- if ("Fonte(s)" %in% names(row)) as.character(row[["Fonte(s)"]][[1]]) else ""
  if (nzchar(sources) && !identical(sources, "Nenhum")) {
    issue <- paste0(issue, "; fonte(s): ", sources)
  }

  issue
}

summarize_snapshot_availability_issues <- function(availability, label, max_examples = 5) {
  if (is.null(availability) || nrow(availability) == 0 || !"Estado" %in% names(availability)) {
    return(character(0))
  }

  issues <- availability %>%
    dplyr::filter(.data$Estado != "Disponível")

  if (nrow(issues) == 0) {
    return(character(0))
  }

  status_summary <- issues %>%
    dplyr::count(.data$Estado, name = "n") %>%
    dplyr::mutate(txt = paste0(.data$Estado, ": ", .data$n)) %>%
    dplyr::pull(.data$txt) %>%
    paste(collapse = "; ")

  shown <- utils::head(seq_len(nrow(issues)), max_examples)
  examples <- vapply(
    shown,
    function(i) format_snapshot_availability_issue(issues[i, , drop = FALSE]),
    character(1)
  )
  if (nrow(issues) > length(shown)) {
    examples <- c(examples, paste0("... e mais ", nrow(issues) - length(shown), " problema(s) de cobertura."))
  }

  c(
    glue::glue("{label}: cobertura RDS incompleta ({status_summary})."),
    paste0(" - ", examples)
  )
}

build_snapshot_request_warnings <- function(
  years,
  areas,
  causes = NULL,
  include_population = TRUE,
  include_deaths = TRUE,
  inventory = get_snapshot_inventory(),
  max_examples = 5
) {
  if (!isTRUE(include_population) && !isTRUE(include_deaths)) {
    return(character(0))
  }

  years <- sort(unique(as.integer(years)))
  years <- years[!is.na(years)]
  areas <- sort(unique(as.character(areas)))
  areas <- areas[!is.na(areas) & nzchar(areas)]

  if (is.null(inventory)) {
    return(glue::glue(
      "Aviso RDS: o inventário de snapshots não foi encontrado em {get_snapshot_inventory_file()}. A importação pode falhar ou recorrer à procura directa de ficheiros."
    ))
  }

  messages <- character(0)

  if (isTRUE(include_population)) {
    population_availability <- build_snapshot_availability_table(
      dataset = "population",
      years = years,
      areas = areas,
      show_missing = TRUE,
      inventory = inventory
    )
    messages <- c(
      messages,
      summarize_snapshot_availability_issues(
        population_availability,
        label = "População",
        max_examples = max_examples
      )
    )
  }

  if (isTRUE(include_deaths)) {
    causes <- sort(unique(as.character(causes)))
    causes <- causes[!is.na(causes) & nzchar(causes)]

    if (length(causes) > 0) {
      death_availability <- build_snapshot_availability_table(
        dataset = "deaths",
        years = years,
        areas = areas,
        causes = causes,
        show_missing = TRUE,
        inventory = inventory
      )
      messages <- c(
        messages,
        summarize_snapshot_availability_issues(
          death_availability,
          label = "Óbitos",
          max_examples = max_examples
        )
      )
    }
  }

  messages
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

summarize_snapshot_chunk <- function(path, dataset, indicator = NA_character_, snapshot_dir = get_snapshot_dir()) {
  chunk <- read_snapshot_object(path)
  cause <- if ("cause" %in% names(chunk)) sort(unique(as.character(chunk$cause))) else NA_character_

  tibble::tibble(
    dataset = dataset,
    indicator = indicator,
    year = if ("year" %in% names(chunk)) as.integer(sort(unique(chunk$year))[[1]]) else NA_integer_,
    cause = if (length(cause) > 0) cause[[1]] else NA_character_,
    path = snapshot_relative_path(path, snapshot_dir),
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
    indicator = "population",
    snapshot_dir = snapshot_dir
  )

  death_inventory <- purrr::map_dfr(death_paths, function(path) {
    relative <- snapshot_relative_path(path, snapshot_dir)
    parts <- strsplit(relative, .Platform$file.sep, fixed = TRUE)[[1]]
    indicator <- if (length(parts) >= 3) parts[[2]] else NA_character_
    summarize_snapshot_chunk(path, dataset = "deaths", indicator = indicator, snapshot_dir = snapshot_dir)
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
      missing_paths <- paths[!snapshot_path_exists(paths)]
      if (length(missing_paths) > 0) {
        stop(
          glue::glue("Missing population snapshot chunk file(s): {paste(missing_paths, collapse = ', ')}."),
          call. = FALSE
        )
      }

      return(collapse_population_snapshot_chunks(paths))
    }
  }

  paths <- snapshot_path_join(get_snapshot_dir(), "population", paste0("year_", years, ".rds"))

  if (!any(snapshot_path_exists(paths))) {
    return(NULL)
  }

  missing_paths <- paths[!snapshot_path_exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      glue::glue("Missing population snapshot chunk(s): {paste(missing_paths, collapse = ', ')}."),
      call. = FALSE
    )
  }

  collapse_population_snapshot_chunks(paths)
}

collapse_population_snapshot_chunks <- function(paths, source_indicator = "RDS population") {
  rows <- lapply(paths, function(path) {
    read_snapshot_object(path) %>%
      ensure_source_indicator_column(source_indicator)
  })

  dplyr::bind_rows(rows)
}

collapse_death_snapshot_chunks <- function(selected) {
  rows <- purrr::pmap(
    list(
      path = selected$path,
      source_priority = selected$source_priority,
      indicator = selected$indicator
    ),
    function(path, source_priority, indicator) {
      chunk <- read_snapshot_object(path)
      if (!"source_priority" %in% names(chunk)) {
        chunk$source_priority <- source_priority
      }
      ensure_source_indicator_column(chunk, indicator)
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
      exists = snapshot_path_exists(.data$path)
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
  if (
    is_url_path(get_snapshot_dir()) &&
      !nzchar(Sys.getenv("MORTALITY_POPULATION_SNAPSHOT_RDS", unset = "")) &&
      !nzchar(Sys.getenv("MORTALITY_SNAPSHOT_RDS", unset = ""))
  ) {
    return(read_chunked_population_snapshot(years, areas = areas))
  }

  tryCatch(
    get_snapshot_dataset("population"),
    error = function(e) read_chunked_population_snapshot(years, areas = areas)
  )
}

get_flat_or_chunked_death_snapshot <- function(years, causes, areas = NULL) {
  if (
    is_url_path(get_snapshot_dir()) &&
      !nzchar(Sys.getenv("MORTALITY_DEATHS_SNAPSHOT_RDS", unset = "")) &&
      !nzchar(Sys.getenv("MORTALITY_SNAPSHOT_RDS", unset = ""))
  ) {
    return(read_chunked_death_snapshot(years, causes, areas = areas))
  }

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
    ensure_source_indicator_column("RDS population") %>%
    dplyr::mutate(
      year = as.integer(year),
      area = as.character(area),
      sex = as.character(sex),
      age_band = as.character(age_band),
      pop = as.numeric(pop),
      source_indicator = as.character(source_indicator)
    ) %>%
    dplyr::filter(
      .data$year %in% .env$year_filter,
      .data$area %in% .env$area_filter
    ) %>%
    dplyr::group_by(year, area, sex, age_band, source_indicator) %>%
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
    ensure_source_indicator_column("RDS deaths") %>%
    dplyr::mutate(
      year = as.integer(year),
      area = as.character(area),
      sex = as.character(sex),
      cause = as.character(cause),
      age_band = as.character(age_band),
      deaths = as.numeric(deaths),
      source_indicator = as.character(source_indicator)
    ) %>%
    dplyr::filter(
      .data$year %in% .env$year_filter,
      .data$area %in% .env$area_filter,
      .data$cause %in% .env$cause_filter
    ) %>%
    dplyr::group_by(year, area, sex, cause, age_band, source_indicator) %>%
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
