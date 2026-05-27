#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(xml2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("build_0008206_snapshot_from_portal.R", mustWork = FALSE)
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

parse_values <- function(value, default_values) {
  value <- value_or_default(value, "ALL")
  if (identical(toupper(value), "ALL")) {
    return(default_values)
  }
  values <- unlist(strsplit(value, "\\|", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  values[nzchar(values)]
}

parse_years <- function(value, default_years) {
  value <- value_or_default(value, "ALL")
  if (identical(toupper(value), "ALL")) {
    return(default_years)
  }

  tokens <- unlist(strsplit(value, "[,|]", perl = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]

  years <- unlist(lapply(tokens, function(token) {
    if (grepl(":", token, fixed = TRUE)) {
      bounds <- suppressWarnings(as.integer(strsplit(token, ":", fixed = TRUE)[[1]]))
      bounds <- bounds[!is.na(bounds)]
      if (length(bounds) == 2) {
        return(seq.int(min(bounds), max(bounds)))
      }
    }

    suppressWarnings(as.integer(token))
  }), use.names = FALSE)

  years[!is.na(years)]
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

normalize_key <- function(x) {
  x <- as.character(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("\\s+", " ", x))
}

prioritize_values <- function(values, preferred) {
  unique(c(preferred[preferred %in% values], values))
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

curl_run <- function(args, label) {
  retries <- as.integer(getOption("ine_curl_retries", env_or_default("CURL_RETRIES", "4")))
  retry_sleep <- as.numeric(getOption("ine_curl_retry_sleep", env_or_default("CURL_RETRY_SLEEP", "15")))
  retries <- ifelse(is.na(retries) || retries < 0L, 0L, retries)
  retry_sleep <- ifelse(is.na(retry_sleep) || retry_sleep < 0, 0, retry_sleep)
  max_attempts <- retries + 1L

  last_status <- NA_integer_
  for (attempt in seq_len(max_attempts)) {
    status <- system2("curl", args = shQuote(args))
    last_status <- status

    if (identical(status, 0L)) {
      return(invisible(TRUE))
    }

    if (attempt < max_attempts) {
      wait_seconds <- retry_sleep * attempt
      message(
        glue(
          "curl failed while {label} (attempt {attempt}/{max_attempts}, exit status {status}); ",
          "retrying in {wait_seconds}s."
        )
      )
      Sys.sleep(wait_seconds)
    }
  }

  stop(
    glue("curl failed while {label} after {max_attempts} attempts. Last exit status: {last_status}."),
    call. = FALSE
  )
}

fetch_portal_form <- function(cookie_jar, timeout_seconds = 120) {
  iframe_url <- "https://www.ine.pt/bddXplorer/htdocs/bddXplorer04.jsp?indOcorrCod=0008206&contexto=bd&userLoadSave=&lang=PT"
  referer <- "https://www.ine.pt/xportal/xmain?xpid=INE&xpgid=ine_indicadores&indOcorrCod=0008206&contexto=bd&selTab=tab2&xlang=pt"
  html_path <- tempfile(fileext = ".html")
  curl_run(
    c(
      "-L", "-sS",
      "-A", "Mozilla/5.0",
      "--compressed",
      "--max-time", as.character(timeout_seconds),
      "-c", cookie_jar,
      "-b", cookie_jar,
      "-e", referer,
      "-o", html_path,
      iframe_url
    ),
    "fetching the INE portal form"
  )

  html <- readChar(html_path, file.info(html_path)$size, useBytes = TRUE)
  unlink(html_path)

  doc <- xml2::read_html(html, encoding = "UTF-8")
  inputs <- xml2::xml_find_all(doc, ".//form[@id='frmIndicador']//input[@name]")
  names <- xml2::xml_attr(inputs, "name")
  values <- xml2::xml_attr(inputs, "value")
  values[is.na(values)] <- ""
  fields <- stats::setNames(as.list(values), names)
  fields <- fields[!startsWith(names(fields), "botao")]

  if (is.null(fields$xmlDocHidden) || !nzchar(fields$xmlDocHidden)) {
    stop("The INE portal form did not expose xmlDocHidden.", call. = FALSE)
  }

  fields
}

encode_form <- function(fields) {
  encode_piece <- function(x) {
    gsub("%20", "+", utils::URLencode(x, reserved = TRUE), fixed = TRUE)
  }
  paste0(
    encode_piece(names(fields)),
    "=",
    vapply(fields, encode_piece, character(1)),
    collapse = "&"
  )
}

submit_csv_request <- function(fields, cookie_jar, timeout_seconds = 180) {
  post_path <- tempfile(fileext = ".txt")
  html_path <- tempfile(fileext = ".html")
  on.exit(unlink(c(post_path, html_path)), add = TRUE)

  fields$botaoDownloadCsv <- "1"
  writeChar(encode_form(fields), post_path, useBytes = TRUE, eos = NULL)
  if (identical(Sys.getenv("DEBUG_PORTAL", unset = "0"), "1")) {
    file.copy(post_path, file.path(getwd(), "ine_portal_last_post.txt"), overwrite = TRUE)
  }

  curl_run(
    c(
      "-L", "-sS",
      "-A", "Mozilla/5.0",
      "--compressed",
      "--max-time", as.character(timeout_seconds),
      "-c", cookie_jar,
      "-b", cookie_jar,
      "-e", "https://www.ine.pt/bddXplorer/htdocs/bddXplorer04.jsp?indOcorrCod=0008206&contexto=bd&userLoadSave=&lang=PT",
      "-H", "Content-Type: application/x-www-form-urlencoded; charset=UTF-8",
      "--data-binary", paste0("@", post_path),
      "-o", html_path,
      "https://www.ine.pt/bddXplorer/htdocs/bddXplorer04.jsp"
    ),
    "submitting the INE CSV export form"
  )

  html <- readChar(html_path, file.info(html_path)$size, useBytes = TRUE)
  csv_links <- unlist(regmatches(html, gregexpr("/clientFiles/[^\"']+\\.csv", html, perl = TRUE)), use.names = FALSE)
  if (length(csv_links) == 0) {
    debug_path <- file.path(getwd(), "ine_portal_export_without_csv.html")
    writeChar(html, debug_path, useBytes = TRUE)
    stop(glue("The INE portal did not return a generated CSV link. Response saved to {debug_path}."), call. = FALSE)
  }

  paste0("https://www.ine.pt", csv_links[[1]])
}

download_generated_csv <- function(csv_url, cookie_jar, timeout_seconds = 180) {
  csv_path <- tempfile(fileext = ".csv")
  curl_run(
    c(
      "-L", "-sS",
      "-A", "Mozilla/5.0",
      "--compressed",
      "--max-time", as.character(timeout_seconds),
      "-b", cookie_jar,
      "-o", csv_path,
      csv_url
    ),
    "downloading the generated INE CSV"
  )
  csv_path
}

copy_failed_export <- function(csv_path, raw_dir, year_batch, batch_area_names, app, attempt, error) {
  if (!file.exists(csv_path)) {
    return(invisible(NULL))
  }

  failed_dir <- file.path(raw_dir, "failed")
  dir.create(failed_dir, recursive = TRUE, showWarnings = FALSE)
  failed_stem <- paste0(
    "failed_0008206_years_",
    paste(year_batch, collapse = "-"),
    "_areas_",
    app$snapshot_file_token(paste(batch_area_names, collapse = "_")),
    "_attempt_",
    attempt,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )

  failed_csv <- file.path(failed_dir, paste0(failed_stem, ".csv"))
  failed_txt <- file.path(failed_dir, paste0(failed_stem, ".txt"))
  file.copy(csv_path, failed_csv, overwrite = TRUE)
  writeLines(conditionMessage(error), failed_txt, useBytes = TRUE)

  message("Saved malformed INE CSV for inspection: ", failed_csv)
  invisible(failed_csv)
}

download_and_parse_export <- function(
  fields,
  cookie_jar,
  timeout_seconds,
  year_batch,
  batch_area_names,
  app,
  raw_dir,
  keep_raw = FALSE
) {
  export_retries <- as.integer(getOption("ine_export_retries", env_or_default("EXPORT_RETRIES", "3")))
  export_retry_sleep <- as.numeric(getOption("ine_export_retry_sleep", env_or_default("EXPORT_RETRY_SLEEP", "30")))
  export_retries <- ifelse(is.na(export_retries) || export_retries < 0L, 0L, export_retries)
  export_retry_sleep <- ifelse(is.na(export_retry_sleep) || export_retry_sleep < 0, 0, export_retry_sleep)
  max_attempts <- export_retries + 1L

  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    csv_path <- NULL
    parsed <- tryCatch({
      csv_url <- submit_csv_request(fields, cookie_jar, timeout_seconds = timeout_seconds)
      csv_path <- download_generated_csv(csv_url, cookie_jar, timeout_seconds = timeout_seconds)

      if (keep_raw) {
        dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
        raw_name <- paste0(
          "0008206_years_",
          paste(year_batch, collapse = "-"),
          "_areas_",
          app$snapshot_file_token(paste(batch_area_names, collapse = "_")),
          "_attempt_",
          attempt,
          ".csv"
        )
        file.copy(csv_path, file.path(raw_dir, raw_name), overwrite = TRUE)
      }

      parse_portal_csv(csv_path)
    }, error = function(e) e)

    if (!inherits(parsed, "error")) {
      if (!is.null(csv_path)) {
        unlink(csv_path)
      }
      return(parsed)
    }

    last_error <- parsed
    if (!is.null(csv_path)) {
      copy_failed_export(csv_path, raw_dir, year_batch, batch_area_names, app, attempt, parsed)
      unlink(csv_path)
    }

    if (attempt < max_attempts) {
      wait_seconds <- export_retry_sleep * attempt
      message(
        glue(
          "INE CSV export could not be parsed for years {paste(year_batch, collapse = ', ')} / ",
          "areas {paste(batch_area_names, collapse = ' | ')} ",
          "(attempt {attempt}/{max_attempts}): {conditionMessage(parsed)}; ",
          "retrying in {wait_seconds}s."
        )
      )
      Sys.sleep(wait_seconds)
    }
  }

  stop(
    glue(
      "INE CSV export could not be parsed after {max_attempts} attempts for years ",
      "{paste(year_batch, collapse = ', ')} / areas {paste(batch_area_names, collapse = ' | ')}: ",
      "{conditionMessage(last_error)}"
    ),
    call. = FALSE
  )
}

find_class_node <- function(root, order) {
  node <- xml2::xml_find_first(
    root,
    glue(".//*[starts-with(local-name(), 'CLASS_') and @order='{order}']")
  )
  if (inherits(node, "xml_missing")) {
    stop(glue("Could not find INE XML class with order {order}."), call. = FALSE)
  }
  node
}

replace_class_categories <- function(class_node, rows) {
  class_id <- xml2::xml_attr(class_node, "Id")
  categories_node <- xml2::xml_find_first(class_node, "./Categorias")
  if (inherits(categories_node, "xml_missing")) {
    stop(glue("Class {class_id} has no Categorias node."), call. = FALSE)
  }

  xml2::xml_remove(xml2::xml_children(categories_node))

  for (i in seq_len(nrow(rows))) {
    xml2::xml_add_child(
      categories_node,
      "Categoria",
      Id = paste0(class_id, "_", rows$cat_id[[i]]),
      cod_categ = rows$cat_id[[i]],
      nivel = as.character(rows$categ_nivel[[i]]),
      order = as.character(rows$categ_ord[[i]]),
      selected = "true",
      label = rows$categ_dsg[[i]]
    )
  }
}

set_existing_category_selection <- function(class_node, requested, dimension_label) {
  category_nodes <- xml2::xml_find_all(class_node, "./Categorias/Categoria")
  labels <- xml2::xml_attr(category_nodes, "label")
  codes <- xml2::xml_attr(category_nodes, "cod_categ")

  requested_key <- normalize_key(requested)
  selected <- normalize_key(labels) %in% requested_key | normalize_key(codes) %in% requested_key

  if (!any(selected)) {
    stop(glue("No {dimension_label} categories matched: {paste(requested, collapse = ', ')}"), call. = FALSE)
  }

  xml2::xml_set_attr(category_nodes, "selected", ifelse(selected, "true", "false"))
  tibble(label = labels[selected], code = codes[selected])
}

prepare_export_xml <- function(
  base_xml,
  years,
  areas,
  causes,
  area_metadata,
  age_metadata
) {
  root <- xml2::read_xml(base_xml)

  cause_class <- find_class_node(root, order = 1)
  period_class <- find_class_node(root, order = 4)
  age_class <- find_class_node(root, order = 3)
  geo_class <- find_class_node(root, order = 5)

  selected_causes <- set_existing_category_selection(cause_class, causes, "cause")
  set_existing_category_selection(period_class, as.character(years), "year")

  replace_class_categories(age_class, age_metadata)
  replace_class_categories(geo_class, area_metadata)

  xml_text <- as.character(root)
  xml_text <- sub("^<\\?xml[^>]+>\\s*", "", xml_text)
  xml_text <- gsub(">\\s+<", "><", xml_text)
  xml_text <- trimws(xml_text)

  list(
    xml = xml_text,
    causes = selected_causes$label
  )
}

prepare_metadata <- function(app, areas) {
  dv <- app$get_dim_values_cached(app$death_indicator_legacy)
  if (!"cat_id" %in% names(dv)) {
    dv$cat_id <- dv$categ_cod
  }
  if (!"categ_ord" %in% names(dv)) {
    dv$categ_ord <- seq_len(nrow(dv))
  }
  if (!"categ_nivel" %in% names(dv)) {
    dv$categ_nivel <- NA_character_
  }

  dv <- dv %>%
    mutate(
      dim_num = as.integer(.data$dim_num),
      cat_id = as.character(.data$cat_id),
      categ_dsg = as.character(.data$categ_dsg),
      categ_ord = as.character(.data$categ_ord),
      categ_nivel = as.character(.data$categ_nivel)
    )

  age_rows <- dv %>%
    filter(.data$dim_num == 4L) %>%
    filter(!.data$categ_dsg %in% c("Total", "Idade ignorada")) %>%
    arrange(suppressWarnings(as.numeric(.data$categ_ord)))

  area_rows_all <- dv %>%
    filter(.data$dim_num == 2L)

  area_key <- normalize_key(areas)
  area_rows <- map_dfr(seq_along(areas), function(i) {
    area <- areas[[i]]
    match <- area_rows_all %>%
      filter(
        normalize_key(.data$categ_dsg) == area_key[[i]] |
          normalize_key(.data$cat_id) == area_key[[i]]
      ) %>%
      slice(1)

    if (nrow(match) == 0) {
      stop(glue("Could not map area '{area}' to indicator 0008206 metadata."), call. = FALSE)
    }

    match
  })

  list(age = age_rows, areas = area_rows)
}

fill_down <- function(x) {
  out <- as.character(x)
  last <- ""
  for (i in seq_along(out)) {
    if (!is.na(out[[i]]) && nzchar(trimws(out[[i]]))) {
      last <- out[[i]]
    } else {
      out[[i]] <- last
    }
  }
  out
}

fill_right <- function(x) {
  out <- as.character(x)
  last <- ""
  for (i in seq_along(out)) {
    if (!is.na(out[[i]]) && nzchar(trimws(out[[i]]))) {
      last <- out[[i]]
    } else {
      out[[i]] <- last
    }
  }
  out
}

clean_cell <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

parse_death_value <- function(x) {
  x <- gsub("\u00a0", " ", as.character(x), fixed = TRUE)
  x <- trimws(x)
  x[grepl("^[-[:space:]]+$", x)] <- NA_character_
  suppressWarnings(readr::parse_number(
    x,
    locale = readr::locale(grouping_mark = " ", decimal_mark = ","),
    na = c("", "x", "X", "-", "...")
  ))
}

parse_portal_csv <- function(csv_path) {
  raw <- utils::read.table(
    csv_path,
    sep = ";",
    header = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    fileEncoding = "latin1",
    colClasses = "character"
  )

  raw <- tibble::as_tibble(raw, .name_repair = "minimal")

  raw <- raw %>% mutate(across(everything(), clean_cell))

  non_empty_cols <- vapply(raw, function(col) any(nzchar(col)), logical(1))
  raw <- raw[, non_empty_cols, drop = FALSE]

  if (ncol(raw) < 4) {
    stop("The INE CSV did not contain value columns.", call. = FALSE)
  }

  value_cols <- seq.int(4L, ncol(raw))
  year_row <- which(apply(raw[, value_cols, drop = FALSE], 1, function(x) any(grepl("^\\d{4}$", x))))[[1]]
  area_row <- which(seq_len(nrow(raw)) > year_row & apply(raw[, value_cols, drop = FALSE], 1, function(x) any(grepl("^.+:\\s*.+$", x))))[[1]]

  if (is.na(year_row) || is.na(area_row)) {
    stop("Could not locate year/area header rows in the INE CSV.", call. = FALSE)
  }

  years_by_col <- fill_right(unlist(raw[year_row, value_cols], use.names = FALSE))
  areas_by_col <- unlist(raw[area_row, value_cols], use.names = FALSE)
  areas_by_col <- trimws(sub("^.*?:\\s*", "", areas_by_col))

  data_start <- area_row + 2L
  data <- raw[data_start:nrow(raw), , drop = FALSE]
  first_col <- data[[1]]
  end_at <- which(
    grepl("^Óbitos \\(", first_col) |
      grepl("^Nota\\(s\\):", first_col) |
      grepl("^Fonte:", first_col) |
      grepl("^Última atualização", first_col)
  )
  if (length(end_at) > 0) {
    data <- data[seq_len(end_at[[1]] - 1L), , drop = FALSE]
  }

  value_col_names_current <- names(data)[value_cols]
  data <- data %>%
    filter(if_any(all_of(value_col_names_current), ~ nzchar(.x)))

  if (nrow(data) == 0) {
    stop("The INE CSV contained no data rows.", call. = FALSE)
  }

  data[[1]] <- fill_down(data[[1]])
  data[[2]] <- fill_down(data[[2]])

  names(data)[1:3] <- c("cause", "sex", "age_band")
  value_names <- paste0("value_", seq_along(value_cols))
  names(data)[value_cols] <- value_names

  column_map <- tibble(
    value_column = value_names,
    year = suppressWarnings(as.integer(years_by_col)),
    area = areas_by_col
  )

  data %>%
    select(cause, sex, age_band, all_of(value_names)) %>%
    pivot_longer(
      cols = all_of(value_names),
      names_to = "value_column",
      values_to = "deaths_raw"
    ) %>%
    left_join(column_map, by = "value_column") %>%
    mutate(
      deaths = parse_death_value(.data$deaths_raw),
      sex = recode(.data$sex, "MF" = "HM", .default = .data$sex),
      age_band = case_when(
        .data$age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ .data$age_band
      )
    ) %>%
    filter(
      !is.na(.data$year),
      nzchar(.data$area),
      !.data$age_band %in% c("Total", "Idade ignorada")
    ) %>%
    group_by(year, area, sex, cause, age_band) %>%
    summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop")
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) {
    stop(glue("Could not move temporary file into {path}."), call. = FALSE)
  }
}

death_path <- function(out_dir, year, cause, app) {
  file.path(
    out_dir,
    "deaths",
    "0008206",
    paste0("year_", year),
    paste0("cause_", app$snapshot_file_token(cause), ".rds")
  )
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

batch_needs_build <- function(out_dir, years, areas, causes, app) {
  for (year in years) {
    for (cause in causes) {
      path <- death_path(out_dir, year, cause, app)
      if (length(missing_areas_for_path(path, areas)) > 0) {
        return(TRUE)
      }
    }
  }
  FALSE
}

save_death_chunks <- function(data, out_dir, app) {
  if (nrow(data) == 0) {
    return(0L)
  }

  age_order <- if (exists("age_levels", envir = app, inherits = FALSE)) app$age_levels else unique(data$age_band)
  chunks <- data %>%
    group_by(year, cause) %>%
    group_split()

  written <- 0L
  for (chunk in chunks) {
    year <- unique(chunk$year)
    cause <- unique(chunk$cause)
    areas <- unique(chunk$area)
    path <- death_path(out_dir, year, cause, app)

    existing <- if (file.exists(path)) readRDS(path) else NULL
    if (!is.null(existing) && nrow(existing) > 0) {
      existing <- existing %>%
        filter(!(.data$year %in% year & .data$cause %in% cause & .data$area %in% areas))
      chunk <- bind_rows(existing, chunk)
    }

    chunk <- chunk %>%
      mutate(age_band = factor(.data$age_band, levels = age_order, ordered = TRUE)) %>%
      arrange(.data$area, .data$sex, .data$age_band) %>%
      mutate(age_band = as.character(.data$age_band)) %>%
      select(year, area, sex, cause, age_band, deaths)

    save_rds_atomic(chunk, path)
    written <- written + 1L
  }

  written
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
app_file <- value_or_default(cli$app, find_default_app_file())
out_dir <- value_or_default(cli$out, file.path(dirname(app_file), "data", "snapshots"))
area_batch_size <- as.integer(value_or_default(cli$area_batch_size, env_or_default("AREA_BATCH_SIZE", "12")))
year_batch_size <- as.integer(value_or_default(cli$year_batch_size, env_or_default("YEAR_BATCH_SIZE", "1")))
max_batches <- parse_limit(value_or_default(cli$max_batches, env_or_default("MAX_BATCHES", "Inf")))
timeout_seconds <- as.integer(value_or_default(cli$timeout, env_or_default("TIMEOUT", "240")))
curl_retries <- as.integer(value_or_default(cli$curl_retries, env_or_default("CURL_RETRIES", "4")))
curl_retry_sleep <- as.numeric(value_or_default(cli$curl_retry_sleep, env_or_default("CURL_RETRY_SLEEP", "15")))
export_retries <- as.integer(value_or_default(cli$export_retries, env_or_default("EXPORT_RETRIES", "3")))
export_retry_sleep <- as.numeric(value_or_default(cli$export_retry_sleep, env_or_default("EXPORT_RETRY_SLEEP", "30")))
keep_raw <- identical(tolower(value_or_default(cli$keep_raw, env_or_default("KEEP_RAW", "0"))), "1")
raw_dir <- value_or_default(cli$raw_dir, file.path(out_dir, "raw", "0008206_portal"))

options(
  ine_curl_retries = ifelse(is.na(curl_retries) || curl_retries < 0L, 0L, curl_retries),
  ine_curl_retry_sleep = ifelse(is.na(curl_retry_sleep) || curl_retry_sleep < 0, 0, curl_retry_sleep),
  ine_export_retries = ifelse(is.na(export_retries) || export_retries < 0L, 0L, export_retries),
  ine_export_retry_sleep = ifelse(is.na(export_retry_sleep) || export_retry_sleep < 0, 0, export_retry_sleep)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

app <- new.env(parent = globalenv())
source(app_file, local = app)

cookie_jar <- tempfile(fileext = ".cookies")
on.exit(unlink(cookie_jar), add = TRUE)

fields <- fetch_portal_form(cookie_jar, timeout_seconds = timeout_seconds)
base_xml <- fields$xmlDocHidden
base_root <- xml2::read_xml(base_xml)
portal_causes <- xml2::xml_attr(
  xml2::xml_find_all(find_class_node(base_root, order = 1), "./Categorias/Categoria"),
  "label"
)
portal_years <- suppressWarnings(as.integer(xml2::xml_attr(
  xml2::xml_find_all(find_class_node(base_root, order = 4), "./Categorias/Categoria"),
  "label"
)))
portal_years <- portal_years[!is.na(portal_years)]

legacy_years <- sort(portal_years)
default_year <- max(legacy_years, na.rm = TRUE)

years <- parse_years(value_or_default(cli$years, env_or_default("YEARS", as.character(default_year))), legacy_years)
areas <- parse_values(value_or_default(cli$areas, env_or_default("AREAS", "Portugal|Norte")), app$local_area)
causes <- parse_values(value_or_default(cli$causes, env_or_default("CAUSES", "ALL")), portal_causes)
priority_areas <- parse_values(value_or_default(cli$priority_areas, env_or_default("PRIORITY_AREAS", "Portugal|Norte")), character())

years <- sort(unique(as.integer(years)), decreasing = TRUE)
years <- years[years %in% legacy_years]
areas <- unique(as.character(areas))
priority_areas <- unique(as.character(priority_areas))
areas <- prioritize_values(areas, priority_areas)
causes <- unique(as.character(causes))

if (length(years) == 0 || length(areas) == 0 || length(causes) == 0) {
  stop("Portal snapshot export requires at least one valid year, area, and cause.", call. = FALSE)
}

metadata <- prepare_metadata(app, areas)
area_rows <- metadata$areas
age_rows <- metadata$age
areas <- area_rows$categ_dsg

year_batches <- chunk_vector(years, year_batch_size)
area_batches <- chunk_vector(seq_len(nrow(area_rows)), area_batch_size)

message("INE portal snapshot export")
message("  Years: ", paste(years, collapse = ", "))
message("  Areas: ", paste(areas, collapse = " | "))
message("  Causes: ", if (identical(sort(causes), sort(portal_causes))) "ALL" else paste(causes, collapse = " | "))
message("  Area batch size: ", area_batch_size)
message("  Curl retries: ", getOption("ine_curl_retries"), " (sleep base ", getOption("ine_curl_retry_sleep"), "s)")
message("  Export parse retries: ", getOption("ine_export_retries"), " (sleep base ", getOption("ine_export_retry_sleep"), "s)")

batches_done <- 0L
chunks_written <- 0L

for (year_batch in year_batches) {
  for (area_batch_i in area_batches) {
    batch_areas <- area_rows[area_batch_i, , drop = FALSE]
    batch_area_names <- batch_areas$categ_dsg

    export_spec <- prepare_export_xml(
      base_xml = base_xml,
      years = year_batch,
      areas = batch_area_names,
      causes = causes,
      area_metadata = batch_areas,
      age_metadata = age_rows
    )

    if (!batch_needs_build(out_dir, year_batch, batch_area_names, export_spec$causes, app)) {
      message("Skipping covered batch: years ", paste(year_batch, collapse = ", "), " / areas ", paste(batch_area_names, collapse = " | "))
      next
    }

    if (batches_done >= max_batches) {
      break
    }

    message("Exporting years ", paste(year_batch, collapse = ", "), " / areas ", paste(batch_area_names, collapse = " | "))
    fields$xmlDocHidden <- export_spec$xml
    fields$lingua_cd <- "PT"

    parsed <- download_and_parse_export(
      fields = fields,
      cookie_jar = cookie_jar,
      timeout_seconds = timeout_seconds,
      year_batch = year_batch,
      batch_area_names = batch_area_names,
      app = app,
      raw_dir = raw_dir,
      keep_raw = keep_raw
    )
    chunks_written <- chunks_written + save_death_chunks(parsed, out_dir, app)
    batches_done <- batches_done + 1L
  }

  if (batches_done >= max_batches) {
    break
  }
}

manifest <- list(
  created_at = Sys.time(),
  source = "INE portal CSV export",
  app_file = app_file,
  years_requested = years,
  areas_requested = areas,
  causes_requested = causes,
  area_batch_size = area_batch_size,
  year_batch_size = year_batch_size,
  batches_built_this_run = batches_done,
  chunks_written_this_run = chunks_written
)
save_rds_atomic(manifest, file.path(out_dir, "chunked_0008206_manifest.rds"))

message("Done. Batches built: ", batches_done, ". RDS chunks written: ", chunks_written, ".")
