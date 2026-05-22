#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
})

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

parse_values <- function(value) {
  value <- value_or_default(value, "")
  values <- unlist(strsplit(value, "\\|", perl = TRUE), use.names = FALSE)
  values <- trimws(values)
  values[nzchar(values)]
}

fill_down <- function(x) {
  out <- as.character(x)
  last <- NA_character_
  for (i in seq_along(out)) {
    value <- out[[i]]
    if (!is.na(value) && nzchar(trimws(value))) {
      last <- value
    } else {
      out[[i]] <- last
    }
  }
  out
}

snapshot_file_token <- function(x) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

parse_death_value <- function(x) {
  x <- gsub("\u00a0", " ", as.character(x), fixed = TRUE)
  x <- trimws(x)
  x[x %in% c("", "x", "X", "-", "...", "*")] <- NA_character_
  suppressWarnings(readr::parse_number(
    x,
    locale = readr::locale(grouping_mark = " ", decimal_mark = ","),
    na = c("", "x", "X", "-", "...", "*")
  ))
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) {
    stop("Could not move temporary file into ", path, call. = FALSE)
  }
  invisible(path)
}

death_path <- function(out_dir, year, cause) {
  file.path(
    out_dir,
    "deaths",
    "0008206",
    paste0("year_", year),
    paste0("cause_", snapshot_file_token(cause), ".rds")
  )
}

read_ine_death_excel <- function(path) {
  raw <- readxl::read_excel(
    path,
    sheet = "Quadro",
    col_names = FALSE,
    .name_repair = "minimal",
    col_types = "text"
  )

  if (ncol(raw) < 7) {
    stop("Workbook does not have the expected INE Quadro layout: ", path, call. = FALSE)
  }

  names(raw) <- paste0("V", seq_len(ncol(raw)))

  parsed <- raw %>%
    transmute(
      year = suppressWarnings(as.integer(fill_down(.data$V1))),
      area = fill_down(.data$V2),
      sex = fill_down(.data$V4),
      cause = fill_down(.data$V5),
      age_band = as.character(.data$V6),
      deaths_raw = as.character(.data$V7),
      deaths = parse_death_value(.data$V7)
    ) %>%
    filter(
      !is.na(.data$year),
      .data$area == "Portugal",
      .data$sex %in% c("HM", "H", "M"),
      !is.na(.data$cause),
      nzchar(.data$cause),
      !is.na(.data$age_band),
      nzchar(.data$age_band)
    )

  year_value_summary <- parsed %>%
    group_by(.data$year) %>%
    summarise(
      numeric_rows = sum(!is.na(.data$deaths)),
      rows = n(),
      .groups = "drop"
    )

  populated_years <- year_value_summary %>%
    filter(.data$numeric_rows > 0) %>%
    pull(.data$year)

  data <- parsed %>%
    filter(.data$year %in% populated_years) %>%
    filter(!.data$age_band %in% c("Total", "Idade ignorada")) %>%
    mutate(
      age_band = case_when(
        .data$age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ .data$age_band
      )
    ) %>%
    group_by(.data$year, .data$area, .data$sex, .data$cause, .data$age_band) %>%
    summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop")

  list(
    path = path,
    data = data,
    year_value_summary = year_value_summary
  )
}

compare_existing_chunks <- function(data, out_dir) {
  comparisons <- data %>%
    distinct(.data$year, .data$cause) %>%
    mutate(path = death_path(out_dir, .data$year, .data$cause)) %>%
    filter(file.exists(.data$path))

  if (nrow(comparisons) == 0) {
    return(tibble())
  }

  bind_rows(lapply(seq_len(nrow(comparisons)), function(i) {
    year <- comparisons$year[[i]]
    cause <- comparisons$cause[[i]]
    path <- comparisons$path[[i]]

    existing <- readRDS(path) %>%
      filter(.data$year == .env$year, .data$area == "Portugal") %>%
      select("sex", "age_band", deaths_existing = "deaths")

    incoming <- data %>%
      filter(.data$year == .env$year, .data$cause == .env$cause) %>%
      select("sex", "age_band", deaths_incoming = "deaths")

    full_join(existing, incoming, by = c("sex", "age_band")) %>%
      mutate(
        year = year,
        cause = cause,
        path = path,
        deaths_existing = coalesce(.data$deaths_existing, NA_real_),
        deaths_incoming = coalesce(.data$deaths_incoming, NA_real_),
        match = coalesce(.data$deaths_existing == .data$deaths_incoming, FALSE)
      ) %>%
      filter(!.data$match)
  }))
}

write_death_chunks <- function(data, out_dir) {
  age_levels <- c(
    "0 - 4 anos", "5 - 9 anos", "10 - 14 anos", "15 - 19 anos",
    "20 - 24 anos", "25 - 29 anos", "30 - 34 anos", "35 - 39 anos",
    "40 - 44 anos", "45 - 49 anos", "50 - 54 anos", "55 - 59 anos",
    "60 - 64 anos", "65 - 69 anos", "70 - 74 anos", "75 - 79 anos",
    "80 - 84 anos", "85 e mais anos"
  )

  chunks <- data %>%
    group_by(.data$year, .data$cause) %>%
    group_split()

  written <- character()
  for (chunk in chunks) {
    year <- unique(chunk$year)
    cause <- unique(chunk$cause)
    path <- death_path(out_dir, year, cause)

    existing <- if (file.exists(path)) readRDS(path) else NULL
    if (!is.null(existing) && nrow(existing) > 0) {
      existing <- existing %>%
        filter(!(.data$year == .env$year & .data$area == "Portugal"))
      chunk <- bind_rows(existing, chunk)
    }

    chunk <- chunk %>%
      mutate(age_band = factor(.data$age_band, levels = age_levels, ordered = TRUE)) %>%
      arrange(.data$area, .data$sex, .data$age_band) %>%
      mutate(age_band = as.character(.data$age_band)) %>%
      select("year", "area", "sex", "cause", "age_band", "deaths")

    save_rds_atomic(chunk, path)
    written <- c(written, path)
  }

  unique(written)
}

cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
input_files <- parse_values(cli$files)
out_dir <- value_or_default(cli$out, file.path(getwd(), "data", "snapshots"))
write_output <- identical(tolower(value_or_default(cli$write, "1")), "1")
strict_existing <- identical(tolower(value_or_default(cli$strict_existing, "1")), "1")

if (length(input_files) == 0) {
  stop("Pass files=/path/a.xls|/path/b.xls", call. = FALSE)
}

imports <- lapply(input_files, read_ine_death_excel)
year_summaries <- bind_rows(lapply(imports, function(x) {
  x$year_value_summary %>%
    mutate(file = basename(x$path), .before = 1)
}))
data <- bind_rows(lapply(imports, `[[`, "data"))

blank_years <- year_summaries %>%
  filter(.data$numeric_rows == 0) %>%
  distinct(.data$file, .data$year)

coverage <- data %>%
  group_by(.data$year) %>%
  summarise(
    causes = n_distinct(.data$cause),
    rows = n(),
    areas = n_distinct(.data$area),
    sexes = paste(sort(unique(.data$sex)), collapse = "|"),
    .groups = "drop"
  )

message("Parsed INE Excel files for indicator 0008206 / Portugal")
message("  Files: ", paste(basename(input_files), collapse = ", "))
message("  Populated years: ", paste(sort(unique(data$year)), collapse = ", "))
if (nrow(blank_years) > 0) {
  message("  Blank year blocks skipped: ", paste(paste0(blank_years$year, " (", blank_years$file, ")"), collapse = ", "))
}
message("  Coverage by year:")
print(coverage, n = Inf)

invalid <- coverage %>%
  filter(.data$causes != 66L | .data$rows != 3564L | .data$areas != 1L | .data$sexes != "H|HM|M")
if (nrow(invalid) > 0) {
  stop("Parsed coverage is not complete for one or more populated years.", call. = FALSE)
}

mismatches <- compare_existing_chunks(data, out_dir)
if (nrow(mismatches) > 0) {
  message("Existing chunk comparison found ", nrow(mismatches), " mismatching Portugal rows.")
  print(head(mismatches, 20), n = 20)
  if (strict_existing) {
    stop("Refusing to write because existing chunks disagree with the Excel import.", call. = FALSE)
  }
}

if (write_output) {
  written <- write_death_chunks(data, out_dir)
  message("Wrote ", length(written), " chunk file(s).")
} else {
  message("Dry run only; no files written.")
}
