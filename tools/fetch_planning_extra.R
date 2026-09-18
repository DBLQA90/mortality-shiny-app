#!/usr/bin/env Rscript
# Fetch the municipal components of the socio-economic, birth and infant
# indicators of the planning tab into data/snapshots/planning_extra.
#
#   Rscript tools/fetch_planning_extra.R [measures=all] [years=ALL] [overwrite=false] [recent=0]
#
# Every measure here is stored as additive municipal components - counts,
# tonnages, value totals - so that any ULS, ARS or region is an exact sum.
# Averages and rates are never stored: a mean pension is kept as pensioners and
# pensioners x mean, purchasing power per capita as the municipality's share of
# the national total and the population weight INE implies.
#
# Each measure is published in up to three INE editions, one per NUTS vintage.
# For every year the newest edition that covers it wins. Areas are matched by
# geography code, never by label: the last four digits of a municipal code are
# its DICO, identical in every vintage, and labels collide ("Lisboa" is also a
# NUTS-2002 region; "Calheta" and "Lagoa" each name two municipalities).
#
# Output, one file per measure and year:
#   data/snapshots/planning_extra/<measure>/year_<year>.rds
#   columns: year, area, category, value, source_indicator
# `category` is the published category of the measure's own dimension (a
# mother's age group, a pregnancy duration, a collection type), or "Total".

suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))

measures_arg <- get_arg("measures", "all")
years_arg <- get_arg("years", "ALL")
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")
# recent=N re-fetches the last N calendar years even when present, so a refresh
# picks up INE's revisions of provisional years; unchanged files stay untouched.
recent_years <- suppressWarnings(as.integer(get_arg("recent", "0")))
if (is.na(recent_years) || recent_years < 0) recent_years <- 0L
recheck_from <- as.integer(format(Sys.Date(), "%Y")) - recent_years

# ---------------------------------------------------------------------------
# Measures
# ---------------------------------------------------------------------------
# `editions`: oldest first. `target`: a regex on the category labels that
# identifies the measure's own dimension (NULL when the measure has none); every
# other dimension is pinned to its total. `value`: what the published number is.
MEASURES <- list(
  rsi_beneficiaries = list(
    editions = c("0004299", "0008251", "0013417"),
    target = NULL
  ),
  pensioners = list(
    # Série 1990-2023, then Série 2017 for the years only it reaches.
    editions = c("0004294", "0010271", "0013395", "0014534"),
    target = NULL
  ),
  pension_mean = list(
    editions = c("0004149", "0010266", "0013398", "0014532"),
    target = NULL
  ),
  waste_collected = list(
    # Tonnes by type of collection.
    editions = c("0000482", "0009612", "0012769"),
    target = "seletiva|selectiva|indiferenciada"
  ),
  purchasing_power_per_capita = list(
    editions = c("0001354", "0008614", "0014580"),
    target = NULL
  ),
  purchasing_power_share = list(
    editions = c("0001355", "0008615", "0014581"),
    target = NULL
  ),
  births_by_mother_age = list(
    editions = c("0005952", "0008092", "0012441"),
    target = "anos$|ignorad"
  ),
  births_by_gestation = list(
    editions = c("0005950", "0008084", "0012434"),
    target = "semanas|^Ignorad"
  ),
  births_by_weight = list(
    # Birth weight bands, for the low-birth-weight share (I36).
    editions = c("0005611", "0008088", "0012438"),
    # Bands read "2 000 - 2 499 g", "Menos de 500 g" and "5 000 g e mais".
    target = "[0-9] g|g e mais|^Ignorad"
  ),
  perinatal_deaths = list(
    # Deaths under 7 days plus stillbirths of 28 or more weeks (I43, I44). The
    # stillbirths are this minus the under-7-day deaths, which come from
    # infant_deaths_by_age.
    editions = c("0003527", "0008173", "0012549"),
    target = NULL
  ),
  # Census series, one value per census year (the app uses 1991 onwards).
  earnings_mean = list(
    # Average monthly earnings of employees (MTSSS/Quadros de Pessoal), by
    # workplace. Not additive: weighted by the number of employees below.
    editions = c("0009047", "0012656"),
    target = NULL
  ),
  employees_by_sector = list(
    editions = c("0010378", "0012648"),
    target = "Agricultura|Indústria|Serviços"
  ),
  census_population = list(
    editions = "0014353",
    target = NULL
  ),
  census_education = list(
    editions = "0014380",
    target = "Primário|Básico|Secundário|Superior"
  ),
  census_population_by_age = list(
    # For the population aged 10 and over, the denominator of the illiteracy rate.
    editions = "0014164",
    # "0 - 4 anos" … "75 e mais anos", plus "Ignorado".
    target = "anos$|^Ignorado$"
  ),
  census_illiteracy_rate = list(
    # INE's published municipal rate: illiterate aged 10+ over population 10+.
    editions = "0014375",
    target = NULL
  ),
  census_literacy = list(
    editions = "0014379",
    target = "Sabe ler e escrever|Não sabe ler e escrever"
  ),
  infant_deaths_by_age = list(
    editions = c("0008181", "0012541"),
    target = "dias|meses|Menos de 1 dia|hora"
  )
)

selected <- if (identical(measures_arg, "all")) names(MEASURES) else strsplit(measures_arg, ",")[[1]]
MEASURES <- MEASURES[intersect(names(MEASURES), selected)]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
dico_lookup <- readRDS("data/nuts_lookup_2024.rds") %>%
  transmute(dico = substr(as.character(municipality_code), 4, 7), municipality)
stopifnot(!anyDuplicated(dico_lookup$dico))

area_for_code <- function(code) {
  code <- as.character(code)
  out <- rep(NA_character_, length(code))
  out[code == "PT"] <- "Portugal"
  out[code == "1"] <- "Continente"
  # Municipal codes are the NUTS III code plus the DICO (7 characters), or, in
  # the historical census series, one character plus the DICO (5).
  municipal <- nchar(code) %in% c(5L, 7L)
  dico <- substr(code[municipal], nchar(code[municipal]) - 3L, nchar(code[municipal]))
  out[municipal] <- dico_lookup$municipality[match(dico, dico_lookup$dico)]
  out
}

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)

# INE answers bursts with 429 and then refuses connections for a while: pause
# after every request and back off hard on failure.
with_retries <- function(fn) {
  for (attempt in 1:5) {
    value <- tryCatch(fn(), error = function(e) NULL)
    if (!is.null(value) && NROW(value) > 0) {
      Sys.sleep(4)
      return(value)
    }
    Sys.sleep(60 * attempt)
  }
  NULL
}

metadata_cache <- list()
edition_years <- function(indicator) {
  if (is.null(metadata_cache[[indicator]])) {
    dv <- with_retries(function() client$get_dim_values(indicator))
    if (is.null(dv)) stop("Metadata unreachable for ", indicator, call. = FALSE)
    metadata_cache[[indicator]] <<- dv
  }
  dv <- metadata_cache[[indicator]]
  periods <- as.character(dv$categ_dsg[as.integer(dv$dim_num) == 1])
  sort(unique(suppressWarnings(as.integer(periods))))
}

parse_years <- function(value, default) {
  if (identical(toupper(value), "ALL")) return(default)
  if (grepl(":", value, fixed = TRUE)) {
    b <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(min(b), max(b)))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

total_labels <- c("Total", "HM", "MF")

# Reduce a raw response to one value per (code, category). The measure's own
# dimension is found by its labels and kept whole; every other dimension must
# carry a total category and is pinned to it, or the numbers would be summed
# across unrequested breakdowns.
tidy_response <- function(raw, target) {
  label_cols <- grep("^dim_[0-9]+_t$", names(raw), value = TRUE)
  label_cols <- setdiff(label_cols, "dim_1_t")
  # Birth indicators pad their labels with trailing spaces ("Total      ").
  # Not only ASCII spaces, so trim every Unicode space.
  for (col in label_cols) raw[[col]] <- gsub("^[\\s\\p{Z}]+|[\\s\\p{Z}]+$", "", as.character(raw[[col]]), perl = TRUE)

  target_col <- NULL
  if (!is.null(target)) {
    hits <- label_cols[vapply(label_cols, function(col) {
      values <- unique(as.character(raw[[col]]))
      mean(grepl(target, values, ignore.case = TRUE) | values %in% total_labels) == 1 &&
        any(grepl(target, values, ignore.case = TRUE))
    }, logical(1))]
    if (length(hits) != 1) stop("Cannot identify the measure dimension (", length(hits), " candidates)", call. = FALSE)
    target_col <- hits
  }

  keep <- rep(TRUE, nrow(raw))
  for (col in setdiff(label_cols, target_col)) {
    values <- as.character(raw[[col]])
    if (!any(values %in% total_labels)) {
      stop("Dimension ", col, " has no total category: ", paste(utils::head(unique(values), 5), collapse = " | "), call. = FALSE)
    }
    keep <- keep & values %in% total_labels
  }

  tibble(
    code = as.character(raw$geocod[keep]),
    category = if (is.null(target_col)) "Total" else as.character(raw[[target_col]][keep]),
    value = suppressWarnings(as.numeric(raw$valor[keep]))
  )
}

trim_label <- function(x) gsub("^[\\s\\p{Z}]+|[\\s\\p{Z}]+$", "", as.character(x), perl = TRUE)

total_pins <- function(dv, target) {
  pins <- list()
  for (dim in sort(unique(as.integer(dv$dim_num)))) {
    if (dim <= 2) next
    rows <- dv[as.integer(dv$dim_num) == dim, , drop = FALSE]
    labels <- trim_label(rows$categ_dsg)
    # The measure's own dimension: every label is one of its categories or a
    # total (the same rule tidy_response() uses to find it).
    is_target <- !is.null(target) && any(grepl(target, labels, ignore.case = TRUE)) &&
      all(grepl(target, labels, ignore.case = TRUE) | labels %in% total_labels)
    if (is_target) next
    code <- rows$categ_cod[labels %in% total_labels]
    if (length(code) >= 1) pins[[paste0("dim", dim)]] <- as.character(code[[1]])
  }
  pins
}

# Writes go through R/data_versions.R: identical content is left alone, a
# revised file is archived under data/archive before being replaced, and
# every write is recorded in data/import_log.csv with its date.
sys.source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]))), "..", "R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "fetch_planning_extra.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------
failed <- character(0)
for (measure in names(MEASURES)) {
  spec <- MEASURES[[measure]]

  # Newest edition covering a year wins.
  plan <- list()
  for (indicator in spec$editions) {
    years <- tryCatch(edition_years(indicator), error = function(e) {
      message("  ! ", conditionMessage(e))
      integer(0)
    })
    for (year in years) plan[[as.character(year)]] <- indicator
  }
  years <- intersect(as.integer(names(plan)), parse_years(years_arg, as.integer(names(plan))))
  if (length(years) == 0) next
  message("\n== ", measure, ": ", min(years), "-", max(years), " (", length(years), " years) ==")

  for (year in sort(years)) {
    path <- file.path("data/snapshots/planning_extra", measure, paste0("year_", year, ".rds"))
    if (file.exists(path) && !overwrite && year < recheck_from) next
    indicator <- plan[[as.character(year)]]

    # Ask only for the total of every dimension that is not the measure's own:
    # an unpinned birth indicator returns 835,000 cells a year, a pinned one
    # about 20,000. tidy_response() still checks the pins took.
    pins <- total_pins(metadata_cache[[indicator]], spec$target)
    raw <- with_retries(function() do.call(client$get_data, c(list(indicator, dim1 = paste0("S7A", year)), pins)))
    if (is.null(raw)) {
      message("  ", year, ": FAILED")
      failed <- c(failed, paste(measure, year))
      next
    }

    chunk <- tryCatch(tidy_response(raw, spec$target), error = function(e) {
      message("  ", year, ": ", indicator, " - ", conditionMessage(e))
      NULL
    })
    if (is.null(chunk)) {
      failed <- c(failed, paste(measure, year))
      next
    }

    chunk <- chunk %>%
      mutate(area = area_for_code(code)) %>%
      filter(!is.na(area)) %>%
      group_by(area, category) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      mutate(year = as.integer(year), source_indicator = indicator) %>%
      select(year, area, category, value, source_indicator) %>%
      arrange(area, category)

    n_mun <- n_distinct(setdiff(chunk$area, c("Portugal", "Continente")))
    national <- sum(chunk$value[chunk$area == "Portugal" & chunk$category %in% c("Total", total_labels)])
    message(sprintf("  %d %s: %d municipalities, %d categories, Portugal total %s",
                    year, indicator, n_mun, n_distinct(chunk$category), format(national, big.mark = " ")))
    if (n_mun < 300) {
      message("    REJECTED: only ", n_mun, " municipalities matched")
      failed <- c(failed, paste(measure, year))
      next
    }
    save_rds_atomic(chunk, path)
  }
}
message("\nDone. Failed: ", length(failed), if (length(failed)) paste0(" (", paste(utils::head(failed, 10), collapse = "; "), ")") else "")
