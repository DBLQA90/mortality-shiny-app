#!/usr/bin/env Rscript
# Fetch INE's own regional death rows into data/snapshots/regional_deaths.
#
#   Rscript tools/fetch_regional_deaths.R [indicator=all] [years=ALL] [overwrite=false]
#
# Why these rows are needed
# -------------------------
# The death archive stores municipalities, and the app builds regions by summing
# them. For cause-specific deaths that sum is systematically short: INE publishes
# complete municipal totals but incomplete municipal age breakdowns, and the
# missing detail is concentrated where counts are small. Lung cancer, share of a
# region's age-banded deaths lost by summing its municipalities:
#
#            Norte  Centro  Alentejo  Açores  Madeira
#   2013     -1.2%  -2.1%    -8.9%   -18.4%   -11.8%
#   2014    -30.7% -50.7%   -64.6%   -71.4%   -83.5%
#   2021      0.0%  -0.2%    -3.4%   -11.6%    -6.0%
#
# INE's regional rows do not have the problem: their age bands add up to their
# total exactly, 2014 included. The app therefore prefers them wherever a row
# exists for the selected region's territory (see R/regional_rows.R).
#
# Why a separate dataset, keyed by code
# -------------------------------------
# Four labels repeat across NUTS levels in these indicators - Algarve, the
# Lisbon metropolitan area and both autonomous regions are named identically at
# NUTS II and III (and the islands at NUTS I too). The municipal archive stores
# rows by label, so those regional rows arrive double- or triple-counted there.
# Here every row is selected by its geography code, one level per territory, and
# written beside the municipal archive rather than into it, so nothing already
# verified is touched.
#
# One request per year returns every region, cause, sex and age band.

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

indicator_arg <- get_arg("indicator", "all")
years_arg <- get_arg("years", "ALL")
out_dir <- get_arg("out", "data/snapshots")
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")

# One code per territory. Continente is NUTS I; the islands are taken at NUTS II,
# which is the same territory as their NUTS I row.
#
# NUTS III subregions are fetched too, because they are what lets a redrawn
# region be composed from INE rows in the years its own row does not exist:
# Oeste, Médio Tejo and Lezíria do Tejo were subregions under both definitions,
# so NUTS-2024 Alentejo before 2022 is four NUTS-2013 subregions, and NUTS-2013
# Centro after 2022 is NUTS-2024 Centro plus Oeste and Médio Tejo. The
# compositions are declared in R/regional_rows.R and verified against the
# municipality lookups.
sources <- list(
  "0008206" = list(
    years = 1991:2022,
    codes = c(
      "1", "11", "16", "17", "18", "15", "20", "30",
      # Centro subregions
      "16B", "16D", "16E", "16F", "16G", "16H", "16I", "16J",
      # Alentejo subregions
      "181", "184", "185", "186", "187",
      # Alto Minho, which is also ULS Alto Minho
      "111"
    )
  ),
  "0013166" = list(
    years = 2022:2024,
    codes = c(
      "1", "11", "19", "1A", "1B", "1C", "1D", "15", "20", "30",
      # Oeste e Vale do Tejo subregions
      "1D1", "1D2", "1D3",
      # Subregions that coincide exactly with a ULS: Alto Minho, Viseu Dão
      # Lafões, and the four Alentejo subregions.
      "111", "194", "1C1", "1C2", "1C3", "1C4"
    )
  )
)

if (!identical(indicator_arg, "all")) {
  sources <- sources[intersect(names(sources), strsplit(indicator_arg, ",")[[1]])]
}

parse_years <- function(value, default) {
  if (identical(toupper(value), "ALL")) return(default)
  if (grepl(":", value, fixed = TRUE)) {
    b <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(min(b), max(b)))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)

fetch <- function(indicator, year, codes) {
  for (attempt in 1:5) {
    result <- tryCatch(
      client$get_data(indicator, dim1 = paste0("S7A", year), dim2 = codes),
      error = function(e) NULL
    )
    if (!is.null(result) && nrow(result) > 0) return(result)
  }
  NULL
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Could not move temporary file into ", path, call. = FALSE)
}

# The response names dimensions by position, and the positions are not assumed:
# each role is found by the labels it carries.
find_dim <- function(raw, pattern) {
  candidates <- grep("^dim_[0-9]+_t$", names(raw), value = TRUE)
  hit <- candidates[vapply(candidates, function(col) any(grepl(pattern, raw[[col]])), logical(1))]
  if (length(hit) != 1) stop("Could not identify dimension matching ", pattern, call. = FALSE)
  hit
}

failed <- character(0)
written <- 0L

for (indicator in names(sources)) {
  spec <- sources[[indicator]]
  years <- intersect(parse_years(years_arg, spec$years), spec$years)
  message("== ", indicator, ": ", length(years), " years, ", length(spec$codes), " territories ==")

  for (year in years) {
    path <- file.path(out_dir, "regional_deaths", indicator, paste0("year_", year, ".rds"))
    if (file.exists(path) && !overwrite) {
      next
    }

    raw <- fetch(indicator, year, spec$codes)
    if (is.null(raw)) {
      message("  ", year, ": FAILED")
      failed <- c(failed, paste(indicator, year))
      next
    }

    sex_col <- find_dim(raw, "^HM$")
    age_col <- find_dim(raw, "anos")
    cause_col <- find_dim(raw, "^Todas as causas de morte$")

    tidy <- raw %>%
      dplyr::transmute(
        year = as.integer(year),
        region_code = as.character(geocod),
        area = as.character(geodsg),
        sex = as.character(.data[[sex_col]]),
        cause = as.character(.data[[cause_col]]),
        age_raw = as.character(.data[[age_col]]),
        deaths = suppressWarnings(as.numeric(valor))
      ) %>%
      dplyr::filter(region_code %in% spec$codes)

    # The whole point of these rows is that their age bands are complete, so
    # check it rather than assume it. A mismatch here would mean the regional
    # rows carry the same defect as the municipal ones.
    check <- tidy %>%
      dplyr::group_by(region_code, sex, cause) %>%
      dplyr::summarise(
        total = sum(deaths[age_raw == "Total"], na.rm = TRUE),
        bands = sum(deaths[age_raw != "Total"], na.rm = TRUE),
        .groups = "drop"
      )
    incomplete <- check %>% dplyr::filter(abs(total - bands) > 0.5)
    if (nrow(incomplete) > 0) {
      message("  ", year, ": ", nrow(incomplete), " region/sex/cause cells where age bands do not sum to the total")
    }

    chunk <- tidy %>%
      dplyr::filter(!age_raw %in% c("Total", "Idade ignorada")) %>%
      dplyr::mutate(age_band = dplyr::if_else(
        age_raw %in% c("Menos de 1 ano", "1 - 4 anos"), "0 - 4 anos", age_raw
      )) %>%
      dplyr::group_by(year, region_code, area, sex, cause, age_band) %>%
      dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(source_indicator = indicator) %>%
      dplyr::arrange(region_code, cause, sex, age_band)

    # NUTS II only: subregions would count their region twice.
    national_regions <- spec$codes[spec$codes != "1" & nchar(spec$codes) == 2]
    allcause <- chunk %>% dplyr::filter(sex == "HM", cause == "Todas as causas de morte")
    save_rds_atomic(chunk, path)
    written <- written + 1L
    message(sprintf(
      "  %d: %d rows | all-cause HM: Continente %s, NUTS II sum %s | incomplete cells %d",
      year, nrow(chunk),
      format(sum(allcause$deaths[allcause$region_code == "1"]), big.mark = " "),
      format(sum(allcause$deaths[allcause$region_code %in% national_regions]), big.mark = " "),
      nrow(incomplete)
    ))
  }
}

message("Done: ", written, " written, ", length(failed), " failed.")
if (length(failed) > 0) message("Re-run to retry: ", paste(failed, collapse = "; "))
