#!/usr/bin/env Rscript
# Fetch complete counts of deaths under 1 year by municipality and sex into
# data/snapshots/infant_totals.
#
#   Rscript tools/fetch_infant_totals.R [years=ALL] [overwrite=false]
#
# Why: data/snapshots/infant_deaths takes the "Menos de 1 ano" age band of the
# by-cause death indicators, and INE's municipal age breakdown of those is
# incomplete. In 2014 the municipalities add up to 112 infant deaths against a
# national 236, so every regional infant rate built from them was too low. INE
# also publishes under-1 deaths as a subject of their own, by municipality of
# residence, sex and age in days/months, whose totals are complete:
#
#   0008181  NUTS-2013  up to 2023
#   0012541  NUTS-2024  up to 2025, used from 2022 (same precedence as deaths)
#
# These carry no cause, so they serve the all-cause infant rate and count; the
# by-cause dataset stays for cause-specific requests.
#
# Areas are matched by geography code, never by label: "Calheta" and "Lagoa"
# each name two municipalities. The last four digits of a municipal code are its
# DICO, identical in every NUTS vintage.

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

years_arg <- get_arg("years", "ALL")
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")

dico_lookup <- readRDS("data/nuts_lookup_2024.rds") %>%
  transmute(dico = substr(as.character(municipality_code), 4, 7), municipality)
stopifnot(!anyDuplicated(dico_lookup$dico))

area_for_code <- function(code) {
  code <- as.character(code)
  out <- rep(NA_character_, length(code))
  out[code == "PT"] <- "Portugal"
  out[code == "1"] <- "Continente"
  municipal <- nchar(code) == 7
  out[municipal] <- dico_lookup$municipality[match(substr(code[municipal], 4, 7), dico_lookup$dico)]
  out
}

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)

with_retries <- function(expr_fn) {
  for (attempt in 1:5) {
    value <- tryCatch(expr_fn(), error = function(e) NULL)
    if (!is.null(value) && NROW(value) > 0) {
      Sys.sleep(5)
      return(value)
    }
    Sys.sleep(60 * attempt)
  }
  NULL
}

indicator_years <- function(indicator) {
  dv <- with_retries(function() client$get_dim_values(indicator))
  if (is.null(dv)) stop("Metadata unreachable for ", indicator, call. = FALSE)
  periods <- dv %>% filter(as.integer(dim_num) == 1)
  list(
    dims = dv,
    years = sort(unique(suppressWarnings(as.integer(as.character(periods$categ_dsg)))))
  )
}

# The age dimension's total category, found by label rather than assumed.
total_code <- function(dv, pattern) {
  hit <- dv %>%
    filter(as.integer(dim_num) > 2) %>%
    group_by(dim_num) %>%
    filter(any(grepl(pattern, categ_dsg, ignore.case = TRUE))) %>%
    ungroup()
  list(dim = unique(as.integer(hit$dim_num)), code = hit$categ_cod[hit$categ_dsg == "Total"])
}

parse_years <- function(value, default) {
  if (identical(toupper(value), "ALL")) return(default)
  if (grepl(":", value, fixed = TRUE)) {
    b <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(min(b), max(b)))
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

save_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Could not write ", path, call. = FALSE)
}

meta <- list("0008181" = indicator_years("0008181"), "0012541" = indicator_years("0012541"))
plan <- c(
  setNames(rep("0008181", length(setdiff(meta[["0008181"]]$years, 2022:2100))), setdiff(meta[["0008181"]]$years, 2022:2100)),
  setNames(rep("0012541", length(intersect(meta[["0012541"]]$years, 2022:2100))), intersect(meta[["0012541"]]$years, 2022:2100))
)
years <- intersect(as.integer(names(plan)), parse_years(years_arg, as.integer(names(plan))))
message("Years: ", min(years), "-", max(years), " (", length(years), ")")

written <- 0L
failed <- character(0)
for (year in sort(years)) {
  path <- file.path("data/snapshots/infant_totals", paste0("year_", year, ".rds"))
  if (file.exists(path) && !overwrite) next

  indicator <- plan[[as.character(year)]]
  dv <- meta[[indicator]]$dims
  age <- total_code(dv, "dias|meses|Menos de")
  if (length(age$dim) != 1 || length(age$code) != 1) stop("Cannot find the age total of ", indicator, call. = FALSE)

  dims <- list(indicator, dim1 = paste0("S7A", year))
  dims[[paste0("dim", age$dim)]] <- age$code
  raw <- with_retries(function() do.call(client$get_data, dims))
  if (is.null(raw)) {
    message("  ", year, ": FAILED")
    failed <- c(failed, as.character(year))
    next
  }

  sex_col <- grep("^dim_[0-9]+_t$", names(raw), value = TRUE)
  sex_col <- sex_col[vapply(sex_col, function(col) any(raw[[col]] %in% c("HM", "H", "M")), logical(1))]
  if (length(sex_col) != 1) stop("Cannot find the sex dimension of ", indicator, call. = FALSE)

  chunk <- tibble(
    code = as.character(raw$geocod),
    sex = as.character(raw[[sex_col]]),
    deaths = suppressWarnings(as.numeric(raw$valor))
  ) %>%
    mutate(area = area_for_code(code)) %>%
    filter(!is.na(area), sex %in% c("HM", "H", "M")) %>%
    group_by(area, sex) %>%
    summarise(deaths = sum(coalesce(deaths, 0)), .groups = "drop") %>%
    mutate(year = as.integer(year), source_indicator = indicator) %>%
    select(year, area, sex, deaths, source_indicator) %>%
    arrange(area, sex)

  n_mun <- n_distinct(setdiff(chunk$area, c("Portugal", "Continente")))
  hm <- chunk %>% filter(sex == "HM")
  national <- sum(hm$deaths[hm$area == "Portugal"])
  municipal <- sum(hm$deaths[!hm$area %in% c("Portugal", "Continente")])
  message(sprintf("  %s %d: %d municipalities | municipal sum %.0f, Portugal %.0f", indicator, year, n_mun, municipal, national))
  if (n_mun < 250 || national <= 0 || municipal > national) {
    message("    REJECTED: implausible coverage")
    failed <- c(failed, as.character(year))
    next
  }
  save_rds_atomic(chunk, path)
  written <- written + 1L
}
message("Done: ", written, " written, ", length(failed), " failed.")
