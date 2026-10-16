#!/usr/bin/env Rscript
# Build calendar-year deaths by five-year age group per NUTS II 2024 region,
# for years the by-cause archive does not reach yet, into
# data/snapshots/calendar_deaths_by_age.
#
#   Rscript tools/build_calendar_deaths_by_age.R [years=AUTO] [validate=false]
#
# Why: deaths by age by region arrive with the causes of death (0013166), in
# March of year n+2. INE's annual deaths for year n+1 (0013331) come earlier but
# only in five broad age groups (<1, 1-14, 15-24, 25-64, 65+), too coarse for
# the European standard population or a life table. Two datasets the app
# already holds give the rest, both observed counts:
#
#   weekly_deaths/0012100   deaths by ISO week, NUTS 2024 region and five-year
#                           age group (to 90+): the age profile
#   parish_vitals/deaths    deaths by parish of residence and calendar year
#                           (0012542): the calendar-year total per region
#
# ISO weeks straddle the new year, so each week is split by the number of its
# seven days that fall in the calendar year. The age profile is then scaled to
# the region's parish total. Portugal keeps the weekly total, which, like INE's
# national row, includes deaths of unknown or foreign residence; Continente is
# the sum of its five regions. Validated on 2022-2024 against the by-cause rows:
# totals within 0.13%, the standardised rate within 0.03%, under 75 within 0.3%
# (validate=true prints the rebuilt totals of the archived years).
#
# Output: data/snapshots/calendar_deaths_by_age/year_<year>.rds
#   columns: year, area, code, sex ("HM" only), age_band, deaths, source_indicator

suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "build_calendar_deaths_by_age.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

REGIONS <- c("Norte" = "11", "Centro" = "19", "Grande Lisboa" = "1A", "Península de Setúbal" = "1B",
             "Oeste e Vale do Tejo" = "1D", "Alentejo" = "1C", "Algarve" = "15",
             "Região Autónoma dos Açores" = "20", "Região Autónoma da Madeira" = "30")
MAINLAND <- c("Norte", "Centro", "Grande Lisboa", "Península de Setúbal", "Oeste e Vale do Tejo", "Alentejo", "Algarve")
BANDS <- c(paste0(seq(0, 80, 5), " - ", seq(4, 84, 5), " anos"), "85 e mais anos")

band_of <- function(label) {
  a <- sub(" anos$", "", trimws(label))
  out <- rep(NA_character_, length(a))
  five <- grepl("^\\d+-\\d+$", a)
  lower <- as.integer(sub("-.*$", "", a[five]))
  out[five] <- ifelse(lower >= 85, "85 e mais anos", paste0(lower, " - ", lower + 4, " anos"))
  out[a == "90 ou mais"] <- "85 e mais anos"
  out
}

iso_monday <- function(year, week) {
  jan4 <- as.Date(sprintf("%d-01-04", year))
  jan4 - (as.integer(format(jan4, "%u")) - 1L) + 7L * (week - 1L)
}

weekly <- readRDS("data/snapshots/weekly_deaths/0012100.rds")
weekly <- weekly[weekly$code %in% c("PT", REGIONS), , drop = FALSE]
weekly$band <- band_of(weekly$age_group)
weeks <- distinct(weekly, year, week)
weeks$monday <- iso_monday(weeks$year, weeks$week)
day_share <- bind_rows(lapply(seq_len(nrow(weeks)), function(i) {
  cal <- as.integer(format(weeks$monday[[i]] + 0:6, "%Y"))
  tibble(year = weeks$year[[i]], week = weeks$week[[i]], calendar = as.integer(names(table(cal))), days = as.integer(table(cal)))
}))

parish <- readRDS("data/snapshots/parish_vitals/deaths.rds")
lookup <- readRDS("data/nuts_lookup_2024.rds")
lookup$dico <- substr(as.character(lookup$municipality_code), 4, 7)

# A calendar year is complete when every day of it is covered by a published week.
covered_days <- function(year) {
  sum(day_share$days[day_share$calendar == year])
}
year_days <- function(year) as.integer(as.Date(sprintf("%d-12-31", year)) - as.Date(sprintf("%d-01-01", year))) + 1L

build_year <- function(year) {
  share <- day_share[day_share$calendar == year, c("year", "week", "days")]
  profile <- weekly %>%
    filter(!is.na(band)) %>%
    inner_join(share, by = c("year", "week")) %>%
    group_by(code, band) %>%
    summarise(deaths = sum(deaths * days / 7), .groups = "drop")
  totals <- parish %>%
    filter(year == !!year) %>%
    inner_join(lookup[, c("dico", "nuts2")], by = "dico") %>%
    group_by(nuts2) %>%
    summarise(total = sum(value), .groups = "drop")
  regional <- bind_rows(lapply(names(REGIONS), function(area) {
    p <- profile[profile$code == REGIONS[[area]], ]
    total <- totals$total[totals$nuts2 == area]
    if (nrow(p) == 0 || length(total) != 1) stop("Missing data for ", area, " in ", year, call. = FALSE)
    tibble(area = area, code = REGIONS[[area]], age_band = p$band, deaths = p$deaths * total / sum(p$deaths))
  }))
  # Portugal keeps the weekly total, deaths without a recorded age included:
  # its known-age profile is scaled to it, as the regions are to their parish
  # totals (11 deaths of unknown age in 2025).
  national <- profile[profile$code == "PT", ]
  national_total <- weekly %>%
    filter(code == "PT", age_group == "Total") %>%
    inner_join(share, by = c("year", "week")) %>%
    summarise(total = sum(deaths * days / 7)) %>%
    pull(total)
  if (length(national_total) == 1 && is.finite(national_total) && national_total > 0) {
    national$deaths <- national$deaths * national_total / sum(national$deaths)
  }
  out <- bind_rows(
    regional,
    regional %>% filter(area %in% MAINLAND) %>% group_by(age_band) %>%
      summarise(deaths = sum(deaths), .groups = "drop") %>% mutate(area = "Continente", code = "1"),
    tibble(area = "Portugal", code = "PT", age_band = national$band, deaths = national$deaths)
  )
  out %>%
    mutate(year = as.integer(year), sex = "HM", source_indicator = "0012100+0012542",
           age_band = factor(age_band, levels = BANDS)) %>%
    arrange(area, age_band) %>%
    mutate(age_band = as.character(age_band)) %>%
    select(year, area, code, sex, age_band, deaths, source_indicator)
}

archive_years <- function() {
  dirs <- list.dirs("data/snapshots/deaths", recursive = TRUE, full.names = FALSE)
  as.integer(sub("^.*year_", "", grep("year_\\d{4}$", dirs, value = TRUE)))
}
candidates <- sort(unique(day_share$calendar))
complete <- candidates[vapply(candidates, function(y) covered_days(y) == year_days(y), logical(1)) &
                       candidates %in% parish$year]
years_arg <- get_arg("years", "AUTO")
years <- if (identical(toupper(years_arg), "AUTO")) setdiff(complete, archive_years()) else as.integer(strsplit(years_arg, ",")[[1]])
if (length(setdiff(years, complete)) > 0) stop("Not a complete calendar year in the weekly data: ", paste(setdiff(years, complete), collapse = ", "), call. = FALSE)

for (year in years) {
  x <- build_year(year)
  message(sprintf("  %d: Portugal %.0f, Continente %.0f deaths", year,
                  sum(x$deaths[x$area == "Portugal"]), sum(x$deaths[x$area == "Continente"])))
  save_rds_atomic(x, file.path("data/snapshots/calendar_deaths_by_age", paste0("year_", year, ".rds")))
}
if (length(years) == 0) message("Nothing to build: every complete calendar year is in the by-cause archive.")

if (tolower(get_arg("validate", "false")) %in% c("true", "1", "yes")) {
  for (year in intersect(complete, archive_years())) {
    x <- build_year(year)
    message(year, ": rebuilt totals by region (compare with the by-cause archive)")
    print(x %>% group_by(area) %>% summarise(deaths = round(sum(deaths)), .groups = "drop"))
  }
}
message("Done.")
