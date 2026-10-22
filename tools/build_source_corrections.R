#!/usr/bin/env Rscript
# Build data/source_corrections.rds: replacements for values INE publishes
# wrongly, applied when a death file is read (apply_source_corrections() in
# R/helpers.R). The snapshots keep INE's data as published; the corrections
# sit beside them, so either can be used and every corrected value can be
# traced to the published one.
#
#   Rscript tools/build_source_corrections.R
#
# Each correction is listed below with its evidence. The script also reports
# every municipality-year whose all-cause deaths disagree with the parish
# deaths by more than 2, so a new defect shows up when the data are refreshed.
#
# The one known defect (found 2026-10-22, checked on INE's live API, last
# updated 2025-10-11): in 0008206 for 2011 the all-cause rows of five Alto
# Alentejo municipalities are wrong as published - Nisa 24 deaths (parish 165),
# Gavião 34 (91), Fronteira 1 (60), Monforte 31 (49), Marvão 49 (59). INE's own
# Alto Alentejo NUTS III row (1,888) equals the sum of its municipalities in
# 2010 and 2012 and exceeds it by exactly the missing 285 in 2011, and four of
# the five show more men than deaths in all. Cause-specific rows look right.
#
# Method ("parish total, neighbouring profile"):
#   - the total is the parish deaths of the year, which equal INE's municipal
#     totals in every other municipality-year of 1996-2024;
#   - where the published men and women add up to that total (Nisa), they are
#     kept and both sexes is their sum;
#   - otherwise the total is split by sex, and each sex by age, in the
#     proportions of the two neighbouring years, rounded to whole deaths that
#     keep the total (largest remainder).
suppressMessages(library(dplyr))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source("R/data_versions.R", envir = environment())

CORRECTIONS <- tibble::tribble(
  ~indicator, ~year, ~cause,                     ~area,
  "0008206",  2011L, "Todas as causas de morte", "Nisa",
  "0008206",  2011L, "Todas as causas de morte", "Gavião",
  "0008206",  2011L, "Todas as causas de morte", "Fronteira",
  "0008206",  2011L, "Todas as causas de morte", "Monforte",
  "0008206",  2011L, "Todas as causas de morte", "Marvão"
)
REASON <- paste(
  "INE 0008206, 2011: linha de todas as causas errada na publicação do INE (soma por sexo incoerente; a linha da",
  "NUTS III Alto Alentejo e os óbitos por freguesia dão o total certo). Total = óbitos por freguesia; repartição por",
  "sexo e idade dos anos vizinhos."
)

lookup <- readRDS("data/nuts_lookup_2024.rds")
lookup$dico <- substr(lookup$municipality_code, nchar(lookup$municipality_code) - 3L, nchar(lookup$municipality_code))
parish <- readRDS("data/snapshots/parish_vitals/deaths.rds") %>%
  mutate(dico = substr(code, 1, 4)) %>%
  group_by(year, dico) %>%
  summarise(parish = sum(value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(lookup %>% select(dico, area = municipality), by = "dico")

cause_file <- function(indicator, year, cause) {
  token <- gsub("[^a-z0-9]+", "_", tolower(iconv(cause, to = "ASCII//TRANSLIT")))
  readRDS(file.path("data/snapshots/deaths", indicator, paste0("year_", year), paste0("cause_", gsub("^_|_$", "", token), ".rds")))
}

# Whole numbers in proportion to `weights`, summing to `total`.
apportion <- function(total, weights) {
  if (sum(weights) <= 0) weights <- rep(1, length(weights))
  exact <- total * weights / sum(weights)
  out <- floor(exact)
  short <- total - sum(out)
  if (short > 0) out[order(exact - out, decreasing = TRUE)[seq_len(short)]] <- out[order(exact - out, decreasing = TRUE)[seq_len(short)]] + 1
  out
}

deaths_rows <- list(); totals_rows <- list()
for (i in seq_len(nrow(CORRECTIONS))) {
  k <- CORRECTIONS[i, ]
  now <- cause_file(k$indicator, k$year, k$cause) %>% filter(area == k$area)
  near <- bind_rows(cause_file(k$indicator, k$year - 1L, k$cause), cause_file(k$indicator, k$year + 1L, k$cause)) %>%
    filter(area == k$area) %>% group_by(sex, age_band) %>% summarise(deaths = sum(deaths), .groups = "drop")
  total <- parish$parish[parish$year == k$year & parish$area == k$area]
  stopifnot(length(total) == 1, total > 0)
  bands <- now$age_band[now$sex == "HM"]
  published <- function(s) setNames(now$deaths[now$sex == s], now$age_band[now$sex == s])[bands]
  profile <- function(s) setNames(near$deaths[near$sex == s], near$age_band[near$sex == s])[bands]
  if (abs(sum(published("H")) + sum(published("M")) - total) <= 1) {
    by_sex <- list(H = published("H"), M = published("M"))
    how <- "homens e mulheres publicados; ambos os sexos = a sua soma"
  } else {
    sex_totals <- apportion(total, c(sum(profile("H")), sum(profile("M"))))
    by_sex <- list(H = setNames(apportion(sex_totals[[1]], profile("H")), bands),
                   M = setNames(apportion(sex_totals[[2]], profile("M")), bands))
    how <- "sexo e idade repartidos como em 2010 e 2012"
  }
  by_sex$HM <- by_sex$H + by_sex$M
  if (abs(sum(by_sex$HM) - total) > 1) stop(k$area, ": corrected total ", sum(by_sex$HM), " against parish ", total)
  for (s in c("H", "M", "HM")) {
    deaths_rows[[length(deaths_rows) + 1]] <- tibble(
      indicator = k$indicator, year = k$year, area = k$area, sex = s, cause = k$cause, age_band = bands,
      deaths = as.numeric(by_sex[[s]]), published = as.numeric(published(s)))
    totals_rows[[length(totals_rows) + 1]] <- tibble(
      indicator = k$indicator, year = k$year, area = k$area, sex = s, cause = k$cause,
      deaths = sum(by_sex[[s]]), published = sum(published(s)))
  }
  message(sprintf("  %s %d %s: published %g (H %g, M %g) -> %g (%s)", k$indicator, k$year, k$area,
                  sum(published("HM")), sum(published("H")), sum(published("M")), sum(by_sex$HM), how))
}

corrections <- list(
  deaths = bind_rows(deaths_rows),
  totals = bind_rows(totals_rows),
  notes = CORRECTIONS %>% mutate(reason = REASON, built = format(Sys.Date()))
)
versioned_save_rds(corrections, "data/source_corrections.rds", tool = "build_source_corrections.R",
                   note = "INE 2011 Alto Alentejo all-cause rows")

# Every municipality-year whose published all-cause total differs from the
# parish deaths by more than 2 (confidential totals excluded), after correction.
check <- list()
for (indicator in c("0008206", "0013166")) {
  for (f in list.files(file.path("data/snapshots/death_totals", indicator), full.names = TRUE)) {
    t <- readRDS(f) %>% filter(sex == "HM", cause == "Todas as causas de morte", area %in% lookup$municipality)
    check[[f]] <- t %>% select(year, area, deaths)
  }
}
check <- bind_rows(check) %>% distinct(year, area, .keep_all = TRUE) %>%
  left_join(corrections$totals %>% filter(sex == "HM") %>% select(year, area, corrected = deaths), by = c("year", "area")) %>%
  mutate(deaths = coalesce(corrected, deaths)) %>%
  inner_join(parish, by = c("year", "area")) %>%
  filter(!is.na(deaths), abs(deaths - parish) > 2)
message("All-cause municipal totals differing from the parish deaths by more than 2, after correction: ", nrow(check))
if (nrow(check) > 0) print(as.data.frame(check %>% select(year, area, deaths, parish)))
