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

fixed_totals <- bind_rows(totals_rows) %>% mutate(method = "freguesias")

# ---- Totals by sex that INE withholds -------------------------------------
# INE publishes a municipality's all-ages total of a cause for both sexes and
# leaves the split by sex blank ("Dado confidencial") where it is small: 1 200 to
# 1 700 municipality x cause cells a year, 6 800 in 2014, when it also blanked
# 4 300 totals of both sexes. The app completes the age bands only up to a
# published total, so a blank left the counts by sex at what the bands hold:
# municipal sums 95-98% of Portugal in a year, 65% in 2014, with no flag.
#
# The NUTS II rows are complete and equal the sum of their municipalities
# wherever none is blank (they hold no deaths of unknown residence: all NUTS II
# together are Portugal to 0.02%). So, never below what the age bands record:
#   - one sex blank, the other and both sexes published: the difference (it is
#     always 0: this is how INE writes a zero);
#   - both sexes blank: the NUTS II residual not yet recorded, shared among the
#     blank municipalities by their expected deaths not yet recorded (their
#     population at Portugal's rates by age);
#   - men and women blank: the municipality's total split between them by
#     iterative proportional fitting, rows = each total, columns = the NUTS II
#     residual by sex, seed = expected deaths not yet recorded by sex.
# Whole deaths throughout (largest remainder). Validated in 2013 by blanking
# both sexes of 30% of the municipalities and their age bands by sex: 4.5%
# error in cells above 100 deaths, under half a death in cells of 1-2
# (forecast_smallarea_2026/40_recover_sex_split.R).
SEX_REASON <- paste(
  "Total por sexo ocultado pelo INE («Dado confidencial»): recuperado do total de ambos os sexos e da linha da NUTS II,",
  "que é completa, nunca abaixo do que as idades registam."
)
lookups <- list("2013" = readRDS("data/nuts_lookup_2013.rds"), "2024" = readRDS("data/nuts_lookup_2024.rds"))

ipf2 <- function(seed, rows, cols, iterations = 50L) {
  x <- seed
  for (k in seq_len(iterations)) {
    x <- x * ifelse(rowSums(x) > 0, rows / rowSums(x), 0)
    cs <- colSums(x)
    x <- sweep(x, 2, ifelse(cs > 0, cols / cs, 0), `*`)
  }
  x * ifelse(rowSums(x) > 0, rows / rowSums(x), 0)
}

recover_sex_year <- function(indicator, year) {
  totals_path <- file.path("data/snapshots/death_totals", indicator, paste0("year_", year, ".rds"))
  regional_path <- file.path("data/snapshots/regional_deaths", indicator, paste0("year_", year, ".rds"))
  bands_dir <- file.path("data/snapshots/deaths", indicator, paste0("year_", year))
  if (!file.exists(totals_path) || !file.exists(regional_path) || !dir.exists(bands_dir)) return(NULL)
  lookup <- lookups[[if (indicator == "0013166") "2024" else "2013"]]
  pop_years <- as.integer(sub("^year_(\\d+)\\.rds$", "\\1", list.files("data/snapshots/population")))
  pop <- readRDS(file.path("data/snapshots/population", paste0("year_", max(pop_years[pop_years <= year]), ".rds")))

  # Published totals, with the corrections above in force.
  totals <- readRDS(totals_path) %>% filter(area %in% lookup$municipality) %>% select(area, sex, cause, deaths) %>%
    left_join(fixed_totals %>% filter(.data$indicator == !!indicator, .data$year == !!year) %>% select(area, sex, cause, fixed = deaths),
              by = c("area", "sex", "cause")) %>%
    mutate(deaths = coalesce(fixed, deaths)) %>% select(-fixed)
  if (!anyNA(totals$deaths)) return(NULL)
  bands <- bind_rows(lapply(list.files(bands_dir, full.names = TRUE), readRDS))
  recorded <- bands %>% filter(area %in% lookup$municipality) %>% group_by(area, cause, sex) %>%
    summarise(recorded = sum(deaths, na.rm = TRUE), .groups = "drop")
  rates <- bands %>% filter(area == "Portugal") %>% select(cause, sex, age_band, d = deaths) %>%
    inner_join(pop %>% filter(area == "Portugal") %>% select(sex, age_band, p = pop), by = c("sex", "age_band")) %>%
    mutate(rate = ifelse(p > 0, d / p, 0))
  expected <- pop %>% filter(area %in% lookup$municipality) %>%
    inner_join(rates %>% select(cause, sex, age_band, rate), by = c("sex", "age_band"), relationship = "many-to-many") %>%
    group_by(area, cause, sex) %>% summarise(expected = sum(pop * rate), .groups = "drop")
  regional <- readRDS(regional_path) %>% filter(nchar(region_code) == 2) %>%
    group_by(region = area, cause, sex) %>% summarise(regional = sum(deaths, na.rm = TRUE), .groups = "drop")
  stopifnot(all(unique(lookup$nuts2) %in% regional$region))

  w <- totals %>% left_join(recorded, by = c("area", "cause", "sex")) %>% left_join(expected, by = c("area", "cause", "sex")) %>%
    mutate(recorded = coalesce(recorded, 0), expected = coalesce(expected, 0), region = lookup$nuts2[match(area, lookup$municipality)]) %>%
    left_join(regional, by = c("region", "cause", "sex")) %>%
    tidyr::pivot_wider(id_cols = c(area, cause, region), names_from = sex, values_from = c(deaths, recorded, expected, regional))
  w <- w %>% filter(!is.na(regional_HM), is.na(deaths_HM) | is.na(deaths_H) | is.na(deaths_M))
  if (nrow(w) == 0) return(NULL)
  w <- w %>% mutate(HM = deaths_HM, H = deaths_H, M = deaths_M)

  out <- list()
  for (g in split(w, list(w$region, w$cause), drop = TRUE)) {
    # Both sexes blank: the NUTS II residual over the blank municipalities.
    all_rows <- totals$cause == g$cause[1] & totals$sex == "HM" & lookup$nuts2[match(totals$area, lookup$municipality)] == g$region[1]
    blank <- is.na(g$HM)
    if (any(blank)) {
      floor_hm <- pmax(g$recorded_HM[blank], coalesce(g$H[blank], g$recorded_H[blank]) + coalesce(g$M[blank], g$recorded_M[blank]))
      free <- max(g$regional_HM[1] - sum(totals$deaths[all_rows], na.rm = TRUE) - sum(floor_hm), 0)
      g$HM[blank] <- floor_hm + apportion(round(free), pmax(g$expected_HM[blank] - floor_hm, 0) + 1e-6)
    }
    # One sex blank: the difference.
    g$H <- ifelse(is.na(g$H) & !is.na(g$M), pmax(g$HM - g$M, g$recorded_H), g$H)
    g$M <- ifelse(is.na(g$M) & !is.na(g$H), pmax(g$HM - g$H, g$recorded_M), g$M)
    # Men and women blank: the total split by IPF against the residual by sex.
    both <- is.na(g$H) & is.na(g$M)
    if (any(both)) {
      lo <- cbind(g$recorded_H[both], g$recorded_M[both])
      rows <- pmax(g$HM[both] - rowSums(lo), 0)
      sex_rows <- function(s) totals$cause == g$cause[1] & totals$sex == s & lookup$nuts2[match(totals$area, lookup$municipality)] == g$region[1]
      known <- c(sum(totals$deaths[sex_rows("H")], na.rm = TRUE) + sum(g$H[!both & is.na(g$deaths_H)]),
                 sum(totals$deaths[sex_rows("M")], na.rm = TRUE) + sum(g$M[!both & is.na(g$deaths_M)]))
      need <- pmax(c(g$regional_H[1], g$regional_M[1]) - known - colSums(lo), 0)
      need <- if (sum(need) > 0) need / sum(need) * sum(rows) else rep(sum(rows) / 2, 2)
      seed <- cbind(pmax(g$expected_H[both] - lo[, 1], 0), pmax(g$expected_M[both] - lo[, 2], 0)) + 1e-6
      fit <- ipf2(seed, rows, need)
      parts <- t(vapply(seq_along(rows), function(i) apportion(rows[i], fit[i, ]), numeric(2)))
      g$H[both] <- lo[, 1] + parts[, 1]
      g$M[both] <- lo[, 2] + parts[, 2]
    }
    out[[length(out) + 1]] <- g
  }
  r <- bind_rows(out)
  method_of <- function(published, s) {
    how <- if (s == "HM") rep("resíduo da NUTS II", nrow(r)) else
      ifelse(is.na(r$deaths_H) & is.na(r$deaths_M), "repartição do total por sexo", "diferença")
    ifelse(is.na(published), how, NA_character_)
  }
  bind_rows(lapply(c("H", "M", "HM"), function(s) {
    tibble(indicator = indicator, year = as.integer(year), area = r$area, sex = s, cause = r$cause,
           deaths = as.numeric(r[[s]]), published = r[[paste0("deaths_", s)]], method = method_of(r[[paste0("deaths_", s)]], s))
  })) %>% filter(!is.na(method))
}

recovered <- list()
for (indicator in c("0008206", "0013166", "0008461")) {
  years <- sort(as.integer(sub("^year_(\\d+)\\.rds$", "\\1", list.files(file.path("data/snapshots/death_totals", indicator)))))
  for (year in years) {
    r <- recover_sex_year(indicator, year)
    if (!is.null(r)) recovered[[paste(indicator, year)]] <- r
  }
  message(sprintf("  %s: %d blank totals recovered", indicator, sum(vapply(recovered[startsWith(names(recovered), indicator)], nrow, integer(1)))))
}
recovered <- bind_rows(recovered)
message("  by method: ", paste(names(table(recovered$method)), table(recovered$method), collapse = ", "))

corrections <- list(
  deaths = bind_rows(deaths_rows),
  totals = bind_rows(fixed_totals, recovered),
  notes = bind_rows(
    CORRECTIONS %>% mutate(reason = REASON, built = format(Sys.Date())),
    tibble(indicator = c("0008206", "0013166", "0008461"), year = NA_integer_, cause = NA_character_, area = NA_character_,
           reason = SEX_REASON, built = format(Sys.Date()))
  )
)
versioned_save_rds(corrections, "data/source_corrections.rds", tool = "build_source_corrections.R",
                   note = "INE 2011 Alto Alentejo all-cause rows; totals by sex withheld by INE")

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
  left_join(corrections$totals %>% filter(sex == "HM", cause == "Todas as causas de morte", indicator != "0008461") %>%
              distinct(year, area, .keep_all = TRUE) %>% select(year, area, corrected = deaths), by = c("year", "area")) %>%
  mutate(deaths = coalesce(corrected, deaths)) %>%
  inner_join(parish, by = c("year", "area")) %>%
  filter(!is.na(deaths), abs(deaths - parish) > 2)
message("All-cause municipal totals differing from the parish deaths by more than 2, after correction: ", nrow(check))
if (nrow(check) > 0) print(as.data.frame(check %>% select(year, area, deaths, parish)))
