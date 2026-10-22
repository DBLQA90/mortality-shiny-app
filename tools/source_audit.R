#!/usr/bin/env Rscript
# Audit of the INE data as stored: internal coherence of each dataset, and
# agreement between datasets that measure the same thing. Reports; changes
# nothing. Every finding is described, with what the app does about it, in
# ESCOLHAS_METODOLOGICAS.md and METHODOLOGY.md.
#
#   Rscript tools/source_audit.R [--report=FILE]
#
# The data are read as INE publishes them (source corrections off), so a
# defect at the source still shows; the corrected values are listed apart.
suppressMessages(library(dplyr))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
report_file <- sub("^--report=", "", grep("^--report=", commandArgs(TRUE), value = TRUE))
out <- character(0)
say <- function(...) { line <- paste0(...); out <<- c(out, line); cat(line, "\n", sep = "") }
table_lines <- function(df) { lines <- utils::capture.output(print(as.data.frame(df), row.names = FALSE)); for (l in lines) say("    ", l) }

lookup <- readRDS("data/nuts_lookup_2024.rds")
lookup$dico <- substr(lookup$municipality_code, nchar(lookup$municipality_code) - 3L, nchar(lookup$municipality_code))
munis <- lookup$municipality
ALL <- "Todas as causas de morte"
death_indicator <- function(y) if (y >= 2022) "0013166" else "0008206"
totals_of <- function(y) { f <- sprintf("data/snapshots/death_totals/%s/year_%d.rds", death_indicator(y), y); if (file.exists(f)) readRDS(f) else NULL }
bands_of <- function(y, cause = "todas_as_causas_de_morte") {
  f <- sprintf("data/snapshots/deaths/%s/year_%d/cause_%s.rds", death_indicator(y), y, cause); if (file.exists(f)) readRDS(f) else NULL
}
years <- 1991:2025
say("== INE source audit, ", format(Sys.time(), "%Y-%m-%d %H:%M"), " ==")

# 1. Deaths: municipal totals against Portugal, confidential cells, sexes.
say("\n[1] Deaths, all causes: municipal totals")
rows <- list()
for (y in years) {
  t <- totals_of(y); b <- bands_of(y); if (is.null(t) || is.null(b)) next
  ta <- t[t$cause == ALL, ]
  hm <- ta[ta$sex == "HM" & ta$area %in% munis, ]
  w <- ta[ta$area %in% munis, ] %>% tidyr::pivot_wider(id_cols = area, names_from = sex, values_from = deaths)
  bs <- b[b$sex == "HM" & b$area %in% munis, ] %>% group_by(area) %>% summarise(bands = sum(deaths), .groups = "drop")
  j <- inner_join(hm %>% select(area, total = deaths), bs, by = "area")
  rows[[length(rows) + 1]] <- tibble(
    year = y,
    portugal = sum(ta$deaths[ta$area == "Portugal" & ta$sex == "HM"]),
    municipal_sum = sum(hm$deaths, na.rm = TRUE),
    unattributed_pct = round((portugal - municipal_sum) / portugal * 100, 2),
    confidential_totals = sum(is.na(t$deaths[t$area %in% munis & t$sex == "HM"])),
    hm_not_h_plus_m = sum(abs(w$HM - (w$H + w$M)) > 2, na.rm = TRUE),
    without_age_pct = round((sum(j$total, na.rm = TRUE) - sum(j$bands)) / sum(j$total, na.rm = TRUE) * 100, 2),
    munis_over_10pct_without_age = sum((j$total - j$bands) / j$total > 0.10, na.rm = TRUE)
  )
}
deaths_tbl <- bind_rows(rows)
table_lines(deaths_tbl)

# 2. Deaths by cause: how much of the municipal age detail is missing.
say("\n[2] Deaths by cause, municipal age bands against the all-ages totals (lung cancer and all causes)")
lung <- sub("^cause_|\\.rds$", "", grep("pulm", list.files("data/snapshots/deaths/0008206/year_2014"), value = TRUE)[[1]])
lung <- sub("\\.rds$", "", lung)
rows <- list()
for (y in years) {
  t <- totals_of(y); b <- bands_of(y, lung); if (is.null(t) || is.null(b)) next
  cause <- unique(b$cause)
  tot <- sum(t$deaths[t$cause == cause & t$sex == "HM" & t$area %in% munis], na.rm = TRUE)
  rows[[length(rows) + 1]] <- tibble(year = y, lung_total = tot, lung_in_bands = sum(b$deaths[b$sex == "HM" & b$area %in% munis]),
                                     lung_without_age_pct = round((tot - lung_in_bands) / tot * 100, 1))
}
table_lines(bind_rows(rows) %>% filter(lung_without_age_pct > 1 | year %in% c(1991, 2002, 2013, 2014, 2015, 2021, 2024)))

# 3. Parish deaths against municipal totals.
say("\n[3] Parish deaths against INE's municipal all-cause totals (|difference| > 2)")
parish <- readRDS("data/snapshots/parish_vitals/deaths.rds") %>% mutate(dico = substr(code, 1, 4)) %>%
  group_by(year, dico) %>% summarise(parish = sum(value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(lookup %>% select(dico, area = municipality), by = "dico")
tot_all <- bind_rows(lapply(years, function(y) { t <- totals_of(y); if (is.null(t)) NULL else t[t$cause == ALL & t$sex == "HM" & t$area %in% munis, c("year", "area", "deaths")] }))
pm <- inner_join(parish, tot_all, by = c("year", "area")) %>% mutate(diff = deaths - parish)
say("    municipality-years compared: ", nrow(pm), " (", min(pm$year), "-", max(pm$year), "); exact: ", sum(pm$diff == 0, na.rm = TRUE),
    "; confidential total: ", sum(is.na(pm$diff)))
table_lines(pm %>% filter(abs(diff) > 2) %>% select(year, area, ine_total = deaths, parish))

# 4. Regional rows: bands against totals, and against the municipal sums.
say("\n[4] INE regional rows (all causes, both sexes)")
rows <- list()
for (y in years) {
  for (ind in unique(c(death_indicator(y), if (y == 2022) "0008206"))) {
    f <- sprintf("data/snapshots/regional_deaths/%s/year_%d.rds", ind, y); if (!file.exists(f)) next
    r <- readRDS(f); r <- r[r$cause == ALL & r$sex == "HM", ]
    rows[[length(rows) + 1]] <- r %>% group_by(year, region_code, area) %>% summarise(deaths = sum(deaths), .groups = "drop") %>% mutate(indicator = ind)
  }
}
regional <- bind_rows(rows)
say("    region-years stored: ", nrow(regional), "; regions: ", n_distinct(regional$region_code))
lk13 <- readRDS("data/nuts_lookup_2013.rds")
comp <- tot_all %>% left_join(lk13 %>% select(area = municipality, nuts3_13 = nuts3), by = "area") %>% group_by(year, nuts3_13) %>%
  summarise(municipal = sum(deaths, na.rm = TRUE), has_na = any(is.na(deaths)), .groups = "drop")
nuts3 <- regional %>% filter(indicator == "0008206") %>% inner_join(comp, by = c("year", "area" = "nuts3_13")) %>%
  mutate(diff = deaths - municipal) %>% filter(!has_na, abs(diff) > 2)
say("    NUTS III rows (NUTS 2013) against the sum of their municipalities, |difference| > 2: ", nrow(nuts3))
table_lines(nuts3 %>% select(year, area, ine_row = deaths, municipal_sum = municipal, diff))

# 5. Population: sexes, municipal sum against Portugal, year-on-year steps.
say("\n[5] Population (31 December), all ages")
rows <- list(); prev <- NULL
for (y in years) {
  f <- sprintf("data/snapshots/population/year_%d.rds", y); if (!file.exists(f)) next
  p <- readRDS(f); m <- p[p$area %in% munis, ]
  w <- m %>% group_by(area, sex) %>% summarise(pop = sum(pop), .groups = "drop") %>% tidyr::pivot_wider(names_from = sex, values_from = pop)
  now <- w %>% select(area, HM)
  steps <- if (is.null(prev)) NA else { s <- inner_join(prev, now, by = "area"); sum(abs(s$HM.y / s$HM.x - 1) > 0.03) }
  rows[[length(rows) + 1]] <- tibble(year = y, municipalities = nrow(w), municipal_sum = round(sum(w$HM)),
    portugal = round(sum(p$pop[p$area == "Portugal" & p$sex == "HM"])), hm_not_h_plus_m = sum(abs(w$HM - (w$H + w$M)) > 1),
    munis_moving_over_3pct = steps)
  prev <- now
}
table_lines(bind_rows(rows))

# 6. Births and infant deaths: municipal sums against Portugal.
say("\n[6] Live births and under-1 deaths: municipal sums against Portugal")
rows <- list()
for (y in years) {
  fb <- sprintf("data/snapshots/births/year_%d.rds", y); fi <- sprintf("data/snapshots/infant_totals/year_%d.rds", y)
  b <- if (file.exists(fb)) readRDS(fb) else NULL; i <- if (file.exists(fi)) readRDS(fi) else NULL
  if (is.null(b) && is.null(i)) next
  rows[[length(rows) + 1]] <- tibble(year = y,
    births_municipal = if (is.null(b)) NA else sum(b$births[b$area %in% munis]),
    births_portugal = if (is.null(b)) NA else sum(b$births[b$area == "Portugal"]),
    infant_municipal = if (is.null(i)) NA else sum(i$deaths[i$area %in% munis & i$sex == "HM"]),
    infant_portugal = if (is.null(i)) NA else sum(i$deaths[i$area == "Portugal" & i$sex == "HM"]))
}
table_lines(bind_rows(rows))

# 7. Weekly deaths against the calendar-year parish deaths, by NUTS II.
say("\n[7] Weekly deaths (ISO weeks, whole weeks summed by year) against parish deaths")
wk <- readRDS("data/snapshots/weekly_deaths/0012100.rds")
wk_year <- wk %>% filter(code == "PT", age_group == "Total") %>% group_by(year) %>% summarise(weekly = sum(deaths), .groups = "drop")
par_year <- parish %>% group_by(year) %>% summarise(parish = sum(parish), .groups = "drop")
table_lines(inner_join(wk_year, par_year, by = "year") %>% mutate(diff_pct = round((weekly / parish - 1) * 100, 2)))

# 8. Datasets published for every municipality: municipalities without a value.
say("\n[8] Other municipal datasets: municipality-years without a value")
rows <- list()
for (measure in list.files("data/snapshots/planning_extra")) {
  for (f in list.files(file.path("data/snapshots/planning_extra", measure), full.names = TRUE)) {
    d <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(d) || !"area" %in% names(d)) next
    if ("category" %in% names(d) && "Total" %in% d$category) d <- d[d$category == "Total", , drop = FALSE]
    vcol <- intersect(c("value", "count", "births", "deaths", "pop", "employees", "amount"), names(d))
    if (length(vcol) == 0) next
    have <- unique(d$area[d$area %in% munis & !is.na(d[[vcol[[1]]]])])
    rows[[length(rows) + 1]] <- tibble(measure = measure, year = as.integer(sub("\\D*(\\d{4}).*", "\\1", basename(f))), missing = length(setdiff(munis, have)))
  }
}
gaps <- bind_rows(rows) %>% filter(missing > 0) %>% group_by(measure) %>%
  summarise(years_with_gaps = n(), worst_year = year[which.max(missing)], max_missing = max(missing), .groups = "drop")
table_lines(gaps)

# 9. Corrections applied on reading (data/source_corrections.rds).
say("\n[9] Source corrections in force")
if (file.exists("data/source_corrections.rds")) {
  corr <- readRDS("data/source_corrections.rds")
  table_lines(corr$totals %>% filter(sex == "HM") %>% select(indicator, year, area, cause, published, corrected = deaths))
} else say("    none")

if (length(report_file) == 1) writeLines(out, report_file)
