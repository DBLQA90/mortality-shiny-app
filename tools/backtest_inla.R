#!/usr/bin/env Rscript
# Compare INLA backtests built by tools/precompute_inla.R origin=YYYY (one
# directory per variant) with the rates observed after the origin: by variant
# and level, the mean relative error of the median (bias), the share of
# observed rates inside the 95% predictive interval as fitted (not widened, so
# the variants compare like for like), and the median absolute error per
# 100 000. Written for the 2021 step (roadmap item 5): with the COVID years
# separated it can only be told from them by the years after 2022.
#
#   Rscript tools/backtest_inla.R dirs=data/inla_backtests/2023_C,data/inla_backtests/2023_Cs2021 [report=FILE]
suppressMessages({library(dplyr); library(tibble)})
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
arg_of <- function(name, default) {
  x <- sub(paste0("^", name, "="), "", grep(paste0("^", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE))
  if (length(x)) x else default
}
dirs <- strsplit(arg_of("dirs", ""), ",", fixed = TRUE)[[1]]
stopifnot(length(dirs) >= 1, all(dir.exists(dirs)))
report <- arg_of("report", "")
source("R/config.R", local = TRUE)  # age_levels, esp2013_df
weights <- esp2013_df$stdpop[match(age_levels, as.character(esp2013_df$age_band))]
names(weights) <- age_levels

# Observed standardised rate of each unit and year, over the bands the fit kept.
observed <- function(extract, members, kept, years) {
  d <- extract %>% filter(year %in% years, band %in% kept)
  units <- bind_rows(tibble(unit = unique(d$area), area = unique(d$area)),
                     distinct(members, unit = aggregate, area))
  # Weighted by the whole standard population, as the forecasts and the app.
  w <- weights[kept] / sum(weights)
  d %>% inner_join(units, by = "area", relationship = "many-to-many") %>%
    group_by(unit, year, band) %>% summarise(deaths = sum(deaths), py = sum(py), .groups = "drop") %>%
    filter(py > 0) %>% group_by(unit, year) %>%
    summarise(obs = sum(w[as.character(band)] * deaths / py) * 1e5, .groups = "drop")
}

rows <- list()
for (dir in dirs) {
  manifest <- readRDS(file.path(dir, "manifest.rds"))
  setup <- readRDS(file.path(dir, "extract", "setup.rds"))
  years <- seq.int(manifest$origin + 1L, manifest$last_observed)
  if (length(years) == 0 || manifest$origin >= manifest$last_observed) stop(dir, ": no observed year after the origin")
  for (f in list.files(dir, "__B[01]C[01]\\.rds$", full.names = TRUE)) {
    x <- readRDS(f)
    key <- sub("__B[01]C[01]\\.rds$", "", basename(f))
    obs <- observed(readRDS(file.path(dir, "extract", paste0(key, ".rds"))), setup$members, x$kept, years)
    rows[[length(rows) + 1]] <- x$forecast %>% filter(year %in% years) %>%
      inner_join(obs, by = c(area = "unit", "year")) %>%
      transmute(variant = basename(dir), key, level, area, year, horizon = year - manifest$origin,
                point = q500, lower = p25, upper = p975, obs)
  }
}
cases <- bind_rows(rows)
levels <- c("municipio", "ULS", "NUTS III", "NUTS II", "Portugal")
summary <- cases %>% mutate(level = factor(level, levels)) %>% group_by(variant, level) %>%
  summarise(n = n(), bias = round(100 * mean((point - obs)[obs > 0] / obs[obs > 0]), 1),
            coverage = round(mean(obs >= lower & obs <= upper), 3), error = round(stats::median(abs(point - obs)), 2),
            .groups = "drop") %>% arrange(level, variant)
out <- c(sprintf("INLA backtest, %s; years %s; %d cases", paste(basename(dirs), collapse = " vs "),
                 paste(range(cases$year), collapse = "-"), nrow(cases)),
         utils::capture.output(print(as.data.frame(summary), row.names = FALSE)))
writeLines(out)
if (nzchar(report)) writeLines(out, report)
