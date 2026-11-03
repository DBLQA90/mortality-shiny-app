#!/usr/bin/env Rscript
# Precompute the small-area forecasts of the INLA model (roadmap item 7): the
# cause-specific mortality of all 308 municipalities, fitted together, and its
# posterior samples summed to the ULS, NUTS III, NUTS II and Portugal - so every
# aggregate is exactly the sum of its municipalities. The app reads the results
# (R/inla_forecasts.R); fitting takes about 6 GB and one to two minutes per
# cause and sex, too slow to run while someone waits.
#
#   Rscript tools/precompute_inla.R [causes=a;b] [combos=B1C1,B0C1,...] [horizon=11] [extract_only=1]
#
# causes= replaces the 8 backtested causes (INLA_CAUSES) by others, separated
# by ";"; extract_only=1 builds the data only (no INLA, little memory).
#
# INLA is needed here only, not by the app: install it from its own repository,
#   install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
# Output in data/inla_forecasts/ (not in git: rebuilt from the snapshots).
# Resumable: a finished cause x sex x combination is skipped; extraction is
# cached in data/inla_forecasts/extract/. Fits run one at a time: two at once
# need more memory than 23 GB.
#
# The model (V3 in ~/Desktop/DGS/Mortality/forecast_smallarea_2026/README.md):
# Poisson deaths by municipality, age band and year, offset log person-years;
# age effect, a national linear trend that may change slope in the last 6
# years, AR1 deviations of the years, BYM2 municipality effect (neighbours
# within the mainland and within each island), a trend of each municipality
# and an iid municipality-year deviation; PC priors. Age bands with fewer than
# 20 deaths in the country over the period are left out. Fitted on 1999 to the
# last year of causes of death; the projected years keep that year's
# population. Two switches (B off and C on by default, R/inla_forecasts.R):
#   B  a level shift where INE's population series breaks (2014, new
#      indicator; 2021, census revision, only without C), so the jump in
#      the denominator is not read as trend;
#   C  an effect for each of 2020, 2021 and 2022, so the COVID shocks keep
#      their level but bend no trend; a projected year never has one.
# Backtested at 9 origins, 2008-2022, on the corrected data (decisions_*_all.txt
# there).
suppressMessages({library(dplyr); library(tibble)})
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA is not installed; see the header of this script.")
suppressMessages(library(INLA))
Sys.setenv(MORTALITY_INSTALL_MISSING_PACKAGES = "false")
options(mortality.data_mode = "corrected")

app <- new.env(parent = globalenv())
assign("app_dir", getwd(), app); assign("get_app_dir", function() getwd(), app)
for (f in c("R/config.R", "R/helpers.R", "R/cache.R", "R/snapshots.R", "R/ine_client.R", "R/regions.R", "R/infant.R",
            "R/regional_rows.R", "R/metadata.R", "R/metrics.R", "R/standardisation.R", "R/avoidable.R",
            "R/planning_parish.R", "R/planning_indicators.R", "R/life_expectancy.R", "R/planning_under75.R",
            "R/planning_standardised.R", "R/planning_sns.R", "R/planning_weekly.R", "R/planning_export.R",
            "R/planning_charts.R", "R/planning_profile.R", "R/planning_maps.R", "R/cause_maps.R",
            "R/data_versions.R", "R/forecast_helpers.R", "R/trend.R", "R/data_access.R", "R/inla_forecasts.R")) {
  sys.source(f, envir = app)
}
attach(app, warn.conflicts = FALSE)

arg_of <- function(name, default) {
  x <- sub(paste0("^", name, "="), "", grep(paste0("^", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE))
  if (length(x)) x else default
}
causes <- arg_of("causes", "default")
causes <- if (identical(causes, "default")) INLA_CAUSES else strsplit(causes, ";", fixed = TRUE)[[1]]
combos <- strsplit(arg_of("combos", paste(INLA_COMBOS, collapse = ",")), ",", fixed = TRUE)[[1]]
stopifnot(all(combos %in% INLA_COMBOS))
HORIZON <- as.integer(arg_of("horizon", "11"))
extract_only <- identical(arg_of("extract_only", "0"), "1")
N_SAMPLES <- 300L
THREADS <- Sys.getenv("INLA_THREADS", "16:1")
QUANTILES <- c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)
out_dir <- inla_forecast_dir()
dir.create(file.path(out_dir, "extract"), recursive = TRUE, showWarnings = FALSE)
lookup <- get_nuts_lookup("2024")
areas <- sort(unique(as.character(lookup$municipality)))
stopifnot(length(areas) == 308)
last_year <- max(planning_cause_years())
years <- INLA_FIRST_YEAR:last_year
weights <- esp2013_df$stdpop[match(age_levels, as.character(esp2013_df$age_band))]
weights <- weights / sum(weights); names(weights) <- age_levels

# 1. Deaths and person-years by municipality, age band and year (the planning
#    pool, as the maps' annual series, so the observed rates are the app's).
extract <- function(cause, sex) {
  path <- file.path(out_dir, "extract", paste0(inla_key(cause, sex), ".rds"))
  if (file.exists(path)) return(readRDS(path))
  measure <- planning_cause_measure(cause)
  rows <- lapply(years, function(y) {
    pool <- planning_standardised_pool(areas, y, sex, lookup, extra = measure, span = 1L)
    if (is.null(pool)) stop("no data for ", cause, " ", sex, " ", y)
    d <- pool$deaths[, , measure]
    tibble(area = rep(rownames(d), times = ncol(d)), band = rep(colnames(d), each = nrow(d)), year = y,
           deaths = as.vector(d), py = as.vector(pool$person_years[rownames(d), colnames(d)]))
  })
  data <- bind_rows(rows)
  saveRDS(data, path)
  data
}

# 2. The neighbours (within the mainland and within each island) and the
#    aggregates the municipalities are summed to.
setup_path <- file.path(out_dir, "extract", "setup.rds")
if (file.exists(setup_path)) {
  setup <- readRDS(setup_path)
} else {
  b <- planning_map_boundaries("municipio", "2024")$polygons
  b <- b[b$id %in% areas, ]
  geom <- lapply(areas, function(a) {
    p <- b[b$id == a, ]
    polys <- lapply(split(p, factor(p$piece, levels = unique(p$piece))), function(pc) {
      lapply(split(pc, factor(pc$hole, levels = unique(pc$hole))), function(r) {
        m <- as.matrix(r[, c("x", "y")]); if (!all(m[1, ] == m[nrow(m), ])) m <- rbind(m, m[1, ]); m })
    })
    sf::st_multipolygon(unname(lapply(polys, unname)))
  })
  g <- sf::st_make_valid(sf::st_sfc(geom, crs = 4326))
  touch <- sf::st_intersects(sf::st_buffer(sf::st_transform(g, 3763), 50))
  adj <- Matrix::sparseMatrix(i = rep(seq_along(touch), lengths(touch)), j = unlist(touch), x = 1, dims = c(length(areas), length(areas)))
  diag(adj) <- 0; adj <- (adj > 0) * 1
  uls <- readRDS("data/uls_lookup.rds"); uls <- uls[uls$kind == "ULS", c("unit", "municipality")]
  members <- bind_rows(
    tibble(level = "ULS", aggregate = uls$unit, area = uls$municipality),
    tibble(level = "NUTS III", aggregate = lookup$nuts3, area = lookup$municipality),
    tibble(level = "NUTS II", aggregate = lookup$nuts2, area = lookup$municipality),
    tibble(level = "Portugal", aggregate = "Portugal", area = areas)
  ) %>% filter(area %in% areas) %>% distinct()
  setup <- list(adj = adj, members = members)
  saveRDS(setup, setup_path)
}
# The Algarve, the Azores and Madeira are a NUTS II and a NUTS III of one name:
# their members are taken once.
stopifnot(!any(setup$members$aggregate %in% areas))
aggregates <- distinct(setup$members, aggregate, area)
units <- c(setNames(as.list(areas), areas), split(aggregates$area, aggregates$aggregate))
level_of <- c(setNames(rep("municipio", length(areas)), areas),
              with(distinct(setup$members, level, aggregate), setNames(level, aggregate)))

h_bym <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)), phi = list(prior = "pc", param = c(0.5, 0.5)))
h1 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
h01 <- list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))
q <- function(m) t(apply(m, 1, stats::quantile, QUANTILES, names = FALSE))

# The standardised rate of a set of rows, per year: the expected deaths
# (credible) and Poisson draws around them (predictive), summed by band.
summarise_unit <- function(r, d, mu_all, noisy_all, w) {
  key_by <- paste(d$year[r], d$band[r])
  mu <- rowsum(mu_all[r, , drop = FALSE], key_by); noisy <- rowsum(noisy_all[r, , drop = FALSE], key_by)
  py <- as.vector(tapply(d$py[r], key_by, sum)[rownames(mu)])
  yr <- as.integer(sub(" .*$", "", rownames(mu))); band <- sub("^[0-9]+ ", "", rownames(mu))
  wr <- w[band] / py * 1e5
  rate <- rowsum(mu * wr, yr); obs <- rowsum(noisy * wr, yr)
  qr <- q(rate); qo <- q(obs)
  colnames(qr) <- paste0("q", QUANTILES * 1000); colnames(qo) <- paste0("p", QUANTILES * 1000)
  tibble(year = as.integer(rownames(rate)), mean = rowMeans(rate), as_tibble(qr), as_tibble(qo))
}

fit_one <- function(cause, sex, combo) {
  path <- file.path(out_dir, paste0(inla_key(cause, sex), "__", combo, ".rds"))
  if (file.exists(path)) return(invisible(NULL))
  t0 <- Sys.time()
  steps <- grepl("B1", combo, fixed = TRUE); covid <- grepl("C1", combo, fixed = TRUE)
  all <- extract(cause, sex)
  totals <- all %>% group_by(band) %>% summarise(deaths = sum(deaths), .groups = "drop")
  kept <- intersect(age_levels, as.character(totals$band[totals$deaths >= 20]))
  origin <- last_year
  future <- all %>% filter(year == origin) %>% select(area, band, py) %>%
    tidyr::crossing(year = origin + seq_len(HORIZON)) %>% mutate(deaths = NA_real_)
  d <- bind_rows(all, future) %>% filter(band %in% kept, py > 0) %>%
    mutate(band = factor(band, levels = kept), y = ifelse(year <= origin, round(deaths), NA_real_),
           area_id = match(area, areas), slope_id = area_id, year_id = year - INLA_FIRST_YEAR + 1L,
           year_c = year - origin, recent_c = pmax(year - (origin - 6L), 0L), st_id = (area_id - 1L) * 100L + year_id,
           step2014 = as.numeric(year >= 2014L), step2021 = as.numeric(year >= 2021L),
           covid2020 = as.numeric(year == 2020L), covid2021 = as.numeric(year == 2021L), covid2022 = as.numeric(year == 2022L))
  formula <- y ~ band + year_c + recent_c +
    f(area_id, model = "bym2", graph = setup$adj, scale.model = TRUE, constr = TRUE, hyper = h_bym) +
    f(year_id, model = "ar1", constr = TRUE, hyper = h1) +
    f(slope_id, year_c, model = "iid", constr = TRUE, hyper = h01) + f(st_id, model = "iid", hyper = h1)
  # With the COVID years separated, the 2021 step cannot be told from them:
  # only 2023 onwards would set it apart, and it took the rebound after
  # COVID as a jump of the population (pneumonia: +59%, the projection
  # falling by half). It is left out then, as in the backtests.
  extra <- c(if (steps) "step2014", if (steps && !covid) "step2021", if (covid) paste0("covid", 2020:2022))
  if (length(extra)) formula <- update(formula, stats::as.formula(paste(". ~ . +", paste(extra, collapse = " + "))))
  set.seed(2026)
  fit <- inla(formula, family = "poisson", data = d, offset = log(d$py), control.compute = list(config = TRUE),
              control.predictor = list(link = 1), num.threads = THREADS, verbose = FALSE)
  samples <- inla.posterior.sample(N_SAMPLES, fit, selection = list(Predictor = seq_len(nrow(d))), num.threads = THREADS)
  mu_all <- exp(do.call(cbind, lapply(samples, function(s) s$latent[, 1]))); rm(samples)
  noisy_all <- matrix(stats::rpois(length(mu_all), mu_all), nrow(mu_all))
  check <- sum(rowMeans(mu_all)[!is.na(d$y)]) / sum(d$y[!is.na(d$y)])
  d$keep <- d$year >= origin - 9L
  rows_of <- split(which(d$keep), d$area[d$keep])
  w <- weights[kept]
  result <- bind_rows(lapply(names(units), function(u) {
    r <- unlist(rows_of[units[[u]]], use.names = FALSE)
    if (length(r) == 0) return(NULL)
    mutate(summarise_unit(r, d, mu_all, noisy_all, w), level = level_of[[u]], area = u, .before = 1)
  }))
  saveRDS(list(cause = cause, sex = sex, combo = combo, origin = origin, horizon = HORIZON, kept = kept, check = check,
               terms = extra, fixed = fit$summary.fixed[intersect(extra, rownames(fit$summary.fixed)), c("mean", "0.025quant", "0.975quant"), drop = FALSE],
               hyper = fit$summary.hyperpar[, c("mean", "0.025quant", "0.975quant")], forecast = result,
               created = Sys.time()), path)
  cat(sprintf("%-60s %s %s: %6d rows, %4.0f s, expected/observed %.3f\n", cause, sex, combo, nrow(d),
              as.numeric(difftime(Sys.time(), t0, units = "secs")), check))
  rm(mu_all, noisy_all, fit); invisible(gc())
}

for (cause in causes) for (sex in c("HM", "H", "M")) {
  extract(cause, sex)
  if (extract_only) { cat(cause, sex, "extracted\n"); next }
  for (combo in combos) {
    res <- tryCatch(fit_one(cause, sex, combo), error = function(e) e)
    if (inherits(res, "error")) cat(cause, sex, combo, "FAILED:", conditionMessage(res), "\n")
  }
}
if (extract_only) quit(save = "no")
saveRDS(list(origin = last_year, first_year = INLA_FIRST_YEAR, horizon = HORIZON, data_mode = "corrected",
             population_basis = default_population_basis(), standard = "ESP 2013", created = Sys.time()),
        file.path(out_dir, "manifest.rds"))
cat("done:", length(list.files(out_dir, pattern = "__B[01]C[01]\\.rds$")), "forecasts in", out_dir, "\n")
