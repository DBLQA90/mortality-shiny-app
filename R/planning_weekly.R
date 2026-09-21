# =========================================================
# Planning indicators: weekly deaths and excess mortality
# =========================================================
# INE's weekly deaths by region of residence and age group
# (tools/fetch_weekly_deaths.R) - the most recent mortality the app has, a few
# weeks behind real time, against the deaths expected from recent years.
#
#   Regions   NUTS III, NUTS II, Continente and Portugal: INE publishes no
#             municipality or ULS. A location reads the smallest region that
#             holds it, which the tab states.
#   Expected  age-specific weekly death rates of the baseline years, applied to
#             the population of the year: an ageing or growing population is
#             not read as excess. Age groups under 65, 65-74, 75-84, 85+;
#             population from the annual estimates (the latest available for
#             years not yet estimated).
#   Baseline  up to five years before the year shown, from 2023: after the
#             COVID-19 excess of 2020-2022 and on the same population series
#             (INE revised the population from 2021, 1.7-5% up). Rates of
#             2018-2019 over the superseded population came out 4-7% above
#             those of 2023-2025 and made every recent year a "deficit"; so
#             2026 has three baseline years, 2025 two, and 2024 none. The
#             weekly counts of 2018-2020 (0010112, NUTS 2013, only for regions
#             whose counts agree between the editions within 0.5% in
#             2021-2024) are kept for the observed series.
#   Band      95% prediction interval: the between-year variability of the
#             baseline rates (at least Poisson), for a new year.
#
# Deaths of unknown age (0.01%) are left out on both sides.

WEEKLY_AGE_GROUPS <- c(
  "Todas as idades" = "all", "Menos de 65 anos" = "lt65", "65-74 anos" = "65_74", "75-84 anos" = "75_84", "85 e mais anos" = "85plus"
)
WEEKLY_EXCLUDED_YEARS <- 2020:2022
WEEKLY_FIRST_BASELINE_YEAR <- 2023L
WEEKLY_BASELINE_YEARS <- 5L

weekly_age_group <- function(band) {
  lower <- suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", band)))
  dplyr::case_when(
    is.na(lower) ~ NA_character_,
    lower < 65 ~ "lt65",
    lower < 75 ~ "65_74",
    lower < 85 ~ "75_84",
    TRUE ~ "85plus"
  )
}

planning_weekly_dir <- function() file.path(infant_snapshot_root(), "weekly_deaths")
planning_weekly_available <- function() file.exists(file.path(planning_weekly_dir(), "0012100.rds"))

# Weekly deaths by region, year, week and broad age group, from both editions
# (the older one only for regions it shares unchanged).
planning_weekly_data <- function() {
  key <- paste(infant_snapshot_root(), "weekly", sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  read <- function(indicator) {
    path <- file.path(planning_weekly_dir(), paste0(indicator, ".rds"))
    if (!file.exists(path)) return(NULL)
    x <- readRDS(path)
    x$region[x$region == "Alto Tâmega"] <- "Alto Tâmega e Barroso"
    x
  }
  current <- read("0012100")
  if (is.null(current)) {
    assign(key, NULL, envir = planning_cache)
    return(NULL)
  }
  older <- read("0010112")
  collapse <- function(x) {
    x$group <- weekly_age_group(x$age_group)
    x <- x[!is.na(x$group), , drop = FALSE]
    stats::aggregate(deaths ~ region + year + week + group, data = x, FUN = sum, na.action = stats::na.pass)
  }
  out <- collapse(current)
  if (!is.null(older)) {
    total <- function(x) {
      t <- x[x$age_group == "Total" & x$year %in% 2021:2024, , drop = FALSE]
      stats::aggregate(deaths ~ region + year, data = t, FUN = sum)
    }
    both <- merge(total(current), total(older), by = c("region", "year"), suffixes = c("_new", "_old"))
    agree <- tapply(abs(both$deaths_new / both$deaths_old - 1), both$region, max)
    same <- names(agree)[agree < 0.005]
    extra <- collapse(older[older$region %in% same & older$year < min(current$year), , drop = FALSE])
    out <- rbind(extra, out)
  }
  out <- out[order(out$region, out$group, out$year, out$week), ]
  assign(key, out, envir = planning_cache)
  out
}

# The weekly region for an area: itself when INE publishes it, otherwise the
# smallest NUTS region holding all its municipalities.
planning_weekly_region <- function(area, lookup = get_nuts_lookup()) {
  data <- planning_weekly_data()
  if (is.null(data)) return(NULL)
  regions <- unique(data$region)
  if (area == PLANNING_PORTUGAL_MUNICIPAL) area <- "Portugal"
  if (area %in% regions) return(list(region = area, note = NULL))
  members <- planning_area_members(area, lookup)
  rows <- lookup[lookup$municipality %in% members, , drop = FALSE]
  for (column in c("nuts3", "nuts2", "nuts1")) {
    value <- unique(rows[[column]])
    if (length(value) == 1 && value %in% regions) {
      return(list(region = value, note = paste0("O INE publica os óbitos semanais por NUTS III: ", area, " mostra ", value, ".")))
    }
  }
  list(region = "Portugal", note = paste0("O INE publica os óbitos semanais por NUTS III: ", area, " atravessa várias regiões e mostra Portugal."))
}

# Population by year and broad age group for a weekly region.
planning_weekly_population <- function(region, years, lookup = get_nuts_lookup()) {
  available <- snapshot_years_for("population")
  dplyr::bind_rows(lapply(years, function(year) {
    use <- max(c(available[available <= year], min(available)))
    pop <- read_year_file("population", use)
    if (is.null(pop)) return(NULL)
    rows <- if (region %in% planning_published_areas && region %in% pop$area) {
      pop[pop$area == region & pop$sex == "HM", , drop = FALSE]
    } else {
      pop[pop$area %in% planning_area_members(region, lookup) & pop$sex == "HM", , drop = FALSE]
    }
    rows$group <- weekly_age_group(as.character(rows$age_band))
    agg <- stats::aggregate(pop ~ group, data = rows, FUN = sum)
    tibble::tibble(year = as.integer(year), group = agg$group, pop = agg$pop, estimate_year = as.integer(use))
  }))
}

# Observed and expected weekly deaths for one region and age selection, for
# the years in `years`. Expected is NA when fewer than two baseline years exist.
planning_weekly_excess <- function(region, years, age = "all", lookup = get_nuts_lookup()) {
  data <- planning_weekly_data()
  if (is.null(data)) return(tibble::tibble())
  groups <- if (identical(age, "all")) c("lt65", "65_74", "75_84", "85plus") else age
  rows <- data[data$region == region & data$group %in% groups, , drop = FALSE]
  if (nrow(rows) == 0) return(tibble::tibble())
  have <- sort(unique(rows$year))
  pop <- planning_weekly_population(region, sort(unique(c(have, years))), lookup)
  rows <- merge(rows, pop[, c("year", "group", "pop")], by = c("year", "group"))
  rows$rate <- rows$deaths / rows$pop

  dplyr::bind_rows(lapply(as.integer(years), function(year) {
    observed <- rows[rows$year == year, , drop = FALSE]
    if (nrow(observed) == 0) return(NULL)
    baseline_years <- utils::tail(setdiff(have[have < year & have >= max(year - 8L, WEEKLY_FIRST_BASELINE_YEAR)], WEEKLY_EXCLUDED_YEARS), WEEKLY_BASELINE_YEARS)
    weeks <- sort(unique(observed$week))
    target_pop <- pop[pop$year == year, c("group", "pop")]
    per_group <- lapply(groups, function(g) {
      obs <- observed[observed$group == g, , drop = FALSE]
      base <- rows[rows$group == g & rows$year %in% baseline_years, , drop = FALSE]
      p <- target_pop$pop[target_pop$group == g]
      vapply(weeks, function(w) {
        # Week 53 exists only in some years: it borrows week 52.
        r <- base$rate[base$week == min(w, 52L)]
        o <- obs$deaths[obs$week == w]
        if (length(r) < 2 || length(p) == 0) return(c(if (length(o)) o else NA_real_, NA_real_, NA_real_))
        expected <- mean(r) * p
        variance <- max(stats::var(r) * p^2, expected) * (1 + 1 / length(r))
        c(if (length(o)) o else NA_real_, expected, variance)
      }, numeric(3))
    })
    parts <- Reduce(`+`, per_group)
    tibble::tibble(
      region = region, age = age, year = as.integer(year), week = weeks,
      observed = parts[1, ], expected = parts[2, ],
      lower = pmax(parts[2, ] - 1.96 * sqrt(parts[3, ]), 0), upper = parts[2, ] + 1.96 * sqrt(parts[3, ]),
      variance = parts[3, ], baseline = paste(baseline_years, collapse = ", "), baseline_n = length(baseline_years)
    )
  }))
}

# Year-to-date totals: observed, expected, excess and its 95% interval, over
# the weeks each year has so far.
planning_weekly_summary <- function(excess) {
  if (nrow(excess) == 0) return(tibble::tibble())
  excess %>%
    dplyr::filter(!is.na(.data$observed), !is.na(.data$expected)) %>%
    dplyr::group_by(.data$region, .data$age, .data$year, .data$baseline) %>%
    dplyr::summarise(
      weeks = dplyr::n(), last_week = max(.data$week),
      observed = sum(.data$observed), expected = sum(.data$expected), sd = sqrt(sum(.data$variance)), .groups = "drop"
    ) %>%
    dplyr::mutate(
      excess = .data$observed - .data$expected,
      excess_lower = .data$excess - 1.96 * .data$sd, excess_upper = .data$excess + 1.96 * .data$sd,
      excess_pct = .data$excess / .data$expected * 100
    )
}
