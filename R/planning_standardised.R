# =========================================================
# Planning indicators: standardised, premature and avoidable mortality
# =========================================================
# Mortality compared on a common age structure, per triennium and for any area:
#
#   SMR                 observed deaths over those expected at the benchmark's
#                       (Portugal's) age-specific rates in the same triennium,
#                       x 100 - indirect standardisation
#   standardised rate   the rate the area would have with the European Standard
#                       Population 2013 - direct standardisation, per 100,000
#   premature           the standardised rate under 75
#   preventable,        under-75 standardised rates of the Eurostat/OECD (2019)
#   treatable           lists, as mapped to INE's shortlist in R/avoidable.R -
#                       a lower bound, see avoidable_unresolved_note()
#   YPLL                years of potential life lost before 70, per 100,000
#                       residents under 70
#
# All of it needs deaths by age and cause per municipality, which INE publishes
# incomplete: in 2023 the municipal age bands hold 97.9% of circulatory deaths,
# in 2014 as little as 65% of suicides, while Portugal's row is complete.
# Summed raw, every region and ULS would read 1-3% below Portugal and the large
# ones would come out "significantly lower" by construction. So, as for life
# expectancy, each municipality's count per cause is brought up to INE's
# complete all-ages total for that cause (death_totals), the missing deaths
# spread over ages with the profile of the deaths missing nationally (Portugal's
# row by age less the municipal sum), and values where more than 2% of the
# deaths were spread carry the flag ‡.
#
# Areas are matrix sums over municipalities, with Portugal and Continente
# reading INE's published rows (or Portugal as the municipal sum); person-years
# are the mean of consecutive end-of-year populations, as elsewhere.

# Cause measures: all causes, the 13 I45 groups, and the two avoidable groups
# (sums of their leaf causes, each completed on its own).
planning_standardised_measures <- function() {
  c(all = planning_all_causes, stats::setNames(PLANNING_CAUSE_GROUPS$cause, PLANNING_CAUSE_GROUPS$code),
    preventable = "preventable", treatable = "treatable")
}

STANDARDISED_INDICATORS <- tibble::tribble(
  ~id,                ~measure,      ~kind,    ~max_age,
  "smr_all",          "all",         "smr",    Inf,
  "dsr_all",          "all",         "dsr",    Inf,
  "dsr_premature",    "all",         "dsr",    75,
  "premature_deaths", "all",         "count",  75,
  "dsr_preventable",  "preventable", "dsr",    75,
  "dsr_treatable",    "treatable",   "dsr",    75,
  "avoidable_deaths", "avoidable",   "count",  75,
  "ypll_rate",        "all",         "ypll",   70
)
standardised_ids <- STANDARDISED_INDICATORS$id

PLANNING_YPLL_CUTOFF <- 70
# Deaths without a published age: flagged above 2%, and withheld above 25%,
# where the age distribution would be more assumption than measurement (1.6%
# of municipal triennia, nearly all before 1999 and around 2014).
PLANNING_SPREAD_FLAG <- 0.02
PLANNING_SPREAD_SUPPRESS <- 0.25

# Fill in the deaths a municipality has without a published age.
#
# Two things are known and both must hold. Each municipality knows how many
# deaths it is missing (its complete all-ages total less what its bands hold),
# and the country knows in which age bands they are missing (the national row
# by age less the sum of the municipal rows). Neither margin alone gives a
# usable answer: each municipality's own recorded profile reproduces neither
# (2012-2014 premature mortality of the municipal sum read 338.8 against
# Portugal's 350.5), and the national gap profile, which INE's suppression of
# small cells skews young, gave Alvito 19 deaths under 5 in 1997 - a rate 40
# times the national one. So the missing deaths start where the municipality's
# own population makes them likely (its population by age at the national
# age-specific rates, less what it records) and are then fitted to both
# margins by a few proportional passes.
planning_complete_by_age <- function(recorded, population, total, labels, iterations = 20L) {
  bands <- colnames(recorded)
  recorded_sum <- rowSums(recorded)
  missing <- pmax(ifelse(is.na(total), recorded_sum, total) - recorded_sum, 0)
  municipal <- !labels %in% planning_published_areas
  national <- recorded["Portugal", ]
  gap <- pmax(national - colSums(recorded[municipal, , drop = FALSE]), 0)
  national_rate <- ifelse(population["Portugal", ] > 0, national / population["Portugal", ], 0)
  fallback <- if (sum(gap) > 0) gap / sum(gap) else if (sum(national) > 0) national / sum(national) else rep(1 / length(bands), length(bands))

  expected <- population * matrix(national_rate, length(labels), length(bands), byrow = TRUE)
  deficit <- pmax(expected - recorded, 0)
  deficit_sum <- rowSums(deficit)
  profile <- deficit / ifelse(deficit_sum > 0, deficit_sum, 1)
  if (any(deficit_sum <= 0)) profile[deficit_sum <= 0, ] <- matrix(fallback, sum(deficit_sum <= 0), length(bands), byrow = TRUE)

  allocation <- profile[municipal, , drop = FALSE] * missing[municipal]
  target_bands <- if (sum(gap) > 0) gap / sum(gap) * sum(missing[municipal]) else NULL
  if (!is.null(target_bands) && sum(allocation) > 0) {
    for (pass in seq_len(iterations)) {
      columns <- colSums(allocation)
      allocation <- allocation * matrix(ifelse(columns > 0, target_bands / columns, 1), nrow(allocation), length(bands), byrow = TRUE)
      rows_now <- rowSums(allocation)
      allocation <- allocation * ifelse(rows_now > 0, missing[municipal] / rows_now, 1)
    }
  }
  completed <- recorded
  completed[municipal, ] <- recorded[municipal, , drop = FALSE] + allocation
  # The published rows are complete by age; top them up proportionally if not.
  published <- which(!municipal)
  own <- recorded[published, , drop = FALSE] / pmax(recorded_sum[published], 1)
  completed[published, ] <- recorded[published, , drop = FALSE] + own * missing[published]
  list(deaths = completed, spread = missing)
}

planning_death_cause_file <- function(year, cause) {
  for (indicator in c("0013166", "0008206")) {
    path <- file.path(infant_snapshot_root(), "deaths", indicator, paste0("year_", year),
                      paste0("cause_", planning_cause_file_token(cause), ".rds"))
    if (file.exists(path)) return(path)
  }
  NULL
}

# Per year and sex: completed deaths [label x band x measure], deaths spread
# [label x measure], end-of-year population [label x band] and infant deaths.
# Cached per year.
planning_cause_age_block <- function(year, municipalities) {
  key <- paste(infant_snapshot_root(), "causeage", year, length(municipalities), sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))

  measures <- planning_standardised_measures()
  leaf <- list(preventable = AVOIDABLE_PREVENTABLE, treatable = AVOIDABLE_TREATABLE)
  causes <- unique(c(measures[!names(measures) %in% names(leaf)], unlist(leaf)))
  pop <- read_year_file("population", year)
  totals <- read_death_totals_year(year)
  files <- lapply(stats::setNames(causes, causes), planning_death_cause_file, year = year)
  if (is.null(pop) || is.null(files[[planning_all_causes]])) {
    assign(key, NULL, envir = planning_cache)
    return(NULL)
  }
  frames <- lapply(files, function(f) if (is.null(f)) NULL else readRDS(f))
  infant <- read_year_file("infant_totals", year)
  if (is.null(infant)) {
    infant <- read_year_file("infant_deaths", year)
    if (!is.null(infant)) infant <- infant[infant$cause == planning_all_causes, , drop = FALSE]
  }

  labels <- unique(c(municipalities, planning_published_areas))
  bands <- age_levels
  # Fill the matrix by index: aggregate() sorts, and this runs once per cause,
  # sex and year.
  matrix_of <- function(frame, value, sex) {
    out <- matrix(0, length(labels), length(bands), dimnames = list(labels, bands))
    if (is.null(frame)) return(out)
    rows <- frame[frame$sex == sex & frame$area %in% labels & frame$age_band %in% bands, , drop = FALSE]
    if (nrow(rows) == 0) return(out)
    index <- (match(as.character(rows$age_band), bands) - 1L) * length(labels) + match(rows$area, labels)
    values <- rows[[value]]
    if (anyDuplicated(index)) {
      summed <- rowsum(values, index, reorder = FALSE)
      out[as.integer(rownames(summed))] <- summed[, 1]
    } else {
      out[index] <- values
    }
    out
  }

  by_sex <- lapply(c("HM", "H", "M"), function(sex) {
    population <- matrix_of(pop, "pop", sex)
    complete_cause <- function(cause) {
      recorded <- matrix_of(frames[[cause]], "deaths", sex)
      total <- stats::setNames(rep(NA_real_, length(labels)), labels)
      if (!is.null(totals)) {
        rows <- totals[totals$sex == sex & totals$cause == cause & totals$area %in% labels, , drop = FALSE]
        total[rows$area] <- rows$deaths
      }
      planning_complete_by_age(recorded, population, total, labels)
    }
    completed <- lapply(stats::setNames(causes, causes), complete_cause)
    deaths <- array(0, c(length(labels), length(bands), length(measures) + 1L),
                    dimnames = list(labels, bands, c(names(measures), "avoidable")))
    spread <- matrix(0, length(labels), length(measures) + 1L, dimnames = list(labels, c(names(measures), "avoidable")))
    for (m in names(measures)) {
      parts <- if (m %in% names(leaf)) leaf[[m]] else measures[[m]]
      for (cause in parts) {
        deaths[, , m] <- deaths[, , m] + completed[[cause]]$deaths
        spread[, m] <- spread[, m] + completed[[cause]]$spread
      }
    }
    deaths[, , "avoidable"] <- deaths[, , "preventable"] + deaths[, , "treatable"]
    spread[, "avoidable"] <- spread[, "preventable"] + spread[, "treatable"]

    infant_deaths <- stats::setNames(rep(0, length(labels)), labels)
    if (!is.null(infant)) {
      rows <- infant[infant$sex == sex & infant$area %in% labels, , drop = FALSE]
      if (nrow(rows) > 0) {
        agg <- tapply(rows$deaths, rows$area, sum)
        infant_deaths[intersect(names(agg), labels)] <- agg[intersect(names(agg), labels)]
      }
    }
    list(deaths = deaths, spread = spread, population = population, infant = infant_deaths,
         has_row = rowSums(matrix_of(frames[[planning_all_causes]], "deaths", sex)) > 0)
  })
  names(by_sex) <- c("HM", "H", "M")
  assign(key, by_sex, envir = planning_cache)
  by_sex
}

planning_standardised_years <- function() life_expectancy_years()

# Pooled triennium sums for `areas` (and the benchmark), by band: deaths per
# measure, person-years, deaths spread, infant deaths. NULL when a year is
# missing.
planning_standardised_pool <- function(areas, end_year, sex, lookup, mode = PLANNING_DEFAULT_SPLIT_MODE) {
  # 0/1 here: planning_band_product() applies the parish weights per age band.
  membership <- planning_membership_matrix(areas, lookup)
  areas <- rownames(membership)
  municipalities <- colnames(membership)
  window <- seq.int(end_year - 2L, end_year)
  blocks <- lapply(window, planning_cause_age_block, municipalities = municipalities)
  if (any(vapply(blocks, is.null, logical(1)))) return(NULL)
  previous <- lapply(window - 1L, planning_cause_age_block, municipalities = municipalities)

  measures <- dimnames(blocks[[1]][[sex]]$deaths)[[3]]
  bands <- age_levels
  deaths <- array(0, c(length(areas), length(bands), length(measures)), dimnames = list(areas, bands, measures))
  spread <- matrix(0, length(areas), length(measures), dimnames = list(areas, measures))
  person_years <- matrix(0, length(areas), length(bands), dimnames = list(areas, bands))
  infant <- stats::setNames(numeric(length(areas)), areas)
  published <- intersect(areas, planning_published_areas)

  for (j in seq_along(window)) {
    block <- blocks[[j]][[sex]]
    before <- if (is.null(previous[[j]])) block else previous[[j]][[sex]]
    for (m in measures) {
      summed <- planning_band_product(membership, block$deaths[municipalities, , m, drop = TRUE], mode,
                                      target = planning_parish_weights(areas, municipalities, "mortality", mode, window[[j]]))
      for (area in published) if (isTRUE(block$has_row[[area]])) summed[area, ] <- block$deaths[area, , m]
      deaths[, , m] <- deaths[, , m] + summed
    }
    s <- membership %*% block$spread[municipalities, , drop = FALSE]
    if (identical(mode, "parish")) {
      s <- apply(block$spread[municipalities, , drop = FALSE], 2, function(v) planning_weighted_sum(membership, v, "mortality", mode, window[[j]]))
      dimnames(s) <- list(rownames(membership), colnames(block$spread))
    }
    for (area in published) if (isTRUE(block$has_row[[area]])) s[area, ] <- block$spread[area, ]
    spread <- spread + s
    mid <- (block$population + before$population) / 2
    py <- planning_band_product(membership, mid, mode)
    for (area in published) if (sum(mid[area, ]) > 0) py[area, ] <- mid[area, ]
    person_years <- person_years + py
    inf <- planning_weighted_sum(membership, block$infant, "female_15_49", mode, window[[j]])
    for (area in published) if (block$infant[[area]] > 0) inf[match(area, areas)] <- block$infant[[area]]
    infant <- infant + inf
  }
  # Deaths up to 1998 sit with the parent of Odivelas, Trofa and Vizela.
  joint <- Reduce(`|`, lapply(window, function(y) planning_joint_split(membership, y, "deaths")))
  list(deaths = deaths, spread = spread, person_years = person_years, infant = infant, joint = stats::setNames(joint, areas))
}

# Directly standardised rate per 100,000 with Dobson's interval, for each row of
# `deaths` / `person_years` (areas x bands), over the bands in `keep`.
planning_dsr <- function(deaths, person_years, keep = rep(TRUE, ncol(deaths)), confidence = 0.95) {
  weights <- esp2013_df$stdpop[match(colnames(deaths), as.character(esp2013_df$age_band))]
  weights[!keep] <- 0
  w <- weights / sum(weights)
  rates <- ifelse(person_years > 0, deaths / person_years, NA_real_)
  usable <- rowSums(person_years[, keep, drop = FALSE] > 0) == sum(keep)
  value <- as.vector(rates[, keep, drop = FALSE] %*% w[keep])
  variance <- as.vector(ifelse(person_years > 0, deaths / person_years^2, 0)[, keep, drop = FALSE] %*% (w[keep]^2))
  observed <- rowSums(deaths[, keep, drop = FALSE])
  ci <- planning_poisson_ci(observed, confidence)
  # Dobson: the Poisson limits of the observed count, scaled by the DSR's
  # standard error per death.
  scale <- ifelse(observed > 0, sqrt(variance / observed), 0)
  lower <- value + scale * (ci$lower - observed)
  upper <- value + scale * (ci$upper - observed)
  value[!usable] <- lower[!usable] <- upper[!usable] <- NA_real_
  list(value = value * 1e5, lower = pmax(lower, 0) * 1e5, upper = upper * 1e5, observed = observed)
}

planning_band_upper <- function(bands) {
  lower <- planning_age_lower(bands)
  ifelse(grepl("e mais", bands), Inf, lower + 5)
}

# The standardised indicators for `areas` over the triennia ending in
# `end_years`, in the shape of planning_indicator_table().
planning_standardised_table <- function(areas, end_years, ids = standardised_ids, lookup = get_nuts_lookup(),
                                        benchmark = "Portugal", sex = "HM", mode = PLANNING_DEFAULT_SPLIT_MODE) {
  wanted <- STANDARDISED_INDICATORS[STANDARDISED_INDICATORS$id %in% ids, , drop = FALSE]
  areas <- unique(as.character(areas))
  all_areas <- unique(c(areas, benchmark))
  available <- planning_standardised_years()
  bands <- age_levels
  upper <- planning_band_upper(bands)
  midpoint <- (planning_age_lower(bands) + pmin(upper, PLANNING_YPLL_CUTOFF)) / 2

  rows <- list()
  for (end_year in as.integer(end_years)) {
    pool <- if (end_year %in% available) planning_standardised_pool(all_areas, end_year, sex, lookup, mode) else NULL
    for (i in seq_len(nrow(wanted))) {
      spec <- wanted[i, ]
      value <- lower <- upper_ci <- numerator <- denominator <- rep(NA_real_, length(areas))
      flag <- rep("", length(areas))
      if (!is.null(pool)) {
        idx <- match(areas, rownames(pool$person_years))
        d <- pool$deaths[, , spec$measure]
        py <- pool$person_years
        keep <- upper <= spec$max_age
        if (spec$kind == "dsr") {
          r <- planning_dsr(d, py, keep)
          value <- r$value[idx]; lower <- r$lower[idx]; upper_ci <- r$upper[idx]
          numerator <- r$observed[idx]; denominator <- rowSums(py[, keep, drop = FALSE])[idx]
        } else if (spec$kind == "count") {
          value <- numerator <- rowSums(d[, keep, drop = FALSE])[idx]
          ci <- planning_poisson_ci(value)
          lower <- ci$lower; upper_ci <- ci$upper
        } else if (spec$kind == "smr") {
          ref <- match(benchmark, rownames(py))
          ref_rates <- ifelse(py[ref, ] > 0, d[ref, ] / py[ref, ], 0)
          expected <- as.vector(py %*% ref_rates)
          observed <- rowSums(d)
          ci <- planning_poisson_ci(observed)
          value <- ifelse(expected > 0, observed / expected * 100, NA_real_)[idx]
          lower <- ifelse(expected > 0, ci$lower / expected * 100, NA_real_)[idx]
          upper_ci <- ifelse(expected > 0, ci$upper / expected * 100, NA_real_)[idx]
          numerator <- observed[idx]; denominator <- expected[idx]
        } else if (spec$kind == "ypll") {
          under <- upper <= PLANNING_YPLL_CUTOFF
          lost <- d[, under, drop = FALSE] %*% (PLANNING_YPLL_CUTOFF - midpoint[under])
          # Infant deaths sit at the bottom of 0-4: weighted at 0.5, not 2.5
          # (see split_infant_age_band()).
          infant <- pmin(pool$infant, d[, 1])
          lost <- as.vector(lost) + infant * ((PLANNING_YPLL_CUTOFF - 0.5) - (PLANNING_YPLL_CUTOFF - midpoint[[1]]))
          pop_under <- rowSums(py[, under, drop = FALSE])
          value <- ifelse(pop_under > 0, lost / pop_under * 1e5, NA_real_)[idx]
          numerator <- lost[idx]; denominator <- pop_under[idx]
        }
        spread_share <- pool$spread[idx, if (spec$measure %in% colnames(pool$spread)) spec$measure else "all"] /
          pmax(rowSums(d)[idx], 1)
        flag <- ifelse(is.finite(value) & spread_share > PLANNING_SPREAD_FLAG, "‡", "")
        withheld <- !is.na(spread_share) & spread_share > PLANNING_SPREAD_SUPPRESS
        value[withheld] <- lower[withheld] <- upper_ci[withheld] <- NA_real_
        flag[withheld] <- ""
        joint <- pool$joint[idx]
        value[joint] <- lower[joint] <- upper_ci[joint] <- NA_real_
      }
      rows[[length(rows) + 1L]] <- tibble::tibble(
        area = areas, year = end_year, indicator = spec$id, value = unname(value), lower = unname(lower), upper = unname(upper_ci),
        numerator = unname(numerator), denominator = unname(denominator), flag = unname(flag)
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Mortality by cause group for one triennium: observed, expected, SMR against
# the benchmark, and the standardised rates at all ages and under 75.
planning_cause_standardised <- function(areas, end_year, lookup = get_nuts_lookup(), benchmark = "Portugal", sex = "HM",
                                        mode = PLANNING_DEFAULT_SPLIT_MODE) {
  areas <- unique(as.character(areas))
  pool <- planning_standardised_pool(unique(c(areas, benchmark)), as.integer(end_year), sex, lookup, mode)
  if (is.null(pool)) return(tibble::tibble())
  measures <- planning_standardised_measures()
  measures <- measures[!names(measures) %in% c("preventable", "treatable")]
  labels <- c(all = "Todas as causas de morte", stats::setNames(PLANNING_CAUSE_GROUPS$label, PLANNING_CAUSE_GROUPS$code))
  upper <- planning_band_upper(age_levels)
  py <- pool$person_years
  ref <- match(benchmark, rownames(py))
  idx <- match(areas, rownames(py))
  dplyr::bind_rows(lapply(names(measures), function(m) {
    d <- pool$deaths[, , m]
    ref_rates <- ifelse(py[ref, ] > 0, d[ref, ] / py[ref, ], 0)
    expected <- as.vector(py %*% ref_rates)
    observed <- rowSums(d)
    ci <- planning_poisson_ci(observed)
    all_ages <- planning_dsr(d, py)
    under75 <- planning_dsr(d, py, upper <= 75)
    spread_share <- pool$spread[, m] / pmax(observed, 1)
    withheld <- !is.na(spread_share[idx]) & spread_share[idx] > PLANNING_SPREAD_SUPPRESS
    out <- tibble::tibble(
      area = areas, code = m, group = unname(labels[[m]]), period = paste0(end_year - 2L, "-", end_year), end_year = as.integer(end_year),
      observed = unname(observed[idx]), expected = unname(expected[idx]),
      smr = ifelse(expected > 0, observed / expected * 100, NA_real_)[idx],
      smr_lower = ifelse(expected > 0, ci$lower / expected * 100, NA_real_)[idx],
      smr_upper = ifelse(expected > 0, ci$upper / expected * 100, NA_real_)[idx],
      dsr = all_ages$value[idx], dsr_lower = all_ages$lower[idx], dsr_upper = all_ages$upper[idx],
      dsr75 = under75$value[idx], dsr75_lower = under75$lower[idx], dsr75_upper = under75$upper[idx],
      flag = ifelse(spread_share[idx] > PLANNING_SPREAD_FLAG & !withheld, "‡", "")
    )
    out[withheld, setdiff(names(out), c("area", "code", "group", "period", "end_year", "flag"))] <- NA_real_
    joint <- pool$joint[idx]
    out[joint, setdiff(names(out), c("area", "code", "group", "period", "end_year", "flag"))] <- NA_real_
    out$significance <- planning_significance(out$smr_lower, out$smr_upper, rep(100, nrow(out)))
    out$significance[out$area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL)] <- NA_character_
    out
  }))
}
