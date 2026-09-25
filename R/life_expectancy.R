# =========================================================
# Life expectancy at birth (planning indicator I10)
# =========================================================
# An abridged period life table (Chiang II) over a triennium, for any area the
# app builds, by sex. It follows the method and variance of
# PHEindicatormethods::phe_life_expectancy() - Chiang's formulas, Silcocks'
# adjustment for the open interval, suppression when the population is 5,000 or
# less or the 95% interval spans more than 20 years - with one difference: the
# first age band is 0-4, because the death archive publishes no separate
# under-1 band with a matching population. The fraction of the interval lived
# by those dying in it (a) is taken from the known split of those deaths into
# under-1 (0.1 of a year lived) and ages 1-4 (2.5 years on average), so the
# band still carries the weight of infant mortality.
#
# Inputs, pooled over the three years:
#   deaths      all causes by five-year band, from the death archive. INE's
#               municipal age breakdown misses some deaths (4% in 2014); each
#               municipality's missing deaths are spread over ages in proportion
#               to its recorded ones, up to its complete total, and the share
#               spread is reported
#   population  mid-year person-years: the mean of the end-of-year estimates of
#               the previous year and the current one, summed over the three
#   infant      deaths under 1 year, for the first band's a
#
# Validation (METHODOLOGY.md): for Portugal the result matches Eurostat's
# standard life tables (2017-2019: 81.9 years here, 82.0 at Eurostat for 2019;
# at 65, 20.7 and 20.6). INE's own published values (Metodologia 2007, Portugal
# and NUTS III) are 0.8-0.9 years lower, a steady offset: across the 26 NUTS III
# of 2021-2023 the correlation with INE is 0.97 and the difference has a
# standard deviation of 0.3 years, also with INE's pre-revision population. The
# app's values are comparable with each other, not with INE's.

LIFE_OPEN_BAND <- "85 e mais anos"
LIFE_WIDTHS <- c(rep(5, length(age_levels) - 1), NA)
LIFE_AGE_STARTS <- as.integer(sub("^(\\d+).*$", "\\1", age_levels))
LIFE_MIN_POPULATION <- 5000
LIFE_REDISTRIBUTED_FLAG <- 0.02

# The six indicators the life table feeds: expectancy at birth and at 65, each
# for both sexes and for men and women.
LIFE_INDICATORS <- tibble::tribble(
  ~id,                        ~sex, ~age,
  "life_expectancy",          "HM", 0L,
  "life_expectancy_men",      "H",  0L,
  "life_expectancy_women",    "M",  0L,
  "life_expectancy_65",       "HM", 65L,
  "life_expectancy_65_men",   "H",  65L,
  "life_expectancy_65_women", "M",  65L
)

life_expectancy_ids <- LIFE_INDICATORS$id

# One life table. `deaths` and `population` are per band, youngest first; `a` is
# the fraction of each closed interval lived by those dying in it; the last band
# is open. Returns e0 and its interval, or NA with a reason.
# The whole table in one pass: life expectancy and its standard error at every
# age band, or a reason why it cannot be built.
abridged_life_table <- function(deaths, population, widths, a) {
  k <- length(deaths)
  out <- function(reason) list(e = rep(NA_real_, k), se = rep(NA_real_, k), reason = reason)
  if (any(is.na(deaths)) || any(is.na(population))) return(out("sem dados"))
  if (any(population <= 0)) return(out("população nula num grupo etário"))
  if (sum(population) <= LIFE_MIN_POPULATION) return(out("população igual ou inferior a 5.000"))
  if (any(deaths > population)) return(out("mais óbitos do que população"))

  m <- deaths / population
  n <- widths
  n[k] <- 2 / m[k]
  q <- ifelse(deaths <= population / n / a, m * n / (1 + m * n * (1 - a)), 1)
  q[k] <- 1
  l <- numeric(k)
  l[1] <- 100000
  for (i in 2:k) l[i] <- l[i - 1] * (1 - q[i - 1])
  d <- l - c(l[-1], 0)
  L <- numeric(k)
  L[-k] <- n[-k] * (l[-1] + a[-k] * d[-k])
  L[k] <- l[k] / m[k]
  T <- rev(cumsum(rev(L)))
  e <- ifelse(l == 0, 0, T / l)

  variance_q <- ifelse(d == 0, 0, q^2 * (1 - q) / deaths)
  variance_q[k] <- 4 / (deaths[k] * m[k]^2)
  weighted <- numeric(k)
  weighted[-k] <- variance_q[-k] * l[-k]^2 * ((1 - a[-k]) * n[-k] + e[-1])^2
  weighted[k] <- (l[k] / 2)^2 * variance_q[k]
  # The variance of e(x) accumulates only the bands from x on.
  se <- sqrt(rev(cumsum(rev(weighted))) / l^2)

  list(e = e, se = se, reason = "")
}

# Life expectancy at one age, with its interval. Suppressed, like PHE, when the
# 95% interval spans more than 20 years.
abridged_life_expectancy <- function(deaths, population, widths, a, confidence = 0.95, at_age = 0) {
  out <- function(value, lower = NA_real_, upper = NA_real_, reason = "") {
    list(value = value, lower = lower, upper = upper, reason = reason)
  }
  table <- abridged_life_table(deaths, population, widths, a)
  if (nzchar(table$reason)) return(out(NA_real_, reason = table$reason))

  index <- match(at_age, LIFE_AGE_STARTS)
  if (is.na(index)) stop("No age band starts at ", at_age, call. = FALSE)
  value <- table$e[[index]]
  se <- table$se[[index]]
  z <- stats::qnorm(confidence + (1 - confidence) / 2)

  if (!is.finite(value)) return(out(NA_real_, reason = "valor não finito"))
  if (stats::qnorm(0.975) * se > 10) return(out(NA_real_, reason = "intervalo de confiança superior a 20 anos"))
  out(value, value - z * se, value + z * se)
}

# ---------------------------------------------------------
# Data
# ---------------------------------------------------------
life_death_file <- function(year) {
  root <- file.path(infant_snapshot_root(), "deaths")
  for (indicator in c("0013166", "0008206")) {
    path <- file.path(root, indicator, paste0("year_", year), "cause_todas_as_causas_de_morte.rds")
    if (file.exists(path)) return(path)
  }
  NULL
}

life_death_years <- function() {
  root <- file.path(infant_snapshot_root(), "deaths")
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  dirs <- dirs[grepl("year_\\d{4}$", dirs) & file.exists(file.path(dirs, "cause_todas_as_causas_de_morte.rds"))]
  sort(unique(as.integer(sub(".*year_(\\d{4})$", "\\1", dirs))))
}

life_expectancy_years <- function() {
  deaths <- life_death_years()
  population <- snapshot_years_for("population")
  complete <- intersect(deaths, population)
  sort(complete[vapply(complete, function(y) all(c(y - 1L, y - 2L) %in% complete), logical(1))])
}

# Per year and sex: deaths by band per area label, with each municipality's
# unrecorded deaths spread over its ages; the deaths spread; end-of-year
# population by band; infant deaths. Cached per year.
life_year_block <- function(year, municipalities) {
  key <- paste(infant_snapshot_root(), "life", year, length(municipalities), sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))

  path <- life_death_file(year)
  pop <- read_year_file("population", year)
  totals <- read_death_totals_year(year)
  if (is.null(path) || is.null(pop)) {
    assign(key, NULL, envir = planning_cache)
    return(NULL)
  }
  deaths <- readRDS(path)
  infant <- read_year_file("infant_totals", year)
  if (is.null(infant)) {
    infant <- read_year_file("infant_deaths", year)
    if (!is.null(infant)) infant <- infant[infant$cause == planning_all_causes, , drop = FALSE]
  }

  labels <- unique(c(municipalities, planning_published_areas))
  bands <- age_levels
  by_sex <- lapply(c("HM", "H", "M"), function(sex) {
    # Fill by index rather than aggregate(), which sorts.
    matrix_of <- function(frame, value) {
      out <- matrix(0, length(labels), length(bands), dimnames = list(labels, bands))
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
    recorded <- matrix_of(deaths, "deaths")
    population <- matrix_of(pop, "pop")

    total <- stats::setNames(rep(NA_real_, length(labels)), labels)
    if (!is.null(totals)) {
      rows <- totals[totals$sex == sex & totals$cause == planning_all_causes & totals$area %in% labels, , drop = FALSE]
      total[rows$area] <- rows$deaths
    }
    recorded_sum <- rowSums(recorded)
    # Deaths without a published age, spread over the ages where the
    # municipality's own population and the country's missing ages both put
    # them; see planning_complete_by_age() in R/planning_standardised.R.
    filled <- planning_complete_by_age(recorded, population, total, labels)
    completed <- filled$deaths
    missing <- filled$spread

    infant_deaths <- stats::setNames(rep(0, length(labels)), labels)
    if (!is.null(infant)) {
      rows <- infant[infant$sex == sex & infant$area %in% labels, , drop = FALSE]
      if (nrow(rows) > 0) {
        agg <- tapply(rows$deaths, rows$area, sum)
        infant_deaths[names(agg)] <- agg
      }
    }
    list(deaths = completed, spread = missing, population = population, infant = infant_deaths,
         has_row = recorded_sum > 0, has_population = rowSums(population) > 0)
  })
  names(by_sex) <- c("HM", "H", "M")
  assign(key, by_sex, envir = planning_cache)
  by_sex
}

# Life expectancy for `areas` over the triennia ending in `end_years`, by sex.
planning_life_expectancy_table <- function(areas, end_years, ids = LIFE_INDICATORS$id, lookup = get_nuts_lookup(),
                                           mode = PLANNING_DEFAULT_SPLIT_MODE, vintage = planning_lookup_vintage(lookup)) {
  wanted <- LIFE_INDICATORS[LIFE_INDICATORS$id %in% ids, , drop = FALSE]
  sexes <- unique(wanted$sex)
  # 0/1 here: the parish weights are applied per age band below, and applying
  # them to the membership as well would count them twice.
  membership <- planning_membership_matrix(areas, lookup)
  areas <- rownames(membership)
  municipalities <- colnames(membership)
  available <- life_expectancy_years()

  rows <- list()
  for (end_year in as.integer(end_years)) {
    window <- seq.int(end_year - 2L, end_year)
    if (!end_year %in% available) {
      for (id in wanted$id) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          area = areas, year = end_year, indicator = id,
          value = NA_real_, lower = NA_real_, upper = NA_real_, numerator = NA_real_, denominator = NA_real_, flag = ""
        )
      }
      next
    }
    blocks <- lapply(window, life_year_block, municipalities = municipalities)
    previous <- lapply(window - 1L, life_year_block, municipalities = municipalities)

    for (sex in sexes) {
      deaths <- matrix(0, length(areas), length(age_levels))
      person_years <- matrix(0, length(areas), length(age_levels))
      spread <- infant <- numeric(length(areas))

      for (j in seq_along(window)) {
        block <- blocks[[j]][[sex]]
        before <- if (is.null(previous[[j]])) block else previous[[j]][[sex]]
        sum_areas <- function(values, own_ok, target = NULL) {
          summed <- planning_band_product(membership, values, mode, target = target)
          for (area in intersect(areas, planning_published_areas)) {
            if (isTRUE(own_ok[[area]])) summed[area, ] <- values[area, ]
          }
          summed
        }
        # Where INE publishes the area's own rows by age, they replace the sum
        # of its municipalities: the bands of a regional row add up to its
        # total, so none of its deaths has to be spread over ages.
        from_rows <- stats::setNames(lapply(areas, function(area) {
          planning_regional_area_deaths(area, window[[j]], sex, block, vintage, lookup, measure = "all")
        }), areas)
        summed <- sum_areas(block$deaths, block$has_row, target = planning_parish_weights(areas, municipalities, "mortality", mode, window[[j]]))
        for (area in areas) if (!is.null(from_rows[[area]])) summed[match(area, areas), ] <- from_rows[[area]]
        deaths <- deaths + summed
        mid <- (block$population + before$population) / 2
        person_years <- person_years + sum_areas(mid, block$has_population & before$has_population)
        spread_vec <- planning_weighted_sum(membership, block$spread, "mortality", mode, window[[j]])
        infant_vec <- planning_weighted_sum(membership, block$infant, "female_15_49", mode, window[[j]])
        for (area in intersect(areas, planning_published_areas)) {
          i <- match(area, areas)
          if (isTRUE(block$has_row[[area]])) spread_vec[i] <- block$spread[[area]]
          if (block$infant[[area]] > 0) infant_vec[i] <- block$infant[[area]]
        }
        for (area in areas) if (!is.null(from_rows[[area]])) spread_vec[match(area, areas)] <- 0
        spread <- spread + spread_vec
        infant <- infant + infant_vec
      }

      # One table per area, read at each age the caller asked for.
      tables <- lapply(seq_along(areas), function(i) {
        first <- deaths[i, 1]
        infant_share <- if (first > 0) min(infant[i], first) / first else 1
        a_first <- (infant_share * 0.1 + (1 - infant_share) * 2.5) / 5
        a <- c(a_first, rep(0.5, length(age_levels) - 1))
        abridged_life_table(deaths[i, ], person_years[i, ], LIFE_WIDTHS, a)
      })
      total_deaths <- rowSums(deaths)
      z <- stats::qnorm(0.975)

      for (age in unique(wanted$age[wanted$sex == sex])) {
        index <- match(age, LIFE_AGE_STARTS)
        value <- vapply(tables, function(t) t$e[[index]], numeric(1))
        se <- vapply(tables, function(t) t$se[[index]], numeric(1))
        reason <- vapply(tables, `[[`, character(1), "reason")
        too_wide <- !is.na(se) & z * se > 10
        reason[too_wide] <- "intervalo de confiança superior a 20 anos"
        # Where most deaths have no published age the life table would be
        # mostly assumption: municipalities before 1999 and around 2014.
        spread_share <- ifelse(total_deaths > 0, spread / total_deaths, 0)
        withheld <- spread_share > PLANNING_SPREAD_SUPPRESS
        reason[withheld] <- "mais de 25% dos óbitos do triénio sem idade publicada por município"
        value[withheld] <- NA_real_
        value[too_wide | !is.finite(value)] <- NA_real_
        # Deaths up to 1998 sit with the parent of Odivelas, Trofa and Vizela.
        joint <- Reduce(`|`, lapply(window, function(y) planning_joint_split(membership, y, "deaths")))
        reason[joint] <- "óbitos registados no município de origem (criado em 1998)"
        value[joint] <- NA_real_

        rows[[length(rows) + 1L]] <- tibble::tibble(
          area = areas, year = end_year,
          indicator = wanted$id[wanted$sex == sex & wanted$age == age],
          value = value,
          lower = ifelse(is.na(value), NA_real_, value - z * se),
          upper = ifelse(is.na(value), NA_real_, value + z * se),
          numerator = unname(total_deaths),
          denominator = unname(rowSums(person_years)),
          flag = unname(ifelse(!is.na(value) & total_deaths > 0 & spread / total_deaths > LIFE_REDISTRIBUTED_FLAG, "‡", "")),
          reason = reason
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}
