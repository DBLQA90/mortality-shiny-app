# =========================================================
# Compute metrics (crude rates + DSR with proper Poisson CIs)
# =========================================================

compute_metrics <- function(df) {
  df_age <- df %>%
    dplyr::group_by(year, sex, cause, age_band) %>%
    dplyr::summarise(
      deaths = sum(deaths, na.rm = TRUE),
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    )

  # Crude rates and Poisson CIs
  crude <- df_age %>%
    dplyr::group_by(year, sex, cause) %>%
    dplyr::summarise(
      deaths_total = sum(deaths),
      pop_total    = sum(pop),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      crude_rate = dplyr::if_else(pop_total > 0, deaths_total / pop_total * 1e5, NA_real_),
      ci = purrr::map2(
        deaths_total, pop_total,
        ~ if (is.na(.y) || .y <= 0) c(NA_real_, NA_real_) else stats::poisson.test(.x)$conf.int / .y * 1e5
      ),
      crude_lower = purrr::map_dbl(ci, 1),
      crude_upper = purrr::map_dbl(ci, 2)
    ) %>%
    dplyr::select(-ci)

  # Age-standardised (DSR). calculate_dsr() normalises by the sum of the
  # supplied standard weights, so the "under 75" scope yields the conventional
  # premature-mortality rate standardised to the ESP-2013 0-74 sub-population.
  # That rate uses a different standard structure from the all-age rate, so the
  # two are not on the same standard; the app labels the under-75 standardised
  # rate accordingly (see get_rate_mapping) rather than rescaling it.
  dsr <- df_age %>%
    dplyr::left_join(esp2013_df, by = "age_band") %>%
    dplyr::group_by(year, sex, cause) %>%
    PHEindicatormethods::calculate_dsr(
      x          = deaths,
      n          = pop,
      stdpop     = stdpop,
      multiplier = 1e5,
      confidence = 0.95
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      dsr       = value,
      dsr_lower = lowercl,
      dsr_upper = uppercl
    )

  dplyr::left_join(crude, dsr, by = c("year", "sex", "cause"))
}

get_age_midpoint <- function(age_band) {
  age_band <- as.character(age_band)

  dplyr::case_when(
    age_band == "0 - 4 anos" ~ 2.5,
    grepl("^[0-9]+\\s*-\\s*[0-9]+", age_band) ~ {
      lower <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", age_band)))
      upper <- suppressWarnings(as.numeric(sub("^[0-9]+\\s*-\\s*([0-9]+).*", "\\1", age_band)))
      (lower + upper) / 2
    },
    TRUE ~ NA_real_
  )
}

compute_ypll <- function(df, cutoff = 70) {
  df %>%
    dplyr::mutate(
      age_midpoint = get_age_midpoint(age_band),
      years_lost = pmax(cutoff - age_midpoint, 0)
    ) %>%
    dplyr::filter(!is.na(years_lost)) %>%
    dplyr::summarise(AVPP = sum(deaths * years_lost, na.rm = TRUE)) %>%
    dplyr::pull(AVPP)
}

compute_ypll_interval <- function(df, cutoff = 70, confidence = 0.95) {
  # Only deaths before the cutoff contribute years of potential life lost, so
  # the interval is driven by the premature-death count.
  ypll_data <- df %>%
    dplyr::mutate(
      age_midpoint = get_age_midpoint(age_band),
      years_lost = pmax(cutoff - age_midpoint, 0)
    ) %>%
    dplyr::filter(!is.na(years_lost), years_lost > 0)

  estimate <- sum(ypll_data$deaths * ypll_data$years_lost, na.rm = TRUE)
  variance <- sum(ypll_data$deaths * ypll_data$years_lost^2, na.rm = TRUE)
  total_deaths <- sum(ypll_data$deaths, na.rm = TRUE)
  count <- round(total_deaths)

  if (!is.finite(estimate) || count <= 0 || !is.finite(variance) || variance <= 0) {
    return(c(
      estimate = if (is.finite(estimate)) estimate else 0,
      lower = 0,
      upper = if (is.finite(estimate)) estimate else 0
    ))
  }

  # Dobson et al. (1991) confidence limits for a weighted sum of Poisson
  # counts: scale the exact Poisson limits of the total premature-death count
  # by sqrt(variance / count). Unlike a plain normal approximation this is
  # asymmetric and stays sensible for sparse local counts.
  pois_ci <- stats::poisson.test(count, conf.level = confidence)$conf.int
  scale <- sqrt(variance / count)

  c(
    estimate = estimate,
    lower = max(estimate + scale * (pois_ci[[1]] - count), 0),
    upper = estimate + scale * (pois_ci[[2]] - count)
  )
}

compute_count_interval <- function(count, confidence = 0.95) {
  if (!is.finite(count) || count < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  stats::poisson.test(round(count), conf.level = confidence)$conf.int
}

compute_proportion_interval <- function(numerator, denominator, confidence = 0.95) {
  if (!is.finite(numerator) || !is.finite(denominator) || denominator <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  numerator <- round(max(min(numerator, denominator), 0))
  denominator <- round(denominator)
  stats::binom.test(numerator, denominator, conf.level = confidence)$conf.int * 100
}
