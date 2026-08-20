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
    # The two halves of `0 - 4 anos`, produced by split_infant_age_band(). Both
    # are stated explicitly because the generic rule below would read
    # "1 - 4 anos" as (1 + 4) / 2 = 2.5, which is the midpoint of the wrong
    # interval: the band covers exact ages 1 up to 5, so its midpoint is 3.
    age_band == "< 1 ano" ~ 0.5,
    age_band == "1 - 4 anos" ~ 3,
    grepl("^[0-9]+\\s*-\\s*[0-9]+", age_band) ~ {
      lower <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", age_band)))
      upper <- suppressWarnings(as.numeric(sub("^[0-9]+\\s*-\\s*([0-9]+).*", "\\1", age_band)))
      (lower + upper) / 2
    },
    TRUE ~ NA_real_
  )
}

# Split `0 - 4 anos` into `< 1 ano` and `1 - 4 anos` using a known count of
# infant deaths.
#
# Years of potential life lost weight each death by how far short of the cutoff
# it fell, taken from the midpoint of its age band. Infant deaths are the
# extreme case: they sit at the very bottom of a five-year band, so the band
# midpoint of 2.5 credits each one with 67.5 lost years against a cutoff of 70
# when the true figure is close to 70. And they are not a minor part of the
# band - for Portugal in 2024, 254 of the 286 deaths in `0 - 4 anos` were
# infants - so the band's weight is wrong for most of the deaths in it.
#
# The under-1 counts fetched for infant mortality make the correction possible:
# subtract them out and the remainder sits at ages one to four, where a midpoint
# of 3 is right. Deaths are unchanged; only their weights move.
#
# Applied only when the count is known for the whole window. With NA infant
# deaths the frame is returned untouched, so AVPP falls back to its previous
# behaviour rather than to a half-corrected figure.
split_infant_age_band <- function(df, infant_deaths, band = "0 - 4 anos") {
  if (!is.finite(infant_deaths) || !(band %in% as.character(df$age_band))) {
    return(df)
  }

  # The new labels are not levels of whatever factor the caller supplied.
  df$age_band <- as.character(df$age_band)

  # More than one row for the band (an unpooled frame still split by year or
  # area) collapses into the two new rows, which is what AVPP sums anyway.
  band_row <- df[df$age_band == band, , drop = FALSE]
  band_deaths <- sum(band_row$deaths, na.rm = TRUE)

  # The two counts come from different indicators, and a vintage revision could
  # in principle leave more infant deaths on record than the band holds. Cap
  # rather than emit a negative count for ages one to four.
  under_one <- max(min(infant_deaths, band_deaths), 0)

  rest <- df[df$age_band != band, , drop = FALSE]
  split_rows <- band_row[c(1, 1), , drop = FALSE]
  split_rows$age_band <- c("< 1 ano", "1 - 4 anos")
  split_rows$deaths <- c(under_one, band_deaths - under_one)

  # Population is a property of the original band and does not divide with the
  # deaths. No rate is computed from a split frame - the correction exists for
  # AVPP, which uses deaths alone - so mark it unknown rather than invent a
  # share of it.
  if ("pop" %in% names(split_rows)) {
    split_rows$pop <- NA_real_
  }

  dplyr::bind_rows(rest, split_rows)
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
