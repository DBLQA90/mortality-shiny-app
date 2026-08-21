# Tests for R/metrics.R (mortality-rate, DSR, YPLL and interval helpers).

as_bands <- function(x) factor(x, levels = age_levels, ordered = TRUE)

test_that("get_age_midpoint uses band midpoints and drops open/total bands", {
  expect_equal(get_age_midpoint("0 - 4 anos"), 2.5)
  expect_equal(get_age_midpoint("30 - 34 anos"), 32)
  expect_equal(get_age_midpoint("70 - 74 anos"), 72)
  expect_true(is.na(get_age_midpoint("85 e mais anos")))
  expect_true(is.na(get_age_midpoint("Total")))
})

test_that("compute_metrics crude rate and Poisson CI match the direct formula", {
  df <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(c("30 - 34 anos", "60 - 64 anos")),
    deaths = c(20, 130), pop = c(50000, 40000)
  )
  m <- compute_metrics(df)

  deaths_total <- 150; pop_total <- 90000
  expect_equal(m$crude_rate, deaths_total / pop_total * 1e5)

  exp_ci <- stats::poisson.test(deaths_total)$conf.int / pop_total * 1e5
  expect_equal(m$crude_lower, exp_ci[[1]])
  expect_equal(m$crude_upper, exp_ci[[2]])
})

test_that("all-age DSR equals the ESP-2013 weighted average of age-specific rates", {
  # Constant 100/100k rate in every band -> DSR must be exactly 100 per 100k.
  df <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(age_levels),
    pop = rep(100000, length(age_levels)),
    deaths = rep(100, length(age_levels))
  )
  m <- compute_metrics(df)
  expect_equal(m$dsr, 100, tolerance = 1e-8)
})

test_that("under-75 DSR keeps the conventional 0-74 renormalisation (not rescaled)", {
  set.seed(42)
  full <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(age_levels),
    deaths = round(runif(length(age_levels), 5, 60)),
    pop = round(runif(length(age_levels), 8000, 20000))
  )
  under <- dplyr::filter(full, !age_band %in% c("75 - 79 anos", "80 - 84 anos", "85 e mais anos"))

  got <- compute_metrics(under)

  renorm <- under %>%
    dplyr::left_join(esp2013_df, by = "age_band") %>%
    dplyr::group_by(year, sex, cause) %>%
    PHEindicatormethods::calculate_dsr(x = deaths, n = pop, stdpop = stdpop, multiplier = 1e5, confidence = 0.95) %>%
    dplyr::ungroup()

  # Conventional premature ASMR: normalised to the 0-74 sub-population.
  expect_equal(got$dsr, renorm$value, tolerance = 1e-8)
  # Guard against re-introducing the full-standard rescaling (x 0.91).
  rescaled <- renorm$value * (sum(esp2013_df$stdpop[1:15]) / sum(esp2013_df$stdpop))
  expect_false(isTRUE(all.equal(got$dsr, rescaled)))
})

test_that("get_age_midpoint places the split infant bands correctly", {
  expect_equal(get_age_midpoint("< 1 ano"), 0.5)
  # Not 2.5: the band covers exact ages 1 up to 5, so the generic
  # (lower + upper) / 2 rule would read the wrong interval.
  expect_equal(get_age_midpoint("1 - 4 anos"), 3)
})

test_that("split_infant_age_band moves weight without moving deaths", {
  df <- tibble(
    age_band = c("0 - 4 anos", "30 - 34 anos"),
    deaths = c(286, 10)
  )
  # Portugal 2024: 254 of the 286 deaths in the band were infants.
  split <- split_infant_age_band(df, 254)

  expect_setequal(split$age_band, c("< 1 ano", "1 - 4 anos", "30 - 34 anos"))
  expect_equal(sum(split$deaths), sum(df$deaths))
  expect_equal(split$deaths[split$age_band == "< 1 ano"], 254)
  expect_equal(split$deaths[split$age_band == "1 - 4 anos"], 32)

  # The point of the split: infants lose nearly the whole cutoff, so AVPP rises.
  expect_equal(compute_ypll(df, cutoff = 70), 286 * 67.5 + 10 * 38)
  expect_equal(compute_ypll(split, cutoff = 70), 254 * 69.5 + 32 * 67 + 10 * 38)
  expect_gt(compute_ypll(split, cutoff = 70), compute_ypll(df, cutoff = 70))
})

test_that("split_infant_age_band declines rather than guessing", {
  df <- tibble(age_band = c("0 - 4 anos", "30 - 34 anos"), deaths = c(286, 10))

  # An unknown count leaves AVPP exactly as it was, rather than half-corrected.
  expect_equal(split_infant_age_band(df, NA_real_), df)

  # Nothing to split.
  no_band <- tibble(age_band = "30 - 34 anos", deaths = 10)
  expect_equal(split_infant_age_band(no_band, 254), no_band)

  # The two counts come from different indicators. More infant deaths than the
  # band holds must cap, never produce a negative count for ages one to four.
  capped <- split_infant_age_band(df, 400)
  expect_equal(capped$deaths[capped$age_band == "< 1 ano"], 286)
  expect_equal(capped$deaths[capped$age_band == "1 - 4 anos"], 0)
  expect_equal(sum(capped$deaths), sum(df$deaths))
})

test_that("split_infant_age_band handles factors and repeated band rows", {
  # An ordered factor of the original levels cannot hold the new labels.
  df <- tibble(age_band = as_bands("0 - 4 anos"), deaths = 100)
  split <- split_infant_age_band(df, 80)
  expect_equal(split$deaths[split$age_band == "< 1 ano"], 80)

  # An unpooled frame still split by year collapses into the two new rows.
  multi <- tibble(age_band = c("0 - 4 anos", "0 - 4 anos"), deaths = c(60, 40))
  collapsed <- split_infant_age_band(multi, 70)
  expect_equal(nrow(collapsed), 2)
  expect_equal(sum(collapsed$deaths), 100)
  expect_equal(collapsed$deaths[collapsed$age_band == "1 - 4 anos"], 30)

  # Population does not divide with the deaths, so it is marked unknown rather
  # than apportioned. No rate is computed from a split frame.
  with_pop <- tibble(age_band = "0 - 4 anos", deaths = 100, pop = 50000)
  expect_true(all(is.na(split_infant_age_band(with_pop, 80)$pop)))
})

test_that("compute_ypll sums deaths x years-lost with a 70-year cutoff", {
  df <- tibble(
    age_band = as_bands(c("0 - 4 anos", "30 - 34 anos", "60 - 64 anos", "80 - 84 anos")),
    deaths = c(3, 10, 25, 40)
  )
  # years lost: 67.5, 38, 8, 0 (>=70 contributes nothing)
  expect_equal(compute_ypll(df, cutoff = 70), 3 * 67.5 + 10 * 38 + 25 * 8)
})

test_that("compute_ypll_interval uses the Dobson weighted-Poisson method", {
  df <- tibble(
    age_band = as_bands(c("0 - 4 anos", "30 - 34 anos", "60 - 64 anos")),
    deaths = c(3, 10, 25)
  )
  ci <- compute_ypll_interval(df, cutoff = 70)

  yl <- c(67.5, 38, 8)
  S <- sum(df$deaths * yl); V <- sum(df$deaths * yl^2); O <- sum(df$deaths)
  pois <- stats::poisson.test(round(O))$conf.int
  scale <- sqrt(V / round(O))
  expect_equal(unname(ci[["estimate"]]), S)
  expect_equal(unname(ci[["lower"]]), max(S + scale * (pois[[1]] - round(O)), 0))
  expect_equal(unname(ci[["upper"]]), S + scale * (pois[[2]] - round(O)))
  # Poisson skew -> asymmetric interval
  expect_gt(ci[["upper"]] - ci[["estimate"]], ci[["estimate"]] - ci[["lower"]])
  expect_gte(ci[["lower"]], 0)
})

test_that("compute_ypll_interval is degenerate with no premature deaths", {
  df <- tibble(age_band = as_bands("80 - 84 anos"), deaths = 25)
  ci <- compute_ypll_interval(df, cutoff = 70)
  expect_equal(unname(ci[["estimate"]]), 0)
  expect_equal(unname(ci[["lower"]]), 0)
  expect_equal(unname(ci[["upper"]]), 0)
})

test_that("count and proportion intervals match the exact reference tests", {
  ci_count <- compute_count_interval(25)
  expect_equal(as.numeric(ci_count), as.numeric(stats::poisson.test(25)$conf.int))

  ci_prop <- compute_proportion_interval(30, 100)
  expect_equal(as.numeric(ci_prop), as.numeric(stats::binom.test(30, 100)$conf.int * 100))

  # numerator is clamped to the denominator
  ci_clamped <- compute_proportion_interval(120, 100)
  expect_equal(as.numeric(ci_clamped), as.numeric(stats::binom.test(100, 100)$conf.int * 100))

  # invalid denominator -> NA bounds
  expect_true(all(is.na(compute_proportion_interval(5, 0))))
})

test_that("population seams are announced when a window crosses them", {
  # Two handovers between population indicators: 0003182 -> 0008273 at 2014,
  # and 0008273 -> the revised 0012918 at 2021.
  expect_setequal(vapply(POPULATION_SEAMS, function(s) s$year, integer(1)), c(2014L, 2021L))

  msg <- population_revision_warning(2019:2022, "crude")
  expect_match(msg, "2020/2021")
  expect_match(msg, "denominador")
  expect_false(grepl("2013/2014", msg))

  # The 2013/2014 seam: measured at 3% on standardised rates, and the reason a
  # breakpoint lands there in Portugal's series.
  seam14 <- population_revision_warning(2010:2018, "dsr")
  expect_match(seam14, "2013/2014")
  expect_false(grepl("2020/2021", seam14))

  # A long series crosses both and both are named.
  both <- population_revision_warning(1991:2024, "dsr")
  expect_match(both, "2013/2014")
  expect_match(both, "2020/2021")
  expect_match(both, "2 mudanças")

  # Wholly on one side of a seam, nothing to say.
  expect_null(population_revision_warning(2015:2019, "crude"))
  expect_null(population_revision_warning(2021:2024, "crude"))
  expect_null(population_revision_warning(2005:2013, "dsr"))

  # Counts do not divide by population, so no seam reaches them.
  for (metric in c("deaths", "ypll", "infant", "infant_deaths", "proportional")) {
    expect_null(population_revision_warning(1991:2024, metric))
  }

  # Unspecified metric errs toward warning.
  expect_match(population_revision_warning(2019:2022), "2020/2021")

  # An empty or non-finite selection is not an error.
  expect_null(population_revision_warning(integer(0), "dsr"))
  expect_null(population_revision_warning(NA_integer_, "dsr"))
})

test_that("the revised population series outranks the older ones", {
  # get_source_year_plan() resolves overlaps by DESCENDING priority, so the
  # largest number wins. Listing the revised indicator first and giving it 1
  # made it lose every overlapping year to the series it was meant to replace -
  # invisible in snapshot mode, where the files had already been overwritten,
  # and wrong only when reading INE live.
  expect_equal(
    names(which.max(population_source_priorities)),
    population_indicator_revised
  )

  # The 2011-2013 overlap keeps its existing winner.
  expect_gt(
    population_source_priorities[[population_indicator_legacy]],
    population_source_priorities[[population_indicator_current]]
  )

  # Distinct, or the tie-break is arbitrary.
  expect_equal(length(unique(population_source_priorities)), 3L)
})
