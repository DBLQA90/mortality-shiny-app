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
