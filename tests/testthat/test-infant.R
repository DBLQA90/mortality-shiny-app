# Tests for R/infant.R (infant mortality rate).

test_that("infant mortality is deaths per 1,000 live births", {
  res <- compute_infant_mortality(deaths = 254, births = 84642)

  expect_equal(res$value, 254 / 84642 * 1000)
  # Portugal 2024, which INE publishes as roughly 3.0 per 1,000.
  expect_equal(round(res$value, 1), 3.0)
  expect_equal(res$deaths, 254)
  expect_equal(res$births, 84642)
  expect_true(is.na(res$note))
})

test_that("the interval comes from the death count, scaled by births", {
  res <- compute_infant_mortality(deaths = 20, births = 5000)
  expected <- stats::poisson.test(20)$conf.int / 5000 * 1000

  expect_equal(res$lower, expected[[1]])
  expect_equal(res$upper, expected[[2]])
  expect_lt(res$lower, res$value)
  expect_gt(res$upper, res$value)

  # Births are treated as a fixed denominator: the same rate from ten times the
  # data must produce a materially tighter interval.
  big <- compute_infant_mortality(deaths = 200, births = 50000)
  expect_equal(big$value, res$value)
  expect_lt(big$upper - big$lower, res$upper - res$lower)
})

test_that("zero infant deaths gives a real zero rate with an upper bound", {
  res <- compute_infant_mortality(deaths = 0, births = 1200)

  expect_equal(res$value, 0)
  expect_equal(res$lower, 0)
  # A municipality with no infant deaths has not proved the risk is zero.
  expect_gt(res$upper, 0)
})

test_that("an absent denominator is refused rather than guessed", {
  no_births <- compute_infant_mortality(deaths = 3, births = 0)
  expect_true(is.na(no_births$value))
  expect_match(no_births$note, "nados-vivos")

  missing <- compute_infant_mortality(deaths = NA_real_, births = 1000)
  expect_true(is.na(missing$value))
  expect_match(missing$note, "Sem dados")

  # Negative or non-finite denominators must not slip through as a rate.
  expect_true(is.na(compute_infant_mortality(deaths = 3, births = -5)$value))
})

# The test runner's working directory is tests/testthat, so point the snapshot
# reader at the repository copy explicitly rather than relying on app_dir.
with_snapshots <- function(code) {
  withr::with_envvar(
    c(MORTALITY_SNAPSHOT_DIR = normalizePath("../../data/snapshots", mustWork = FALSE)),
    code
  )
}

test_that("the coverage gap is reported by year", {
  skip_if_not(dir.exists("../../data/snapshots/births"), "infant snapshots not built")

  with_snapshots({
  covered <- infant_mortality_years()
  expect_null(infant_gap_message(covered[[1]]))

  # A year outside the covered range must be named as a source limitation
  # rather than silently returning nothing.
  gap_year <- setdiff(c(1990L, 2030L), covered)
  if (length(gap_year) > 0) {
    msg <- infant_gap_message(gap_year[[1]])
    expect_match(msg, "nados-vivos")
    expect_match(msg, as.character(gap_year[[1]]))
  }
  })
})

test_that("the count covers years the rate cannot", {
  skip_if_not(dir.exists("../../data/snapshots/births"), "infant snapshots not built")

  with_snapshots({
    rate_years <- infant_metric_years("infant")
    count_years <- infant_metric_years("infant_deaths")

    # The count needs only the numerator, so it can never cover less.
    expect_true(all(rate_years %in% count_years))
    expect_gt(length(count_years), length(rate_years))

    # 1991-1994 are the years this buys: deaths are on record, births are not.
    early <- setdiff(count_years, rate_years)
    expect_true(1991L %in% early)

    # And each metric is refused only outside its own coverage.
    expect_null(infant_gap_message(1991L, "infant_deaths"))
    expect_match(infant_gap_message(1991L, "infant"), "nados-vivos")
  })
})

test_that("a year published only as a total is refused at finer detail", {
  skip_if_not(dir.exists("../../data/snapshots/infant_deaths"), "infant snapshots not built")

  with_snapshots({
    latest <- max(infant_death_years())

    # The by-cause indicators end a year before the total-only one, so the last
    # year holds a single all-cause, both-sexes figure per area.
    expect_true(infant_year_has_detail(latest, "Todas as causas de morte", "HM"))

    detailed <- infant_detail_years(latest, "Doenças do aparelho circulatório", "HM")
    if (length(detailed) == 0) {
      # Asking it for a cause must be refused, not answered with zero - which
      # is what an unguarded filter on no matching rows would produce.
      msg <- infant_detail_message(latest, "Doenças do aparelho circulatório", "HM")
      expect_match(msg, as.character(latest))
      expect_match(msg, "apenas no total")
      expect_true(is.na(infant_deaths_total(latest, "Portugal",
                                            cause = "Doenças do aparelho circulatório")))

      # And the same request one year earlier still works.
      expect_null(infant_detail_message(latest - 1L, "Doenças do aparelho circulatório", "HM"))
    }

    # An all-cause request is never blocked.
    expect_null(infant_detail_message(infant_death_years(), "Todas as causas de morte", "HM"))
  })
})

test_that("a thin denominator is marked, not suppressed", {
  # The threshold is about resolution: below it a single death moves the rate
  # by more than one whole unit per 1,000.
  expect_true(infant_rate_is_unstable(9))
  expect_true(infant_rate_is_unstable(999))
  expect_false(infant_rate_is_unstable(1000))
  expect_false(infant_rate_is_unstable(84642))

  # An unknown denominator is not an unstable one; it has its own message.
  expect_false(infant_rate_is_unstable(NA_real_))

  # The value itself must still be computed - marking is not suppression.
  thin <- compute_infant_mortality(deaths = 0, births = 9)
  expect_equal(thin$value, 0)
  expect_gt(thin$upper, 100)

  expect_match(infant_instability_footnote(), "^\\*")
})

test_that("partial coverage of a window reports NA rather than a partial sum", {
  skip_if_not(dir.exists("../../data/snapshots/infant_deaths"), "infant snapshots not built")

  with_snapshots({
    covered <- infant_death_years()
    full <- infant_deaths_total(max(covered), "Portugal")
    expect_true(is.finite(full))
    expect_gt(full, 0)

    # A window straying outside the archive must not return the covered part as
    # if it were the whole: AVPP subtracts this from a wider band, so a partial
    # sum would silently move infant deaths into ages one to four.
    expect_true(is.na(infant_deaths_total(c(max(covered), max(covered) + 1L), "Portugal")))
    expect_true(is.na(infant_deaths_total(integer(0), "Portugal")))
  })
})

test_that("the committed snapshots reconcile with published national figures", {
  skip_if_not(dir.exists("../../data/snapshots/births"), "infant snapshots not built")

  with_snapshots({
    covered <- infant_mortality_years()
    expect_gt(length(covered), 15)

    # Portugal's published infant mortality has run between about 2 and 8 per
    # 1,000 across the covered period, falling over time. Anything outside that
    # means the numerator or denominator has been assembled wrongly.
    for (year in utils::tail(covered, 3)) {
      deaths <- sum(get_infant_death_data(year, "Portugal")$deaths, na.rm = TRUE)
      births <- sum(get_births_data(year, "Portugal")$births, na.rm = TRUE)
      rate <- compute_infant_mortality(deaths, births)$value

      expect_gt(births, 50000)
      expect_gt(rate, 1)
      expect_lt(rate, 10)
    }

    # Anchored against INE's published figures, so a future refetch that
    # silently changes the pinning is caught rather than merely looking sane.
    imr_2024 <- compute_infant_mortality(
      sum(get_infant_death_data(2024L, "Portugal")$deaths, na.rm = TRUE),
      sum(get_births_data(2024L, "Portugal")$births, na.rm = TRUE)
    )$value
    expect_equal(round(imr_2024, 1), 3.0)

    imr_2000 <- compute_infant_mortality(
      sum(get_infant_death_data(2000L, "Portugal")$deaths, na.rm = TRUE),
      sum(get_births_data(2000L, "Portugal")$births, na.rm = TRUE)
    )$value
    expect_equal(round(imr_2000, 1), 5.5)
  })
})
