# Tests for the annual percentage change (R/trend.R).

test_that("a constant relative fall gives its APC", {
  years <- 2000:2019
  values <- 100 * 0.98^(years - 2000)
  trend <- trend_apc(years, values, values * 0.97, values * 1.03)
  expect_equal(trend$aapc$apc, -2, tolerance = 1e-6)
  expect_equal(nrow(trend$segments), 1)
  expect_equal(trend_reading(trend$aapc), "descida significativa")
})

test_that("a bend is found and the AAPC weights the segments by length", {
  set.seed(2)
  years <- 2000:2024
  slope <- ifelse(years <= 2012, -0.04 * (years - 2000), -0.04 * 12 + 0.01 * (years - 2012))
  values <- 100 * exp(slope + stats::rnorm(length(years), 0, 0.003))
  trend <- trend_apc(years, values, values * 0.98, values * 1.02, max_joinpoints = 2)
  expect_equal(trend$joinpoints, 2012)
  expect_equal(trend$segments$apc, 100 * (exp(c(-0.04, 0.01)) - 1), tolerance = 0.02)
  expect_equal(trend$aapc$apc, 100 * (exp((12 * -0.04 + 12 * 0.01) / 24) - 1), tolerance = 0.02)
})

test_that("noise around a flat line is not a trend, and short series give none", {
  set.seed(3)
  values <- 50 * exp(stats::rnorm(12, 0, 0.05))
  trend <- trend_apc(2010:2021, values, values * 0.8, values * 1.25)
  expect_lt(trend$aapc$lower, 0)
  expect_gt(trend$aapc$upper, 0)
  expect_equal(trend_reading(trend$aapc), "sem tendência significativa")
  expect_null(trend_apc(2020:2023, c(1, 2, 3, 4)))
})

test_that("an imprecise year counts less", {
  years <- 2010:2019
  values <- 100 * 0.97^(years - 2010)
  values[10] <- 200
  precise <- trend_apc(years, values, values * 0.98, values * 1.02)
  # The outlier's interval is wide, so it barely moves the slope.
  lower <- values * 0.98; upper <- values * 1.02
  lower[10] <- 20; upper[10] <- 2000
  weighted <- trend_apc(years, values, lower, upper)
  expect_lt(abs(weighted$aapc$apc + 3), abs(precise$aapc$apc + 3))
})
