# Tests for R/forecast_helpers.R (shared forecast-metric + validation helpers).

test_that("compute_validation_test_size clamps to the requested fraction and training floor", {
  expect_equal(compute_validation_test_size(21, 0.25), 5L)   # round(5.25)
  expect_equal(compute_validation_test_size(21, 0.40), 8L)
  expect_equal(compute_validation_test_size(4, 0.25), 1L)    # min 1 test year
  expect_equal(compute_validation_test_size(3, 0.25), 0L)    # too short -> fallback
  expect_equal(compute_validation_test_size(5, 0.90), 2L)    # clamped to n - min_train
  expect_equal(compute_validation_test_size(NA, 0.25), 0L)
  expect_equal(compute_validation_test_size(20, NA), 0L)
})

test_that("compute_forecast_error_metrics computes the standard accuracy measures", {
  actual <- c(100, 200, 300)
  predicted <- c(110, 190, 330)   # errors: +10, -10, +30
  m <- compute_forecast_error_metrics(actual, predicted, scale_denom = 50)

  expect_equal(m$ME, mean(c(10, -10, 30)))
  expect_equal(m$RMSE, sqrt(mean(c(10, 10, 30)^2)))
  expect_equal(m$MAE, mean(c(10, 10, 30)))
  expect_equal(m$MAPE, mean(c(10 / 100, 10 / 200, 30 / 300)) * 100)
  expect_equal(m$MASE, mean(c(10, 10, 30)) / 50)
})

test_that("compute_forecast_error_metrics guards inputs and degenerate cases", {
  expect_error(compute_forecast_error_metrics(1:3, 1:2), "mesmo comprimento")
  # all-zero actuals -> MAPE undefined
  expect_true(is.na(compute_forecast_error_metrics(c(0, 0), c(1, 2))$MAPE))
  # no valid scale -> MASE undefined
  expect_true(is.na(compute_forecast_error_metrics(c(1, 2), c(1, 2), scale_denom = NA_real_)$MASE))
})

test_that("build_accuracy_table reports metrics on the original (back-transformed) scale", {
  skip_if_not_installed("forecast")

  set.seed(7)
  years <- 2000:2019
  values <- 500 - 8 * (years - 2000) + rnorm(20, 0, 10)  # hundreds-scale rates
  offset <- min(values[values > 0]) / 2
  ts_log <- stats::ts(log(values + offset), start = min(years))
  inv <- function(x) pmax(exp(as.numeric(x)) - offset, 0)

  fit <- forecast::ets(ts_log)
  obs <- tibble::tibble(year = years, value = values)
  acc <- build_accuracy_table(list(ets = fit), obs = obs, inverse = inv)

  # Manual original-scale RMSE from back-transformed fitted values.
  fitted_orig <- inv(as.numeric(stats::fitted(fit)))
  expected_rmse <- sqrt(mean((fitted_orig - values)^2))
  expect_equal(acc$RMSE, expected_rmse, tolerance = 1e-8)

  # Sanity: original-scale RMSE is on the rates' scale (tens), not log-units (<1).
  expect_gt(acc$RMSE, 1)

  # Contrast: a naive log-scale accuracy would be tiny -> confirms the fix matters.
  log_rmse <- sqrt(mean((as.numeric(stats::fitted(fit)) - as.numeric(ts_log))^2))
  expect_lt(log_rmse, 1)
})

test_that("build_accuracy_table handles empty and NA-observation cases", {
  empty <- build_accuracy_table(list())
  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("Model", "ME", "RMSE", "MAE", "MAPE", "MASE") %in% names(empty)))
})

test_that("build_fitted_values_df back-transforms and aligns fitted values to years", {
  skip_if_not_installed("forecast")

  set.seed(9)
  years <- 2005:2019
  values <- 300 - 5 * (years - 2005) + rnorm(15, 0, 6)
  offset <- min(values[values > 0]) / 2
  inv <- function(x) pmax(exp(as.numeric(x)) - offset, 0)
  fit <- forecast::ets(stats::ts(log(values + offset), start = min(years)))
  obs <- tibble::tibble(year = years, value = values)

  fdf <- build_fitted_values_df(obs, list(ets = fit), inverse = inv)

  expect_equal(nrow(fdf), length(years))
  expect_equal(fdf$year, years)
  expect_equal(fdf$model, rep("ets", length(years)))
  # values back-transformed to the rate scale (not log-units)
  expect_gt(mean(fdf$fitted), 50)
  expect_equal(fdf$fitted, inv(as.numeric(stats::fitted(fit)))[seq_len(length(years))], tolerance = 1e-8)
})
