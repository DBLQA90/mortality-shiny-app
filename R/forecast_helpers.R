# =========================================================
# Forecast metric + validation helpers
# =========================================================
# Pure helpers shared by the forecasting paths. Kept out of the Shiny server so
# they can be unit tested directly (see tests/testthat/test-forecast-helpers.R).

# Minimum training years required before an out-of-sample test set can be held
# out; below this the app falls back to in-sample selection.
MIN_VALIDATION_TRAIN <- 3L

# Forecast accuracy metrics for a pair of equal-length series, on whatever scale
# the caller supplies. `scale_denom` is the in-sample naive MAE used for MASE.
compute_forecast_error_metrics <- function(actual, predicted, scale_denom = NA_real_) {
  if (length(actual) != length(predicted)) {
    stop("As séries observada e prevista têm de ter o mesmo comprimento.", call. = FALSE)
  }

  err <- predicted - actual
  nonzero_actual <- abs(actual) > 1e-8

  tibble::tibble(
    ME = mean(err, na.rm = TRUE),
    RMSE = sqrt(mean(err^2, na.rm = TRUE)),
    MAE = mean(abs(err), na.rm = TRUE),
    MAPE = if (any(nonzero_actual)) mean(abs(err[nonzero_actual] / actual[nonzero_actual]), na.rm = TRUE) * 100 else NA_real_,
    MASE = if (is.finite(scale_denom) && scale_denom > 0) mean(abs(err), na.rm = TRUE) / scale_denom else NA_real_
  )
}

# Size of the out-of-sample test window (in years) for a series of length `n`,
# given a requested test fraction. Returns 0 when the series is too short to
# leave `min_train` training years, signalling an in-sample fallback.
compute_validation_test_size <- function(n, test_fraction, min_train = MIN_VALIDATION_TRAIN) {
  if (!is.finite(n) || !is.finite(test_fraction) || n < min_train + 1L) {
    return(0L)
  }
  k <- max(1L, as.integer(round(n * test_fraction)))
  k <- min(k, as.integer(n - min_train))
  if (k < 1L) 0L else k
}

# In-sample accuracy table for a set of fitted models, on the original rate
# scale. Models are fitted on the (possibly log-transformed) modelling scale, so
# fitted values are back-transformed via `inverse` before comparison. This keeps
# ME/RMSE/MAE/MAPE/MASE on the same per-100,000 scale as the forecasts and
# holdout metrics, and makes the recommended-model choice reflect fit quality on
# the scale the user actually sees.
build_accuracy_table <- function(fits, obs = NULL, inverse = function(x) as.numeric(x)) {
  metrics <- c("ME", "RMSE", "MAE", "MAPE", "MASE")

  na_row <- function(m) {
    df <- as.data.frame(stats::setNames(as.list(rep(NA_real_, length(metrics))), metrics))
    df$Model <- m
    df[, c("Model", metrics), drop = FALSE]
  }

  if (length(fits) == 0) {
    return(
      tibble::tibble(
        Model = character(0),
        ME = numeric(0),
        RMSE = numeric(0),
        MAE = numeric(0),
        MAPE = numeric(0),
        MASE = numeric(0)
      )
    )
  }

  actual <- if (!is.null(obs)) as.numeric(obs$value) else NULL
  scale_denom <- if (!is.null(actual) && sum(is.finite(actual)) >= 2) {
    mean(abs(diff(actual)), na.rm = TRUE)
  } else {
    NA_real_
  }

  dplyr::bind_rows(lapply(names(fits), function(m) {
    tryCatch({
      if (is.null(actual)) {
        stop("Observed series not supplied for accuracy calculation.", call. = FALSE)
      }

      fitted_vals <- inverse(as.numeric(stats::fitted(fits[[m]])))
      n <- min(length(fitted_vals), length(actual))
      a <- actual[seq_len(n)]
      p <- fitted_vals[seq_len(n)]
      finite <- is.finite(a) & is.finite(p)

      if (!any(finite)) {
        stop("No overlapping finite fitted values for accuracy calculation.", call. = FALSE)
      }

      df <- compute_forecast_error_metrics(
        actual = a[finite],
        predicted = p[finite],
        scale_denom = scale_denom
      )
      df$Model <- m
      df[, c("Model", metrics), drop = FALSE]
    }, error = function(e) {
      na_row(m)
    })
  }), .id = NULL)
}

# Back-transformed fitted values per model, aligned to the observed years, so a
# comparison plot draws fitted lines on the same rate scale as the observations.
build_fitted_values_df <- function(obs, fits, inverse = function(x) as.numeric(x)) {
  if (length(fits) == 0) {
    return(tibble::tibble(year = numeric(0), fitted = numeric(0), model = character(0)))
  }

  dplyr::bind_rows(lapply(names(fits), function(model_id) {
    fitted_vals <- inverse(as.numeric(stats::fitted(fits[[model_id]])))
    year_count <- min(length(fitted_vals), nrow(obs))

    tibble::tibble(
      year = obs$year[seq_len(year_count)],
      fitted = fitted_vals[seq_len(year_count)],
      model = model_id
    )
  }))
}

# ---------------------------------------------------------
# Infant mortality over several years
# ---------------------------------------------------------
# Infant mortality is usually reported over three or five years, because a
# region has a few dozen infant deaths a year. Two ways to project it:
#   "annual"  fit the models to the yearly rates and add up the projected years
#             into each window afterwards (infant_window_forecast());
#   "rolling" fit them to the series of window rates itself.
# The rolling series shares four of its five years with its neighbour, so it is
# smooth by construction: the models read that smoothness as signal, choose
# trends too readily, and give intervals far narrower than the yearly
# variation warrants. The annual fit sees the real year-to-year noise.

# Under-1 deaths and live births of a selection, by year.
infant_annual_series <- function(years, areas) {
  rows <- lapply(as.integer(years), function(year) {
    tibble::tibble(year = year, deaths = infant_deaths_total(year, areas), births = get_births_total(year, areas))
  })
  out <- dplyr::bind_rows(rows)
  out[is.finite(out$deaths) & is.finite(out$births) & out$births > 0, , drop = FALSE]
}

# Rate per 1,000 live births over the `window` years ending in each year, with
# an exact Poisson interval; only years whose whole window is present.
infant_window_rates <- function(annual, window = 1L) {
  annual <- annual[order(annual$year), , drop = FALSE]
  window <- as.integer(window)
  rows <- lapply(annual$year, function(year) {
    span <- seq.int(year - window + 1L, year)
    if (!all(span %in% annual$year)) return(NULL)
    take <- annual$year %in% span
    r <- compute_infant_mortality(sum(annual$deaths[take]), sum(annual$births[take]))
    tibble::tibble(year = year, value = r$value, lower = r$lower, upper = r$upper,
                   deaths = r$deaths, births = r$births)
  })
  dplyr::bind_rows(rows)
}

# Window rates from forecasts of the yearly rate (`fc`: year, model, mean,
# lower, upper, per 1,000). Each window adds the observed years it still holds
# to the projected ones; projected years count with the last observed births,
# so the rate is a births-weighted mean. The interval adds the yearly limits as
# if the forecast errors were perfectly correlated across years - conservative,
# and close for the trend and random-walk models, whose errors are.
infant_window_forecast <- function(annual, fc, window) {
  window <- as.integer(window)
  annual <- annual[order(annual$year), , drop = FALSE]
  last_births <- utils::tail(annual$births, 1)
  rows <- lapply(split(fc, fc$model), function(f) {
    f <- f[order(f$year), , drop = FALSE]
    dplyr::bind_rows(lapply(f$year, function(year) {
      span <- seq.int(year - window + 1L, year)
      observed <- annual[annual$year %in% span, , drop = FALSE]
      projected <- f[f$year %in% span, , drop = FALSE]
      if (nrow(observed) + nrow(projected) < window) return(NULL)
      births <- sum(observed$births) + last_births * nrow(projected)
      combine <- function(rate) (sum(observed$deaths) + sum(rate) / 1000 * last_births) / births * 1000
      tibble::tibble(year = year, model = f$model[[1]], mean = combine(projected$mean),
                     lower = combine(projected$lower), upper = combine(projected$upper))
    }))
  })
  dplyr::bind_rows(rows)
}
