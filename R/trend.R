# =========================================================
# Trend: annual percentage change
# =========================================================
#
# The method of the US National Cancer Institute's Joinpoint program: a
# log-linear regression of the rate on the year, so the slope b is a constant
# relative change and APC = 100 * (exp(b) - 1) per year. Each year is weighted
# by the inverse of the sampling variance of its log rate, taken from its own
# 95% interval, so a year with few events counts less; the residual variance
# is still estimated, which absorbs year-to-year variation beyond chance.
#
# With joinpoints the line may bend at up to `max_joinpoints` years and stays
# continuous there (hinge terms). Each candidate year is tried; the number of
# joinpoints is chosen by BIC (one of Joinpoint's own selection methods), and
# each segment has at least `min_segment` years. The average annual
# percentage change (AAPC) over the whole period weights each segment's slope
# by its length; its interval comes from the same covariance matrix, so it is
# exact for the fitted model rather than an average of intervals.
#
# Not a Mann-Kendall test: that says whether a series goes up or down, not how
# fast, and a planner needs the rate.

trend_log_weights <- function(values, lower, upper, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  se <- (log(upper) - log(lower)) / (2 * z)
  ok <- is.finite(se) & se > 0
  if (!any(ok)) return(rep(1, length(values)))
  se[!ok] <- stats::median(se[ok])
  1 / se^2
}

trend_design <- function(years, joins) {
  x <- cbind(1, years)
  for (j in joins) x <- cbind(x, pmax(years - j, 0))
  x
}

# The fit with the given joinpoints, or NULL when it cannot be made.
trend_fit_joins <- function(years, y, w, joins) {
  x <- trend_design(years, joins)
  if (nrow(x) <= ncol(x)) return(NULL)
  fit <- stats::lm.wfit(x, y, w)
  if (fit$rank < ncol(x)) return(NULL)
  residual <- y - as.vector(x %*% fit$coefficients)
  df <- nrow(x) - ncol(x)
  sigma2 <- sum(w * residual^2) / df
  xtwx <- crossprod(x * sqrt(w))
  vcov <- tryCatch(sigma2 * solve(xtwx), error = function(e) NULL)
  if (is.null(vcov)) return(NULL)
  n <- length(y)
  list(joins = joins, coef = fit$coefficients, vcov = vcov, df = df,
       fitted = as.vector(x %*% fit$coefficients),
       # Joinpoint counts each joinpoint as a parameter, as well as its slope.
       bic = n * log(sum(w * residual^2) / n) + (2 + 2 * length(joins)) * log(n))
}

trend_combination <- function(fit, weights, confidence) {
  slope <- sum(weights * fit$coef)
  se <- sqrt(as.numeric(t(weights) %*% fit$vcov %*% weights))
  q <- stats::qt(1 - (1 - confidence) / 2, fit$df)
  tibble::tibble(
    apc = 100 * (exp(slope) - 1),
    lower = 100 * (exp(slope - q * se) - 1),
    upper = 100 * (exp(slope + q * se) - 1),
    p_value = 2 * stats::pt(-abs(slope / se), fit$df)
  )
}

# APC per segment and AAPC over the period for one series. `lower` and `upper`
# are the 95% interval of each value (NA: equal weights). NULL with fewer than
# `min_points` usable years.
trend_apc <- function(years, values, lower = NA, upper = NA, max_joinpoints = 0L, min_segment = 4L,
                      confidence = 0.95, min_points = 5L) {
  lower <- rep_len(lower, length(values))
  upper <- rep_len(upper, length(values))
  ok <- is.finite(years) & is.finite(values) & values > 0
  years <- as.numeric(years[ok]); values <- values[ok]; lower <- lower[ok]; upper <- upper[ok]
  o <- order(years)
  years <- years[o]; y <- log(values[o]); w <- trend_log_weights(values[o], lower[o], upper[o])
  n <- length(years)
  if (n < min_points) return(NULL)

  # Joinpoints can sit on any year that leaves `min_segment` years on each side.
  candidates <- years[years >= years[min_segment] & years <= years[n - min_segment + 1L]]
  fits <- list(trend_fit_joins(years, y, w, numeric(0)))
  for (k in seq_len(max(0L, as.integer(max_joinpoints)))) {
    if (length(candidates) < k) break
    sets <- utils::combn(candidates, k, simplify = FALSE)
    sets <- Filter(function(s) all(diff(c(years[1], s, years[n])) >= min_segment - 1L), sets)
    best <- NULL
    for (s in sets) {
      f <- trend_fit_joins(years, y, w, s)
      if (!is.null(f) && (is.null(best) || f$bic < best$bic)) best <- f
    }
    if (!is.null(best)) fits[[length(fits) + 1L]] <- best
  }
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) return(NULL)
  fit <- fits[[which.min(vapply(fits, `[[`, numeric(1), "bic"))]]

  bounds <- c(years[1], fit$joins, years[n])
  p <- length(fit$coef)
  segments <- dplyr::bind_rows(lapply(seq_len(length(bounds) - 1L), function(m) {
    weights <- c(0, 1, rep(1, m - 1L), rep(0, length(fit$joins) - m + 1L))[seq_len(p)]
    dplyr::bind_cols(tibble::tibble(from = bounds[m], to = bounds[m + 1L]), trend_combination(fit, weights, confidence))
  }))
  lengths <- diff(bounds)
  # The average slope over the period: each segment's slope weighted by its length.
  slope_weights <- vapply(seq_len(length(fit$joins)), function(j) sum(lengths[(j + 1L):length(lengths)]) / sum(lengths), numeric(1))
  aapc <- dplyr::bind_cols(tibble::tibble(from = years[1], to = years[n]),
                           trend_combination(fit, c(0, 1, slope_weights), confidence))
  list(segments = segments, aapc = aapc, joinpoints = fit$joins,
       fitted = tibble::tibble(year = years, fitted = exp(fit$fitted)), n = n)
}

# One line of text: "-1,8% por ano (-2,4; -1,2)" and whether the interval excludes 0.
trend_label <- function(row, digits = 1) {
  fmt <- function(x) formatC(x, format = "f", digits = digits, decimal.mark = ",", flag = "+")
  paste0(fmt(row$apc), "% por ano (", fmt(row$lower), "; ", fmt(row$upper), ")")
}

trend_reading <- function(row) {
  if (!is.finite(row$lower) || !is.finite(row$upper)) return("")
  if (row$upper < 0) "descida significativa" else if (row$lower > 0) "subida significativa" else "sem tendência significativa"
}

# A warning when the period crosses 2020/2021, where two things move every
# rate at once: COVID-19 (2020-2022) and INE's revised population series from
# 2021, which steps the denominator by about 1.7% (METHODOLOGY.md). A straight
# line through them absorbs both; a joinpoint near those years may be them
# rather than a change of trend.
trend_break_note <- function(years, joinpoints = numeric(0)) {
  years <- as.integer(years)
  if (length(years) == 0 || min(years) > 2020L || max(years) < 2021L) return(NULL)
  near <- joinpoints[joinpoints >= 2019 & joinpoints <= 2022]
  paste0(
    "O período atravessa 2020-2021: a COVID-19 (2020-2022) e a série de população revista do INE a partir de 2021 ",
    "(um degrau de cerca de 1,7% no denominador) afectam a tendência. ",
    if (length(near) > 0) {
      paste0("A mudança de tendência em ", paste(near, collapse = " e "), " pode reflectir isso e não uma mudança real. ")
    } else {
      ""
    },
    "Compare com o período anterior a 2020 antes de concluir."
  )
}
