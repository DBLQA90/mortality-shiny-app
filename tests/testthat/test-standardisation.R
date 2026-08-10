# Tests for R/standardisation.R (indirect standardisation and multi-year pooling).

as_bands <- function(x) factor(x, levels = age_levels, ordered = TRUE)

# Reference schedule with deliberately different risk by age, so that an area
# with a different age structure would get the wrong answer under a crude
# comparison but the right one under indirect standardisation.
ref_bands <- function(year = 2020L, sex = "HM", cause = "C") {
  tibble(
    year = year, sex = sex, cause = cause,
    age_band = as_bands(c("30 - 34 anos", "60 - 64 anos", "80 - 84 anos")),
    deaths = c(50, 500, 2000),
    pop = c(100000, 100000, 50000)
  )
}

test_that("pooling sums deaths and person-years, and is a no-op at window 1", {
  df <- tibble(
    year = 2000:2004, area = "A", sex = "HM", cause = "C",
    age_band = as_bands("60 - 64 anos"),
    deaths = c(10, 20, 30, 40, 50), pop = c(1000, 1000, 1000, 1000, 1000)
  )

  single <- pool_age_data(df, window = 1)
  expect_equal(single$deaths, df$deaths)
  expect_equal(single$pop, df$pop)
  expect_equal(single$n_years, rep(1L, 5))
  expect_equal(single$period_label, as.character(2000:2004))

  pooled <- pool_age_data(df, window = 3) %>% dplyr::arrange(year)
  # Centre year 2002 pools 2001-2003: deaths 20+30+40, person-years 3000.
  mid <- pooled %>% dplyr::filter(year == 2002)
  expect_equal(mid$deaths, 90)
  expect_equal(mid$pop, 3000)
  expect_equal(mid$n_years, 3L)
  expect_equal(mid$period_label, "2001-2003")

  # The pooled rate is total deaths over total person-years, not the mean of
  # the annual rates - those differ whenever population varies across the window.
  expect_equal(mid$deaths / mid$pop * 1e5, 90 / 3000 * 1e5)
})

test_that("pooling truncates at the series edges and reports the real span", {
  df <- tibble(
    year = 2000:2004, area = "A", sex = "HM", cause = "C",
    age_band = as_bands("60 - 64 anos"),
    deaths = rep(10, 5), pop = rep(1000, 5)
  )

  pooled <- pool_age_data(df, window = 5) %>% dplyr::arrange(year)

  first <- pooled %>% dplyr::filter(year == 2000)
  expect_equal(first$n_years, 3L)
  expect_equal(first$period_label, "2000-2002")
  expect_equal(first$deaths, 30)

  middle <- pooled %>% dplyr::filter(year == 2002)
  expect_equal(middle$n_years, 5L)
  expect_equal(middle$period_label, "2000-2004")
  expect_equal(middle$deaths, 50)

  last <- pooled %>% dplyr::filter(year == 2004)
  expect_equal(last$n_years, 3L)
  expect_equal(last$period_label, "2002-2004")
})

test_that("pooling keeps groups separate and rejects even windows", {
  df <- tibble(
    year = rep(2000:2002, each = 2), area = rep(c("A", "B"), 3),
    sex = "HM", cause = "C", age_band = as_bands("60 - 64 anos"),
    deaths = c(1, 100, 2, 200, 3, 300), pop = 1000
  )

  pooled <- pool_age_data(df, window = 3) %>% dplyr::filter(year == 2001)
  expect_equal(nrow(pooled), 2)
  expect_equal(pooled$deaths[pooled$area == "A"], 6)
  expect_equal(pooled$deaths[pooled$area == "B"], 600)

  expect_error(pool_age_data(df, window = 2), "ímpar")
})

test_that("expected deaths apply reference rates to the local age structure", {
  ref <- ref_bands()
  # Area with the same size but a much older structure than the reference.
  area <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(c("30 - 34 anos", "60 - 64 anos", "80 - 84 anos")),
    deaths = c(1, 20, 300), pop = c(1000, 2000, 5000)
  )

  parts <- compute_expected_deaths(area, ref)
  # 1000*(50/100000) + 2000*(500/100000) + 5000*(2000/50000)
  expect_equal(parts$expected, 0.5 + 10 + 200)
  expect_equal(parts$observed, 321)
  expect_equal(parts$bands_used, 3L)
  expect_length(parts$missing_bands, 0)
})

test_that("SMR of the reference against itself is exactly the reference value", {
  ref <- ref_bands()
  res <- compute_smr(ref, ref, refvalue = 100)

  expect_equal(res$smr, 100)
  expect_equal(res$observed, sum(ref$deaths))
  expect_equal(res$expected, sum(ref$deaths))
  # The indirectly standardised rate collapses to the reference crude rate.
  expect_equal(res$isr, sum(ref$deaths) / sum(ref$pop) * 1e5)
  expect_equal(classify_smr(res$smr_lower, res$smr_upper), "Sem diferença significativa")
})

test_that("SMR equals observed over expected and brackets it with an interval", {
  ref <- ref_bands()
  area <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(c("30 - 34 anos", "60 - 64 anos", "80 - 84 anos")),
    deaths = c(2, 30, 400), pop = c(1000, 2000, 5000)
  )

  res <- compute_smr(area, ref, refvalue = 100)
  expect_equal(res$smr, res$observed / res$expected * 100)
  expect_lt(res$smr_lower, res$smr)
  expect_gt(res$smr_upper, res$smr)
  expect_equal(classify_smr(res$smr_lower, res$smr_upper), "Acima da referência")

  # ISR is the SMR rescaled by the reference crude rate, so the ratio of the
  # two carries no extra information - only the scale differs.
  expect_equal(res$isr, res$smr / 100 * res$ref_rate)
})

test_that("SMR reports rather than hides unusable references", {
  ref <- ref_bands()
  zero_ref <- ref %>% dplyr::mutate(deaths = 0)
  area <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands("60 - 64 anos"), deaths = 5, pop = 1000
  )

  res <- compute_smr(area, zero_ref, refvalue = 100)
  expect_true(is.na(res$smr))
  expect_equal(res$expected, 0)
  expect_match(res$note, "esperados nulos")

  # An age band the reference cannot cover is excluded and named, not treated
  # as zero risk.
  partial_ref <- ref %>% dplyr::filter(age_band != "80 - 84 anos")
  older_area <- tibble(
    year = 2020L, sex = "HM", cause = "C",
    age_band = as_bands(c("60 - 64 anos", "80 - 84 anos")),
    deaths = c(10, 90), pop = c(1000, 2000)
  )
  res2 <- compute_smr(older_area, partial_ref, refvalue = 100)
  expect_equal(res2$observed, 10)
  expect_equal(res2$bands_used, 1L)
  expect_match(res2$note, "80 - 84 anos")

  expect_true(is.na(compute_smr(area[0, ], ref)$smr))
})

test_that("classify_smr only calls a difference when the interval excludes the reference", {
  expect_equal(classify_smr(120, 150), "Acima da referência")
  expect_equal(classify_smr(50, 90), "Abaixo da referência")
  expect_equal(classify_smr(90, 150), "Sem diferença significativa")
  expect_true(is.na(classify_smr(NA_real_, 150)))
})

test_that("compute_smr_series standardises each area against the matching reference", {
  ref <- dplyr::bind_rows(ref_bands(2020L), ref_bands(2021L)) %>%
    dplyr::mutate(area = "Portugal")

  areas <- dplyr::bind_rows(
    ref %>% dplyr::mutate(area = "Igual à referência"),
    ref %>% dplyr::mutate(area = "Dobro do risco", deaths = deaths * 2)
  )

  res <- compute_smr_series(areas, ref, refvalue = 100) %>%
    dplyr::arrange(area, year)

  expect_equal(nrow(res), 4)
  expect_equal(res$smr[res$area == "Igual à referência"], c(100, 100))
  expect_equal(res$smr[res$area == "Dobro do risco"], c(200, 200))
  expect_equal(
    unique(res$smr_class[res$area == "Dobro do risco"]),
    "Acima da referência"
  )
})

test_that("pooling then standardising narrows the interval for a sparse area", {
  ref <- dplyr::bind_rows(lapply(2018:2022, function(y) ref_bands(y) %>% dplyr::mutate(area = "Portugal")))
  sparse <- dplyr::bind_rows(lapply(2018:2022, function(y) {
    tibble(
      year = y, area = "Pequeno", sex = "HM", cause = "C",
      age_band = as_bands(c("60 - 64 anos", "80 - 84 anos")),
      deaths = c(1, 3), pop = c(400, 300)
    )
  }))

  single <- compute_smr(
    sparse %>% dplyr::filter(year == 2020),
    ref %>% dplyr::filter(year == 2020)
  )
  pooled_area <- pool_age_data(sparse, window = 5) %>% dplyr::filter(year == 2020)
  pooled_ref <- pool_age_data(ref, window = 5) %>% dplyr::filter(year == 2020)
  pooled <- compute_smr(pooled_area, pooled_ref)

  expect_equal(pooled_area$n_years[[1]], 5L)
  # Same underlying risk, five times the data: the point estimate is unchanged
  # but the interval is materially tighter.
  expect_equal(pooled$smr, single$smr, tolerance = 1e-8)
  expect_lt(pooled$smr_upper - pooled$smr_lower, single$smr_upper - single$smr_lower)
})
