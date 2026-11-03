# Tests for R/inla_forecasts.R: reading the precomputed INLA forecasts, on a
# small file written here (the real ones are built by tools/precompute_inla.R
# and are not in git).

with_inla_dir <- function(code) {
  dir <- tempfile("inla_")
  dir.create(dir)
  old <- Sys.getenv("MORTALITY_INLA_DIR", unset = NA)
  Sys.setenv(MORTALITY_INLA_DIR = dir)
  on.exit(if (is.na(old)) Sys.unsetenv("MORTALITY_INLA_DIR") else Sys.setenv(MORTALITY_INLA_DIR = old), add = TRUE)
  forecast <- tibble::tibble(
    level = rep(c("municipio", "NUTS II"), each = 3), area = rep(c("Beja", "Alentejo"), each = 3),
    year = rep(2024:2026, 2), mean = 30, q500 = c(30, 31, 32, 40, 41, 42),
    p25 = c(20, 21, 0, 35, 36, 37), p975 = c(45, 46, 47, 46, 47, 48)
  )
  saveRDS(list(origin = 2024L, horizon = 2L, forecast = forecast,
               fixed = data.frame(mean = c(0.1, -0.2), `0.025quant` = c(0, -0.3), `0.975quant` = c(0.2, -0.1),
                                  row.names = c("step2014", "covid2020"), check.names = FALSE)),
          file.path(dir, paste0(inla_key("Diabetes mellitus", "HM"), "__", inla_combo(), ".rds")))
  saveRDS(list(origin = 2024L), file.path(dir, "manifest.rds"))
  force(code)
}

test_that("the INLA forecast serves one area, the standardised rate and the corrected data", {
  with_inla_dir({
    reason <- function(...) inla_unavailable_reason(cause = "Diabetes mellitus", sex = "HM", mode = "corrected", ...)
    expect_null(reason(areas = "Beja"))
    expect_null(reason(areas = "Alentejo", last_year = 2024))
    expect_match(reason(areas = c("Beja", "Moura")), "um local")
    expect_match(reason(areas = "Beja", rate_type = "crude"), "padronizada")
    expect_match(reason(areas = "Beja", population = "Menos de 75 anos"), "total")
    expect_match(inla_unavailable_reason("Beja", "Diabetes mellitus", "HM", mode = "original"), "corrigidos")
    expect_match(inla_unavailable_reason("Beja", "Pneumonia", "HM", mode = "corrected"), "calculada")
    expect_match(inla_unavailable_reason("Beja", "Todas as causas de morte", "HM", mode = "corrected"), "8 causas")
    expect_match(reason(areas = "Lisboa"), "municípios")
    expect_match(reason(areas = "Alentejo", vintage = "2013"), "2024")
    expect_null(reason(areas = "Beja", vintage = "2013"))
    expect_match(reason(areas = "Beja", last_year = 2020), "2024")
    expect_equal(inla_origin_year(), 2024L)
  })
})

test_that("the forecast starts after the origin, keeps the 95% interval and scales other levels", {
  with_inla_dir({
    f <- inla_forecast("Beja", "Diabetes mellitus", "HM", widening = c(municipio = 1))
    expect_equal(f$year, 2025:2026)
    expect_equal(f$model, rep(INLA_MODEL_ID, 2))
    expect_equal(f$mean, c(31, 32))
    expect_equal(f$lower, c(21, 0))
    expect_equal(f$upper, c(46, 47))
    expect_equal(nrow(inla_forecast("Beja", "Diabetes mellitus", "HM", horizon = 1)), 1)
    f80 <- inla_forecast("Beja", "Diabetes mellitus", "HM", conf_level = 80, widening = c(municipio = 1))
    expect_true(all(f80$lower >= f$lower & f80$upper <= f$upper))
    expect_true(all(f80$lower <= f80$mean & f80$mean <= f80$upper))
    # Widening keeps the median and moves both limits out, on the log scale.
    w <- inla_forecast("Beja", "Diabetes mellitus", "HM", widening = c(municipio = 1.5))
    expect_equal(w$mean, f$mean)
    expect_equal(w$lower[1], 31 * (21 / 31)^1.5)
    expect_equal(w$upper[1], 31 * (46 / 31)^1.5)
    expect_equal(w$lower[2], 0)
    expect_null(inla_forecast("Beja", "Diabetes mellitus", "HM", combo = setdiff(INLA_COMBOS, inla_combo())[[1]]))
  })
})

test_that("the effects are rate ratios with their intervals", {
  with_inla_dir({
    e <- inla_effects("Diabetes mellitus", "HM")
    expect_equal(e$Termo, c("Quebra da série de população, 2014", "COVID-19, 2020"))
    expect_equal(e$`Razão de taxas`, exp(c(0.1, -0.2)))
    expect_true(all(e$`IC 95% inf.` <= e$`Razão de taxas` & e$`Razão de taxas` <= e$`IC 95% sup.`))
  })
})

test_that("the map takes INLA's areas, measured against the level it estimates for the last year", {
  with_inla_dir({
    rows <- inla_map_rows(c("Portugal", "Beja", "Alentejo"), "Diabetes mellitus", "HM", 2026, mode = "corrected")
    expect_equal(rows$area, c("Beja", "Alentejo"))
    expect_equal(rows$last_value, c(30, 40))
    expect_equal(rows$value, c(32, 42))
    expect_equal(rows$method, rep(INLA_MODEL_ID, 2))
    expect_equal(inla_map_rows(c("Beja", "Alentejo"), "Diabetes mellitus", "HM", 2026, mode = "corrected", vintage = "2013")$area, "Beja")
    expect_null(inla_map_rows("Beja", "Diabetes mellitus", "HM", 2027, mode = "corrected"))
    expect_null(inla_map_rows("Beja", "Diabetes mellitus", "HM", 2026, mode = "corrected", under75 = TRUE))
    expect_null(inla_map_rows("Beja", "Diabetes mellitus", "HM", 2026, mode = "original"))
    hover <- cause_map_projection_hover(dplyr::mutate(rows, index = value / last_value * 100, class = NA, benchmark_value = 33), 2026)
    expect_match(hover[1], "2024 \\(nível estimado\\)")
  })
})

test_that("a forecast past the tested years says its intervals are untested there", {
  with_inla_dir({
    expect_equal(inla_horizon_note(2028), "")
    expect_match(inla_horizon_note(2029), "testados até 4 anos à frente \\(2028\\)")
    expect_equal(inla_horizon_note(-Inf), "")
    expect_equal(inla_horizon_note(2035, origin = NA_integer_), "")
  })
})
