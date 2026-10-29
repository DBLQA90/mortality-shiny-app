# Tests for the data mode (R/helpers.R): with "original", everything the app
# adds to or replaces in INE's data is switched off, and nothing computed in one
# mode is served in the other. DADOS_INE_LIMITACOES.md lists what it covers.

lookup <- get_nuts_lookup("2024")

with_snapshots <- function(code) {
  withr::with_envvar(c(MORTALITY_SNAPSHOT_DIR = normalizePath("../../data/snapshots", mustWork = FALSE)), code)
}

with_mode <- function(mode, code) withr::with_options(list(mortality.data_mode = mode), code)

test_that("the data mode defaults to corrected and accepts only its two values", {
  withr::with_options(list(mortality.data_mode = NULL), expect_equal(data_mode(), "corrected"))
  expect_equal(normalize_data_mode("original"), "original")
  expect_equal(normalize_data_mode("anything"), "corrected")
  expect_true(with_mode("corrected", data_corrections_on()))
  expect_false(with_mode("original", data_corrections_on()))
})

test_that("INE's 2011 Alto Alentejo rows are corrected only in the corrected mode", with_snapshots({
  path <- planning_death_cause_file(2011L, "Todas as causas de morte")
  nisa <- function(frame) sum(frame$deaths[frame$area == "Nisa" & frame$sex == "HM"])
  expect_equal(with_mode("corrected", nisa(read_death_file(path))), 165)
  expect_equal(with_mode("original", nisa(read_death_file(path))), 24)
}))

test_that("deaths without an age are left out in the original mode", with_snapshots({
  deaths <- readRDS(planning_death_cause_file(2024L, "Todas as causas de morte"))
  deaths <- deaths[deaths$area == "Barrancos" & deaths$sex == "HM", ]
  expect_equal(sum(with_mode("corrected", fill_death_ages(deaths))$deaths), 32, tolerance = 1e-6)
  original <- with_mode("original", fill_death_ages(deaths))
  expect_equal(sum(original$deaths), 19)
  expect_equal(sum(original$deaths_filled), 0)
}))

test_that("confidential totals stay missing in the original mode", with_snapshots({
  vimioso <- function(t) t$deaths[t$area == "Vimioso" & t$sex == "HM" & t$cause == "Todas as causas de morte"]
  expect_equal(with_mode("corrected", vimioso(repair_death_totals(death_totals_raw(2024L), 2024L))), 87)
  expect_true(is.na(with_mode("original", vimioso(repair_death_totals(death_totals_raw(2024L), 2024L)))))
}))

test_that("the rebuilt latest year is not offered in the original mode", with_snapshots({
  expect_true(length(with_mode("corrected", planning_calendar_years())) > 0)
  expect_length(with_mode("original", planning_calendar_years()), 0)
  expect_true(with_mode("corrected", calendar_selection_supported("Norte")))
  expect_false(with_mode("original", calendar_selection_supported("Norte")))
}))

test_that("each mode has its own cache", {
  a <- with_mode("corrected", planning_cache)
  b <- with_mode("original", planning_cache)
  expect_false(identical(a, b))
  with_mode("original", assign("data-mode-test", 1, envir = planning_cache))
  expect_false(with_mode("corrected", exists("data-mode-test", envir = planning_cache, inherits = FALSE)))
  with_mode("original", rm("data-mode-test", envir = planning_cache))
})

test_that("suppressed employees are recovered only in the corrected mode", with_snapshots({
  azores <- "Região Autónoma dos Açores"
  value <- function() {
    t <- planning_indicator_table(azores, 2013L, ids = "employees", lookup = lookup)
    t$value[t$area == azores]
  }
  expect_true(is.finite(with_mode("corrected", value())))
  expect_true(is.na(with_mode("original", value())))
}))
