# Tests for the population denominator option (population_basis_choices).

test_that("the mean needs the previous year's estimate as well", {
  expect_equal(population_years_to_load(c(2020L, 2021L), "mean"), 2019:2021)
  expect_equal(population_years_to_load(c(2020L, 2021L), "end"), c(2020L, 2021L))
  expect_equal(normalize_population_basis("nonsense"), "mean")
})

test_that("each year's estimate becomes the mean with the previous year's, or stays its own", {
  pop <- tibble::tibble(
    year = c(2020L, 2021L, 2022L), area = "Alfa", sex = "HM", age_band = "0 - 4 anos",
    pop = c(100, 120, 150), source_indicator = "x"
  )
  mean_basis <- apply_population_basis(pop, 2021:2022, "mean")
  expect_equal(mean_basis$year, 2021:2022)
  expect_equal(mean_basis$pop, c(110, 135))
  # The first year of the series has no predecessor and keeps its own estimate.
  expect_equal(apply_population_basis(pop, 2020L, "mean")$pop, 100)
  end_basis <- apply_population_basis(pop, 2021:2022, "end")
  expect_equal(end_basis$pop, c(120, 150))
})

test_that("a previous year the snapshot does not have is skipped, not an error", {
  loader <- function(year) {
    if (year == 1990L) stop("Missing population snapshot chunk(s) for year(s): 1990.")
    tibble::tibble(year = year, area = "Alfa", sex = "HM", age_band = "0 - 4 anos", pop = 100, source_indicator = "x")
  }
  # 1991 needs 1990, which does not exist; 1992 needs 1991, already requested.
  prev <- load_previous_population(c(1991L, 1992L), "mean", loader)
  expect_equal(nrow(prev), 0)
  expect_equal(load_previous_population(2000L, "mean", loader)$year, 1999L)
  expect_equal(nrow(load_previous_population(2000L, "end", loader)), 0)
})

test_that("the planning rates per inhabitant follow the basis; stocks do not exist here", {
  components <- tibble::tibble(area = "Alfa", year = c(2021L, 2022L), deaths = c(10, 12), pop_total = c(1000, 1200))
  mean_rate <- planning_compute_indicators(components, "death_rate", basis = "mean")
  end_rate <- planning_compute_indicators(components, "death_rate", basis = "end")
  expect_equal(mean_rate$value[mean_rate$year == 2022], 12 / 1100 * 1000)
  expect_equal(end_rate$value[end_rate$year == 2022], 12 / 1200 * 1000)
})
