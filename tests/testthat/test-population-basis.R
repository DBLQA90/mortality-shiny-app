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

test_that("the planning rates per inhabitant follow the basis; stocks do not exist here", {
  components <- tibble::tibble(area = "Alfa", year = c(2021L, 2022L), deaths = c(10, 12), pop_total = c(1000, 1200))
  mean_rate <- planning_compute_indicators(components, "death_rate", basis = "mean")
  end_rate <- planning_compute_indicators(components, "death_rate", basis = "end")
  expect_equal(mean_rate$value[mean_rate$year == 2022], 12 / 1100 * 1000)
  expect_equal(end_rate$value[end_rate$year == 2022], 12 / 1200 * 1000)
})
