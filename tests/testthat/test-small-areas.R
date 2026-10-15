# Tests for the empirical Bayes smoothing of the funnel's units.

poisson <- list(model = "poisson", multiplier = 1000)

test_that("a small unit moves towards the mean, a large one barely", {
  units <- data.frame(area = c("grande", "pequena", "media", "outra"),
                      value = c(12, 30, 9, 10), denominator = c(500000, 800, 40000, 60000))
  smoothed <- planning_eb_smooth(units, poisson)
  mean_all <- attr(smoothed, "smoothing_mean")
  expect_equal(mean_all, sum(units$value * units$denominator) / sum(units$denominator))
  expect_gt(smoothed$shrinkage_weight[1], 0.9)
  expect_lt(smoothed$shrinkage_weight[2], smoothed$shrinkage_weight[3])
  # Each value ends between where it was and the mean.
  expect_true(all(abs(smoothed$smoothed - mean_all) <= abs(units$value - mean_all) + 1e-9))
})

test_that("with no spread beyond chance every unit takes the mean", {
  units <- data.frame(area = c("a", "b", "c"), value = c(10, 10, 10), denominator = c(1000, 5000, 20000))
  smoothed <- planning_eb_smooth(units, poisson)
  expect_equal(smoothed$shrinkage_weight, c(0, 0, 0))
  expect_equal(smoothed$smoothed, c(10, 10, 10))
})

test_that("proportions use the binomial variance", {
  units <- data.frame(area = c("a", "b", "c"), value = c(8, 12, 30), denominator = c(2000, 3000, 50))
  binomial <- planning_eb_smooth(units, list(model = "binomial", multiplier = 100))
  expect_true(all(binomial$shrinkage_weight >= 0 & binomial$shrinkage_weight <= 1))
  expect_lt(binomial$shrinkage_weight[3], binomial$shrinkage_weight[2])
})
