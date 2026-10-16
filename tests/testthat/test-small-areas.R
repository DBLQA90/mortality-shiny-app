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

test_that("every map level has a boundary for every area the app names", {
  lookup <- get_nuts_lookup("2024")
  for (level in c("municipio", "nuts3", "nuts2")) {
    b <- planning_map_boundaries(level, "2024")
    expect_true(all(planning_map_units(level, lookup) %in% b$polygons$id), info = level)
  }
  b <- planning_map_boundaries("uls")
  uls <- readRDS(file.path(app_dir, "data", "uls_lookup.rds"))
  expect_true(all(unique(uls$unit[uls$kind %in% c("ULS", "ULS (partilhada)")]) %in% b$polygons$id))
})

test_that("a click finds the area under it", {
  b <- planning_map_boundaries("municipio")
  expect_equal(planning_map_hit(b, -9.15, 38.73), "Lisboa")
  expect_equal(planning_map_hit(b, -8.68, 41.2), "Matosinhos")
  expect_null(planning_map_hit(b, -12, 36))
})

test_that("maps colour by meaning: red is worse whichever way the indicator runs", {
  expect_equal(planning_indicator_direction(c("dsr_premature", "life_expectancy_65", "ageing_index")), c("worse", "better", "neutral"))
  b <- planning_map_boundaries("nuts2", "2024")
  values <- tibble::tibble(
    indicator = c("dsr_premature", "life_expectancy_65"), direction = c("worse", "better"),
    panel = c("a", "b"), area = "Norte", value = c(1, 1), observed = c(1, 1), flag = "", digits = 1,
    benchmark_value = c(1, 1), index = c(150, 150), quintile = 5L
  )
  chart <- planning_map_chart(values, b, "index")
  expect_s3_class(chart, "planning_map")
  expect_length(chart$plots, 2)
  fill_of <- function(p) {
    built <- ggplot2::ggplot_build(p)$data[[1]]
    unique(built$fill[built$group %in% which(levels(factor(p$data$piece)) %in% p$data$piece[p$data$id == "Norte"])])
  }
  # 150 on a "worse" indicator is the red end; on a "better" one, the blue end.
  expect_equal(toupper(fill_of(chart$plots[[1]])), toupper(PLANNING_MAP_DIVERGING[["high"]]))
  expect_equal(toupper(fill_of(chart$plots[[2]])), toupper(PLANNING_MAP_DIVERGING[["low"]]))
})
