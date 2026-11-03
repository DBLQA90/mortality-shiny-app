# Tests for R/parish_maps.R: the parish SMR on the snapshots (INE's parish
# deaths and the 2021 census) and the parish boundaries.
root <- normalizePath(file.path("..", ".."))
have_parish <- file.exists(file.path(root, "data", "snapshots", "parish_vitals", "deaths.rds")) &&
  file.exists(file.path(root, "data", "snapshots", "census_parish", "year_2021.rds"))

test_that("the parish map serves all causes, both sexes, the SMR of one triennium", {
  expect_null(parish_map_unavailable_reason(planning_all_causes, "HM", "smr", "period"))
  expect_match(parish_map_unavailable_reason("Diabetes mellitus", "HM", "smr", "period"), "todas as causas")
  expect_match(parish_map_unavailable_reason(planning_all_causes, "M", "smr", "period"), "Ambos")
  expect_match(parish_map_unavailable_reason(planning_all_causes, "HM", "dsr", "period"), "SMR")
  expect_match(parish_map_unavailable_reason(planning_all_causes, "HM", "smr", "trend"), "Censos 2021")
})

test_that("the parish SMR adds up to Portugal and every parish has a boundary", {
  skip_if_not(have_parish, "parish snapshots not present")
  withr::local_envvar(MORTALITY_SNAPSHOT_DIR = file.path(root, "data", "snapshots"))
  lookup <- get_nuts_lookup("2024")
  years <- parish_map_end_years()
  skip_if(length(years) == 0, "no triennium with national deaths by age")
  v <- parish_map_values(max(years), lookup, smoothed = FALSE)
  expect_equal(nrow(v), 3092)
  expect_false(anyNA(v$municipality))
  expect_equal(sum(v$expected), sum(v$observed))
  expect_true(all(v$lower <= v$value & v$value <= v$upper))
  # Smoothing keeps the observed SMR and pulls the parishes towards the mean.
  s <- parish_map_values(max(years), lookup, smoothed = TRUE)
  expect_equal(s$observed_value, v$value)
  expect_lt(stats::sd(s$value), stats::sd(v$value))
  expect_equal(parish_map_values(min(PARISH_MAP_END_YEARS) - 1L, lookup), tibble::tibble())
  # The outline of a municipality is its parishes.
  porto <- parish_map_highlight("Porto", lookup)
  expect_equal(length(porto), 7)
  b <- file.path(root, "data", "boundaries", "freguesia.rds")
  skip_if_not(file.exists(b), "parish boundaries not built")
  expect_setequal(unique(stats::na.omit(readRDS(b)$polygons$id)), v$area)
})
