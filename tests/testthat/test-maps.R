root <- normalizePath(file.path("..", ".."))
have_sns <- file.exists(file.path(root, "data", "snapshots", "sns", "diabetes.rds"))
have_aces <- file.exists(file.path(root, "data", "boundaries", "aces.rds"))

test_that("every SNS label before 2024 is one ACES, and each ACES one label per dataset and month", {
  skip_if_not(have_sns, "SNS snapshots not present")
  crosswalk <- planning_sns_aces_crosswalk()
  expect_false(anyDuplicated(crosswalk$portal_label) > 0)
  expect_equal(length(unique(crosswalk$aces)), 55)
  for (f in list.files(file.path(root, "data", "snapshots", "sns"), full.names = TRUE)) {
    x <- readRDS(f)
    x <- x[x$period < SNS_FIRST_PERIOD, , drop = FALSE]
    if (nrow(x) == 0) next
    aces <- crosswalk$aces[match(x$unit, crosswalk$portal_label)]
    expect_false(any(is.na(aces)), info = basename(f))
    labels <- unique(data.frame(period = x$period, unit = x$unit, aces = aces))
    expect_false(anyDuplicated(labels[, c("period", "aces")]) > 0, info = basename(f))
  }
})

test_that("the ACES of the crosswalk are the boundaries' units", {
  skip_if_not(have_aces, "ACES boundaries not built")
  b <- planning_map_boundaries("aces")
  ids <- unique(stats::na.omit(b$polygons$id))
  expect_setequal(ids, unique(planning_sns_aces_crosswalk()$aces))
})

test_that("repeated SNS records count once", {
  skip_if_not(have_sns, "SNS snapshots not present")
  # The portal repeats ULS São José's records of October 2025 to January 2026;
  # summed, its HbA1c read 124.6%.
  withr::local_envvar(MORTALITY_SNAPSHOT_DIR = file.path(root, "data", "snapshots"))
  map <- planning_sns_map_values("sns_diabetes_hba1c", 2025)
  skip_if(is.null(map), "no 2025 data")
  expect_true(all(map$values$value <= 100, na.rm = TRUE))
  aces <- planning_sns_map_values("sns_diabetes_hba1c", 2019)
  expect_equal(aces$kind, "aces")
  expect_equal(nrow(aces$values), 55)
  expect_true(all(aces$values$value <= 100, na.rm = TRUE))
})

test_that("map classes and quintiles", {
  expect_equal(map_quintile(c(1:10, NA)), c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, NA))
  expect_true(all(is.na(map_quintile(c(1, 2, NA)))))
  worse <- map_class_colours(MAP_TREND_LABELS)
  better <- map_class_colours(MAP_TREND_LABELS, "better")
  expect_equal(names(worse), MAP_TREND_LABELS)
  # Higher is red for mortality and blue where higher is better.
  expect_equal(unname(worse[[3]]), PLANNING_MAP_DIVERGING[["high"]])
  expect_equal(unname(better[[3]]), PLANNING_MAP_DIVERGING[["low"]])
  expect_equal(map_significance_class(c("Superior", NA, "Inferior")), unname(MAP_SIGNIFICANCE_LABELS[c("Superior", NA, "Inferior")]))
})

test_that("the choropleth draws every scale, and the mainland alone where islands have no unit", {
  b <- list(
    polygons = data.frame(id = c(rep("A", 4), rep("B", 4), rep(NA, 4)), name = c(rep("A", 4), rep("B", 4), rep("Ilha", 4)),
                          piece = rep(c("1.1", "2.1", "3.1"), each = 4), hole = rep(c("1.1.1", "2.1.1", "3.1.1"), each = 4),
                          x = c(0, 1, 1, 0, 1, 2, 2, 1, 5, 6, 6, 5), y = c(0, 0, 1, 1, 0, 0, 1, 1, 5, 5, 6, 6)),
    frames = data.frame(group = "1", x = c(4, 7, 7, 4), y = c(4, 4, 7, 7)),
    meta = list(attribution = "teste")
  )
  v <- data.frame(area = c("A", "B"), index = c(80, 120), quintile = c(1L, 5L), z = c(-2, 3),
                  class = c(MAP_TREND_LABELS[[1]], NA), hover = c("a", "b"))
  for (scale in c("index", "quintile", "linear", "classes")) {
    p <- interactive_choropleth(v, b, scale = scale, classes = map_class_colours(MAP_TREND_LABELS), limits = c(-5, 5), ticks = c(-5, 0, 5))
    built <- plotly::plotly_build(p)$x$data
    types <- vapply(built, function(t) t$type, character(1))
    expect_true("choropleth" %in% types, info = scale)
    shown <- as.character(unlist(lapply(built[types == "choropleth"], function(t) t$locations)))
    # The island has no id: neither drawn nor framed.
    expect_setequal(shown, c("A", "B"))
    expect_false("scattergeo" %in% types)
  }
  # B has no class: drawn grey, as without a value.
  p <- plotly::plotly_build(interactive_choropleth(v, b, scale = "classes", classes = map_class_colours(MAP_TREND_LABELS)))$x$data
  grey <- Filter(function(t) identical(t$name, "Sem valor"), p)
  expect_equal(as.character(unlist(grey[[1]]$locations)), "B")
})
