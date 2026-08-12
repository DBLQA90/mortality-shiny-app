# Tests for R/regions.R (NUTS vintage handling and municipal region rebuilds).

test_lookup <- function() {
  tibble(
    municipality = c("Beja", "Serpa", "Santarém", "Almeirim", "Braga", "Fafe", "Calheta (R.A.M.)"),
    municipality_code = c("1C1401", "1C1408", "1D1418", "1D1402", "111303", "111305", "3003"),
    nuts2 = c("Alentejo", "Alentejo", "Oeste e Vale do Tejo", "Oeste e Vale do Tejo",
              "Norte", "Norte", "Região Autónoma da Madeira"),
    nuts3 = c("Baixo Alentejo", "Baixo Alentejo", "Lezíria do Tejo", "Lezíria do Tejo",
              "Cávado", "Ave", "Região Autónoma da Madeira")
  )
}

test_that("regions are recognised at both NUTS levels", {
  lookup <- test_lookup()

  expect_true(is_region_label("Alentejo", lookup))
  expect_true(is_region_label("Lezíria do Tejo", lookup))
  expect_false(is_region_label("Beja", lookup))
  expect_false(is_region_label("Portugal", lookup))

  expect_equal(region_municipalities("Alentejo", lookup), c("Beja", "Serpa"))
  expect_equal(region_municipalities("Cávado", lookup), "Braga")
})

test_that("Lezíria do Tejo sits outside Alentejo under NUTS-2024", {
  lookup <- test_lookup()

  # This is the whole reason the region mode exists: under NUTS-2013 Santarém
  # would have been part of Alentejo, and mixing the vintages is what produced
  # the ~30% understatement.
  expect_false("Santarém" %in% region_municipalities("Alentejo", lookup))
  expect_true("Santarém" %in% region_municipalities("Oeste e Vale do Tejo", lookup))
})

test_that("municipal mode expands regions and leaves other areas alone", {
  lookup <- test_lookup()

  res <- resolve_region_areas(c("Alentejo", "Braga"), region_mode = "municipal", lookup = lookup)
  expect_equal(res$areas, c("Beja", "Braga", "Serpa"))
  expect_equal(res$expanded, "Alentejo")
  expect_length(res$warnings, 0)

  # Original mode is a pass-through: INE's own regional rows are used as-is.
  keep <- resolve_region_areas(c("Alentejo", "Braga"), region_mode = "original", lookup = lookup)
  expect_equal(keep$areas, c("Alentejo", "Braga"))
  expect_length(keep$expanded, 0)

  # An unknown mode must not silently expand anything.
  expect_equal(
    resolve_region_areas("Alentejo", region_mode = "lixo", lookup = lookup)$areas,
    "Alentejo"
  )
})

test_that("incomplete municipal coverage is reported, not silently undercounted", {
  lookup <- test_lookup()
  available <- c("Beja", "Braga", "Fafe", "Santarém", "Almeirim")

  cov <- region_coverage("Alentejo", lookup, available_areas = available)
  expect_equal(cov$expected, 2)
  expect_equal(cov$found, 1)
  expect_equal(cov$missing, "Serpa")
  expect_false(cov$complete)
  expect_match(region_coverage_warning(cov), "Serpa")
  expect_match(region_coverage_warning(cov), "subestimam")

  res <- resolve_region_areas("Alentejo", region_mode = "municipal", lookup = lookup, available_areas = available)
  expect_equal(res$areas, "Beja")
  expect_length(res$warnings, 1)

  # A complete region produces no warning at all.
  full <- region_coverage("Norte", lookup, available_areas = available)
  expect_true(full$complete)
  expect_null(region_coverage_warning(full))
})

test_that("the vintage warning fires only for affected regions spanning 2022", {
  # Alentejo changed definition in 2022; Norte and Portugal did not.
  expect_match(region_vintage_warning("Alentejo", years = 2019:2023), "quebra")
  expect_null(region_vintage_warning("Norte", years = 2019:2023))
  expect_null(region_vintage_warning("Portugal", years = 2019:2023))

  # A window entirely on one side of the break has no discontinuity to warn about.
  expect_null(region_vintage_warning("Alentejo", years = 2015:2021))
  expect_null(region_vintage_warning("Alentejo", years = 2022:2024))

  # Rebuilding from municipalities removes the problem, so no warning.
  expect_null(region_vintage_warning("Alentejo", years = 2019:2023, region_mode = "municipal"))
})

test_that("the committed lookup covers the mainland regions completely", {
  skip_if_not(file.exists("../../data/nuts_lookup.rds"), "lookup not built")
  lookup <- readRDS("../../data/nuts_lookup.rds")

  expect_gt(nrow(lookup), 300)
  expect_true(all(c("Norte", "Centro", "Alentejo", "Algarve") %in% lookup$nuts2))
  # Every municipality resolves to a NUTS II region.
  expect_equal(sum(is.na(lookup$nuts2)), 0)
  # The Lezíria move is present in the real data, not just the fixture.
  expect_equal(lookup$nuts2[lookup$municipality == "Santarém"], "Oeste e Vale do Tejo")
})

test_that("overlapping area selections are flagged, disjoint ones are not", {
  lookup <- test_lookup()

  # Two distinct NUTS II regions are disjoint: summing them is legitimate.
  expect_null(overlapping_selection_warning(c("Alentejo", "Norte"), lookup))
  expect_null(overlapping_selection_warning(c("Beja", "Braga"), lookup))
  expect_null(overlapping_selection_warning("Portugal", lookup))

  # Portugal already contains everything else.
  expect_match(overlapping_selection_warning(c("Portugal", "Beja"), lookup), "Portugal já inclui")

  # A region alongside one of its own municipalities double-counts it.
  msg <- overlapping_selection_warning(c("Alentejo", "Beja"), lookup)
  expect_match(msg, "Alentejo já inclui Beja")
  expect_match(msg, "duas vezes")

  # A region alongside a municipality that is not part of it is fine.
  expect_null(overlapping_selection_warning(c("Alentejo", "Braga"), lookup))
})
