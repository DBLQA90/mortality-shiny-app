# Tests for R/regions.R (NUTS vintage handling and municipal region rebuilds).

test_lookup <- function() {
  tibble(
    municipality = c("Beja", "Serpa", "Santarém", "Almeirim", "Braga", "Fafe", "Calheta (R.A.M.)"),
    municipality_code = c("1C1401", "1C1408", "1D1418", "1D1402", "111303", "111305", "3003"),
    nuts1 = c(rep("Continente", 6), "Região Autónoma da Madeira"),
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

lookup_file <- function(vintage) sprintf("../../data/nuts_lookup_%s.rds", vintage)

test_that("the committed lookup covers the mainland regions completely", {
  skip_if_not(file.exists(lookup_file("2024")), "lookup not built")
  lookup <- readRDS(lookup_file("2024"))

  expect_gt(nrow(lookup), 300)
  expect_true(all(c("Norte", "Centro", "Alentejo", "Algarve") %in% lookup$nuts2))
  # Every municipality resolves to a NUTS II region.
  expect_equal(sum(is.na(lookup$nuts2)), 0)
  # The Lezíria move is present in the real data, not just the fixture.
  expect_equal(lookup$nuts2[lookup$municipality == "Santarém"], "Oeste e Vale do Tejo")
})

test_that("both NUTS vintages cover the same municipalities under the same labels", {
  skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")
  skip_if_not(file.exists(lookup_file("2024")), "2024 lookup not built")

  old <- readRDS(lookup_file("2013"))
  new <- readRDS(lookup_file("2024"))

  # This is what makes the vintage a regrouping rather than a change of
  # coverage: switching it must never add or drop a place, only move it.
  expect_equal(sort(old$municipality), sort(new$municipality))
  expect_equal(nrow(old), 308)
  expect_equal(sum(is.na(old$nuts2)), 0)

  # NUTS-2013 has seven NUTS II regions, NUTS-2024 has nine.
  expect_equal(length(unique(old$nuts2)), 7)
  expect_equal(length(unique(new$nuts2)), 9)
  expect_true("Área Metropolitana de Lisboa" %in% old$nuts2)
  expect_false("Área Metropolitana de Lisboa" %in% new$nuts2)
  expect_false("Oeste e Vale do Tejo" %in% old$nuts2)

  # The 2024 reform, checked against the real files rather than asserted: the
  # Lezíria leaves Alentejo, Oeste and Médio Tejo leave Centro, and the Lisbon
  # metropolitan area splits in two.
  moved <- merge(
    data.frame(municipality = old$municipality, old = old$nuts2),
    data.frame(municipality = new$municipality, new = new$nuts2),
    by = "municipality"
  )
  expect_equal(sum(moved$old == "Alentejo" & moved$new == "Oeste e Vale do Tejo"), 11)
  expect_equal(sum(moved$old == "Centro" & moved$new == "Oeste e Vale do Tejo"), 23)
  expect_equal(sum(moved$old == "Área Metropolitana de Lisboa" & moved$new == "Grande Lisboa"), 9)
  expect_equal(sum(moved$old == "Área Metropolitana de Lisboa" & moved$new == "Península de Setúbal"), 9)

  # Norte is untouched by the reform, which is why it is safe as the app's
  # fixed second column.
  expect_equal(sum(moved$old == "Norte"), sum(moved$new == "Norte"))
})

test_that("the vintage selects which lookup is read", {
  skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")

  withr::with_envvar(c(MORTALITY_NUTS_VINTAGE = ""), {
    expect_equal(default_nuts_vintage(), "2024")
  })
  withr::with_envvar(c(MORTALITY_NUTS_VINTAGE = "2013"), {
    expect_equal(default_nuts_vintage(), "2013")
  })
  # An unrecognised vintage must fall back to the current definition rather
  # than silently reading no lookup at all.
  expect_equal(normalize_nuts_vintage("1998"), "2024")
  expect_equal(normalize_nuts_vintage(NA), "2024")

  old <- readRDS(lookup_file("2013"))
  new <- readRDS(lookup_file("2024"))

  # Santarém is the clearest case: Alentejo under 2013, Oeste e Vale do Tejo
  # under 2024.
  expect_equal(region_municipalities("Alentejo", old), sort(old$municipality[old$nuts2 == "Alentejo"]))
  expect_true("Santarém" %in% region_municipalities("Alentejo", old))
  expect_false("Santarém" %in% region_municipalities("Alentejo", new))

  # Six names exist in both vintages and mean different things in each.
  expect_length(region_municipalities("Centro", old), 100)
  expect_length(region_municipalities("Centro", new), 77)
})

test_that("region and area choices follow the vintage", {
  skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")

  old <- readRDS(lookup_file("2013"))
  new <- readRDS(lookup_file("2024"))

  # 7 and 9 NUTS II regions, plus Continente from NUTS I. The islands are named
  # identically at both levels, so they are not counted twice.
  expect_length(region_choices_for("2013", old), 8)
  expect_length(region_choices_for("2024", new), 10)
  # Coarsest level first, then INE's own ordering from the geography codes.
  expect_equal(region_choices_for("2024", new)[[1]], "Continente")
  expect_equal(region_choices_for("2024", new)[[4]], "Norte")

  areas_new <- area_choices_for("2024", new)
  expect_equal(areas_new[[1]], "Portugal")
  expect_equal(length(areas_new), 1 + 10 + 308)
  expect_true(all(new$municipality %in% areas_new))

  # A region dropped by the reform is still a known name somewhere, so a
  # selection carried over from the other vintage can be recognised and
  # reported rather than treated as an unknown area.
  expect_true("Área Metropolitana de Lisboa" %in% all_region_choices())
  expect_true("Grande Lisboa" %in% all_region_choices())
})

test_that("NUTS I regions are rebuilt like any other", {
  skip_if_not(file.exists(lookup_file("2024")), "lookup not built")

  for (vintage in c("2013", "2024")) {
    lookup <- readRDS(lookup_file(vintage))
    skip_if_not("nuts1" %in% names(lookup), "lookup has no NUTS I column")

    # Three NUTS I units: the mainland and the two autonomous regions.
    expect_setequal(
      unique(lookup$nuts1),
      c("Continente", "Região Autónoma dos Açores", "Região Autónoma da Madeira")
    )

    # Continente is the mainland: everything except the 30 island municipalities.
    mainland <- region_municipalities("Continente", lookup)
    expect_length(mainland, 278)
    expect_false("Funchal" %in% mainland)
    expect_false("Ponta Delgada" %in% mainland)
    expect_true("Lisboa" %in% mainland)

    # The three NUTS I units partition the country exactly.
    islands <- c(
      region_municipalities("Região Autónoma dos Açores", lookup),
      region_municipalities("Região Autónoma da Madeira", lookup)
    )
    expect_length(islands, 30)
    expect_setequal(c(mainland, islands), lookup$municipality)

    # The islands are the same territory at NUTS I, II and III, so they resolve
    # to the same municipalities however they are named - and appear once in the
    # selector rather than three times.
    expect_equal(sum(region_choices_for(vintage, lookup) == "Região Autónoma dos Açores"), 1)
  }

  # Continente is the only label NUTS I adds to the selector.
  new <- readRDS(lookup_file("2024"))
  skip_if_not("nuts1" %in% names(new), "lookup has no NUTS I column")
  expect_equal(
    setdiff(region_choices_for("2024", new), region_choices_for("2024", new, columns = "nuts2")),
    "Continente"
  )
  expect_length(region_choices_for("2024", new), 10)
  expect_length(region_choices_for("2013", readRDS(lookup_file("2013"))), 8)
})

test_that("containment between two selected regions is flagged", {
  lookup <- test_lookup()

  # The case that adding NUTS I makes reachable. Comparing a region's
  # municipalities against the other *selected areas* never caught this,
  # because Norte is a region label and so is in nobody's municipality list.
  msg <- overlapping_selection_warning(c("Continente", "Norte"), lookup)
  expect_match(msg, "Continente já inclui Norte")
  expect_match(msg, "absorvida")

  # A NUTS II region alongside one of its own NUTS III children.
  expect_match(
    overlapping_selection_warning(c("Alentejo", "Baixo Alentejo"), lookup),
    "Alentejo já inclui Baixo Alentejo"
  )

  # Region plus a municipality of a different region is still fine, and two
  # disjoint regions at the same level are still fine.
  expect_null(overlapping_selection_warning(c("Continente", "Região Autónoma da Madeira"), lookup))
  expect_null(overlapping_selection_warning(c("Alentejo", "Norte"), lookup))
  expect_null(overlapping_selection_warning(c("Norte", "Beja"), lookup))

  # Mixed: a contained region and a contained municipality reported together.
  both <- overlapping_selection_warning(c("Continente", "Norte", "Beja"), lookup)
  expect_match(both, "Beja")
  expect_match(both, "Norte")
})

test_that("switching vintage keeps municipalities and drops absent regions", {
  skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")

  old <- readRDS(lookup_file("2013"))
  new <- readRDS(lookup_file("2024"))

  # 2024 -> 2013. Municipalities are identical across vintages, so they survive;
  # Grande Lisboa has no NUTS-2013 counterpart and is reported.
  to_old <- reconcile_area_selection(c("Grande Lisboa", "Beja", "Braga"), "2013", old)
  expect_equal(to_old$selected, c("Beja", "Braga"))
  expect_equal(to_old$dropped, "Grande Lisboa")
  expect_true("Área Metropolitana de Lisboa" %in% to_old$choices)
  expect_false("Grande Lisboa" %in% to_old$choices)

  # 2013 -> 2024, the other direction.
  to_new <- reconcile_area_selection(c("Área Metropolitana de Lisboa", "Beja"), "2024", new)
  expect_equal(to_new$selected, "Beja")
  expect_equal(to_new$dropped, "Área Metropolitana de Lisboa")

  # A name meaning something in both vintages is kept, because it still selects
  # a real region - a different one, which is what the switch is for.
  both <- reconcile_area_selection(c("Centro", "Alentejo"), "2013", old)
  expect_equal(both$selected, c("Centro", "Alentejo"))
  expect_length(both$dropped, 0)

  # The annual tab's third column excludes the two fixed comparator columns.
  annual <- reconcile_area_selection("Beja", "2024", new, exclude = c("Portugal", "Norte"))
  expect_false("Portugal" %in% annual$choices)
  expect_false("Norte" %in% annual$choices)
  expect_equal(annual$selected, "Beja")

  # An empty selection is not an error.
  expect_length(reconcile_area_selection(NULL, "2024", new)$selected, 0)
})

test_that("a region from the other vintage is refused, not passed through", {
  skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")

  old <- readRDS(lookup_file("2013"))
  new <- readRDS(lookup_file("2024"))

  # Passed through untouched, such a name reaches the loader as an ordinary
  # place and fails with a storage error about a file that was never missing.
  expect_equal(unknown_region_areas("Oeste e Vale do Tejo", old), "Oeste e Vale do Tejo")
  expect_equal(unknown_region_areas("Área Metropolitana de Lisboa", new),
               "Área Metropolitana de Lisboa")

  # Names valid under the active vintage, and ordinary municipalities, pass.
  expect_length(unknown_region_areas(c("Alentejo", "Norte", "Beja", "Portugal"), old), 0)
  expect_length(unknown_region_areas("Oeste e Vale do Tejo", new), 0)

  msg <- vintage_mismatch_message("Grande Lisboa", "2013", old)
  expect_match(msg, "Grande Lisboa")
  expect_match(msg, "NUTS 2013")
  expect_null(vintage_mismatch_message("Alentejo", "2013", old))

  # The "original" escape hatch passes labels straight to INE, which still
  # publishes the other vintage's regional rows, so it must not be narrowed.
  expect_null(vintage_mismatch_message("Grande Lisboa", "2013", old, region_mode = "original"))
})

test_that("overlapping area selections are flagged, disjoint ones are not", {
  lookup <- test_lookup()

  # Two distinct NUTS II regions are disjoint: summing them is legitimate.
  expect_null(overlapping_selection_warning(c("Alentejo", "Norte"), lookup))
  expect_null(overlapping_selection_warning(c("Beja", "Braga"), lookup))
  expect_null(overlapping_selection_warning("Portugal", lookup))

  # Portugal is never expanded into municipalities, so its row is summed whole
  # and the overlap really is counted twice.
  pt <- overlapping_selection_warning(c("Portugal", "Beja"), lookup)
  expect_match(pt, "Portugal já inclui")
  expect_match(pt, "duas vezes")

  # A region is expanded, and the expansion is a unique set, so the contained
  # area is absorbed rather than double-counted. Saying "counted twice" here
  # would describe the wrong failure.
  msg <- overlapping_selection_warning(c("Alentejo", "Beja"), lookup)
  expect_match(msg, "Alentejo já inclui Beja")
  expect_match(msg, "absorvida")
  expect_no_match(msg, "seriam contados duas vezes")

  # A region alongside a municipality that is not part of it is fine.
  expect_null(overlapping_selection_warning(c("Alentejo", "Braga"), lookup))
})

test_that("the default region mode aggregates from municipalities", {
  withr::with_envvar(c(MORTALITY_REGION_MODE = ""), {
    expect_equal(default_region_mode(), "municipal")
  })
  # Documented escape hatch for reproducing a published INE regional figure.
  withr::with_envvar(c(MORTALITY_REGION_MODE = "original"), {
    expect_equal(default_region_mode(), "original")
  })
  # An unrecognised value must not silently disable the aggregation.
  withr::with_envvar(c(MORTALITY_REGION_MODE = "lixo"), {
    expect_equal(default_region_mode(), "municipal")
  })
})

test_that("municipality membership is stable across the whole series", {
  skip_if_not(file.exists(lookup_file("2024")), "lookup not built")
  skip_if_not(dir.exists("../../data/snapshots/population"), "snapshots not present")

  lookup <- readRDS(lookup_file("2024"))
  counts <- vapply(c(1991L, 2005L, 2013L, 2014L, 2023L), function(y) {
    path <- sprintf("../../data/snapshots/population/year_%d.rds", y)
    if (!file.exists(path)) return(NA_integer_)
    length(intersect(unique(readRDS(path)$area), lookup$municipality))
  }, integer(1))

  # Rebuilding a region from one fixed municipality list is only valid because
  # the municipalities themselves never change across the archive.
  expect_true(all(stats::na.omit(counts) == nrow(lookup)))
})
