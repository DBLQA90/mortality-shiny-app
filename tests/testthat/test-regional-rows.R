# Tests for R/regional_rows.R (INE regional rows versus municipal sums).

rr_lookup <- function() {
  tibble(
    municipality = c("Braga", "Fafe", "Beja", "Serpa", "Coimbra"),
    municipality_code = c("111303", "111305", "1C1401", "1C1408", "191306"),
    nuts1 = "Continente",
    nuts2 = c("Norte", "Norte", "Alentejo", "Alentejo", "Centro"),
    nuts3 = c("Cávado", "Ave", "Baixo Alentejo", "Baixo Alentejo", "Região de Coimbra")
  )
}

# Municipal frame as the loaders return it: one row per (year, area, sex, cause,
# band), population repeated per cause.
rr_frame <- function(years = 2014L, causes = "C1", deaths = 1) {
  tidyr::expand_grid(
    year = years, area = c("Braga", "Fafe", "Beja", "Serpa", "Coimbra"),
    sex = "HM", cause = causes, age_band = c("0 - 4 anos", "60 - 64 anos")
  ) %>%
    dplyr::mutate(
      deaths = deaths, pop = 1000,
      population_source = "RDS population", death_source = "0008206"
    )
}

# Regional rows written where read_regional_rows() looks for them.
with_regional_rows <- function(rows_by_file, code) {
  root <- tempfile("rr")
  for (key in names(rows_by_file)) {
    parts <- strsplit(key, "/")[[1]]
    dir <- file.path(root, "regional_deaths", parts[[1]])
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(rows_by_file[[key]], file.path(dir, paste0("year_", parts[[2]], ".rds")))
  }
  rm(list = ls(.regional_rows_cache), envir = .regional_rows_cache)
  withr::with_envvar(c(MORTALITY_SNAPSHOT_DIR = root), code)
}

regional <- function(year, code, cause = "C1", deaths = 10, sex = "HM") {
  tibble(
    year = year, region_code = code, area = "x", sex = sex, cause = cause,
    age_band = c("0 - 4 anos", "60 - 64 anos"), deaths = deaths,
    source_indicator = "i"
  )
}

test_that("stable territories have a row in every year, under either vintage", {
  for (vintage in c("2013", "2024")) {
    plan <- regional_row_plan("Norte", vintage, 1991:2024)
    expect_equal(plan$year, 1991:2024)
    # Follows the death archive's own precedence at 2022.
    expect_equal(plan$indicator[plan$year == 2021], "0008206")
    expect_equal(plan$indicator[plan$year == 2022], "0013166")
  }
  expect_equal(nrow(regional_row_plan("Região Autónoma dos Açores", "2024", 1991:2024)), 34)
  expect_equal(nrow(regional_row_plan("Continente", "2013", 1991:2024)), 34)
})

test_that("redrawn territories are composed from subregions outside their own years", {
  # NUTS 2013: own row to 2022, composed from NUTS-2024 rows after.
  old <- regional_row_plan("Alentejo", "2013", 1991:2024)
  expect_equal(old$year, 1991:2024)
  expect_equal(old$codes[old$year == 2022], "18")
  expect_equal(old$codes[old$year == 2023], "1C,1D3")

  expect_equal(regional_row_plan("Centro", "2013", 2024)$codes, "19,1D1,1D2")
  expect_equal(regional_row_plan("Área Metropolitana de Lisboa", "2013", 2024)$codes, "1A,1B")

  # NUTS 2024: own row from 2022, composed from NUTS-2013 subregions before.
  new <- regional_row_plan("Alentejo", "2024", 1991:2024)
  expect_equal(new$year, 1991:2024)
  expect_equal(new$codes[new$year == 2021], "181,184,186,187")
  expect_equal(new$codes[new$year == 2022], "1C")

  # The two compositions a municipality move keeps from being exact.
  centro <- regional_row_plan("Centro", "2024", 2014)
  expect_equal(centro$plus, "Sertã,Vila de Rei")
  ovt <- regional_row_plan("Oeste e Vale do Tejo", "2024", 2014)
  expect_equal(ovt$codes, "16B,16I,185")
  expect_equal(ovt$minus, "Sertã,Vila de Rei")

  # The Lisbon split has no NUTS-2013 equivalent before 2022.
  expect_equal(regional_row_plan("Grande Lisboa", "2024", 1991:2024)$year, 2022:2024)
  expect_equal(nrow(regional_row_plan("Grande Lisboa", "2013", 1991:2024)), 0)
  expect_equal(nrow(regional_row_plan("Área Metropolitana de Lisboa", "2024", 1991:2024)), 0)
  expect_equal(nrow(regional_row_plan("Beja", "2024", 1991:2024)), 0)
})

test_that("the declared compositions match the municipality lookups", {
  skip_if_not(file.exists("../../data/nuts_lookup_2013.rds"), "lookups not built")
  old <- readRDS("../../data/nuts_lookup_2013.rds")
  new <- readRDS("../../data/nuts_lookup_2024.rds")
  members <- function(lk, level, units) sort(unique(lk$municipality[lk[[level]] %in% units]))

  code_to_unit_2013 <- c("16B" = "Oeste", "16D" = "Região de Aveiro", "16E" = "Região de Coimbra",
    "16F" = "Região de Leiria", "16G" = "Viseu Dão Lafões", "16H" = "Beira Baixa",
    "16I" = "Médio Tejo", "16J" = "Beiras e Serra da Estrela", "181" = "Alentejo Litoral",
    "184" = "Baixo Alentejo", "185" = "Lezíria do Tejo", "186" = "Alto Alentejo", "187" = "Alentejo Central")

  compose <- function(region) {
    p <- regional_row_plan(region, "2024", 2014)
    base <- members(old, "nuts3", code_to_unit_2013[split_list(p$codes)])
    sort(setdiff(union(base, split_list(p$plus)), split_list(p$minus)))
  }

  # Each NUTS-2024 region before 2022 must be exactly its current membership.
  for (region in c("Centro", "Alentejo", "Oeste e Vale do Tejo")) {
    expect_equal(compose(region), members(new, "nuts2", region), info = region)
  }

  # NUTS-2013 regions after 2022, from NUTS-2024 units.
  expect_equal(
    members(old, "nuts2", "Centro"),
    sort(c(members(new, "nuts2", "Centro"), members(new, "nuts3", c("Oeste", "Médio Tejo"))))
  )
  expect_equal(
    members(old, "nuts2", "Alentejo"),
    sort(c(members(new, "nuts2", "Alentejo"), members(new, "nuts3", "Lezíria do Tejo")))
  )
})

test_that("a region's municipal deaths are replaced by its INE row", {
  with_regional_rows(list("0008206/2014" = regional(2014L, "11", deaths = 10)), {
    out <- substitute_regional_deaths(rr_frame(), "Norte", "2024", rr_lookup())

    norte <- out[out$area == "Norte", ]
    expect_equal(nrow(norte), 2)
    expect_equal(sum(norte$deaths), 20)
    # Population is the sum of the members, untouched by the substitution.
    expect_equal(norte$pop, c(2000, 2000))
    expect_match(unique(norte$death_source), "linha regional")

    # The members are gone for that year; everything else is as it was.
    expect_false(any(out$area %in% c("Braga", "Fafe")))
    expect_equal(sum(out$deaths[out$area %in% c("Beja", "Serpa", "Coimbra")]), 6)
    expect_equal(nrow(attr(out, "regional_substitutions")), 1)
  })
})

test_that("a composition sums subregions and applies municipal corrections", {
  frame <- rr_frame() %>% dplyr::mutate(area = dplyr::recode(area, Braga = "Sertã", Fafe = "Coimbra2"))
  lk <- rr_lookup() %>% dplyr::mutate(
    municipality = dplyr::recode(municipality, Braga = "Sertã", Fafe = "Coimbra2"),
    nuts2 = dplyr::if_else(municipality %in% c("Sertã", "Coimbra2", "Coimbra"), "Centro", nuts2)
  )
  rows <- dplyr::bind_rows(
    regional(2014L, "16D", deaths = 3), regional(2014L, "16E", deaths = 4),
    regional(2014L, "16F", deaths = 0), regional(2014L, "16G", deaths = 0),
    regional(2014L, "16H", deaths = 0), regional(2014L, "16J", deaths = 0)
  )
  with_regional_rows(list("0008206/2014" = rows), {
    out <- substitute_regional_deaths(frame, "Centro", "2024", lk)
    centro <- out[out$area == "Centro", ]
    # (3 + 4) from the subregions, plus Sertã's own row of 1, per band.
    expect_equal(centro$deaths, c(8, 8))
    expect_match(unique(centro$death_source), "compostas")
  })
})

test_that("a subtraction uses the loader and falls back without it", {
  lk <- rr_lookup() %>% dplyr::mutate(nuts2 = dplyr::if_else(municipality == "Coimbra", "Oeste e Vale do Tejo", nuts2))
  rows <- dplyr::bind_rows(regional(2014L, "16B", deaths = 5), regional(2014L, "16I", deaths = 5), regional(2014L, "185", deaths = 5))
  loader <- function(areas, years, causes) {
    tidyr::expand_grid(year = years, area = areas, sex = "HM", cause = causes, age_band = c("0 - 4 anos", "60 - 64 anos")) %>%
      dplyr::mutate(deaths = 2)
  }
  with_regional_rows(list("0008206/2014" = rows), {
    out <- substitute_regional_deaths(rr_frame(), "Oeste e Vale do Tejo", "2024", lk, load_municipal_deaths = loader)
    ovt <- out[out$area == "Oeste e Vale do Tejo", ]
    # 15 from three subregions, minus Sertã and Vila de Rei at 2 each.
    expect_equal(ovt$deaths, c(11, 11))

    # Without a loader the composition cannot be completed, so it is skipped
    # rather than applied without its correction.
    plain <- substitute_regional_deaths(rr_frame(), "Oeste e Vale do Tejo", "2024", lk)
    expect_false(any(plain$area == "Oeste e Vale do Tejo"))
  })
})

test_that("a composition missing one of its subregions is not applied", {
  lk <- rr_lookup() %>% dplyr::mutate(nuts2 = dplyr::if_else(nuts2 == "Alentejo", "Alentejo", nuts2))
  partial <- dplyr::bind_rows(regional(2014L, "181"), regional(2014L, "184"), regional(2014L, "186"))
  with_regional_rows(list("0008206/2014" = partial), {
    out <- substitute_regional_deaths(rr_frame(), "Alentejo", "2024", lk)
    expect_false(any(out$area == "Alentejo"))
  })
})

test_that("population is deduplicated across the causes in the frame", {
  frame <- rr_frame(causes = c("C1", "C2"))
  rows <- dplyr::bind_rows(regional(2014L, "11", "C1", 10), regional(2014L, "11", "C2", 4))
  with_regional_rows(list("0008206/2014" = rows), {
    out <- substitute_regional_deaths(frame, "Norte", "2013", rr_lookup())
    norte <- out[out$area == "Norte", ]
    # Two municipalities of 1000, not four rows of 1000.
    expect_true(all(norte$pop == 2000))
    expect_equal(sum(norte$deaths[norte$cause == "C2"]), 8)
  })
})

test_that("years whose rows are missing keep the municipal sum", {
  with_regional_rows(list("0008206/2021" = regional(2021L, "18")), {
    frame <- rr_frame(years = c(2021L, 2023L))
    out <- substitute_regional_deaths(frame, "Alentejo", "2013", rr_lookup())
    # 2021 from the regional row; the 2023 file is absent, so Beja and Serpa stay.
    expect_true(any(out$area == "Alentejo" & out$year == 2021))
    expect_true(all(c("Beja", "Serpa") %in% out$area[out$year == 2023]))
    expect_false(any(out$area == "Alentejo" & out$year == 2023))
  })
})

test_that("a missing file falls back to the municipal sum rather than failing", {
  with_regional_rows(list(), {
    frame <- rr_frame()
    out <- substitute_regional_deaths(frame, "Norte", "2024", rr_lookup())
    expect_equal(sum(out$deaths), sum(frame$deaths))
    expect_equal(nrow(attr(out, "regional_substitutions")), 0)
  })
})

test_that("the municipal-sum source leaves the frame alone", {
  with_regional_rows(list("0008206/2014" = regional(2014L, "11")), {
    frame <- rr_frame()
    out <- substitute_regional_deaths(frame, "Norte", "2024", rr_lookup(), source = "municipal_sum")
    expect_equal(sum(out$deaths), sum(frame$deaths))
    expect_false(any(out$area == "Norte"))
  })
})

test_that("a nested selection substitutes the containing region only", {
  rows <- dplyr::bind_rows(regional(2014L, "1", deaths = 50), regional(2014L, "11", deaths = 10))
  with_regional_rows(list("0008206/2014" = rows), {
    out <- substitute_regional_deaths(rr_frame(), c("Norte", "Continente"), "2024", rr_lookup())
    # Continente covers Norte's municipalities, so Norte must not be added again.
    expect_equal(unique(out$area), "Continente")
    expect_equal(sum(out$deaths), 100)
    expect_equal(attr(out, "regional_substitutions")$region, "Continente")
  })
})

test_that("rows are matched by code, not by label", {
  # A row carrying the right label but the wrong code must not be used - the
  # reason these rows are fetched by code at all.
  impostor <- regional(2014L, "99", deaths = 999) %>% dplyr::mutate(area = "Norte")
  with_regional_rows(list("0008206/2014" = impostor), {
    out <- substitute_regional_deaths(rr_frame(), "Norte", "2024", rr_lookup())
    expect_false(any(out$deaths == 999))
  })
})

test_that("the seam warning is left only for the Lisbon split", {
  msg <- region_source_seam_warning("Grande Lisboa", "2024", 2019:2024)
  expect_match(msg, "Grande Lisboa \\(2022-2024\\)")

  # Composed regions are covered in every year, so they have no seam.
  for (region in c("Norte", "Centro", "Alentejo", "Oeste e Vale do Tejo")) {
    expect_null(region_source_seam_warning(region, "2024", 1991:2024), info = region)
  }
  expect_null(region_source_seam_warning("Área Metropolitana de Lisboa", "2013", 1991:2024))
  expect_null(region_source_seam_warning("Grande Lisboa", "2024", 2019:2024, source = "municipal_sum"))
  expect_null(region_source_seam_warning(character(0), "2024", 2019:2024))
})

test_that("standalone municipalities are flagged, regions are not", {
  lk <- rr_lookup()
  expect_match(municipal_age_detail_warning("Beja", lk, "C1", 2019), "subestimados")
  expect_match(municipal_age_detail_warning("Beja", lk, "Todas as causas de morte", 2014), "2014")
  expect_null(municipal_age_detail_warning("Beja", lk, "Todas as causas de morte", 2019))
  expect_null(municipal_age_detail_warning(c("Norte", "Portugal"), lk, "C1", 2014))
})

# --- INE's regional rows inside the standardised module -------------------

# A snapshot root holding both the regional rows and one municipal death file,
# which is what tells planning_death_indicator() which edition a year is on.
with_regional_snapshot <- function(rows, year, indicator, code) {
  root <- tempfile("rrs")
  dir.create(file.path(root, "regional_deaths", indicator), recursive = TRUE, showWarnings = FALSE)
  saveRDS(rows, file.path(root, "regional_deaths", indicator, paste0("year_", year, ".rds")))
  dir.create(file.path(root, "deaths", indicator, paste0("year_", year)), recursive = TRUE, showWarnings = FALSE)
  saveRDS(rr_frame(year, planning_all_causes),
          file.path(root, "deaths", indicator, paste0("year_", year), "cause_todas_as_causas_de_morte.rds"))
  rm(list = ls(.regional_rows_cache), envir = .regional_rows_cache)
  rm(list = ls(planning_cache), envir = planning_cache)
  on.exit(rm(list = ls(planning_cache), envir = planning_cache), add = TRUE)
  withr::with_envvar(c(MORTALITY_SNAPSHOT_DIR = root), code)
}

# Deaths by band and measure for a handful of municipalities, in the shape
# planning_cause_age_block() returns.
rr_block <- function(municipalities, per_band = 1) {
  measures <- c(names(planning_standardised_measures()), "avoidable")
  deaths <- array(0, c(length(municipalities), length(age_levels), length(measures)),
                  dimnames = list(municipalities, age_levels, measures))
  deaths[, c("0 - 4 anos", "60 - 64 anos"), "all"] <- per_band
  list(deaths = deaths)
}

test_that("the NUTS vintage is recovered from the lookup itself", {
  skip_if_not(file.exists("../../data/nuts_lookup_2013.rds"), "lookups not built")
  expect_equal(planning_lookup_vintage(get_nuts_lookup("2024")), "2024")
  expect_equal(planning_lookup_vintage(get_nuts_lookup("2013")), "2013")
})

test_that("a region's own rows replace the sum of its municipalities", {
  year <- 2014L
  plan <- regional_row_plan("Centro", "2024", year)
  ovt <- regional_row_plan("Oeste e Vale do Tejo", "2024", year)
  codes <- split_list(plan$codes)
  rows <- dplyr::bind_rows(lapply(union(codes, split_list(ovt$codes)), function(code)
    regional(year, code, cause = planning_all_causes, deaths = 10)))

  with_regional_snapshot(rows, year, plan$indicator, {
    block <- rr_block(c("Sertã", "Vila de Rei"), per_band = 3)
    total <- planning_regional_area_deaths("Centro", year, "HM", block, "2024")
    # Every code in the plan, plus the two municipalities it has to add.
    expect_equal(unname(total["0 - 4 anos", "all"]), 10 * length(codes) + 6)
    expect_equal(unname(total["30 - 34 anos", "all"]), 0)
    # Measures with no rows of their own stay at zero rather than going missing.
    expect_equal(unname(total["0 - 4 anos", "preventable"]), 0)

    # Oeste e Vale do Tejo takes the same two municipalities away.
    expect_equal(
      unname(planning_regional_area_deaths("Oeste e Vale do Tejo", year, "HM", block, "2024")["0 - 4 anos", "all"]),
      10 * length(split_list(ovt$codes)) - 6
    )

    # Life expectancy asks for the all-cause vector from its own matrix.
    life_block <- list(deaths = rr_block(c("Sertã", "Vila de Rei"), per_band = 3)$deaths[, , "all"])
    expect_equal(unname(planning_regional_area_deaths("Centro", year, "HM", life_block, "2024", measure = "all")[["0 - 4 anos"]]),
                 10 * length(codes) + 6)

    # Portugal and Continente keep their published rows; a municipality has none.
    expect_null(planning_regional_area_deaths("Portugal", year, "HM", block, "2024"))
    expect_null(planning_regional_area_deaths("Beja", year, "HM", block, "2024"))
  })
})

test_that("rows from another edition than the year's deaths are not used", {
  year <- 2014L
  plan <- regional_row_plan("Norte", "2024", year)
  rows <- regional(year, "11", cause = planning_all_causes, deaths = 10)
  # The year's municipal deaths are on 0013166, the plan on 0008206: no substitution.
  with_regional_snapshot(rows, year, "0013166", {
    expect_equal(plan$indicator, "0008206")
    expect_null(planning_regional_area_deaths("Norte", year, "HM", rr_block("Braga"), "2024"))
  })
})
