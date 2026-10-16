# Tests for the gaps INE leaves in the data, filled the same way in every tab:
# deaths without a published age, values suppressed for a municipality, and the
# latest year, which has all-cause deaths by age for regions only.
# They read the committed snapshots.

lookup <- get_nuts_lookup("2024")

# The runner's working directory is tests/testthat: point the snapshot readers
# at the repository copy, as test-infant.R does.
with_snapshots <- function(code) {
  withr::with_envvar(c(MORTALITY_SNAPSHOT_DIR = normalizePath("../../data/snapshots", mustWork = FALSE)), code)
}

test_that("deaths without a published age are filled to the municipality's total", with_snapshots({
  deaths <- readRDS(planning_death_cause_file(2024L, "Todas as causas de morte"))
  deaths <- deaths[deaths$area %in% c("Barrancos", "Portugal", "Lisboa"), ]
  filled <- fill_death_ages(deaths)
  total <- function(frame, area) sum(frame$deaths[frame$area == area & frame$sex == "HM"])
  # Barrancos 2024: 19 deaths in the age bands, 32 in INE's all-ages total.
  expect_equal(total(deaths, "Barrancos"), 19)
  expect_equal(total(filled, "Barrancos"), 32, tolerance = 1e-6)
  expect_equal(total(filled, "Portugal"), 118396, tolerance = 1e-6)
  expect_equal(sum(filled$deaths_filled[filled$area == "Barrancos" & filled$sex == "HM"]), 13, tolerance = 1e-6)
  # Nothing is invented where the bands already add up.
  expect_equal(total(filled, "Lisboa"), 6570, tolerance = 1e-6)
}))

test_that("employees INE suppressed come back from the published totals", with_snapshots({
  t <- planning_indicator_table(c("Norte", "São João da Madeira", "Vale de Cambra", "Região Autónoma dos Açores"),
                                c(2013L, 2015L), ids = "employees", lookup = lookup)
  value <- function(area, year) t$value[t$area == area & t$year == year]
  flag <- function(area, year) t$flag[t$area == area & t$year == year]
  # 2015: São João da Madeira alone was suppressed, so it gets the whole residual.
  expect_equal(value("São João da Madeira", 2015L), 9820)
  expect_equal(flag("São João da Madeira", 2015L), PLANNING_ESTIMATED_FLAG)
  # 2013: it and Vale de Cambra share one; the Norte holds both, neither alone.
  expect_true(is.finite(value("Norte", 2013L)))
  expect_true(is.na(value("São João da Madeira", 2013L)))
  expect_true(is.na(value("Vale de Cambra", 2013L)))
  expect_equal(value("Região Autónoma dos Açores", 2013L), 33786)
}))

test_that("pensioners are not recovered: INE's totals leave the suppressed municipalities out", with_snapshots({
  t <- planning_indicator_table("Norte", 2005L, ids = "pensioners", lookup = lookup)
  expect_true(is.na(t$value))
}))

test_that("the latest year serves whole regions and all causes only", with_snapshots({
  lvt <- lookup$municipality[lookup$nuts2 %in% c("Grande Lisboa", "Península de Setúbal", "Oeste e Vale do Tejo")]
  expect_true(calendar_selection_supported(lookup$municipality[lookup$nuts2 == "Norte"]))
  expect_true(calendar_selection_supported(lvt))
  expect_true(calendar_selection_supported("Portugal"))
  expect_true(calendar_selection_supported("Alentejo"))
  expect_false(calendar_selection_supported("Lisboa"))
  expect_false(calendar_selection_supported(lvt[-1]))
  expect_false(calendar_selection_supported("Portugal", causes = "Tumores malignos"))
  expect_false(calendar_selection_supported("Portugal", sex = "H"))
}))

test_that("the latest year's deaths shared among municipalities add back to each region", with_snapshots({
  year <- planning_calendar_years()
  skip_if(length(year) == 0, "no rebuilt year in the snapshots")
  year <- max(year)
  calendar <- planning_calendar_deaths(year)
  municipal <- calendar_municipal_deaths(year)
  for (region in c("Norte", "Alentejo", "Região Autónoma dos Açores")) {
    members <- lookup$municipality[lookup$nuts2 == region]
    expect_equal(sum(municipal$deaths[municipal$area %in% members]), sum(calendar$deaths[calendar$area == region]), tolerance = 1e-6)
  }
  # Only the latest year is asked for, so the by-cause archive is not read.
  loaded <- load_snapshot_deaths(year, c("Portugal", lookup$municipality[lookup$nuts2 == "Algarve"]), "Todas as causas de morte")
  expect_equal(sum(loaded$deaths[loaded$area == "Portugal"]), sum(calendar$deaths[calendar$area == "Portugal"]), tolerance = 1e-6)
  expect_equal(sum(loaded$deaths[loaded$area != "Portugal"]), sum(calendar$deaths[calendar$area == "Algarve"]), tolerance = 1e-6)
  expect_error(load_snapshot_deaths(year, "Portugal", "Tumores malignos"), "óbitos por causa")
}))

test_that("the planning indicators reach the latest year where the data allow", with_snapshots({
  year <- planning_calendar_years()
  skip_if(length(year) == 0, "no rebuilt year in the snapshots")
  year <- max(year)
  t <- planning_indicator_table(c("Portugal", "Norte", "Lisboa"), year,
                                ids = c("life_expectancy_65", "dsr_all", "dsr_preventable", "death_rate"), lookup = lookup)
  value <- function(area, id) t$value[t$area == area & t$indicator == id]
  expect_true(is.finite(value("Norte", "life_expectancy_65")))
  expect_true(is.finite(value("Portugal", "dsr_all")))
  # A municipality has no deaths by age; by cause nothing exists yet.
  expect_true(is.na(value("Lisboa", "dsr_all")))
  expect_true(is.na(value("Norte", "dsr_preventable")))
  # All-cause deaths by parish give every area its crude rate.
  expect_true(is.finite(value("Lisboa", "death_rate")))
  expect_equal(sum(planning_parish_death_totals(2024L)$deaths[planning_parish_death_totals(2024L)$area == "Portugal"]), 118396)
}))

test_that("views of small areas and by cause leave out the rebuilt latest year", with_snapshots({
  calendar <- planning_calendar_years()
  skip_if(length(calendar) == 0, "no rebuilt year in the snapshots")
  expect_true(all(calendar %in% planning_indicator_years("dsr_all")))
  expect_false(any(calendar %in% planning_small_area_years("dsr_all")))
  expect_false(any(calendar %in% planning_cause_years()))
  # The crude death rate comes from parish deaths and reaches every area.
  expect_true(all(calendar %in% planning_small_area_years("death_rate")))
  expect_equal(planning_map_year("dsr_all", 2030L, "nuts2", "2024"), max(calendar))
  expect_lt(planning_map_year("dsr_all", 2030L, "nuts2", "2013"), max(calendar))
  expect_lt(planning_map_year("dsr_all", 2030L, "municipio", "2024"), max(calendar))
}))

test_that("a live frame that differs from the snapshot is not filled", with_snapshots({
  deaths <- readRDS(planning_death_cause_file(2024L, "Todas as causas de morte"))
  deaths <- deaths[deaths$area %in% c("Barrancos", "Lisboa"), ]
  # As if INE had revised Barrancos since the snapshot.
  deaths$deaths[deaths$area == "Barrancos" & deaths$age_band == "85 e mais anos" & deaths$sex == "HM"] <-
    deaths$deaths[deaths$area == "Barrancos" & deaths$age_band == "85 e mais anos" & deaths$sex == "HM"] + 5
  filled <- suppressWarnings(fill_death_ages(deaths, check_recorded = TRUE))
  expect_equal(sum(filled$deaths_filled[filled$area == "Barrancos"]), 0)
  expect_warning(fill_death_ages(deaths, check_recorded = TRUE), "differ from the snapshot")
  # Without the check (the snapshot path) it is filled as before.
  expect_gt(sum(fill_death_ages(deaths)$deaths_filled[deaths$area == "Barrancos"]), 0)
}))

test_that("Portugal's rebuilt latest year keeps the deaths without a recorded age", with_snapshots({
  year <- planning_calendar_years()
  skip_if(length(year) == 0, "no rebuilt year in the snapshots")
  calendar <- planning_calendar_deaths(max(year))
  weekly <- readRDS(file.path(Sys.getenv("MORTALITY_SNAPSHOT_DIR"), "weekly_deaths", "0012100.rds"))
  known <- sum(weekly$deaths[weekly$code == "PT" & weekly$year == max(year) & !weekly$age_group %in% c("Total", "Ignorado")])
  # More than its known-age weekly deaths (ISO weeks straddle the year, so
  # the comparison is loose), and at least the five mainland regions.
  expect_gt(sum(calendar$deaths[calendar$area == "Portugal"]), sum(calendar$deaths[calendar$area == "Continente"]))
  expect_gt(abs(sum(calendar$deaths[calendar$area == "Portugal"]) - known), 0)
}))
