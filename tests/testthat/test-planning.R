# Tests for R/planning_indicators.R (the planning-indicator tab).

# A two-municipality country in a temporary snapshot root. Every dataset is
# written in the layout the fetchers produce, so the engine is exercised through
# its real readers.
planning_fixture <- function(years = 2020:2022) {
  root <- tempfile("planning-snap-")
  lookup <- tibble::tibble(
    municipality = c("Alfa", "Beta"),
    municipality_code = c("1110001", "1110002"),
    nuts1 = "Continente", nuts2 = "Norte", nuts3 = "Sub"
  )

  bands <- age_levels
  lower <- as.integer(sub("^(\\d+).*$", "\\1", bands))

  write_year <- function(dataset, year, frame, sub = NULL) {
    dir <- file.path(root, dataset, sub %||% "")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(frame, file.path(dir, paste0("year_", year, ".rds")))
  }

  for (year in years) {
    # Alfa: 100 per band per sex; Beta: 50. Young bands 0-14 = 3 bands.
    pop <- tidyr::expand_grid(area = c("Alfa", "Beta"), sex = c("H", "M"), age_band = bands) %>%
      dplyr::mutate(pop = ifelse(area == "Alfa", 100, 50))
    pop <- dplyr::bind_rows(
      pop,
      pop %>% dplyr::group_by(area, age_band) %>% dplyr::summarise(pop = sum(pop), .groups = "drop") %>% dplyr::mutate(sex = "HM"),
      pop %>% dplyr::group_by(sex, age_band) %>% dplyr::summarise(pop = sum(pop), .groups = "drop") %>% dplyr::mutate(area = "Portugal"),
      pop %>% dplyr::group_by(age_band) %>% dplyr::summarise(pop = sum(pop), .groups = "drop") %>% dplyr::mutate(area = "Portugal", sex = "HM")
    ) %>% dplyr::mutate(year = year)
    write_year("population", year, pop)

    write_year("births", year, tibble::tibble(
      year = year, area = c("Alfa", "Beta", "Portugal"), births = c(300, 100, 401), source_indicator = "x"
    ))

    causes <- c("Todas as causas de morte", "Tumores (neoplasmas) malignos", "Doenças do aparelho circulatório")
    totals <- tidyr::expand_grid(area = c("Alfa", "Beta"), sex = "HM", cause = causes) %>%
      dplyr::mutate(deaths = dplyr::case_when(
        cause == "Todas as causas de morte" ~ ifelse(area == "Alfa", 40, 20),
        cause == "Tumores (neoplasmas) malignos" ~ 10,
        TRUE ~ 5
      ))
    # Portugal carries one death of unknown residence.
    totals <- dplyr::bind_rows(
      totals,
      tibble::tibble(area = "Portugal", sex = "HM", cause = causes, deaths = c(61, 20, 10))
    ) %>% dplyr::mutate(year = year, source_indicator = "0008206")
    write_year("death_totals", year, totals, sub = "0008206")

    write_year("infant_deaths", year, tibble::tibble(
      year = year, area = c("Alfa", "Beta", "Portugal"), sex = "HM",
      cause = "Todas as causas de morte", deaths = c(2, 1, 3), source_indicator = "x"
    ))
  }

  list(root = root, lookup = lookup)
}

with_planning_fixture <- function(code) {
  fx <- planning_fixture()
  old <- Sys.getenv("MORTALITY_SNAPSHOT_DIR", unset = NA)
  Sys.setenv(MORTALITY_SNAPSHOT_DIR = fx$root)
  planning_clear_cache()
  on.exit({
    if (is.na(old)) Sys.unsetenv("MORTALITY_SNAPSHOT_DIR") else Sys.setenv(MORTALITY_SNAPSHOT_DIR = old)
    planning_clear_cache()
    unlink(fx$root, recursive = TRUE)
  })
  code(fx$lookup)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("indicator years follow the files and triennia need all three years", {
  with_planning_fixture(function(lookup) {
    expect_equal(planning_indicator_years("pop_total"), 2020:2022)
    expect_equal(planning_indicator_years("deaths"), 2020:2022)
    expect_equal(planning_indicator_years("infant_rate"), 2022L)
    expect_equal(planning_proportional_years(), 2022L)
  })
})

test_that("areas are ratios of summed components, not averages of rates", {
  with_planning_fixture(function(lookup) {
    tab <- planning_indicator_table(c("Norte", "Alfa"), 2022L, lookup = lookup)
    get <- function(area, id) tab$value[tab$area == area & tab$indicator == id]

    # 18 bands x 2 sexes x (100 + 50)
    expect_equal(get("Norte", "pop_total"), 18 * 2 * 150)
    expect_equal(get("Alfa", "pop_total"), 18 * 2 * 100)

    # 3 young bands, 10 working-age bands (15-64), 5 old bands (65+), 2 of 75+... plus 85+
    expect_equal(get("Norte", "ageing_index"), 5 / 3 * 100)
    expect_equal(get("Norte", "youth_dependency"), 3 / 10 * 100)
    expect_equal(get("Norte", "old_dependency"), 5 / 10 * 100)
    expect_equal(get("Norte", "pct_75_plus"), 3 / 18 * 100)

    # Births and deaths: sums, and rates of sums.
    expect_equal(get("Norte", "births"), 400)
    expect_equal(get("Norte", "birth_rate"), 400 / 5400 * 1000)
    expect_equal(get("Norte", "deaths"), 60)
    expect_equal(get("Alfa", "death_rate"), 40 / 3600 * 1000)

    # Triennial infant rate pools three years of deaths over three of births.
    expect_equal(get("Norte", "infant_rate"), (3 * 3) / (3 * 400) * 1000)
    # 1,200 pooled births clears the 1,000 threshold; Alfa's 900 does not.
    expect_equal(tab$flag[tab$area == "Norte" & tab$indicator == "infant_rate"], "")
    expect_equal(tab$flag[tab$area == "Alfa" & tab$indicator == "infant_rate"], "*")

    # The interval is the exact Poisson one on the pooled count.
    row <- tab[tab$area == "Norte" & tab$indicator == "death_rate", ]
    ci <- stats::poisson.test(60)$conf.int / 5400 * 1000
    expect_equal(c(row$lower, row$upper), as.numeric(ci))

    # Population indices carry no interval.
    expect_true(is.na(tab$lower[tab$area == "Norte" & tab$indicator == "ageing_index"]))
  })
})

test_that("Portugal uses the published row, which includes unknown residence", {
  with_planning_fixture(function(lookup) {
    tab <- planning_indicator_table("Portugal", 2022L, ids = c("deaths", "births"), lookup = lookup)
    expect_equal(tab$value[tab$indicator == "deaths"], 61)
    expect_equal(tab$value[tab$indicator == "births"], 401)
  })
})

test_that("a year outside the files is missing, not zero", {
  with_planning_fixture(function(lookup) {
    tab <- planning_indicator_table("Norte", 2023L, ids = c("deaths", "pop_total"), lookup = lookup)
    expect_true(all(is.na(tab$value)))
    # A triennium needs all three years.
    tab <- planning_indicator_table("Norte", 2021L, ids = "infant_rate", lookup = lookup)
    expect_true(is.na(tab$value))
  })
})

test_that("proportional mortality sums to 100% with the remainder shown", {
  with_planning_fixture(function(lookup) {
    prop <- planning_proportional("Norte", 2022L, lookup = lookup)
    expect_equal(prop$deaths[prop$code == "C00"], 180)
    expect_equal(prop$deaths[prop$code == "C07"], 60)
    expect_equal(prop$deaths[prop$code == "C33"], 30)
    expect_equal(sum(prop$share[prop$code != "C00"]), 100)
    expect_equal(prop$deaths[prop$code == "Outras"], 90)
    expect_equal(prop$period[[1]], "2020-2022")
    ci <- stats::binom.test(60, 180)$conf.int * 100
    expect_equal(c(prop$lower[prop$code == "C07"], prop$upper[prop$code == "C07"]), as.numeric(ci))
  })
})

test_that("the pyramid splits by sex and sums to 100%", {
  with_planning_fixture(function(lookup) {
    pyr <- planning_pyramid("Norte", 2022L, lookup = lookup)
    expect_equal(nrow(pyr), 36)
    expect_equal(sum(pyr$share), 100)
    expect_setequal(unique(pyr$sex), c("H", "M"))
  })
})

test_that("the profile table has one row per indicator and one column per area", {
  with_planning_fixture(function(lookup) {
    tab <- planning_indicator_table(c("Alfa", "Beta"), 2022L, lookup = lookup)
    wide <- planning_profile_wide(tab)
    expect_equal(nrow(wide), nrow(PLANNING_INDICATORS))
    expect_true(all(c("Alfa", "Beta") %in% names(wide)))
    expect_equal(wide$Período[grepl("infantil", wide$Indicador)], "2020-2022")
  })
})

test_that("the cause groups are chapter-level rubrics INE publishes", {
  expect_equal(nrow(PLANNING_CAUSE_GROUPS), 13)
  expect_equal(anyDuplicated(PLANNING_CAUSE_GROUPS$cause), 0L)
  # Real archive labels (the shortlist spells "infeciosas" without the c).
  expect_true("Algumas doenças infeciosas e parasitárias" %in% PLANNING_CAUSE_GROUPS$cause)
  expect_true("Tumores (neoplasmas) malignos" %in% PLANNING_CAUSE_GROUPS$cause)
})

test_that("health units for the ranking are the ULS and groups, not ARS", {
  units <- planning_uls_units()
  expect_equal(length(units), 36)
  expect_false(any(grepl("^ARS ", units)))
  # Together they cover every municipality of the mainland exactly once; the
  # islands have their own regional health services and no ULS.
  members <- unlist(lapply(units, health_unit_members))
  expect_equal(length(members), 278)
  expect_equal(anyDuplicated(members), 0L)
})

test_that("complete under-1 counts replace the band-derived ones for all causes", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "infant_totals"), showWarnings = FALSE)
    # The band-derived file says 2 + 1; INE's own count says 4 + 1.
    saveRDS(
      tibble::tibble(year = 2022L, area = c("Alfa", "Beta", "Portugal"), sex = "HM",
                     deaths = c(4, 1, 6), source_indicator = "0012541"),
      file.path(root, "infant_totals", "year_2022.rds")
    )
    planning_clear_cache()

    expect_equal(infant_deaths_total(2022L, c("Alfa", "Beta")), 5)
    # The AVPP correction keeps the count consistent with the age bands it splits.
    expect_equal(infant_deaths_total(2022L, c("Alfa", "Beta"), complete = FALSE), 3)
    # A year without a complete count falls back.
    expect_equal(infant_deaths_total(2021L, c("Alfa", "Beta")), 3)
    # A specific cause has no complete count to use.
    expect_equal(nrow(get_infant_death_data(2022L, "Alfa", cause = "Pneumonia")), 0)

    tab <- planning_indicator_table("Norte", 2022L, ids = "infant_rate", lookup = lookup)
    expect_equal(tab$numerator, 3 + 3 + 5)
  })
})

test_that("years whose municipal under-1 counts fall short of Portugal are flagged", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    # 2020: municipalities 2 + 1 = 3 of a national 10 -> incomplete.
    saveRDS(
      tibble::tibble(year = 2020L, area = c("Alfa", "Beta", "Portugal"), sex = "HM",
                     cause = "Todas as causas de morte", deaths = c(2, 1, 10), source_indicator = "x"),
      file.path(root, "infant_deaths", "year_2020.rds")
    )
    planning_clear_cache()

    expect_equal(infant_undercount_years(2020:2022, municipalities = lookup$municipality), 2020L)
    expect_match(infant_undercount_message(2020:2022, "Norte", municipalities = lookup$municipality), "2020")
    # Portugal reads its own row.
    expect_null(infant_undercount_message(2020:2022, "Portugal", municipalities = lookup$municipality))

    tab <- planning_indicator_table(c("Portugal", "Norte"), 2022L, ids = "infant_rate", lookup = lookup)
    expect_equal(tab$flag[tab$area == "Norte"], "†")
    expect_equal(tab$flag[tab$area == "Portugal"], "")

    # A complete count for the year clears it.
    dir.create(file.path(root, "infant_totals"), showWarnings = FALSE)
    saveRDS(
      tibble::tibble(year = 2020L, area = c("Alfa", "Beta", "Portugal"), sex = "HM",
                     deaths = c(6, 4, 10), source_indicator = "0008181"),
      file.path(root, "infant_totals", "year_2020.rds")
    )
    expect_length(infant_undercount_years(2020:2022, municipalities = lookup$municipality), 0)
  })
})

test_that("socio-economic, birth and neonatal indicators aggregate exactly", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    put <- function(measure, frame) {
      dir.create(file.path(root, "planning_extra", measure), recursive = TRUE, showWarnings = FALSE)
      saveRDS(dplyr::mutate(frame, year = 2022L, source_indicator = "x"),
              file.path(root, "planning_extra", measure, "year_2022.rds"))
    }
    total <- function(values) tibble::tibble(area = c("Alfa", "Beta"), category = "Total", value = values)

    put("rsi_beneficiaries", total(c(30, 12)))
    put("pensioners", total(c(100, 300)))
    put("pension_mean", total(c(9000, 5000)))
    # Alfa holds 3% of national purchasing power at index 150, Beta 1% at 50.
    put("purchasing_power_share", total(c(3, 1)))
    put("purchasing_power_per_capita", total(c(150, 50)))
    put("waste_collected", dplyr::bind_rows(
      total(c(2000, 700)),
      tibble::tibble(area = c("Alfa", "Beta"), category = "Recolha selectiva", value = c(500, 100))
    ))
    # Mother's age: single years and a 15-49 group sit beside the five-year
    # groups and must not be counted.
    put("births_by_mother_age", tibble::tibble(
      area = "Alfa",
      # INE publishes "50 - 54", "50 e mais" and "55 e mais" together: only
      # "50 e mais" may be read.
      category = c("Total", "10 - 14 anos", "15 - 19 anos", "17 anos", "15 - 49 anos", "20 - 24 anos",
                   "25 - 29 anos", "30 - 34 anos", "35 - 39 anos", "40 - 44 anos", "45 - 49 anos", "50 e mais anos",
                   "50 - 54 anos", "55 e mais anos"),
      value = c(100, 1, 9, 4, 99, 20, 30, 20, 12, 6, 1, 1, 1, 1)
    ))
    put("births_by_gestation", tibble::tibble(
      area = "Alfa",
      category = c("Total", "Menos de 22 semanas", "22 - 27 semanas", "28 - 31 semanas", "32 - 36 semanas", "37 - 41 semanas", "Ignorada"),
      value = c(100, 0, 1, 2, 5, 82, 10)
    ))
    put("infant_deaths_by_age", tibble::tibble(
      area = "Alfa",
      category = c("Total", "Menos de 28 dias", "Menos de 7 dias", "28 - 364 dias", "1 - 6 dias"),
      value = c(2, 1, 1, 1, 1)
    ))
    planning_clear_cache()

    ids <- c("rsi_rate", "pensioners_rate", "pension_mean", "purchasing_power", "waste_per_capita",
             "waste_selective_per_capita", "fertility_index", "rsi_beneficiaries", "pensioners")
    tab <- planning_indicator_table(c("Norte", "Alfa"), 2022L, ids = ids, lookup = lookup)
    get <- function(area, id) tab$value[tab$area == area & tab$indicator == id]

    # 15+ = 15 of 18 bands, both sexes: Alfa 3,000, Norte 4,500.
    expect_equal(get("Norte", "rsi_rate"), 42 / 4500 * 1000)
    expect_equal(get("Norte", "pensioners"), 400)
    expect_equal(get("Norte", "pensioners_rate"), 400 / 4500 * 1000)
    # Weighted by pensioners, not the mean of means (7,000).
    expect_equal(get("Norte", "pension_mean"), (100 * 9000 + 300 * 5000) / 400)
    # Sum of shares over implied population shares: 4 / (2 + 2) = 100, not the
    # mean of the indices (100 by coincidence would hide a bug, so check Alfa too).
    expect_equal(get("Norte", "purchasing_power"), 4 / (3 / 150 + 1 / 50))
    expect_equal(get("Alfa", "purchasing_power"), 150)
    expect_equal(get("Norte", "waste_per_capita"), 2700 * 1000 / 5400)
    expect_equal(get("Norte", "waste_selective_per_capita"), 600 * 1000 / 5400)

    # Fertility: each group's births over 100 women (Alfa 100 per band per sex,
    # the same at both ends of the year, so the mid-year mean is 100 too),
    # under-15 folded into 15-19 and 50+ into 45-49; single years ignored.
    expect_equal(get("Alfa", "fertility_index"), (10 + 20 + 30 + 20 + 12 + 6 + 2) / 100 * 5)
    ages <- planning_indicator_table("Alfa", 2022L, ids = "fertility_index", lookup = lookup)
    comps <- planning_components_year("Alfa", 2022L, lookup)
    expect_equal(comps$births_mother_ge35, 12 + 6 + 1 + 1)
    expect_equal(comps$births_mother_lt20, 1 + 9)
  })
})

test_that("triennial birth and neonatal proportions need all three years", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    for (year in 2020:2022) {
      put <- function(measure, frame) {
        dir.create(file.path(root, "planning_extra", measure), recursive = TRUE, showWarnings = FALSE)
        saveRDS(dplyr::mutate(frame, year = year, source_indicator = "x"),
                file.path(root, "planning_extra", measure, paste0("year_", year, ".rds")))
      }
      put("births_by_mother_age", tibble::tibble(
        area = c("Alfa", "Alfa", "Alfa", "Alfa"),
        category = c("Total", "15 - 19 anos", "35 - 39 anos", "50 e mais anos"),
        value = c(100, 5, 10, 1)
      ))
      put("births_by_gestation", tibble::tibble(
        area = "Alfa", category = c("Total", "32 - 36 semanas", "Ignorada"), value = c(100, 8, 20)
      ))
      put("infant_deaths_by_age", tibble::tibble(
        area = "Alfa", category = c("Menos de 28 dias", "Menos de 7 dias", "28 - 364 dias"), value = c(1, 1, 1)
      ))
    }
    planning_clear_cache()

    tab <- planning_indicator_table("Alfa", 2022L,
                                    ids = c("teen_births_pct", "older_births_pct", "preterm_pct", "neonatal_rate", "postneonatal_rate"),
                                    lookup = lookup)
    get <- function(id) tab[tab$indicator == id, ]
    expect_equal(get("teen_births_pct")$value, 15 / 300 * 100)
    expect_equal(get("older_births_pct")$value, 33 / 300 * 100)
    # Preterm over births of known duration.
    expect_equal(get("preterm_pct")$value, 24 / 240 * 100)
    expect_equal(c(get("preterm_pct")$lower, get("preterm_pct")$upper), as.numeric(stats::binom.test(24, 240)$conf.int * 100))
    # Per 1,000 live births from the births file (Alfa: 300 a year).
    expect_equal(get("neonatal_rate")$value, 3 / 900 * 1000)
    expect_equal(get("postneonatal_rate")$value, 3 / 900 * 1000)
    expect_equal(planning_indicator_years("preterm_pct"), 2022L)
  })
})

test_that("every indicator has a computation and a year rule", {
  for (id in PLANNING_INDICATORS$id) {
    expect_type(planning_indicator_years(id), "integer")
  }
  choices <- planning_indicator_choices()
  expect_setequal(unlist(choices, use.names = FALSE), PLANNING_INDICATORS$id)
})

test_that("the fertility index counts women at mid-year", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    # Double the 2021 population: the 2022 mid-year denominator is 150 women
    # per band for Alfa, not the 100 at the end of 2022.
    pop <- readRDS(file.path(root, "population", "year_2021.rds"))
    pop$pop <- pop$pop * 2
    saveRDS(pop, file.path(root, "population", "year_2021.rds"))
    dir.create(file.path(root, "planning_extra", "births_by_mother_age"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      tibble::tibble(year = 2022L, area = "Alfa", category = paste0(seq(15, 45, 5), " - ", seq(19, 49, 5), " anos"),
                     value = 15, source_indicator = "x"),
      file.path(root, "planning_extra", "births_by_mother_age", "year_2022.rds")
    )
    planning_clear_cache()
    tab <- planning_indicator_table("Alfa", 2022L, ids = "fertility_index", lookup = lookup)
    expect_equal(tab$value, 7 * 15 / 150 * 5)
  })
})
