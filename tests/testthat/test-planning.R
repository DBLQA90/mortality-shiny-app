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

test_that("the cause groups are chapter-level rubrics INE publishes", {
  expect_equal(nrow(PLANNING_CAUSE_GROUPS), 13)
  expect_equal(anyDuplicated(PLANNING_CAUSE_GROUPS$cause), 0L)
  # Real archive labels (the shortlist spells "infeciosas" without the c).
  expect_true("Algumas doenças infeciosas e parasitárias" %in% PLANNING_CAUSE_GROUPS$cause)
  expect_true("Tumores (neoplasmas) malignos" %in% PLANNING_CAUSE_GROUPS$cause)
})

test_that("health units for the ranking are the ULS and groups, not ARS", {
  units <- planning_uls_units()
  # 33 ULS of whole municipalities plus the two exact groups; the six ULS that
  # share a municipality are listed individually only in the "units" view.
  expect_equal(length(units), 35)
  expect_equal(length(planning_uls_units("units")), 39)
  expect_false(any(grepl("\\+", planning_uls_units("units"))))
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

test_that("closed-form intervals equal poisson.test and binom.test", {
  for (n in c(0, 1, 7, 60, 1234)) {
    ci <- planning_poisson_ci(n)
    expect_equal(c(ci$lower, ci$upper), as.numeric(stats::poisson.test(n)$conf.int), tolerance = 1e-8)
  }
  for (pair in list(c(0, 10), c(3, 10), c(10, 10), c(60, 180))) {
    ci <- planning_binomial_ci(pair[[1]], pair[[2]])
    expect_equal(c(ci$lower, ci$upper), as.numeric(stats::binom.test(pair[[1]], pair[[2]])$conf.int * 100), tolerance = 1e-8)
  }
  expect_true(all(is.na(unlist(planning_binomial_ci(3, 0)))))
})

test_that("comparators are the containing areas, one per level, nearest first", {
  lookup <- tibble::tibble(
    municipality = c("Alfa", "Beta", "Gama"),
    municipality_code = c("1110001", "1110002", "1120001"),
    nuts1 = "Continente", nuts2 = "Norte", nuts3 = c("Sub", "Sub", "Outra")
  )
  local_health <- tibble::tibble(unit = character(0), kind = character(0), municipality = character(0), ars = character(0))

  comps <- planning_comparators("Alfa", lookup, local_health)
  expect_equal(comps$level, c("NUTS III", "NUTS II", "Portugal"))
  expect_equal(comps$area, c("Sub", "Norte", "Portugal"))
  # Continente has the same municipalities as Norte here: not repeated.
  expect_false("Continente" %in% comps$area)

  expect_equal(planning_comparators("Sub", lookup, local_health)$area, c("Norte", "Portugal"))
  expect_equal(nrow(planning_comparators("Portugal", lookup, local_health)), 0)

  levels <- planning_area_levels(lookup, local_health)
  expect_equal(levels$level[levels$area == "Norte"], "NUTS II")
  expect_equal(levels$level[levels$area == "Alfa"], "Município")
})

test_that("every indicator is marked comparable or not, and counts are not", {
  expect_true(is.logical(PLANNING_INDICATORS$comparable))
  counts <- PLANNING_INDICATORS$id[PLANNING_INDICATORS$unit == "N.º"]
  expect_true(all(!PLANNING_INDICATORS$comparable[PLANNING_INDICATORS$id %in% counts]))
  expect_true(all(PLANNING_INDICATORS$comparable[!PLANNING_INDICATORS$id %in% counts]))
  expect_setequal(names(PLANNING_SHEET_NAMES), PLANNING_INDICATORS$id)
  expect_true(all(nchar(PLANNING_SHEET_NAMES) <= 31))
  for (forbidden in c("[", "]", ":", "*", "?", "/", "\\")) {
    expect_false(any(grepl(forbidden, PLANNING_SHEET_NAMES, fixed = TRUE)))
  }
})

test_that("the many-area proportional table matches the single-area one", {
  with_planning_fixture(function(lookup) {
    single <- planning_proportional("Norte", 2022L, lookup = lookup)
    many <- planning_proportional_table(c("Alfa", "Norte", "Portugal"), 2022L, lookup = lookup)
    norte <- many[many$area == "Norte", ]
    expect_equal(norte$deaths, single$deaths)
    expect_equal(norte$lower, single$lower)
    # Portugal reads its published row (61 all-cause deaths a year).
    expect_equal(many$deaths[many$area == "Portugal" & many$code == "C00"], 183)
  })
})

test_that("both Excel workbooks are written with a sheet per indicator", {
  skip_if_not_installed("openxlsx")
  with_planning_fixture(function(lookup) {
    path <- tempfile(fileext = ".xlsx")
    areas <- tibble::tibble(area = c("Alfa", "Norte", "Portugal"), level = c("Local", "NUTS II", "Portugal"))
    write_planning_workbook(path, areas, lookup = lookup, data_date = "2026-09-17", focus = "Alfa")
    sheets <- openxlsx::getSheetNames(path)
    expect_true(all(c("Leia-me", "Resumo", "Dados", "I4 Envelhecimento", "I37 Óbitos") %in% sheets))
    readme <- openxlsx::read.xlsx(path, "Leia-me", colNames = FALSE)[[1]]
    expect_true(any(grepl("2026-09-17", readme)))

    ageing <- openxlsx::read.xlsx(path, "I4 Envelhecimento", startRow = 4)
    expect_equal(ageing$Local, c("Alfa", "Norte", "Portugal"))
    expect_equal(ageing[["2022"]][ageing$Local == "Norte"], 5 / 3 * 100)

    summary <- openxlsx::read.xlsx(path, "Resumo")
    deaths <- summary[grepl("^Óbitos", summary$Indicador), ]
    # A count has no comparators.
    expect_true(is.na(deaths[["Norte.(NUTS.II)"]]))
    expect_false(is.na(deaths[["Alfa.(Local)"]]))

    full <- tempfile(fileext = ".xlsx")
    write_planning_workbook(full, areas, lookup = lookup, include_long = FALSE)
    expect_false("Dados" %in% openxlsx::getSheetNames(full))
    expect_false("Resumo" %in% openxlsx::getSheetNames(full))
  })
})

test_that("charts build for rates, counts, rankings, pyramids and causes", {
  with_planning_fixture(function(lookup) {
    areas <- tibble::tibble(area = c("Alfa", "Norte", "Portugal"), level = c("Local", "NUTS II", "Portugal"))
    rate <- planning_indicator_table(areas$area, 2020:2022, ids = "death_rate", lookup = lookup)
    expect_s3_class(planning_trend_chart(rate, planning_indicator_spec("death_rate"), areas), "plotly")
    count <- planning_indicator_table("Alfa", 2020:2022, ids = "deaths", lookup = lookup)
    expect_s3_class(planning_count_chart(count, planning_indicator_spec("deaths"), "Alfa"), "plotly")
    ranking <- planning_indicator_table(c("Portugal", "Alfa", "Beta"), 2022L, ids = "death_rate", lookup = lookup)
    expect_s3_class(planning_ranking_chart(ranking, planning_indicator_spec("death_rate"), 2022L, "Alfa"), "plotly")
    expect_s3_class(planning_pyramid_chart(planning_pyramid("Alfa", 2022L, lookup), planning_pyramid("Portugal", 2022L, lookup), "Alfa", "Portugal"), "plotly")
    expect_s3_class(planning_proportional_chart(planning_proportional_table(areas$area, 2022L, lookup = lookup), areas), "plotly")

    shown <- planning_series_display(rate, planning_indicator_spec("death_rate"), areas)
    expect_equal(names(shown), c("Período", "Alfa", "Norte (NUTS II)", "Portugal"))
    expect_equal(shown$Período, c("2022", "2021", "2020"))

    summary <- planning_summary_display(planning_indicator_table(areas$area, 2022L, lookup = lookup), areas)
    expect_equal(summary$`Norte (NUTS II)`[summary$Indicador == "Óbitos [I37]"], "")
  })
})

test_that("the life table reproduces PHEindicatormethods on its own age structure", {
  skip_if_not_installed("PHEindicatormethods")
  ages <- c(0L, 1L, seq(5L, 90L, 5L))
  population <- c(5000, 21000, rep(26000, 8), rep(30000, 6), 22000, 15000, 9000, 6000)
  rates <- c(0.003, 0.0002, 0.0001, 0.0001, 0.0003, 0.0005, 0.0006, 0.0008, 0.001, 0.0015,
             0.0025, 0.004, 0.006, 0.009, 0.014, 0.022, 0.035, 0.06, 0.1, 0.2)
  deaths <- round(population * rates)
  phe <- PHEindicatormethods::phe_life_expectancy(
    tibble::tibble(age = ages, pop = population, deaths = deaths),
    deaths, pop, age, age_contents = ages, le_age = 0
  )
  ours <- abridged_life_expectancy(deaths, population, c(1, 4, rep(5, 17), NA), c(0.1, rep(0.5, 19)))
  expect_equal(ours$value, phe$value, tolerance = 1e-10)
  expect_equal(c(ours$lower, ours$upper), c(phe$lowercl, phe$uppercl), tolerance = 1e-10)
})

test_that("life tables are refused where PHE refuses them", {
  n <- length(age_levels)
  widths <- c(rep(5, n - 1), NA)
  a <- rep(0.5, n)
  expect_match(abridged_life_expectancy(rep(1, n), rep(200, n), widths, a)$reason, "5.000")
  expect_match(abridged_life_expectancy(c(0, rep(1, n - 1)), c(0, rep(1000, n - 1)), widths, a)$reason, "nula")
  expect_true(is.na(abridged_life_expectancy(rep(NA, n), rep(1000, n), widths, a)$value))
})

test_that("life expectancy for areas pools three years and spreads unrecorded ages", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    rates <- c(0.001, 0.0001, 0.0001, 0.0002, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0015,
               0.0025, 0.004, 0.006, 0.01, 0.016, 0.028, 0.05, 0.14)
    for (year in 2020:2022) {
      # Ten times the fixture's population, so every area and sex clears the
      # 5,000 person-year minimum.
      pop <- readRDS(file.path(root, "population", paste0("year_", year, ".rds")))
      pop$pop <- pop$pop * 10
      saveRDS(pop, file.path(root, "population", paste0("year_", year, ".rds")))
      pop <- pop[pop$area %in% c("Alfa", "Beta", "Portugal"), ]
      deaths <- pop %>%
        dplyr::mutate(cause = "Todas as causas de morte", deaths = pop * 2 * rates[match(age_band, age_levels)]) %>%
        dplyr::select(year, area, sex, cause, age_band, deaths)
      # Beta publishes half its 2021 deaths without an age (a sixth of the
      # triennium's, under the quarter that withholds the value).
      recorded <- deaths
      recorded$deaths[recorded$area == "Beta" & recorded$year == 2021] <-
        recorded$deaths[recorded$area == "Beta" & recorded$year == 2021] / 2
      dir <- file.path(root, "deaths", "0008206", paste0("year_", year))
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(recorded, file.path(dir, "cause_todas_as_causas_de_morte.rds"))

      totals <- deaths %>%
        dplyr::group_by(year, area, sex, cause) %>%
        dplyr::summarise(deaths = sum(deaths), .groups = "drop") %>%
        dplyr::mutate(source_indicator = "0008206")
      saveRDS(totals, file.path(root, "death_totals", "0008206", paste0("year_", year, ".rds")))
    }
    planning_clear_cache()

    expect_equal(life_expectancy_years(), 2022L)
    table <- planning_indicator_table(c("Alfa", "Beta", "Norte"), 2022L, ids = life_expectancy_ids, lookup = lookup)
    expect_equal(nrow(table), 3 * nrow(LIFE_INDICATORS))
    expect_true(all(is.finite(table$value)))
    # Same age-specific rates everywhere, so the same life expectancy.
    hm <- table[table$indicator == "life_expectancy", ]
    expect_equal(hm$value[hm$area == "Alfa"], hm$value[hm$area == "Norte"], tolerance = 1e-6)
    # Beta's unrecorded 2021 ages were spread back over the ages its own
    # population and the national gap imply: its value matches too, and is flagged.
    expect_equal(hm$value[hm$area == "Beta"], hm$value[hm$area == "Alfa"], tolerance = 1e-6)
    expect_equal(hm$flag[hm$area == "Beta"], "‡")
    expect_equal(hm$flag[hm$area == "Alfa"], "")
    # Life expectancy at 65 is shorter than at birth, and shares the flag.
    at65 <- table[table$indicator == "life_expectancy_65", ]
    expect_true(all(at65$value < hm$value))
    expect_equal(at65$flag[at65$area == "Beta"], "\u2021")

    # The smaller area has the wider interval.
    expect_gt(hm$upper[hm$area == "Beta"] - hm$lower[hm$area == "Beta"], hm$upper[hm$area == "Alfa"] - hm$lower[hm$area == "Alfa"])
  })
})

test_that("birth weight, late fetal and perinatal rates follow INE's definitions", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    put <- function(measure, year, frame) {
      dir.create(file.path(root, "planning_extra", measure), recursive = TRUE, showWarnings = FALSE)
      saveRDS(dplyr::mutate(frame, year = year, source_indicator = "x"),
              file.path(root, "planning_extra", measure, paste0("year_", year, ".rds")))
    }
    for (year in 2020:2022) {
      # 300 births: 12 under 2,500 g, 10 of unknown weight.
      put("births_by_weight", year, tibble::tibble(
        area = "Alfa",
        category = c("Total", "Menos de 500 g", "500 - 999 g", "2 000 - 2 499 g", "3 000 - 3 499 g", "5 000 g e mais", "Ignorada"),
        value = c(300, 1, 3, 8, 275, 3, 10)
      ))
      put("perinatal_deaths", year, tibble::tibble(area = "Alfa", category = "Total", value = 5))
      put("infant_deaths_by_age", year, tibble::tibble(
        area = "Alfa", category = c("Menos de 28 dias", "Menos de 7 dias", "28 - 364 dias"), value = c(3, 2, 1)
      ))
    }
    planning_clear_cache()

    tab <- planning_indicator_table("Alfa", 2022L,
                                    ids = c("low_birth_weight_pct", "late_fetal_rate", "perinatal_rate"),
                                    lookup = lookup)
    get <- function(id) tab[tab$indicator == id, ]

    # Under 2,500 g over births of known weight, pooled over three years.
    expect_equal(get("low_birth_weight_pct")$value, 36 / 870 * 100)
    expect_equal(
      c(get("low_birth_weight_pct")$lower, get("low_birth_weight_pct")$upper),
      as.numeric(stats::binom.test(36, 870)$conf.int * 100)
    )

    # Stillbirths are perinatal deaths less deaths under 7 days: 15 - 6 = 9.
    # Both rates divide by live births plus stillbirths (Alfa: 900 + 9).
    expect_equal(get("late_fetal_rate")$value, 9 / 909 * 1000)
    expect_equal(get("late_fetal_rate")$numerator, 9)
    expect_equal(get("perinatal_rate")$value, 15 / 909 * 1000)
    expect_equal(get("perinatal_rate")$denominator, 909)
    expect_equal(planning_indicator_years("perinatal_rate"), 2022L)
  })
})

test_that("birth-weight bands parse, including the open top band", {
  expect_equal(planning_weight_lower(c("Menos de 500 g", "500 - 999 g", "2 000 - 2 499 g", "5 000 g e mais", "Total", "Ignorada")),
               c(0, 500, 2000, 5000, NA, NA))
})

test_that("census indicators use the census years and their own denominators", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    put <- function(measure, year, frame) {
      dir.create(file.path(root, "planning_extra", measure), recursive = TRUE, showWarnings = FALSE)
      saveRDS(dplyr::mutate(frame, year = year, source_indicator = "x"),
              file.path(root, "planning_extra", measure, paste0("year_", year, ".rds")))
    }
    both <- function(area, value_2011, value_2021) list("2011" = value_2011, "2021" = value_2021)
    for (year in c(2011L, 2021L)) {
      scale <- if (year == 2011L) 1 else 1.1
      put("census_population", year, tibble::tibble(area = c("Alfa", "Beta"), category = "Total", value = c(1000, 500) * scale))
      put("census_population_by_age", year, tibble::tibble(
        area = rep(c("Alfa", "Beta"), each = 3),
        category = rep(c("0 - 4 anos", "5 - 9 anos", "10 - 14 anos"), 2),
        value = c(100, 100, 800, 50, 50, 400) * scale
      ))
      put("census_illiteracy_rate", year, tibble::tibble(area = c("Alfa", "Beta"), category = "Total", value = c(5, 10)))
      put("census_education", year, tibble::tibble(
        area = rep(c("Alfa", "Beta"), each = 4),
        category = rep(c("Total", "Primário/Básico", "Secundário", "Superior"), 2),
        value = c(800, 500, 200, 100, 400, 250, 100, 50) * scale
      ))
    }
    planning_clear_cache()

    ids <- c("census_population", "census_population_change", "pct_education_none", "pct_education_higher", "illiteracy_rate")
    tab <- planning_indicator_table("Norte", c(2011L, 2021L), ids = ids, lookup = lookup)
    get <- function(id, year) tab$value[tab$indicator == id & tab$year == year]

    expect_equal(get("census_population", 2011), 1500)
    # Ten years apart, not one.
    expect_true(is.na(get("census_population_change", 2011)))
    expect_equal(get("census_population_change", 2021), 10)
    # Those with no level are the census population less those with one.
    expect_equal(get("pct_education_none", 2011), (1500 - 1200) / 1500 * 100)
    expect_equal(get("pct_education_higher", 2011), 150 / 1500 * 100)
    # The rate is weighted by each municipality's population aged 10 and over.
    expect_equal(get("illiteracy_rate", 2011), (0.05 * 800 + 0.10 * 400) / 1200 * 100)
    expect_equal(planning_indicator_years("illiteracy_rate"), c(2011L, 2021L))
  })
})

test_that("the life table gives every age, and e65 matches a table read at 65", {
  skip_if_not_installed("PHEindicatormethods")
  ages <- c(0L, 1L, seq(5L, 90L, 5L))
  population <- c(5000, 21000, rep(26000, 8), rep(30000, 6), 22000, 15000, 9000, 6000)
  rates <- c(0.003, 0.0002, 0.0001, 0.0001, 0.0003, 0.0005, 0.0006, 0.0008, 0.001, 0.0015,
             0.0025, 0.004, 0.006, 0.009, 0.014, 0.022, 0.035, 0.06, 0.1, 0.2)
  deaths <- round(population * rates)
  phe <- PHEindicatormethods::phe_life_expectancy(
    tibble::tibble(age = ages, pop = population, deaths = deaths),
    deaths, pop, age, age_contents = ages, le_age = 65
  )
  table <- abridged_life_table(deaths, population, c(1, 4, rep(5, 17), NA), c(0.1, rep(0.5, 19)))
  index <- match(65L, ages)
  expect_equal(table$e[[index]], phe$value, tolerance = 1e-10)
  z <- stats::qnorm(0.975)
  expect_equal(table$e[[index]] - z * table$se[[index]], phe$lowercl, tolerance = 1e-10)
  expect_equal(nchar(table$reason), 0L)
})

test_that("indicator sheets carry the interval bounds below the values", {
  skip_if_not_installed("openxlsx")
  with_planning_fixture(function(lookup) {
    path <- tempfile(fileext = ".xlsx")
    areas <- tibble::tibble(area = c("Alfa", "Norte"), level = c("Local", "NUTS II"))
    write_planning_workbook(path, areas, lookup = lookup, include_long = FALSE)

    sheet <- openxlsx::read.xlsx(path, "I38 Mortalidade", startRow = 4, colNames = TRUE)
    labels <- sheet[[1]]
    expect_true(any(grepl("Limite inferior", labels)))
    expect_true(any(grepl("Limite superior", labels)))

    values <- openxlsx::read.xlsx(path, "I38 Mortalidade", startRow = 4, rows = 4:6, colNames = TRUE)
    rate <- planning_indicator_table("Alfa", 2022L, ids = "death_rate", lookup = lookup)
    expect_equal(values[["2022"]][values$Local == "Alfa"], round(rate$value, 6), tolerance = 1e-5)

    # A count has no interval, so no bound blocks.
    counts <- openxlsx::read.xlsx(path, "I1 População", startRow = 4, colNames = TRUE)
    expect_false(any(grepl("Limite", counts[[1]])))
  })
})

test_that("proportional mortality under 75 prefers INE's rows and flags 2014", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    causes <- c("Todas as causas de morte", "Tumores (neoplasmas) malignos", "Doenças do aparelho circulatório")
    bands <- age_levels
    young <- bands[1:15]

    for (year in 2012:2022) {
      dir <- file.path(root, "deaths", "0008206", paste0("year_", year))
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      for (cause in causes) {
        per_band <- switch(cause, "Todas as causas de morte" = 10, "Tumores (neoplasmas) malignos" = 4, 1)
        frame <- tidyr::expand_grid(area = c("Alfa", "Beta", "Portugal"), sex = "HM", age_band = bands) %>%
          dplyr::mutate(year = year, cause = cause,
                        deaths = per_band * ifelse(area == "Portugal", 3, ifelse(area == "Alfa", 2, 1)))
        saveRDS(frame, file.path(dir, paste0("cause_", planning_cause_file_token(cause), ".rds")))
      }
      # The death totals decide which triennia exist.
      totals <- readRDS(file.path(root, "death_totals", "0008206", "year_2022.rds"))
      saveRDS(dplyr::mutate(totals, year = year), file.path(root, "death_totals", "0008206", paste0("year_", year, ".rds")))
    }
    planning_clear_cache()

    table <- planning_under75_table(c("Alfa", "Norte"), c(2022L, 2014L), lookup = lookup, vintage = "2024")
    alfa <- table[table$area == "Alfa" & table$end_year == 2022, ]
    # 15 bands under 75, three years, two per band for Alfa.
    expect_equal(alfa$deaths[alfa$code == "C00"], 15 * 10 * 2 * 3)
    expect_equal(alfa$share[alfa$code == "C07"], 4 / 10 * 100)
    # Municipal sums: the triennium containing 2014 is flagged, 2020-2022 is not.
    expect_equal(unique(table$flag[table$end_year == 2014]), "§")
    expect_equal(unique(table$flag[table$end_year == 2022]), "")

    # With an INE regional row for Norte (code 11), that row wins and the flag goes.
    dir.create(file.path(root, "regional_deaths", "0008206"), recursive = TRUE, showWarnings = FALSE)
    for (year in 2012:2022) {
      rows <- tidyr::expand_grid(cause = causes, age_band = bands) %>%
        dplyr::mutate(year = year, region_code = "11", area = "Norte", sex = "HM",
                      deaths = ifelse(cause == "Todas as causas de morte", 40, ifelse(cause == "Tumores (neoplasmas) malignos", 20, 4)),
                      source_indicator = "0008206")
      saveRDS(rows, file.path(root, "regional_deaths", "0008206", paste0("year_", year, ".rds")))
    }
    planning_clear_cache()

    # The 2019-2021 triennium: every year of it reads the same INE edition as
    # the fixture's rows (a 2022 window would expect the newer one).
    with_rows <- planning_under75_table("Norte", c(2021L, 2014L), lookup = lookup, vintage = "2024")
    norte <- with_rows[with_rows$end_year == 2021, ]
    expect_equal(norte$deaths[norte$code == "C00"], 15 * 40 * 3)
    expect_equal(norte$share[norte$code == "C07"], 50)
    expect_equal(unique(with_rows$flag), "")
    expect_match(planning_under75_note(table), "2014")
    expect_null(planning_under75_note(with_rows))
  })
})

test_that("earnings are weighted by employees, and sector shares sum to the total", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    put <- function(measure, frame) {
      dir.create(file.path(root, "planning_extra", measure), recursive = TRUE, showWarnings = FALSE)
      saveRDS(dplyr::mutate(frame, year = 2022L, source_indicator = "x"),
              file.path(root, "planning_extra", measure, "year_2022.rds"))
    }
    put("employees_by_sector", tibble::tibble(
      area = rep(c("Alfa", "Beta"), each = 4),
      category = rep(c("Total", "Agricultura, produção animal, caça, floresta e pesca",
                       "Indústria, construção, energia e água", "Serviços"), 2),
      value = c(1000, 100, 300, 600, 500, 50, 250, 200)
    ))
    put("earnings_mean", tibble::tibble(area = c("Alfa", "Beta"), category = "Total", value = c(1200, 900)))
    planning_clear_cache()

    ids <- c("earnings_mean", "employees", "pct_employees_primary", "pct_employees_tertiary")
    tab <- planning_indicator_table(c("Norte", "Alfa"), 2022L, ids = ids, lookup = lookup)
    get <- function(area, id) tab$value[tab$area == area & tab$indicator == id]

    # Weighted by employees, not the mean of the two means (1,050).
    expect_equal(get("Norte", "earnings_mean"), (1200 * 1000 + 900 * 500) / 1500)
    expect_equal(get("Alfa", "earnings_mean"), 1200)
    expect_equal(get("Norte", "employees"), 1500)
    expect_equal(get("Norte", "pct_employees_primary"), 150 / 1500 * 100)
    expect_equal(get("Norte", "pct_employees_tertiary"), 800 / 1500 * 100)

    # The three sectors account for the total.
    shares <- planning_indicator_table("Norte", 2022L,
                                       ids = c("pct_employees_primary", "pct_employees_secondary", "pct_employees_tertiary"),
                                       lookup = lookup)
    expect_equal(sum(shares$value), 100)

    # A count has no comparators; a mean and a share do.
    expect_false(planning_indicator_spec("employees")$comparable)
    expect_true(planning_indicator_spec("earnings_mean")$comparable)
  })
})

test_that("crude rates, RSI and waste use the mean population; pensioners the year-end one", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    # Alfa doubles in 2022: 3,600 at the end of 2021, 7,200 at the end of 2022.
    path <- file.path(root, "population", "year_2022.rds")
    pop <- readRDS(path)
    pop$pop[pop$area == "Alfa"] <- pop$pop[pop$area == "Alfa"] * 2
    saveRDS(pop, path)
    dir.create(file.path(root, "planning_extra", "pensioners"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(root, "planning_extra", "rsi_beneficiaries"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(90, 10), source_indicator = "x"),
            file.path(root, "planning_extra", "pensioners", "year_2022.rds"))
    dir.create(file.path(root, "planning_extra", "pension_mean"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(6000, 6000), source_indicator = "x"),
            file.path(root, "planning_extra", "pension_mean", "year_2022.rds"))
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(45, 5), source_indicator = "x"),
            file.path(root, "planning_extra", "rsi_beneficiaries", "year_2022.rds"))
    planning_clear_cache()

    tab <- planning_indicator_table("Alfa", 2020:2022, ids = c("birth_rate", "death_rate", "rsi_rate", "pensioners_rate"), lookup = lookup)
    get <- function(id, year) tab$value[tab$indicator == id & tab$year == year]
    expect_equal(get("birth_rate", 2022L), 300 / 5400 * 1000)
    expect_equal(get("death_rate", 2022L), 40 / 5400 * 1000)
    # 15+ is 15 of 18 bands: 3,000 before, 6,000 after.
    expect_equal(get("rsi_rate", 2022L), 45 / 4500 * 1000)
    expect_equal(get("pensioners_rate", 2022L), 90 / 6000 * 1000)
    # The first year of the series has no previous estimate.
    expect_equal(get("birth_rate", 2020L), 300 / 3600 * 1000)
  })
})

test_that("a blank all-ages death total is rebuilt from its age bands", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    path <- file.path(root, "death_totals", "0008206", "year_2022.rds")
    totals <- readRDS(path)
    # INE left Alfa's all-causes cell blank; older fetches stored it as 0.
    totals$deaths[totals$area == "Alfa" & totals$cause == "Todas as causas de morte"] <- 0
    saveRDS(totals, path)
    dir.create(file.path(root, "deaths", "0008206", "year_2022"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tibble::tibble(year = 2022L, area = "Alfa", sex = "HM", cause = "Todas as causas de morte",
                           age_band = age_levels[1:4], deaths = c(1, 2, 3, 30)),
            file.path(root, "deaths", "0008206", "year_2022", "cause_todas_as_causas_de_morte.rds"))
    planning_clear_cache()

    repaired <- read_death_totals_year(2022L)
    row <- repaired$area == "Alfa" & repaired$cause == "Todas as causas de morte"
    expect_equal(repaired$deaths[row], 36)
    expect_true(repaired$repaired[row])
    # A total above its (incomplete) bands is left alone.
    expect_false(any(repaired$repaired[repaired$area == "Beta"]))
    tab <- planning_indicator_table("Norte", 2022L, ids = "deaths", lookup = lookup)
    expect_equal(tab$value, 56)
  })
})

test_that("a municipality missing from a complete block makes its areas missing", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "planning_extra", "waste_collected"), recursive = TRUE, showWarnings = FALSE)
    # Beta was not reported: INE's blank, stored as 0.
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(2000, 0), source_indicator = "x"),
            file.path(root, "planning_extra", "waste_collected", "year_2022.rds"))
    planning_clear_cache()
    tab <- planning_indicator_table(c("Alfa", "Beta", "Norte"), 2022L, ids = "waste_per_capita", lookup = lookup)
    expect_equal(tab$value[tab$area == "Alfa"], 2000 * 1000 / 3600)
    expect_true(is.na(tab$value[tab$area == "Beta"]))
    expect_true(is.na(tab$value[tab$area == "Norte"]))
  })
})

test_that("a municipality filed under another has values only together with it", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    env <- environment(planning_joint_split)
    old <- get("PLANNING_JOINT_REPORTING", envir = env)
    on.exit(assign("PLANNING_JOINT_REPORTING", old, envir = env), add = TRUE)
    # Beta was split from Alfa: joint for everything to 2020, and for waste
    # whenever Beta is blank.
    assign("PLANNING_JOINT_REPORTING", tibble::tribble(
      ~municipality, ~holder, ~until, ~blocks,
      "Beta",        "Alfa",  2020L,  list("waste")
    ), envir = env)
    dir.create(file.path(root, "planning_extra", "waste_collected"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(2700, NA), source_indicator = "x"),
            file.path(root, "planning_extra", "waste_collected", "year_2022.rds"))
    planning_clear_cache()

    tab <- planning_indicator_table(c("Alfa", "Beta", "Norte"), 2020:2022, ids = c("birth_rate", "waste_per_capita"), lookup = lookup)
    get <- function(area, id, year) tab$value[tab$area == area & tab$indicator == id & tab$year == year]
    expect_true(is.na(get("Alfa", "birth_rate", 2020L)))
    expect_true(is.na(get("Beta", "birth_rate", 2020L)))
    expect_equal(get("Norte", "birth_rate", 2020L), 400 / 5400 * 1000)
    expect_equal(get("Beta", "birth_rate", 2021L), 100 / 1800 * 1000)
    expect_true(is.na(get("Alfa", "waste_per_capita", 2022L)))
    expect_true(is.na(get("Beta", "waste_per_capita", 2022L)))
    expect_equal(get("Norte", "waste_per_capita", 2022L), 2700 * 1000 / 5400)

    # Proportional mortality: the triennium holding 2020 is withheld.
    prop <- planning_proportional_table(c("Alfa", "Norte"), c(2022L), lookup = lookup)
    expect_true(all(is.na(prop$share[prop$area == "Alfa"])))
    expect_false(any(is.na(prop$share[prop$area == "Norte"])))
    expect_equal(nrow(planning_proportional("Alfa", 2022L, lookup = lookup)), 0)
  })
})

test_that("suppressed employment sectors are estimated from the rest of the NUTS III and flagged", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "planning_extra", "employees_by_sector"), recursive = TRUE, showWarnings = FALSE)
    sectors <- c("Total", "Agricultura, produção animal, caça, floresta e pesca", "Indústria, construção, energia e água", "Serviços")
    # Alfa: 10 / 40 / 50. Beta publishes only services; INE hid the other two.
    saveRDS(tibble::tibble(year = 2022L, area = rep(c("Alfa", "Beta"), each = 4), category = rep(sectors, 2),
                           value = c(100, 10, 40, 50, 100, 0, 0, 30), source_indicator = "x"),
            file.path(root, "planning_extra", "employees_by_sector", "year_2022.rds"))
    planning_clear_cache()
    ids <- c("pct_employees_primary", "pct_employees_secondary", "pct_employees_tertiary")
    tab <- planning_indicator_table(c("Beta", "Norte"), 2022L, ids = ids, lookup = lookup)
    get <- function(area, id) tab$value[tab$area == area & tab$indicator == id]
    # The hidden 70 split 10:40, as in Alfa.
    expect_equal(get("Beta", "pct_employees_primary"), 14)
    expect_equal(get("Beta", "pct_employees_secondary"), 56)
    expect_equal(get("Beta", "pct_employees_tertiary"), 30)
    expect_equal(get("Norte", "pct_employees_secondary"), 48)
    expect_true(all(tab$flag == "≈"))
  })
})

test_that("Portugal as the sum of its municipalities leaves out unknown residence", {
  with_planning_fixture(function(lookup) {
    tab <- planning_indicator_table(c("Portugal", PLANNING_PORTUGAL_MUNICIPAL), 2022L, ids = c("deaths", "births"), lookup = lookup)
    get <- function(area, id) tab$value[tab$area == area & tab$indicator == id]
    # The fixture's Portugal row carries one death and one birth of unknown residence.
    expect_equal(get("Portugal", "deaths"), 61)
    expect_equal(get(PLANNING_PORTUGAL_MUNICIPAL, "deaths"), 60)
    expect_equal(get(PLANNING_PORTUGAL_MUNICIPAL, "births"), 400)
    expect_equal(planning_portugal_area("municipal"), PLANNING_PORTUGAL_MUNICIPAL)
    expect_equal(planning_portugal_area("published"), "Portugal")
  })
})

test_that("significance is the benchmark against the whole interval, for comparable indicators only", {
  expect_equal(planning_significance(c(1, 3, 0.5, NA), c(2, 4, 5, 1), c(2.5, 2.5, 2.5, 2.5)),
               c("Inferior", "Superior", "Semelhante", NA))
  table <- tibble::tibble(
    area = c("Portugal", "A", "Portugal", "A"), indicator = c("death_rate", "death_rate", "deaths", "deaths"),
    year = 2022L, value = c(10, 12, 1000, 50), lower = c(9.9, 11, 990, 40), upper = c(10.1, 13, 1010, 60)
  )
  out <- planning_add_significance(table, "Portugal")
  expect_equal(out$significance[out$area == "A" & out$indicator == "death_rate"], "Superior")
  # A count is never compared, however far apart.
  expect_true(is.na(out$significance[out$area == "A" & out$indicator == "deaths"]))
  expect_true(all(is.na(out$significance[out$area == "Portugal"])))
  expect_equal(planning_significance_mark(c("Superior", NA)), c(" ▲", ""))
})

test_that("funnel limits narrow with size and bracket the benchmark", {
  n <- c(100, 1000, 100000)
  p <- planning_funnel_limits(n, 3, "poisson", 1000, 0.95)
  expect_true(all(p$lower < 3 & p$upper > 3))
  expect_true(all(diff(p$upper - p$lower) < 0))
  # Close to the normal approximation for a large expected count (300).
  expect_equal(p$upper[[3]], 3 + 1.96 * sqrt(300) / 100, tolerance = 0.01)
  # A zero in a small unit is never significantly low.
  expect_lte(planning_funnel_limits(500, 3, "poisson", 1000, 0.998)$lower, 0)
  b <- planning_funnel_limits(n, 10, "binomial", 100, 0.998)
  expect_true(all(b$lower >= 0 & b$upper <= 100 & b$lower < 10 & b$upper > 10))
  units <- tibble::tibble(area = c("Portugal", "Big", "Small"), value = c(3, 4.5, 4.5), lower = NA, upper = NA,
                          flag = "", denominator = c(1e6, 1e5, 500), indicator = "infant_rate")
  data <- planning_funnel_data(units, planning_indicator_spec("infant_rate"), "Portugal")
  expect_equal(data$position[data$area == "Big"], "Acima do limite de 99,8%")
  expect_equal(data$position[data$area == "Small"], "Dentro dos limites")
  expect_null(planning_funnel_data(units, planning_indicator_spec("ageing_index"), "Portugal"))
})

test_that("education can be restricted to an age, from the census by age group", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "planning_extra", "census_education_by_age"), recursive = TRUE, showWarnings = FALSE)
    # Alfa, two groups: under 15 (100, all without a level) and 25-29 (50: 5 none, 45 higher).
    saveRDS(tibble::tibble(
      year = 2021L, area = "Alfa", age = c(0L, 0L, 25L, 25L, 25L),
      category = c("Total", "Nenhum", "Total", "Nenhum", "Superior"), value = c(100, 100, 50, 5, 45), source_indicator = "x"
    ), file.path(root, "planning_extra", "census_education_by_age", "year_2021.rds"))
    planning_clear_cache()
    tab <- planning_indicator_table("Alfa", 2021L, ids = c("pct_education_none", "pct_education_higher"), lookup = lookup, education_min_age = 15L)
    expect_equal(tab$value[tab$indicator == "pct_education_none"], 10)
    expect_equal(tab$value[tab$indicator == "pct_education_higher"], 90)
    all_ages <- planning_indicator_table("Alfa", 2021L, ids = "pct_education_none", lookup = lookup, education_min_age = 25L)
    expect_equal(all_ages$value, 10)
    expect_match(planning_indicator_label("pct_education_none", 25L), "25 e mais anos")
  })
})

test_that("an area holding both of a joint pair keeps its value when nothing else is split", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    env <- environment(planning_joint_split)
    old <- get("PLANNING_JOINT_REPORTING", envir = env)
    on.exit(assign("PLANNING_JOINT_REPORTING", old, envir = env), add = TRUE)
    assign("PLANNING_JOINT_REPORTING", tibble::tribble(~municipality, ~holder, ~until, ~blocks, "Beta", "Alfa", 2000L, list("waste")), envir = env)
    dir.create(file.path(root, "planning_extra", "waste_collected"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(tibble::tibble(year = 2022L, area = c("Alfa", "Beta"), category = "Total", value = c(2700, NA), source_indicator = "x"),
            file.path(root, "planning_extra", "waste_collected", "year_2022.rds"))
    planning_clear_cache()
    # Only Norte is asked for: no split area in the request.
    tab <- planning_indicator_table("Norte", 2022L, ids = "waste_per_capita", lookup = lookup)
    expect_equal(tab$value, 2700 * 1000 / 5400)
  })
})

test_that("the location profile is a Word document", {
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  with_planning_fixture(function(lookup) {
    path <- tempfile(fileext = ".docx")
    on.exit(unlink(path), add = TRUE)
    areas <- tibble::tibble(area = c("Alfa", "Norte"), level = c("Local", "NUTS II"))
    write_planning_profile(path, areas, lookup = lookup, vintage = "2024", data_date = "2026-09-21")
    expect_true(file.exists(path) && file.size(path) > 10000)
    text <- officer::docx_summary(officer::read_docx(path))$text
    expect_true(any(grepl("Perfil do local: Alfa", text, fixed = TRUE)))
    expect_true(any(grepl("Esperança de vida", text, fixed = TRUE)))
  })
})

# Deaths by age for all causes: Alfa complete (10 at 60-64, 30 at 80-84 = its
# total of 40); Beta records only 15 of its 20, all at 80-84; Portugal's row
# holds the 5 missing ones at 30-34 plus one death of unknown residence.
with_cause_age_fixture <- function(code) {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    for (year in 2020:2022) {
      dir <- file.path(root, "deaths", "0008206", paste0("year_", year))
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      rows <- tibble::tribble(
        ~area,      ~age_band,       ~deaths,
        "Alfa",     "60 - 64 anos",  10,
        "Alfa",     "80 - 84 anos",  30,
        "Beta",     "80 - 84 anos",  15,
        "Portugal", "30 - 34 anos",  5,
        "Portugal", "60 - 64 anos",  10,
        "Portugal", "80 - 84 anos",  46
      )
      saveRDS(dplyr::mutate(rows, year = year, sex = "HM", cause = "Todas as causas de morte"),
              file.path(dir, "cause_todas_as_causas_de_morte.rds"))
    }
    planning_clear_cache()
    code(lookup)
  })
}

test_that("missing deaths by age follow the national gap, so the municipal sum matches Portugal by age", {
  with_cause_age_fixture(function(lookup) {
    pool <- planning_standardised_pool(c("Beta", "Portugal", PLANNING_PORTUGAL_MUNICIPAL), 2022L, "HM", lookup)
    beta <- pool$deaths["Beta", , "all"]
    # Gap: 5 at 30-34 and 1 at 80-84 (Portugal 46 against 30 + 15), so Beta's
    # 5 missing go 5/6 and 1/6.
    expect_equal(unname(beta[["30 - 34 anos"]]), 3 * 5 * 5 / 6)
    expect_equal(unname(beta[["80 - 84 anos"]]), 3 * (15 + 5 / 6))
    expect_equal(sum(beta), 3 * 20)
    tab <- planning_indicator_table(c("Alfa", "Beta"), 2022L, ids = c("premature_deaths", "dsr_all"), lookup = lookup)
    expect_equal(tab$value[tab$area == "Beta" & tab$indicator == "premature_deaths"], 3 * 5 * 5 / 6)
    # A quarter of Beta's deaths were spread: flagged.
    expect_equal(tab$flag[tab$area == "Beta" & tab$indicator == "dsr_all"], "‡")
    expect_equal(tab$flag[tab$area == "Alfa" & tab$indicator == "dsr_all"], "")
  })
})

test_that("standardised rates match PHEindicatormethods and the SMR is 100 for the benchmark", {
  skip_if_not_installed("PHEindicatormethods")
  with_cause_age_fixture(function(lookup) {
    tab <- planning_indicator_table(c("Alfa", "Portugal"), 2022L, ids = c("dsr_all", "smr_all", "dsr_premature"), lookup = lookup)
    get <- function(area, id, column = "value") tab[[column]][tab$area == area & tab$indicator == id]
    deaths <- stats::setNames(rep(0, length(age_levels)), age_levels)
    deaths[c("60 - 64 anos", "80 - 84 anos")] <- c(30, 90)
    phe <- PHEindicatormethods::calculate_dsr(
      tibble::tibble(x = deaths, n = rep(600, length(age_levels)), stdpop = esp2013_df$stdpop),
      x = x, n = n, stdpop = stdpop
    )
    expect_equal(get("Alfa", "dsr_all"), phe$value, tolerance = 1e-9)
    # Dobson's interval; the Poisson limits are exact here and Byar's
    # approximation in PHEindicatormethods, which agree to 0.002% at 120 deaths.
    expect_equal(get("Alfa", "dsr_all", "lower"), phe$lowercl, tolerance = 1e-4)
    expect_equal(get("Alfa", "dsr_all", "upper"), phe$uppercl, tolerance = 1e-4)
    expect_equal(get("Portugal", "smr_all"), 100)
    # Alfa's expected deaths at Portugal's rates (Portugal: 300 per band a year).
    expected <- 600 * (15 + 30 + 138) / 900
    expect_equal(get("Alfa", "smr_all"), 120 / expected * 100)
    expect_equal(get("Alfa", "smr_all", "denominator"), expected)
  })
})

test_that("years of potential life lost weight each death to 70", {
  with_cause_age_fixture(function(lookup) {
    tab <- planning_indicator_table("Alfa", 2022L, ids = "ypll_rate", lookup = lookup)
    # Alfa: 10 deaths a year at 60-64 (midpoint 62.5): 7.5 years each, over
    # 14 bands under 70 of 200 people, three years.
    expect_equal(tab$value, 3 * 10 * 7.5 / (3 * 14 * 200) * 1e5)
  })
})

test_that("SNS units map onto the app's ULS, one by one", {
  expect_equal(planning_sns_unit("Área dos CSP da ULS Gaia / Espinho"), "ULS Vila Nova de Gaia/Espinho")
  expect_equal(planning_sns_unit("CSP da ULS Póvoa Varzim / Vila Conde"), "ULS Póvoa de Varzim/Vila do Conde")
  # The portal reports every ULS separately, including those that share a
  # municipality: its units follow the parishes and do not overlap.
  expect_equal(planning_sns_unit("CSP da ULS São José"), "ULS São José")
  expect_equal(planning_sns_unit("CSP da ULS Lisboa Ocidental"), "ULS Lisboa Ocidental")
  expect_equal(planning_sns_unit("Área dos CSP da ULS Guarda"), "ULS Guarda")
})

test_that("SNS proportions are rebuilt as numerators and denominators, and compared at the end of their cycle", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "sns"), recursive = TRUE, showWarnings = FALSE)
    units <- c("Área dos CSP da ULS Santa Maria", "Área dos CSP da ULS São José", "Área dos CSP da ULS Lisboa Ocidental")
    rows <- tidyr::expand_grid(period = c("2024-11", "2024-12", "2025-01"), unit = units)
    # Screened women: 50, 30 and 20 of 100, 100 and 50 eligible.
    counts <- c(50, 30, 20); eligible <- c(100, 100, 50)
    rows$count <- rep(counts, 3); rows$prop <- rep(counts / eligible * 100, 3)
    sns <- dplyr::bind_rows(
      dplyr::transmute(rows, period, unit, region = "LVT", field = "contagem_de_mulheres_com_registo_de_mamografia_nos_ultimos_dois_anos", value = count),
      dplyr::transmute(rows, period, unit, region = "LVT", field = "proporcao_mulheres_50_70_a_c_mamogr_2_anos", value = prop)
    ) %>% dplyr::mutate(dataset = "rastreios-oncologicos")
    saveRDS(sns, file.path(root, "sns", "rastreios-oncologicos.rds"))
    planning_clear_cache()
    group <- c("ULS Santa Maria", "ULS São José", "ULS Lisboa Ocidental")
    # The other screenings' fields are absent here: skipped, with a warning.
    expect_warning(t <- planning_sns_table(list(grupo = group), "sns_mammography"), "skipped")
    expect_equal(t$denominator, rep(250, 3))
    expect_equal(t$value, rep(100 / 250 * 100, 3))
    # Mammography accumulates over the calendar year: only December compares.
    expect_equal(t$complete, c(FALSE, TRUE, FALSE))
    expect_equal(t$provisional, c(FALSE, FALSE, TRUE))
  })
})

test_that("weekly expected deaths come from baseline rates, leaving out the pandemic years", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    dir.create(file.path(root, "weekly_deaths"), recursive = TRUE, showWarnings = FALSE)
    weekly <- tidyr::expand_grid(year = c(2021L, 2023L, 2024L, 2025L), week = 1:3, age_group = c("40-44 anos", "85-89 anos", "Total")) %>%
      dplyr::mutate(
        deaths = dplyr::case_when(
          age_group == "Total" ~ NA_real_,
          year == 2021L ~ 1000,                                   # pandemic: must not count
          year == 2023L ~ ifelse(age_group == "85-89 anos", 10, 2),
          year == 2024L ~ ifelse(age_group == "85-89 anos", 14, 2),
          TRUE ~ ifelse(age_group == "85-89 anos", 20, 2)
        ),
        code = "11", region = "Norte", source_indicator = "0012100"
      )
    weekly$deaths[weekly$age_group == "Total"] <- NA
    saveRDS(weekly, file.path(root, "weekly_deaths", "0012100.rds"))
    planning_clear_cache()
    expect_equal(planning_weekly_region("Alfa", lookup)$region, "Norte")
    x <- planning_weekly_excess("Norte", 2025L, "85plus", lookup)
    # Population is constant, so the expected count is the mean of 2023-2024.
    expect_equal(x$expected, rep(12, 3))
    expect_equal(unique(x$baseline), "2023, 2024")
    s <- planning_weekly_summary(x)
    expect_equal(s$excess, 3 * (20 - 12))
    # 2024 has only one baseline year: no expected value.
    expect_true(all(is.na(planning_weekly_excess("Norte", 2024L, "85plus", lookup)$expected)))
  })
})

test_that("the six ULS that share a municipality read whole or by parish weights", {
  skip_if_not(file.exists("../../data/uls_parish.rds"), "parish lookup not built")
  skip_if_not(file.exists("../../data/snapshots/census_parish/year_2021.rds"), "census parishes not fetched")
  # The engine reads snapshots relative to the app directory; point it at the
  # repository's own data for this test.
  old <- Sys.getenv("MORTALITY_SNAPSHOT_DIR", unset = NA)
  Sys.setenv(MORTALITY_SNAPSHOT_DIR = normalizePath("../../data/snapshots"))
  planning_clear_cache()
  on.exit({
    if (is.na(old)) Sys.unsetenv("MORTALITY_SNAPSHOT_DIR") else Sys.setenv(MORTALITY_SNAPSHOT_DIR = old)
    planning_clear_cache()
  }, add = TRUE)

  lookup <- get_nuts_lookup("2024")
  units <- planning_parish_units()
  expect_setequal(units, c("ULS Lisboa Ocidental", "ULS Loures/Odivelas", "ULS Santa Maria",
                           "ULS Santo António", "ULS São João", "ULS São José"))
  shares <- planning_parish_shares()
  # Every basis splits each shared municipality into shares that add to one.
  totals <- shares %>% dplyr::group_by(.data$municipality, .data$basis) %>% dplyr::summarise(total = sum(.data$weight), .groups = "drop")
  expect_true(all(abs(totals$total - 1) < 1e-9))
  expect_setequal(unique(shares$municipality), c("Lisboa", "Loures", "Porto"))

  year <- max(intersect(planning_indicator_years("pop_total"), planning_indicator_years("deaths")))
  all_uls <- planning_uls_units("units")
  whole <- planning_indicator_table(c(all_uls, "Continente"), year, ids = c("pop_total", "deaths"), lookup = lookup, mode = "whole")
  parish <- planning_indicator_table(c(all_uls, "Continente"), year, ids = c("pop_total", "deaths"), lookup = lookup, mode = "parish")
  sum_of <- function(table, id) sum(table$value[table$indicator == id & table$area != "Continente"], na.rm = TRUE)
  continente <- function(table, id) table$value[table$indicator == id & table$area == "Continente"]
  # Whole municipalities overlap; parish weights add up to the country.
  expect_gt(sum_of(whole, "pop_total"), continente(whole, "pop_total"))
  expect_equal(sum_of(parish, "pop_total"), continente(parish, "pop_total"))
  expect_equal(sum_of(parish, "deaths"), continente(parish, "deaths"))

  # A ULS that shares nothing is untouched by the choice.
  matosinhos <- lapply(c("whole", "parish"), function(m) {
    planning_indicator_table("ULS Matosinhos", year, ids = "pop_total", lookup = lookup, mode = m)$value
  })
  expect_equal(matosinhos[[1]], matosinhos[[2]])
  # The exact group does not depend on it either.
  group <- lapply(c("whole", "parish"), function(m) {
    planning_indicator_table("ULS Santo António + São João", year, ids = "pop_total", lookup = lookup, mode = m)$value
  })
  expect_equal(group[[1]], group[[2]])
  # Porto's two ULS split it, so each is smaller than under the whole reading.
  porto <- planning_indicator_table(c("ULS Santo António", "ULS São João"), year, ids = "pop_total", lookup = lookup, mode = "parish")
  expect_equal(sum(porto$value), group[[1]])

  # Births and deaths are not estimated where INE publishes them by parish:
  # each ULS gets its parishes' own count, plus its whole municipalities.
  registers <- planning_parish_actual("deaths")
  expect_true(year %in% registers$year)
  totals <- registers %>% dplyr::group_by(.data$year, .data$municipality) %>% dplyr::summarise(total = sum(.data$weight), .groups = "drop")
  expect_true(all(abs(totals$total - 1) < 1e-9))

  parishes <- readRDS("../../data/uls_parish.rds")
  deaths <- readRDS("../../data/snapshots/parish_vitals/deaths.rds")
  deaths <- deaths[deaths$year == year, , drop = FALSE]
  for (unit in c("ULS São José", "ULS Santo António")) {
    own <- sum(deaths$value[deaths$code %in% parishes$code[parishes$unit == unit]])
    rest <- setdiff(planning_area_members(unit, lookup), unique(parishes$municipality))
    whole <- if (length(rest) == 0) 0 else sum(planning_indicator_table(rest, year, ids = "deaths", lookup = lookup)$value)
    got <- planning_indicator_table(unit, year, ids = "deaths", lookup = lookup, mode = "parish")$value
    expect_equal(got, own + whole)
  }
})

test_that("life expectancy and standardised rates are withheld where most deaths have no age", {
  with_planning_fixture(function(lookup) {
    root <- Sys.getenv("MORTALITY_SNAPSHOT_DIR")
    rates <- c(0.001, 0.0001, 0.0001, 0.0002, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0015,
               0.0025, 0.004, 0.006, 0.01, 0.016, 0.028, 0.05, 0.14)
    for (year in 2020:2022) {
      pop <- readRDS(file.path(root, "population", paste0("year_", year, ".rds")))
      pop$pop <- pop$pop * 10
      saveRDS(pop, file.path(root, "population", paste0("year_", year, ".rds")))
      pop <- pop[pop$area %in% c("Alfa", "Beta", "Portugal"), ]
      deaths <- pop %>%
        dplyr::mutate(cause = "Todas as causas de morte", deaths = pop * 2 * rates[match(age_band, age_levels)]) %>%
        dplyr::select(year, area, sex, cause, age_band, deaths)
      # Beta publishes no ages at all in 2021 and 2022: two thirds of the
      # triennium, far beyond the quarter the rule allows.
      recorded <- deaths
      recorded$deaths[recorded$area == "Beta" & recorded$year %in% c(2021, 2022)] <- 0
      dir <- file.path(root, "deaths", "0008206", paste0("year_", year))
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(recorded, file.path(dir, "cause_todas_as_causas_de_morte.rds"))
      totals <- deaths %>%
        dplyr::group_by(year, area, sex, cause) %>%
        dplyr::summarise(deaths = sum(deaths), .groups = "drop") %>%
        dplyr::mutate(source_indicator = "0008206")
      saveRDS(totals, file.path(root, "death_totals", "0008206", paste0("year_", year, ".rds")))
    }
    planning_clear_cache()

    table <- planning_indicator_table(c("Alfa", "Beta"), 2022L, ids = c("life_expectancy", "dsr_all", "deaths"), lookup = lookup)
    got <- function(area, id) table$value[table$area == area & table$indicator == id]
    expect_true(is.finite(got("Alfa", "life_expectancy")))
    expect_true(is.na(got("Beta", "life_expectancy")))
    expect_true(is.na(got("Beta", "dsr_all")))
    # The count itself does not depend on the ages, so it stays.
    expect_true(is.finite(got("Beta", "deaths")))
    expect_equal(table$flag[table$area == "Beta" & table$indicator == "life_expectancy"], "")
  })
})
