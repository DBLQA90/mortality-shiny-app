# =========================================================
# Planning indicators (the "Indicadores de Planeamento" tab)
# =========================================================
# Demographic and mortality indicators that local health plans report for every
# ULS, ARS, NUTS region and municipality: the population structure, dependency
# and ageing indices, births, deaths, crude rates, infant mortality and
# proportional mortality by large cause groups. They mirror the indicators of
# the DRS/PNS2030 support workbook (I1, I3-I9, I13-I17, I28, I32, I33, I35,
# I37-I42, I45, I64, I65), recomputed from the
# app's own INE snapshots.
#
# Every area is built the same way: the components (population by broad age
# group, births, deaths, infant deaths) are summed over the area's
# municipalities, and each indicator is a ratio of those sums. Averaging
# municipal rates would weight a village the same as a city. Portugal and
# Continente are read from INE's published rows where the dataset carries them,
# because those include events whose municipality of residence is unknown.
#
# Components come from datasets that are complete at municipal level:
#
#   population     data/snapshots/population     age bands sum to the total
#   births         data/snapshots/births         live births, mother's residence
#   deaths         data/snapshots/death_totals   INE's all-ages "Total" row
#   infant deaths  data/snapshots/infant_totals  deaths under 1 year (infant_deaths
#                                                where no complete count exists)
#
# Deaths deliberately avoid the main archive: INE's municipal breakdown of
# deaths by age is incomplete (2014 worst), while the all-ages total is not.
# A count, a crude rate or an all-ages proportion needs no age, so it reads the
# total and is exact at every level.

PLANNING_INDICATORS <- tibble::tribble(
  ~id,                   ~theme,                ~label,                                                    ~ref,  ~unit,                      ~window, ~digits,
  "pop_total",           "Demografia",          "População residente (estimativa)",                        "I1",  "N.º",                      1L,      0L,
  "pct_0_14",            "Demografia",          "Proporção de jovens (0-14 anos)",                         "I1",  "%",                        1L,      1L,
  "pct_65_plus",         "Demografia",          "Proporção de idosos (65 e mais anos)",                    "I1",  "%",                        1L,      1L,
  "pct_75_plus",         "Demografia",          "Proporção de 75 e mais anos",                             "I1",  "%",                        1L,      1L,
  "ageing_index",        "Demografia",          "Índice de envelhecimento",                                "I4",  "por 100 jovens",           1L,      1L,
  "youth_dependency",    "Demografia",          "Índice de dependência de jovens",                         "I5",  "por 100 em idade activa",  1L,      1L,
  "old_dependency",      "Demografia",          "Índice de dependência de idosos",                         "I6",  "por 100 em idade activa",  1L,      1L,
  "births",              "Natalidade",          "Nados-vivos",                                             "I7",  "N.º",                      1L,      0L,
  "birth_rate",          "Natalidade",          "Taxa bruta de natalidade",                                "I8",  "‰",                        1L,      1L,
  "fertility_index",     "Natalidade",          "Índice sintético de fecundidade",                         "I9",  "filhos por mulher",        1L,      2L,
  "teen_births_pct",     "Natalidade",          "Nascimentos em mães com menos de 20 anos (triénio)",      "I32", "%",                        3L,      1L,
  "older_births_pct",    "Natalidade",          "Nascimentos em mães com 35 e mais anos (triénio)",        "I33", "%",                        3L,      1L,
  "preterm_pct",         "Natalidade",          "Nascimentos pré-termo (triénio)",                         "I35", "%",                        3L,      1L,
  "rsi_beneficiaries",   "Contexto social",     "Beneficiários do rendimento social de inserção",          "I13", "N.º",                      1L,      0L,
  "rsi_rate",            "Contexto social",     "Beneficiários do RSI por 1.000 habitantes com 15+ anos",  "I14", "‰",                        1L,      1L,
  "pensioners",          "Contexto social",     "Pensionistas da segurança social",                        "I15", "N.º",                      1L,      0L,
  "pensioners_rate",     "Contexto social",     "Pensionistas por 1.000 habitantes com 15+ anos",          "I16", "‰",                        1L,      1L,
  "pension_mean",        "Contexto social",     "Valor médio anual das pensões",                           "I17", "€",                        1L,      0L,
  "purchasing_power",    "Contexto social",     "Poder de compra per capita",                              "I28", "Portugal = 100",           1L,      1L,
  "waste_per_capita",    "Ambiente",            "Resíduos urbanos recolhidos por habitante",               "I64", "kg/hab.",                  1L,      0L,
  "waste_selective_per_capita", "Ambiente",     "Resíduos recolhidos selectivamente por habitante",        "I65", "kg/hab.",                  1L,      0L,
  "deaths",              "Mortalidade",         "Óbitos",                                                  "I37", "N.º",                      1L,      0L,
  "death_rate",          "Mortalidade",         "Taxa bruta de mortalidade",                               "I38", "‰",                        1L,      1L,
  "infant_rate",         "Mortalidade",         "Taxa de mortalidade infantil (triénio)",                  "I39", "‰ nados-vivos",            3L,      1L,
  "neonatal_rate",       "Mortalidade",         "Taxa de mortalidade neonatal (triénio)",                  "I40", "‰ nados-vivos",            3L,      1L,
  "early_neonatal_rate", "Mortalidade",         "Taxa de mortalidade neonatal precoce (triénio)",          "I41", "‰ nados-vivos",            3L,      1L,
  "postneonatal_rate",   "Mortalidade",         "Taxa de mortalidade pós-neonatal (triénio)",              "I42", "‰ nados-vivos",            3L,      1L
)

# Where a source changes definition, so the series is not continuous across the
# year named. Shown on the evolution chart and in the notes.
PLANNING_SERIES_BREAKS <- tibble::tribble(
  ~indicator,        ~year, ~note,
  "pensioners",      2017L, "Série 2017 da segurança social substitui a Série 1990-2023 (cerca de -5,5% de pensionistas).",
  "pensioners_rate", 2017L, "Série 2017 da segurança social substitui a Série 1990-2023 (cerca de -5,5% de pensionistas).",
  "pension_mean",    2017L, "Série 2017 da segurança social substitui a Série 1990-2023."
)

# Grouped by theme, which selectInput renders as option groups.
planning_indicator_choices <- function() {
  themes <- unique(PLANNING_INDICATORS$theme)
  stats::setNames(lapply(themes, function(theme) {
    rows <- PLANNING_INDICATORS[PLANNING_INDICATORS$theme == theme, , drop = FALSE]
    stats::setNames(as.list(rows$id), paste0(rows$label, " (", rows$ref, ")"))
  }), themes)
}

planning_indicator_spec <- function(id) {
  spec <- PLANNING_INDICATORS[PLANNING_INDICATORS$id == id, , drop = FALSE]
  if (nrow(spec) != 1) stop("Unknown planning indicator: ", id, call. = FALSE)
  as.list(spec)
}

# The large cause groups of the workbook's I45/I46, in INE's shortlist wording.
# Each is a chapter-level rubric, so none contains another and the shares can be
# shown side by side; what the thirteen leave out (mental disorders, skin,
# pregnancy, congenital malformations) is reported as the remainder, so the
# column always sums to 100%.
PLANNING_CAUSE_GROUPS <- tibble::tribble(
  ~code, ~cause,                                                                                                 ~label,
  "C01", "Algumas doenças infeciosas e parasitárias",                                                           "Doenças infecciosas e parasitárias",
  "C07", "Tumores (neoplasmas) malignos",                                                                       "Tumores malignos",
  "C25", "Doenças do sangue e dos órgãos hematopoéticos e alguns transtornos imunitários",                      "Doenças do sangue e órgãos hematopoéticos",
  "C26", "Doenças endócrinas, nutricionais e metabólicas",                                                      "Doenças endócrinas, nutricionais e metabólicas",
  "C31", "Doenças do sistema nervoso e dos órgãos dos sentidos",                                                "Doenças do sistema nervoso e órgãos dos sentidos",
  "C33", "Doenças do aparelho circulatório",                                                                    "Doenças do aparelho circulatório",
  "C37", "Doenças do aparelho respiratório",                                                                    "Doenças do aparelho respiratório",
  "C42", "Doenças do aparelho digestivo",                                                                       "Doenças do aparelho digestivo",
  "C46", "Doenças do sistema osteomuscular/ tecido conjuntivo",                                                 "Doenças do sistema osteomuscular e tecido conjuntivo",
  "C48", "Doenças do aparelho geniturinário",                                                                   "Doenças do aparelho geniturinário",
  "C51", "Algumas afecções originadas no período perinatal",                                                    "Afecções originadas no período perinatal",
  "C55", "Sintomas, sinais e achados anormais de exames clínicos e de laboratório não classificados em outra parte", "Sintomas, sinais e achados anormais (mal definidas)",
  "C58", "Causas externas de lesão e envenenamento",                                                            "Causas externas"
)

planning_all_causes <- "Todas as causas de morte"
planning_published_areas <- c("Portugal", "Continente")

# ---------------------------------------------------------
# Reading the snapshots
# ---------------------------------------------------------
# Year files are read once per session. The key carries the snapshot root, so a
# test pointing MORTALITY_SNAPSHOT_DIR elsewhere never sees another root's rows.
planning_cache <- new.env(parent = emptyenv())

planning_clear_cache <- function() {
  rm(list = ls(planning_cache, all.names = TRUE), envir = planning_cache)
  invisible(NULL)
}

planning_read <- function(dataset, year) {
  key <- paste(infant_snapshot_root(), dataset, year, sep = "|")
  if (!exists(key, envir = planning_cache, inherits = FALSE)) {
    value <- if (identical(dataset, "death_totals")) {
      read_death_totals_year(year)
    } else {
      read_year_file(dataset, year)
    }
    assign(key, value, envir = planning_cache)
  }
  get(key, envir = planning_cache, inherits = FALSE)
}

# Death totals are stored per source indicator. 0013166 (NUTS-2024) wins from
# 2022, the same precedence as the main death archive.
death_totals_dir <- function() file.path(infant_snapshot_root(), "death_totals")

death_totals_years <- function() {
  dirs <- list.dirs(death_totals_dir(), recursive = FALSE)
  files <- unlist(lapply(dirs, list.files, pattern = "^year_\\d+\\.rds$"), use.names = FALSE)
  sort(unique(as.integer(sub("^year_(\\d+)\\.rds$", "\\1", files))))
}

read_death_totals_year <- function(year) {
  for (indicator in c("0013166", "0008206")) {
    path <- file.path(death_totals_dir(), indicator, paste0("year_", year, ".rds"))
    if (file.exists(path)) return(readRDS(path))
  }
  NULL
}

planning_dataset_years <- function(dataset) {
  if (startsWith(dataset, "extra:")) return(planning_extra_years(sub("^extra:", "", dataset)))
  switch(
    dataset,
    death_totals = death_totals_years(),
    infant = sort(union(snapshot_years_for("infant_totals"), snapshot_years_for("infant_deaths"))),
    snapshot_years_for(dataset)
  )
}

# Years each indicator can be computed for, from the files actually present.
planning_indicator_years <- function(id) {
  needs <- switch(
    id,
    births = "births",
    birth_rate = c("births", "population"),
    deaths = "death_totals",
    death_rate = c("death_totals", "population"),
    infant_rate = c("births", "infant"),
    fertility_index = c("extra:births_by_mother_age", "population"),
    teen_births_pct = , older_births_pct = "extra:births_by_mother_age",
    preterm_pct = "extra:births_by_gestation",
    rsi_beneficiaries = "extra:rsi_beneficiaries",
    rsi_rate = c("extra:rsi_beneficiaries", "population"),
    pensioners = "extra:pensioners",
    pensioners_rate = c("extra:pensioners", "population"),
    pension_mean = c("extra:pensioners", "extra:pension_mean"),
    purchasing_power = c("extra:purchasing_power_share", "extra:purchasing_power_per_capita"),
    waste_per_capita = , waste_selective_per_capita = c("extra:waste_collected", "population"),
    neonatal_rate = , early_neonatal_rate = , postneonatal_rate = c("extra:infant_deaths_by_age", "births"),
    "population"
  )
  years <- Reduce(intersect, lapply(needs, planning_dataset_years))
  window <- planning_indicator_spec(id)$window
  if (window > 1) {
    # A triennial value is labelled by its last year and needs all three.
    years <- years[vapply(years, function(y) all(seq.int(y - window + 1L, y) %in% years), logical(1))]
  }
  sort(as.integer(years))
}

planning_proportional_years <- function(window = 3L) {
  years <- death_totals_years()
  sort(years[vapply(years, function(y) all(seq.int(y - window + 1L, y) %in% years), logical(1))])
}

# ---------------------------------------------------------
# Areas
# ---------------------------------------------------------
# Municipalities an area is built from. Portugal is every municipality; a
# region or health unit uses its lookup membership; anything else is taken to
# be a municipality.
planning_area_members <- function(area, lookup = get_nuts_lookup()) {
  municipalities <- as.character(lookup$municipality)
  if (identical(area, "Portugal")) return(sort(unique(municipalities)))
  if (area %in% municipalities) return(area)
  sort(unique(region_members(area, lookup)))
}

# One compact row per area label present in a year's files: population by broad
# age group, births, deaths and infant deaths, all sexes combined. Built once per
# year and cached, so an area of any size is a sum over at most 310 rows rather
# than a filter over the 70,000-row death file.
planning_component_columns <- c(
  "pop_total", "pop_0_14", "pop_15_64", "pop_15_plus", "pop_65_plus", "pop_75_plus",
  paste0("pop_f_", seq(15, 45, by = 5)),
  "births", "deaths", "infant_deaths",
  "rsi", "pensioners", "pension_value",
  "pp_share", "pp_weight",
  "waste_total", "waste_selective",
  "births_mother_total", "births_mother_lt20", "births_mother_ge35",
  paste0("births_mage_", seq(15, 45, by = 5)),
  "births_gest_total", "births_gest_known", "births_preterm",
  "neonatal_deaths", "early_neonatal_deaths", "postneonatal_deaths"
)

planning_extra_dir <- function(measure) file.path(infant_snapshot_root(), "planning_extra", measure)

planning_extra_years <- function(measure) {
  dir <- planning_extra_dir(measure)
  if (!dir.exists(dir)) return(integer(0))
  files <- list.files(dir, pattern = "^year_\\d+\\.rds$")
  sort(as.integer(sub("^year_(\\d+)\\.rds$", "\\1", files)))
}

read_planning_extra <- function(measure, year) {
  path <- file.path(planning_extra_dir(measure), paste0("year_", year, ".rds"))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

# Lower and upper bound of a five-year age label ("15 - 19 anos"), NA otherwise.
# Single years and wider groups ("15 - 49 anos") are rejected so a category is
# never counted twice.
planning_five_year_lower <- function(label) {
  label <- as.character(label)
  m <- regmatches(label, regexec("^(\\d+) - (\\d+) anos$", label))
  vapply(m, function(x) {
    if (length(x) == 3 && as.integer(x[3]) - as.integer(x[2]) == 4) as.integer(x[2]) else NA_integer_
  }, integer(1))
}

planning_open_lower <- function(label) {
  suppressWarnings(as.integer(sub("^(\\d+) (e|ou) mais anos$", "\\1", as.character(label))))
}

# Components of every area label present in a year's files, as one row per area.
# Each block returns a tibble keyed by area, or NULL when its source is absent
# for the year; the names of the non-NULL blocks record which datasets exist.
planning_component_blocks <- function(year) {
  blocks <- list()

  pop <- read_year_file("population", year)
  if (!is.null(pop)) {
    lower <- planning_age_lower(pop$age_band)
    hm <- pop$sex == "HM"
    female <- pop$sex == "M"
    frame <- tibble::tibble(area = pop$area, pop = pop$pop, lower = lower, hm = hm, female = female)
    blocks$pop <- frame %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        pop_total = sum(pop[hm], na.rm = TRUE),
        pop_0_14 = sum(pop[hm & lower < 15], na.rm = TRUE),
        pop_15_64 = sum(pop[hm & lower >= 15 & lower < 65], na.rm = TRUE),
        pop_15_plus = sum(pop[hm & lower >= 15], na.rm = TRUE),
        pop_65_plus = sum(pop[hm & lower >= 65], na.rm = TRUE),
        pop_75_plus = sum(pop[hm & lower >= 75], na.rm = TRUE),
        pop_f_15 = sum(pop[female & lower == 15], na.rm = TRUE),
        pop_f_20 = sum(pop[female & lower == 20], na.rm = TRUE),
        pop_f_25 = sum(pop[female & lower == 25], na.rm = TRUE),
        pop_f_30 = sum(pop[female & lower == 30], na.rm = TRUE),
        pop_f_35 = sum(pop[female & lower == 35], na.rm = TRUE),
        pop_f_40 = sum(pop[female & lower == 40], na.rm = TRUE),
        pop_f_45 = sum(pop[female & lower == 45], na.rm = TRUE),
        .groups = "drop"
      )
  }

  births <- read_year_file("births", year)
  if (!is.null(births)) {
    blocks$births <- births %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(births = sum(births, na.rm = TRUE), .groups = "drop")
  }

  totals <- read_death_totals_year(year)
  if (!is.null(totals)) {
    blocks$deaths <- totals %>%
      dplyr::filter(.data$sex == "HM", .data$cause == planning_all_causes) %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  }

  # The complete under-1 counts where the year has them; see
  # get_infant_death_data() for why the band-derived ones undercount.
  infant <- read_year_file("infant_totals", year)
  if (is.null(infant)) {
    infant <- read_year_file("infant_deaths", year)
    if (!is.null(infant)) infant <- infant[infant$cause == planning_all_causes, , drop = FALSE]
  }
  if (!is.null(infant)) {
    blocks$infant <- infant %>%
      dplyr::filter(.data$sex == "HM") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(infant_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  }

  total_of <- function(measure, column) {
    x <- read_planning_extra(measure, year)
    if (is.null(x)) return(NULL)
    x %>%
      dplyr::filter(.data$category == "Total") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(!!column := sum(value, na.rm = TRUE), .groups = "drop")
  }

  blocks$rsi <- total_of("rsi_beneficiaries", "rsi")

  # A mean is not additive: keep pensioners and pensioners x mean.
  pensioners <- total_of("pensioners", "pensioners")
  mean_pension <- total_of("pension_mean", "pension_mean")
  if (!is.null(pensioners) && !is.null(mean_pension)) {
    blocks$pensions <- dplyr::inner_join(pensioners, mean_pension, by = "area") %>%
      dplyr::transmute(area, pensioners, pension_value = pensioners * pension_mean)
  }

  # Purchasing power per capita is an index (Portugal = 100). Each
  # municipality's share of the national total divided by its index is its
  # share of the population INE used, so an area's index is sum(share) /
  # sum(share / index) x 100 - exact, with no population estimate needed.
  share <- total_of("purchasing_power_share", "pp_share")
  index <- total_of("purchasing_power_per_capita", "pp_index")
  if (!is.null(share) && !is.null(index)) {
    blocks$purchasing_power <- dplyr::inner_join(share, index, by = "area") %>%
      dplyr::filter(.data$pp_index > 0) %>%
      dplyr::transmute(area, pp_share, pp_weight = pp_share / pp_index * 100)
  }

  waste <- read_planning_extra("waste_collected", year)
  if (!is.null(waste)) {
    blocks$waste <- waste %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        waste_total = sum(value[category == "Total"], na.rm = TRUE),
        waste_selective = sum(value[grepl("selec?tiva", category, ignore.case = TRUE)], na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Births by mother's age. The indicator carries overlapping categories side
  # by side - single years, five-year groups, a 15-49 group, and both "50 - 54",
  # "50 e mais" and "55 e mais" - so it reads the five-year groups below the
  # lowest open group, plus that open group, and nothing else. Mothers under 15
  # are counted with 15-19 and mothers of 50 and over with 45-49 for the
  # fertility index, as INE does.
  mother <- read_planning_extra("births_by_mother_age", year)
  if (!is.null(mother)) {
    lower <- planning_five_year_lower(mother$category)
    open <- planning_open_lower(mother$category)
    open_min <- suppressWarnings(min(open, na.rm = TRUE))
    lower[!is.na(lower) & lower >= open_min] <- NA_integer_
    open[!is.na(open) & open > open_min] <- NA_integer_
    grouped <- dplyr::coalesce(lower, open)
    tfr_group <- pmin(pmax(grouped, 15L), 45L)
    frame <- tibble::tibble(area = mother$area, value = mother$value, category = mother$category,
                            grouped = grouped, tfr_group = tfr_group)
    blocks$mother <- frame %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        births_mother_total = sum(value[category == "Total"], na.rm = TRUE),
        births_mother_lt20 = sum(value[!is.na(grouped) & grouped < 20], na.rm = TRUE),
        births_mother_ge35 = sum(value[!is.na(grouped) & grouped >= 35], na.rm = TRUE),
        births_mage_15 = sum(value[!is.na(tfr_group) & tfr_group == 15], na.rm = TRUE),
        births_mage_20 = sum(value[!is.na(tfr_group) & tfr_group == 20], na.rm = TRUE),
        births_mage_25 = sum(value[!is.na(tfr_group) & tfr_group == 25], na.rm = TRUE),
        births_mage_30 = sum(value[!is.na(tfr_group) & tfr_group == 30], na.rm = TRUE),
        births_mage_35 = sum(value[!is.na(tfr_group) & tfr_group == 35], na.rm = TRUE),
        births_mage_40 = sum(value[!is.na(tfr_group) & tfr_group == 40], na.rm = TRUE),
        births_mage_45 = sum(value[!is.na(tfr_group) & tfr_group == 45], na.rm = TRUE),
        .groups = "drop"
      )
  }

  gestation <- read_planning_extra("births_by_gestation", year)
  if (!is.null(gestation)) {
    blocks$gestation <- gestation %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        births_gest_total = sum(value[category == "Total"], na.rm = TRUE),
        births_gest_known = sum(value[category == "Total"], na.rm = TRUE) -
          sum(value[grepl("ignorad", category, ignore.case = TRUE)], na.rm = TRUE),
        births_preterm = sum(value[category %in% c("Menos de 22 semanas", "22 - 27 semanas", "28 - 31 semanas", "32 - 36 semanas")], na.rm = TRUE),
        .groups = "drop"
      )
  }

  infant_age <- read_planning_extra("infant_deaths_by_age", year)
  if (!is.null(infant_age)) {
    blocks$infant_age <- infant_age %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        neonatal_deaths = sum(value[category == "Menos de 28 dias"], na.rm = TRUE),
        early_neonatal_deaths = sum(value[category == "Menos de 7 dias"], na.rm = TRUE),
        postneonatal_deaths = sum(value[category == "28 - 364 dias"], na.rm = TRUE),
        .groups = "drop"
      )
  }

  blocks[!vapply(blocks, is.null, logical(1))]
}

# Which block each component column comes from, so a column whose source is
# absent for the year stays NA instead of summing to zero.
planning_column_block <- function(column) {
  if (column %in% c("pop_total", "pop_0_14", "pop_15_64", "pop_15_plus", "pop_65_plus", "pop_75_plus") ||
      startsWith(column, "pop_f_")) return("pop")
  if (startsWith(column, "births_mage_") || startsWith(column, "births_mother_")) return("mother")
  switch(
    column,
    births = "births", deaths = "deaths", infant_deaths = "infant", rsi = "rsi",
    pensioners = , pension_value = "pensions",
    pp_share = , pp_weight = "purchasing_power",
    waste_total = , waste_selective = "waste",
    births_gest_total = , births_gest_known = , births_preterm = "gestation",
    neonatal_deaths = , early_neonatal_deaths = , postneonatal_deaths = "infant_age",
    stop("Unknown planning component: ", column, call. = FALSE)
  )
}

# One compact row per area label present in a year's files. Built once per year
# and cached, so an area of any size is a sum over at most 310 rows rather than
# a filter over the 70,000-row death file.
planning_year_components <- function(year) {
  key <- paste(infant_snapshot_root(), "components", year, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) {
    return(get(key, envir = planning_cache, inherits = FALSE))
  }

  blocks <- planning_component_blocks(year)
  table <- if (length(blocks) == 0) {
    tibble::tibble(area = character(0))
  } else {
    Reduce(function(x, y) dplyr::full_join(x, y, by = "area"), blocks)
  }
  for (column in planning_component_columns) {
    if (!column %in% names(table)) table[[column]] <- NA_real_
  }
  # An area absent from a present file counts as zero (INE omits empty cells);
  # a block absent for the year leaves its columns NA for every area.
  attr(table, "present") <- names(blocks)

  assign(key, table, envir = planning_cache)
  table
}

planning_age_lower <- function(age_band) {
  suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", as.character(age_band))))
}

# Components of one area in one year. Portugal and Continente use INE's
# published row where the year's file has one; everything else is the sum of
# its members. `members_found` reports how many members had a population row,
# so partial coverage is visible rather than read as a small area.
planning_components_year <- function(area, year, lookup = get_nuts_lookup()) {
  members <- planning_area_members(area, lookup)
  table <- planning_year_components(year)
  present <- attr(table, "present")

  published <- area %in% planning_published_areas && area %in% table$area
  value_of <- function(column) {
    if (!planning_column_block(column) %in% present) return(NA_real_)
    values <- table[[column]]
    if (published) {
      own <- values[table$area == area][[1]]
      if (!is.na(own)) return(own)
    }
    sum(values[table$area %in% members], na.rm = TRUE)
  }

  values <- lapply(planning_component_columns, value_of)
  names(values) <- planning_component_columns

  dplyr::bind_cols(
    tibble::tibble(
      area = area,
      year = as.integer(year),
      members = length(members),
      members_found = sum(members %in% table$area[!is.na(table$pop_total)])
    ),
    tibble::as_tibble(values)
  )
}

planning_components <- function(areas, years, lookup = get_nuts_lookup()) {
  grid <- tidyr::expand_grid(area = as.character(areas), year = as.integer(years))
  dplyr::bind_rows(purrr::map2(grid$area, grid$year, planning_components_year, lookup = lookup))
}

# ---------------------------------------------------------
# Indicators
# ---------------------------------------------------------
planning_poisson_rate <- function(events, denominator, multiplier) {
  if (!is.finite(events) || !is.finite(denominator) || denominator <= 0) {
    return(c(value = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  ci <- stats::poisson.test(round(events))$conf.int
  c(
    value = events / denominator * multiplier,
    lower = ci[[1]] / denominator * multiplier,
    upper = ci[[2]] / denominator * multiplier
  )
}

# One indicator for one area, ending in `year`. Takes the component table so a
# ranking over many areas reads the snapshots once.
#
# Population-based indices carry no interval: the population estimates are not
# a sample of events, and INE publishes no error for them. Event rates get an
# exact Poisson interval on the count, the denominator treated as fixed.
compute_planning_indicator <- function(components, id, area, year, undercount_years = integer(0)) {
  spec <- planning_indicator_spec(id)
  window_years <- seq.int(as.integer(year) - spec$window + 1L, as.integer(year))
  rows <- components[components$area == area & components$year %in% window_years, , drop = FALSE]

  empty <- tibble::tibble(
    area = area, year = as.integer(year), indicator = id, value = NA_real_,
    lower = NA_real_, upper = NA_real_, numerator = NA_real_, denominator = NA_real_, flag = ""
  )
  if (nrow(rows) < length(window_years)) return(empty)

  total <- function(col) {
    v <- rows[[col]]
    if (any(is.na(v))) NA_real_ else sum(v)
  }
  ratio <- function(num, den, multiplier) {
    if (!is.finite(num) || !is.finite(den) || den <= 0) NA_real_ else num / den * multiplier
  }

  out <- empty
  set <- function(value, numerator = NA_real_, denominator = NA_real_, lower = NA_real_, upper = NA_real_, flag = "") {
    out$value <<- value
    out$numerator <<- numerator
    out$denominator <<- denominator
    out$lower <<- lower
    out$upper <<- upper
    out$flag <<- flag
  }

  switch(
    id,
    pop_total = set(total("pop_total")),
    pct_0_14 = set(ratio(total("pop_0_14"), total("pop_total"), 100), total("pop_0_14"), total("pop_total")),
    pct_65_plus = set(ratio(total("pop_65_plus"), total("pop_total"), 100), total("pop_65_plus"), total("pop_total")),
    pct_75_plus = set(ratio(total("pop_75_plus"), total("pop_total"), 100), total("pop_75_plus"), total("pop_total")),
    ageing_index = set(ratio(total("pop_65_plus"), total("pop_0_14"), 100), total("pop_65_plus"), total("pop_0_14")),
    youth_dependency = set(ratio(total("pop_0_14"), total("pop_15_64"), 100), total("pop_0_14"), total("pop_15_64")),
    old_dependency = set(ratio(total("pop_65_plus"), total("pop_15_64"), 100), total("pop_65_plus"), total("pop_15_64")),
    births = {
      n <- total("births")
      ci <- if (is.finite(n)) stats::poisson.test(round(n))$conf.int else c(NA_real_, NA_real_)
      set(n, lower = ci[[1]], upper = ci[[2]])
    },
    birth_rate = {
      r <- planning_poisson_rate(total("births"), total("pop_total"), 1000)
      set(r[["value"]], total("births"), total("pop_total"), r[["lower"]], r[["upper"]])
    },
    deaths = {
      n <- total("deaths")
      ci <- if (is.finite(n)) stats::poisson.test(round(n))$conf.int else c(NA_real_, NA_real_)
      set(n, lower = ci[[1]], upper = ci[[2]])
    },
    death_rate = {
      r <- planning_poisson_rate(total("deaths"), total("pop_total"), 1000)
      set(r[["value"]], total("deaths"), total("pop_total"), r[["lower"]], r[["upper"]])
    },
    infant_rate = {
      births <- total("births")
      r <- planning_poisson_rate(total("infant_deaths"), births, 1000)
      flag <- if (is.finite(births) && infant_rate_is_unstable(births)) "*" else ""
      # Municipal under-1 counts incomplete in the source (1995-2001).
      if (!identical(area, "Portugal") && any(window_years %in% undercount_years)) {
        flag <- paste0(flag, "\u2020")
      }
      set(r[["value"]], total("infant_deaths"), births, r[["lower"]], r[["upper"]], flag)
    },
    neonatal_rate = , early_neonatal_rate = , postneonatal_rate = {
      column <- switch(id, neonatal_rate = "neonatal_deaths", early_neonatal_rate = "early_neonatal_deaths",
                       postneonatal_rate = "postneonatal_deaths")
      births <- total("births")
      r <- planning_poisson_rate(total(column), births, 1000)
      flag <- if (is.finite(births) && infant_rate_is_unstable(births)) "*" else ""
      set(r[["value"]], total(column), births, r[["lower"]], r[["upper"]], flag)
    },
    fertility_index = {
      # Sum of age-specific fertility rates over the seven five-year groups
      # 15-49, times five: births to mothers of each group per woman of that
      # group. Women are counted at mid-year, as the mean of the estimates at the
      # end of the previous year and of this one - which reproduces INE's
      # published index (Portugal 2019 1.43, 2023 1.32, 2024 1.27); the
      # end-of-year estimate alone reads up to 0.03 low.
      previous <- components[components$area == area & components$year == as.integer(year) - 1L, , drop = FALSE]
      women <- function(g) {
        column <- paste0("pop_f_", g)
        now <- total(column)
        before <- if (nrow(previous) == 1) previous[[column]] else NA_real_
        if (is.finite(before)) (now + before) / 2 else now
      }
      groups <- seq(15, 45, by = 5)
      rates <- vapply(groups, function(g) ratio(total(paste0("births_mage_", g)), women(g), 1), numeric(1))
      set(if (any(is.na(rates))) NA_real_ else sum(rates) * 5)
    },
    teen_births_pct = , older_births_pct = , preterm_pct = {
      num_col <- switch(id, teen_births_pct = "births_mother_lt20", older_births_pct = "births_mother_ge35",
                        preterm_pct = "births_preterm")
      den_col <- if (identical(id, "preterm_pct")) "births_gest_known" else "births_mother_total"
      num <- total(num_col)
      den <- total(den_col)
      ci <- if (is.finite(num) && is.finite(den) && den > 0) {
        compute_proportion_interval_safe(num, den)
      } else {
        c(NA_real_, NA_real_)
      }
      set(ratio(num, den, 100), num, den, ci[[1]], ci[[2]])
    },
    rsi_beneficiaries = set(total("rsi")),
    pensioners = set(total("pensioners")),
    rsi_rate = set(ratio(total("rsi"), total("pop_15_plus"), 1000), total("rsi"), total("pop_15_plus")),
    pensioners_rate = set(ratio(total("pensioners"), total("pop_15_plus"), 1000), total("pensioners"), total("pop_15_plus")),
    pension_mean = set(ratio(total("pension_value"), total("pensioners"), 1), total("pension_value"), total("pensioners")),
    purchasing_power = set(ratio(total("pp_share"), total("pp_weight"), 100), total("pp_share"), total("pp_weight")),
    waste_per_capita = set(ratio(total("waste_total"), total("pop_total"), 1000), total("waste_total"), total("pop_total")),
    waste_selective_per_capita = set(ratio(total("waste_selective"), total("pop_total"), 1000), total("waste_selective"), total("pop_total"))
  )

  out
}

planning_indicator_table <- function(areas, years, ids = PLANNING_INDICATORS$id, lookup = get_nuts_lookup()) {
  years <- as.integer(years)
  max_window <- max(PLANNING_INDICATORS$window[PLANNING_INDICATORS$id %in% ids])
  # One year earlier than the widest window: the fertility index needs the
  # previous year's population for its mid-year denominator.
  component_years <- seq.int(min(years) - max_window, max(years))
  components <- planning_components(areas, component_years, lookup)

  undercount <- if ("infant_rate" %in% ids) {
    infant_undercount_years(component_years, municipalities = lookup$municipality)
  } else {
    integer(0)
  }

  grid <- tidyr::expand_grid(area = as.character(areas), year = years, indicator = ids)
  dplyr::bind_rows(purrr::pmap(grid, function(area, year, indicator) {
    compute_planning_indicator(components, indicator, area, year, undercount_years = undercount)
  }))
}

planning_period_label <- function(id, year) {
  window <- planning_indicator_spec(id)$window
  if (window > 1) paste0(year - window + 1L, "-", year) else as.character(year)
}

planning_format_value <- function(value, digits) {
  digits <- rep_len(as.integer(digits), length(value))
  vapply(seq_along(value), function(i) {
    if (is.na(value[[i]])) return("—")
    formatC(value[[i]], format = "f", digits = digits[[i]], big.mark = ".", decimal.mark = ",")
  }, character(1))
}

# Wide profile for display: one row per indicator, one column per area, each
# cell "value (lower-upper)" with the thin-denominator star where it applies.
planning_profile_wide <- function(table) {
  if (nrow(table) == 0) return(tibble::tibble())
  specs <- PLANNING_INDICATORS
  table %>%
    dplyr::left_join(specs, by = c("indicator" = "id")) %>%
    dplyr::mutate(
      cell = paste0(
        planning_format_value(.data$value, .data$digits),
        .data$flag,
        ifelse(
          is.na(.data$lower),
          "",
          paste0(" (", planning_format_value(.data$lower, .data$digits), "-", planning_format_value(.data$upper, .data$digits), ")")
        )
      ),
      Período = purrr::map2_chr(.data$indicator, .data$year, planning_period_label),
      Tema = .data$theme,
      Indicador = paste0(.data$label, " [", .data$ref, "]"),
      Unidade = .data$unit,
      order = match(.data$indicator, specs$id)
    ) %>%
    dplyr::select(order, Tema, Indicador, Unidade, Período, area, cell) %>%
    tidyr::pivot_wider(names_from = area, values_from = cell) %>%
    dplyr::arrange(order) %>%
    dplyr::select(-order)
}

# ---------------------------------------------------------
# Proportional mortality by large cause groups (I45)
# ---------------------------------------------------------
# All ages, both sexes, pooled over a triennium: deaths of each group over all
# deaths. Read from the death totals, so it is exact for every area.
planning_year_causes <- function(year, sex = "HM") {
  key <- paste(infant_snapshot_root(), "causes", year, sex, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) {
    return(get(key, envir = planning_cache, inherits = FALSE))
  }
  totals <- read_death_totals_year(year)
  table <- if (is.null(totals)) {
    NULL
  } else {
    totals %>%
      dplyr::filter(.data$sex == .env$sex, .data$cause %in% c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)) %>%
      dplyr::group_by(area, cause) %>%
      dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  }
  assign(key, table, envir = planning_cache)
  table
}

planning_proportional <- function(area, end_year, window = 3L, sex = "HM", lookup = get_nuts_lookup()) {
  years <- seq.int(as.integer(end_year) - window + 1L, as.integer(end_year))
  members <- planning_area_members(area, lookup)
  causes <- c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)

  frames <- lapply(years, planning_year_causes, sex = sex)
  if (any(vapply(frames, is.null, logical(1)))) return(tibble::tibble())

  pooled <- dplyr::bind_rows(lapply(frames, function(frame) {
    rows <- if (area %in% planning_published_areas && area %in% frame$area) {
      frame[frame$area == area, , drop = FALSE]
    } else {
      frame[frame$area %in% members, , drop = FALSE]
    }
    rows
  }))
  counts <- vapply(causes, function(cause) sum(pooled$deaths[pooled$cause == cause]), numeric(1))

  total <- counts[[1]]
  groups <- counts[-1]
  other <- max(total - sum(groups), 0)

  rows <- tibble::tibble(
    code = c(PLANNING_CAUSE_GROUPS$code, "Outras"),
    group = c(PLANNING_CAUSE_GROUPS$label, "Restantes causas"),
    deaths = c(unname(groups), other)
  )

  intervals <- lapply(rows$deaths, compute_proportion_interval_safe, total = total)

  dplyr::bind_rows(
    tibble::tibble(code = "C00", group = "Todas as causas de morte", deaths = total),
    rows
  ) %>%
    dplyr::mutate(
      area = area,
      period = paste0(min(years), "-", max(years)),
      share = if (total > 0) .data$deaths / total * 100 else NA_real_,
      lower = c(100, vapply(intervals, `[[`, numeric(1), 1)),
      upper = c(100, vapply(intervals, `[[`, numeric(1), 2))
    )
}

compute_proportion_interval_safe <- function(deaths, total) {
  if (!is.finite(total) || total <= 0) return(c(NA_real_, NA_real_))
  ci <- stats::binom.test(min(round(deaths), round(total)), round(total))$conf.int * 100
  c(ci[[1]], ci[[2]])
}

# ---------------------------------------------------------
# Population pyramid (I3)
# ---------------------------------------------------------
planning_pyramid <- function(area, year, lookup = get_nuts_lookup()) {
  pop <- planning_read("population", year)
  if (is.null(pop)) return(tibble::tibble())
  members <- planning_area_members(area, lookup)
  use_published <- area %in% planning_published_areas && area %in% pop$area
  rows <- if (use_published) pop[pop$area == area, , drop = FALSE] else pop[pop$area %in% members, , drop = FALSE]

  rows %>%
    dplyr::filter(.data$sex %in% c("H", "M")) %>%
    dplyr::group_by(age_band, sex) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      area = area,
      year = as.integer(year),
      lower = planning_age_lower(.data$age_band),
      share = .data$pop / sum(.data$pop) * 100
    ) %>%
    dplyr::arrange(lower)
}

# ---------------------------------------------------------
# Ranking across health units
# ---------------------------------------------------------
# Every ULS (and the grouped units covering the split municipalities) for one
# indicator and year, so a unit can be read against the others and against
# Portugal. ARS are left out: they are sums of the same units.
planning_uls_units <- function() {
  lookup <- get_health_lookup()
  if (is.null(lookup) || nrow(lookup) == 0) return(character(0))
  sort(unique(as.character(lookup$unit[lookup$kind != "ARS"])))
}
