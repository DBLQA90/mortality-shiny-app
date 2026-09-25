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

# `comparable`: whether the indicator does not depend on the size of the area -
# a rate, share, index or per-capita value. Only those are set beside
# comparators; a count of people or events is not, since a municipality's count
# says nothing next to its region's.
PLANNING_INDICATORS <- tibble::tribble(
  ~id,                   ~theme,                ~label,                                                    ~ref,  ~unit,                      ~window, ~digits, ~comparable,
  "pop_total",           "Demografia",          "População residente (estimativa)",                        "I1",  "N.º",                      1L,      0L, FALSE,
  "census_population",   "Demografia",          "População residente nos Censos",                          "I2",  "N.º",                      1L,      0L, FALSE,
  "census_population_change", "Demografia",     "Variação da população desde o censo anterior",            "I2",  "%",                        1L,      1L, TRUE,
  "pct_0_14",            "Demografia",          "Proporção de jovens (0-14 anos)",                         "I1",  "%",                        1L,      1L, TRUE,
  "pct_65_plus",         "Demografia",          "Proporção de idosos (65 e mais anos)",                    "I1",  "%",                        1L,      1L, TRUE,
  "pct_75_plus",         "Demografia",          "Proporção de 75 e mais anos",                             "I1",  "%",                        1L,      1L, TRUE,
  "ageing_index",        "Demografia",          "Índice de envelhecimento",                                "I4",  "por 100 jovens",           1L,      1L, TRUE,
  "youth_dependency",    "Demografia",          "Índice de dependência de jovens",                         "I5",  "por 100 em idade activa",  1L,      1L, TRUE,
  "old_dependency",      "Demografia",          "Índice de dependência de idosos",                         "I6",  "por 100 em idade activa",  1L,      1L, TRUE,
  "pct_education_none",  "Educação",            "População sem nível de escolaridade completo (Censos)",   "I24", "%",                        1L,      1L, TRUE,
  "pct_education_basic", "Educação",            "População com o ensino básico (Censos)",                  "I24", "%",                        1L,      1L, TRUE,
  "pct_education_secondary", "Educação",        "População com o ensino secundário (Censos)",              "I24", "%",                        1L,      1L, TRUE,
  "pct_education_higher", "Educação",           "População com o ensino superior (Censos)",                "I24", "%",                        1L,      1L, TRUE,
  "illiteracy_rate",     "Educação",            "Taxa de analfabetismo (Censos)",                          "I26", "%",                        1L,      1L, TRUE,
  "births",              "Natalidade",          "Nados-vivos",                                             "I7",  "N.º",                      1L,      0L, FALSE,
  "birth_rate",          "Natalidade",          "Taxa bruta de natalidade",                                "I8",  "‰",                        1L,      1L, TRUE,
  "fertility_index",     "Natalidade",          "Índice sintético de fecundidade",                         "I9",  "filhos por mulher",        1L,      2L, TRUE,
  "teen_births_pct",     "Natalidade",          "Nascimentos em mães com menos de 20 anos (triénio)",      "I32", "%",                        3L,      1L, TRUE,
  "older_births_pct",    "Natalidade",          "Nascimentos em mães com 35 e mais anos (triénio)",        "I33", "%",                        3L,      1L, TRUE,
  "preterm_pct",         "Natalidade",          "Nascimentos pré-termo (triénio)",                         "I35", "%",                        3L,      1L, TRUE,
  "rsi_beneficiaries",   "Contexto social",     "Beneficiários do rendimento social de inserção",          "I13", "N.º",                      1L,      0L, FALSE,
  "rsi_rate",            "Contexto social",     "Beneficiários do RSI por 1.000 habitantes com 15+ anos",  "I14", "‰",                        1L,      1L, TRUE,
  "pensioners",          "Contexto social",     "Pensionistas da segurança social",                        "I15", "N.º",                      1L,      0L, FALSE,
  "pensioners_rate",     "Contexto social",     "Pensionistas por 1.000 habitantes com 15+ anos",          "I16", "‰",                        1L,      1L, TRUE,
  "pension_mean",        "Contexto social",     "Valor médio anual das pensões",                           "I17", "€",                        1L,      0L, TRUE,
  "purchasing_power",    "Contexto social",     "Poder de compra per capita",                              "I28", "Portugal = 100",           1L,      1L, TRUE,
  "earnings_mean",       "Contexto social",     "Ganho médio mensal dos trabalhadores por conta de outrem", "I27", "€",                        1L,      0L, TRUE,
  "employees",           "Contexto social",     "Trabalhadores por conta de outrem",                        "I12", "N.º",                      1L,      0L, FALSE,
  "pct_employees_primary", "Contexto social",   "Trabalhadores no sector primário",                        "I12", "%",                        1L,      1L, TRUE,
  "pct_employees_secondary", "Contexto social", "Trabalhadores no sector secundário",                      "I12", "%",                        1L,      1L, TRUE,
  "pct_employees_tertiary", "Contexto social",  "Trabalhadores no sector terciário",                       "I12", "%",                        1L,      1L, TRUE,
  "waste_per_capita",    "Ambiente",            "Resíduos urbanos recolhidos por habitante",               "I64", "kg/hab.",                  1L,      0L, TRUE,
  "waste_selective_per_capita", "Ambiente",     "Resíduos recolhidos selectivamente por habitante",        "I65", "kg/hab.",                  1L,      0L, TRUE,
  "life_expectancy",     "Mortalidade",         "Esperança de vida à nascença (triénio)",                  "I10", "anos",                     3L,      1L, TRUE,
  "life_expectancy_men", "Mortalidade",         "Esperança de vida à nascença, homens (triénio)",          "I10", "anos",                     3L,      1L, TRUE,
  "life_expectancy_women", "Mortalidade",       "Esperança de vida à nascença, mulheres (triénio)",        "I10", "anos",                     3L,      1L, TRUE,
  "life_expectancy_65",  "Mortalidade",         "Esperança de vida aos 65 anos (triénio)",                 "I10", "anos",                     3L,      1L, TRUE,
  "life_expectancy_65_men", "Mortalidade",      "Esperança de vida aos 65 anos, homens (triénio)",         "I10", "anos",                     3L,      1L, TRUE,
  "life_expectancy_65_women", "Mortalidade",    "Esperança de vida aos 65 anos, mulheres (triénio)",       "I10", "anos",                     3L,      1L, TRUE,
  "low_birth_weight_pct", "Natalidade",         "Nascimentos com baixo peso, menos de 2.500 g (triénio)", "I36", "%",                       3L,      1L, TRUE,
  "deaths",              "Mortalidade",         "Óbitos",                                                  "I37", "N.º",                      1L,      0L, FALSE,
  "death_rate",          "Mortalidade",         "Taxa bruta de mortalidade",                               "I38", "‰",                        1L,      1L, TRUE,
  "infant_rate",         "Mortalidade",         "Taxa de mortalidade infantil (triénio)",                  "I39", "‰ nados-vivos",            3L,      1L, TRUE,
  "neonatal_rate",       "Mortalidade",         "Taxa de mortalidade neonatal (triénio)",                  "I40", "‰ nados-vivos",            3L,      1L, TRUE,
  "early_neonatal_rate", "Mortalidade",         "Taxa de mortalidade neonatal precoce (triénio)",          "I41", "‰ nados-vivos",            3L,      1L, TRUE,
  "postneonatal_rate",   "Mortalidade",         "Taxa de mortalidade pós-neonatal (triénio)",              "I42", "‰ nados-vivos",            3L,      1L, TRUE,
  "late_fetal_rate",     "Mortalidade",         "Taxa de mortalidade fetal tardia (triénio)",              "I43", "‰ nascimentos",            3L,      1L, TRUE,
  "perinatal_rate",      "Mortalidade",         "Taxa de mortalidade perinatal (triénio)",                 "I44", "‰ nascimentos",            3L,      1L, TRUE,
  # Not in the workbook: standardised, premature and avoidable mortality (R/planning_standardised.R).
  "smr_all",             "Mortalidade padronizada", "Razão padronizada de mortalidade, todas as causas (triénio)", "novo", "Portugal = 100",     3L,      1L, TRUE,
  "dsr_all",             "Mortalidade padronizada", "Taxa de mortalidade padronizada pela idade (triénio)",  "novo", "por 100.000 hab.",       3L,      1L, TRUE,
  "dsr_premature",       "Mortalidade padronizada", "Taxa de mortalidade prematura padronizada, < 75 anos (triénio)", "novo", "por 100.000 hab.", 3L,     1L, TRUE,
  "premature_deaths",    "Mortalidade padronizada", "Óbitos prematuros, < 75 anos (triénio)",               "novo", "N.º",                    3L,      0L, FALSE,
  "dsr_preventable",     "Mortalidade padronizada", "Mortalidade evitável por prevenção, < 75 anos, padronizada (triénio)", "novo", "por 100.000 hab.", 3L, 1L, TRUE,
  "dsr_treatable",       "Mortalidade padronizada", "Mortalidade evitável por cuidados de saúde (tratável), < 75 anos, padronizada (triénio)", "novo", "por 100.000 hab.", 3L, 1L, TRUE,
  "avoidable_deaths",    "Mortalidade padronizada", "Óbitos evitáveis, < 75 anos (triénio)",                "novo", "N.º",                    3L,      0L, FALSE,
  "ypll_rate",           "Mortalidade padronizada", "Anos potenciais de vida perdidos antes dos 70 anos (triénio)", "novo", "por 100.000 hab. < 70", 3L, 0L, TRUE
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
# Portugal as the sum of its 308 municipalities, beside INE's published row.
# The published row also counts events of residents whose municipality is
# unknown (0.3-0.9% of deaths), which no region, ULS or municipality can hold;
# the municipal sum compares like with like. Not in planning_published_areas,
# so every path treats it as an ordinary municipal sum.
PLANNING_PORTUGAL_MUNICIPAL <- "Portugal (soma dos municípios)"

# Which Portugal a view compares against: "published" or "municipal".
planning_portugal_area <- function(mode = "published") {
  if (identical(mode, "municipal")) PLANNING_PORTUGAL_MUNICIPAL else "Portugal"
}

# Significance against a benchmark, as PHE Fingertips reports it: an area is
# higher (lower) when its whole 95% interval lies above (below) the benchmark's
# value, similar when the interval contains it. Only values with an interval
# (event rates, proportions, life expectancy) are tested; the benchmark itself
# is not.
PLANNING_SIGNIFICANCE_MARKS <- c(Superior = "\u25b2", Inferior = "\u25bc", Semelhante = "=")

planning_significance <- function(lower, upper, benchmark) {
  out <- rep(NA_character_, length(lower))
  ok <- !is.na(lower) & !is.na(upper) & !is.na(benchmark)
  out[ok] <- ifelse(lower[ok] > benchmark[ok], "Superior", ifelse(upper[ok] < benchmark[ok], "Inferior", "Semelhante"))
  out
}

planning_add_significance <- function(table, benchmark_area) {
  if (nrow(table) == 0) return(dplyr::mutate(table, benchmark = numeric(0), significance = character(0)))
  reference <- table[table$area == benchmark_area, c("indicator", "year", "value"), drop = FALSE]
  names(reference)[3] <- "benchmark"
  out <- dplyr::left_join(dplyr::select(table, -dplyr::any_of(c("benchmark", "significance"))), reference, by = c("indicator", "year"))
  out$significance <- planning_significance(out$lower, out$upper, out$benchmark)
  out$significance[out$area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL)] <- NA_character_
  # Absolute counts depend on the size of the area: never compared.
  out$significance[!out$indicator %in% PLANNING_INDICATORS$id[PLANNING_INDICATORS$comparable]] <- NA_character_
  out
}

planning_significance_mark <- function(significance) {
  mark <- unname(PLANNING_SIGNIFICANCE_MARKS[significance])
  ifelse(is.na(mark), "", paste0(" ", mark))
}

# Version of the method, as in the methodological note. Part of the cache key of
# the all-areas export: raise it whenever a change alters published values.
PLANNING_METHOD_VERSION <- "1.7"

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
  key <- paste(infant_snapshot_root(), "totals", year, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  out <- NULL
  for (indicator in c("0013166", "0008206")) {
    path <- file.path(death_totals_dir(), indicator, paste0("year_", year, ".rds"))
    if (file.exists(path)) {
      out <- repair_death_totals(readRDS(path), year)
      break
    }
  }
  assign(key, out, envir = planning_cache)
  out
}

# INE leaves the all-ages cell blank for a few municipality-years (Vimioso
# 2024, Alfândega da Fé and Miranda do Douro 2015, some in the 1990s) while its
# age bands are filled in; older snapshots stored the blank as 0. A total can
# never be below the sum of its own bands, so for all causes and the I45 groups
# take the larger of the two. Repaired cells are marked in `repaired`.
repair_death_totals <- function(totals, year) {
  totals$repaired <- FALSE
  causes <- intersect(unique(c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)), unique(totals$cause))
  for (indicator in c("0013166", "0008206")) {
    dir <- file.path(infant_snapshot_root(), "deaths", indicator, paste0("year_", year))
    if (dir.exists(dir)) break
  }
  if (!dir.exists(dir)) return(totals)
  for (cause in causes) {
    file <- file.path(dir, paste0("cause_", planning_cause_file_token(cause), ".rds"))
    if (!file.exists(file)) next
    banded <- readRDS(file)
    sums <- stats::aggregate(banded$deaths, by = list(area = banded$area, sex = banded$sex), FUN = sum, na.rm = TRUE)
    rows <- which(totals$cause == cause)
    banded_sum <- sums$x[match(paste(totals$area[rows], totals$sex[rows]), paste(sums$area, sums$sex))]
    short <- !is.na(banded_sum) & (is.na(totals$deaths[rows]) | totals$deaths[rows] < banded_sum)
    totals$deaths[rows[short]] <- banded_sum[short]
    totals$repaired[rows[short]] <- TRUE
  }
  totals
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
  if (id %in% life_expectancy_ids) return(life_expectancy_years())
  if (id %in% standardised_ids) return(planning_standardised_years())
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
    low_birth_weight_pct = "extra:births_by_weight",
    earnings_mean = c("extra:earnings_mean", "extra:employees_by_sector"),
    employees = , pct_employees_primary = , pct_employees_secondary = , pct_employees_tertiary =
      "extra:employees_by_sector",
    census_population = "extra:census_population",
    census_population_change = "extra:census_population",
    pct_education_none = , pct_education_basic = , pct_education_secondary = , pct_education_higher =
      c("extra:census_education", "extra:census_population"),
    illiteracy_rate = c("extra:census_illiteracy_rate", "extra:census_population_by_age"),
    late_fetal_rate = , perinatal_rate = c("extra:perinatal_deaths", "extra:infant_deaths_by_age", "births"),
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
  if (area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL)) return(sort(unique(municipalities)))
  # A ULS that shares a municipality serves every municipality it touches; the
  # parish weights decide how much of each it takes (R/planning_parish.R).
  if (area %in% planning_parish_units()) {
    return(sort(unique(c(planning_parish_members(area), intersect(region_members(area, lookup), municipalities)))))
  }
  if (area %in% municipalities) return(area)
  sort(unique(region_members(area, lookup)))
}

# One compact row per area label present in a year's files: population by broad
# age group, births, deaths and infant deaths, all sexes combined. Built once per
# year and cached, so an area of any size is a sum over at most 310 rows rather
# than a filter over the 70,000-row death file.
# Blocks INE reports for every municipality, by their key column: a blank is
# "not reported", never zero. See planning_year_components().
PLANNING_COMPLETE_BLOCKS <- c(waste = "waste_total", pensions = "pensioners", employees = "employees_total", earnings = "earnings_value")

# Municipalities whose values INE files under another one. An area then has a
# value only if it holds both or neither; otherwise it is missing.
# - Odivelas, Trofa and Vizela were created in 1998, out of Loures, Santo Tirso
#   and (mostly) Guimarães. INE's population estimates are back-cast to today's
#   boundaries, but births, deaths and the 1991 census up to 1998 sit with the
#   parent: Loures alone would read its births over a population without
#   Odivelas, and Odivelas zero. Every block but the population is joint to
#   `until`. (Vizela also took a few parishes of Felgueiras and Lousada, which
#   this treats as negligible.)
# - With `blocks`, only when the municipality has no value that year: the
#   joint Loures-Odivelas waste service (SIMAR) is recorded wholly under
#   Loures, with Odivelas blank; purchasing power has no row for the three
#   until 2000.
PLANNING_JOINT_REPORTING <- tibble::tribble(
  ~municipality, ~holder,       ~until, ~blocks,
  "Odivelas",    "Loures",      1998L,  list(c("waste", "purchasing_power")),
  "Trofa",       "Santo Tirso", 1998L,  list("purchasing_power"),
  "Vizela",      "Guimarães",   1998L,  list("purchasing_power")
)

# Areas that hold one of a joint pair but not the other in `year`, for a block
# (FALSE for every area when nothing is joint). `unreported` is the
# municipality's own value missing that year, for the data-driven pairs.
planning_joint_split <- function(membership, year, block, unreported = NULL) {
  municipalities <- colnames(membership)
  split <- rep(FALSE, nrow(membership))
  if (identical(block, "pop")) return(split)
  for (i in seq_len(nrow(PLANNING_JOINT_REPORTING))) {
    joint <- PLANNING_JOINT_REPORTING[i, ]
    m <- match(joint$municipality, municipalities); h <- match(joint$holder, municipalities)
    if (is.na(m) || is.na(h)) next
    by_year <- year <= joint$until
    by_data <- block %in% unlist(joint$blocks) && !is.null(unreported) && isTRUE(unreported[[joint$municipality]])
    if (by_year || by_data) split <- split | xor(membership[, m] > 0, membership[, h] > 0)
  }
  split
}

# Education by age group (2011 and 2021 censuses only): one column per level
# and five-year group, 0 standing for "under 15" and 75 for "75 and over", so
# any "N and over" is a sum of columns.
PLANNING_EDUCATION_AGES <- c(0L, seq.int(15L, 75L, by = 5L))
PLANNING_EDUCATION_LEVELS <- c(total = "Total", none = "Nenhum", basic = "Básico", secondary = "Secundário", higher = "Superior")
planning_education_age_columns <- function() {
  as.vector(outer(names(PLANNING_EDUCATION_LEVELS), PLANNING_EDUCATION_AGES, function(l, a) paste0("census_edu_", l, "_", a)))
}
PLANNING_EDUCATION_IDS <- c("pct_education_none", "pct_education_basic", "pct_education_secondary", "pct_education_higher")

# The minimum ages the education indicators can be restricted to.
planning_education_age_choices <- function() {
  c("Toda a população" = 0L, stats::setNames(PLANNING_EDUCATION_AGES[-1], paste0(PLANNING_EDUCATION_AGES[-1], " e mais anos")))
}

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
  "births_weight_known", "births_low_weight",
  "earnings_value", "employees_total", "employees_primary", "employees_secondary", "employees_tertiary", "employees_estimated",
  "census_pop", "census_pop_10plus", "census_illiterate",
  "census_education_total", "census_education_basic", "census_education_secondary", "census_education_higher",
  planning_education_age_columns(),
  "neonatal_deaths", "early_neonatal_deaths", "postneonatal_deaths", "perinatal_deaths"
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

# Lower bound in grams of a birth-weight band, NA for "Total" and "Ignorada".
planning_weight_lower <- function(label) {
  label <- gsub("[[:space:]\u00a0\u202f]", "", as.character(label))
  out <- rep(NA_real_, length(label))
  under <- grepl("^Menosde[0-9]+g", label)
  out[under] <- 0
  band <- grepl("^[0-9]+-[0-9]+g$", label)
  out[band] <- as.numeric(sub("^([0-9]+)-.*$", "\\1", label[band]))
  open <- grepl("^[0-9]+gemais$", label)
  out[open] <- as.numeric(sub("^([0-9]+)g.*$", "\\1", label[open]))
  out
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

  # Birth weight: the bands read "2 000 - 2 499 g" (with non-breaking spaces),
  # "Menos de 500 g" and "5 000 g e mais". Low birth weight is under 2,500 g,
  # over the births whose weight is known.
  weight <- read_planning_extra("births_by_weight", year)
  if (!is.null(weight)) {
    lower <- planning_weight_lower(weight$category)
    blocks$weight <- tibble::tibble(area = weight$area, category = weight$category, value = weight$value, lower = lower) %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        births_weight_known = sum(value[category == "Total"], na.rm = TRUE) -
          sum(value[grepl("ignorad", category, ignore.case = TRUE)], na.rm = TRUE),
        births_low_weight = sum(value[!is.na(lower) & lower < 2500], na.rm = TRUE),
        .groups = "drop"
      )
  }

  perinatal <- read_planning_extra("perinatal_deaths", year)
  if (!is.null(perinatal)) {
    blocks$perinatal <- perinatal %>%
      dplyr::filter(.data$category == "Total") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(perinatal_deaths = sum(value, na.rm = TRUE), .groups = "drop")
  }

  # Employees by sector, and their average earnings. The mean is not additive,
  # so it is kept as earnings x employees and divided back after summing. The
  # source counts employees where they work, not where they live.
  employees <- read_planning_extra("employees_by_sector", year)
  if (!is.null(employees)) {
    blocks$employees <- employees %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        employees_total = sum(value[category == "Total"], na.rm = TRUE),
        employees_primary = sum(value[grepl("^Agricultura", category)], na.rm = TRUE),
        employees_secondary = sum(value[grepl("^Indústria", category)], na.rm = TRUE),
        employees_tertiary = sum(value[grepl("^Serviços", category)], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      planning_estimate_suppressed_sectors()
  }
  earnings <- read_planning_extra("earnings_mean", year)
  if (!is.null(earnings) && !is.null(blocks$employees)) {
    blocks$earnings <- earnings %>%
      dplyr::filter(.data$category == "Total") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(mean_earnings = sum(value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::inner_join(blocks$employees[, c("area", "employees_total")], by = "area") %>%
      dplyr::transmute(area, earnings_value = .data$mean_earnings * .data$employees_total)
  }

  # Census measures: only the census years have files.
  census_population <- read_planning_extra("census_population", year)
  census_age <- read_planning_extra("census_population_by_age", year)
  census_rate <- read_planning_extra("census_illiteracy_rate", year)
  if (!is.null(census_population)) {
    blocks$census_population <- census_population %>%
      dplyr::filter(.data$category == "Total") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(census_pop = sum(value, na.rm = TRUE), .groups = "drop")
  }
  if (!is.null(census_age)) {
    lower <- planning_five_year_lower(census_age$category)
    open <- planning_open_lower(census_age$category)
    blocks$census_age <- tibble::tibble(area = census_age$area, value = census_age$value,
                                        lower = dplyr::coalesce(lower, open)) %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(census_pop_10plus = sum(value[!is.na(lower) & lower >= 10], na.rm = TRUE), .groups = "drop")
  }
  # INE publishes the illiteracy rate by municipality, not the count. The count
  # each rate implies is additive, so an area's rate is the population-weighted
  # mean of its municipalities'.
  if (!is.null(census_rate) && !is.null(blocks$census_age)) {
    blocks$census_illiteracy <- census_rate %>%
      dplyr::filter(.data$category == "Total") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(rate = sum(value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::inner_join(blocks$census_age, by = "area") %>%
      dplyr::transmute(area, census_illiterate = .data$rate / 100 * .data$census_pop_10plus)
  }
  census_education <- read_planning_extra("census_education", year)
  if (!is.null(census_education)) {
    blocks$census_education <- census_education %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        census_education_total = sum(value[category == "Total"], na.rm = TRUE),
        census_education_basic = sum(value[grepl("Básico", category)], na.rm = TRUE),
        census_education_secondary = sum(value[grepl("Secundário", category)], na.rm = TRUE),
        census_education_higher = sum(value[grepl("Superior", category)], na.rm = TRUE),
        .groups = "drop"
      )
  }

  education_age <- read_planning_extra("census_education_by_age", year)
  if (!is.null(education_age)) {
    education_age$column <- paste0("census_edu_", names(PLANNING_EDUCATION_LEVELS)[match(education_age$category, PLANNING_EDUCATION_LEVELS)], "_", education_age$age)
    education_age <- education_age[education_age$column %in% planning_education_age_columns(), , drop = FALSE]
    wide <- tidyr::pivot_wider(
      dplyr::summarise(dplyr::group_by(education_age, area, column), value = sum(value, na.rm = TRUE), .groups = "drop"),
      names_from = column, values_from = value, values_fill = 0
    )
    for (column in setdiff(planning_education_age_columns(), names(wide))) wide[[column]] <- 0
    blocks$census_education_age <- wide
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
  if (startsWith(column, "census_edu_")) return("census_education_age")
  switch(
    column,
    births = "births", deaths = "deaths", infant_deaths = "infant", rsi = "rsi",
    pensioners = , pension_value = "pensions",
    pp_share = , pp_weight = "purchasing_power",
    waste_total = , waste_selective = "waste",
    births_gest_total = , births_gest_known = , births_preterm = "gestation",
    births_weight_known = , births_low_weight = "weight",
    employees_total = , employees_primary = , employees_secondary = , employees_tertiary = , employees_estimated = "employees",
    earnings_value = "earnings",
    census_pop = "census_population",
    census_pop_10plus = "census_age",
    census_illiterate = "census_illiteracy",
    census_education_total = , census_education_basic = , census_education_secondary = ,
    census_education_higher = "census_education",
    perinatal_deaths = "perinatal",
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
  # a block absent for the year leaves its columns NA for every area. The
  # blocks INE reports for every municipality are the exception: there a blank
  # (stored as 0 by older fetches) means not reported, and stays NA.
  published <- table$area %in% planning_published_areas
  for (block in intersect(names(PLANNING_COMPLETE_BLOCKS), names(blocks))) {
    reported <- table[[PLANNING_COMPLETE_BLOCKS[[block]]]]
    blank <- !published & (is.na(reported) | reported <= 0)
    for (column in planning_component_columns[vapply(planning_component_columns, planning_column_block, character(1)) == block]) {
      table[[column]][blank] <- NA_real_
    }
  }
  attr(table, "present") <- names(blocks)

  assign(key, table, envir = planning_cache)
  table
}

planning_age_lower <- function(age_band) {
  suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", as.character(age_band))))
}

# INE suppresses sector cells for confidentiality - always two of the three,
# so neither can be recovered by subtraction (6 municipalities in 2024, 41-49
# in 2013-2016; Vizela, Marinha Grande...). The blanks arrive as 0, which would
# put all of the municipality's employees in its one published sector. The
# hidden remainder (total minus the published sectors) is split between the
# hidden sectors in the proportions they have in the rest of the NUTS III that
# year (Portugal's where that is empty), and counted in `employees_estimated`
# so the shares that include it can be flagged.
PLANNING_ESTIMATED_FLAG <- "≈"
PLANNING_ESTIMATED_SHARE <- 0.01

planning_estimate_suppressed_sectors <- function(employees, lookup = get_nuts_lookup()) {
  sectors <- c("employees_primary", "employees_secondary", "employees_tertiary")
  employees$employees_estimated <- 0
  municipal <- !employees$area %in% planning_published_areas
  published <- as.matrix(employees[, sectors])
  gap <- employees$employees_total - rowSums(published)
  hidden <- municipal & gap > 0.5 & rowSums(published <= 0) >= 1
  if (!any(hidden)) return(employees)

  nuts3 <- lookup$nuts3[match(employees$area, lookup$municipality)]
  complete <- municipal & !hidden
  reference <- function(group) {
    rows <- complete & !is.na(nuts3) & nuts3 == group
    if (!any(rows)) rows <- complete
    colSums(published[rows, , drop = FALSE])
  }
  for (i in which(hidden)) {
    missing <- published[i, ] <= 0
    weights <- reference(nuts3[[i]])[missing]
    if (sum(weights) <= 0) weights[] <- 1
    published[i, missing] <- gap[[i]] * weights / sum(weights)
    employees$employees_estimated[[i]] <- gap[[i]]
  }
  employees[, sectors] <- published
  employees
}

# Membership of each area in the municipalities of the lookup, as a 0/1 matrix
# (areas x municipalities). Summing components over any set of areas is then one
# matrix product per year, which is what lets an export cover every area.
planning_membership_matrix <- function(areas, lookup = get_nuts_lookup(), mode = PLANNING_DEFAULT_SPLIT_MODE, basis = "total", year = NULL) {
  areas <- unique(as.character(areas))
  municipalities <- sort(unique(as.character(lookup$municipality)))
  membership <- matrix(0, length(areas), length(municipalities), dimnames = list(areas, municipalities))
  for (i in seq_along(areas)) {
    members <- intersect(planning_area_members(areas[[i]], lookup), municipalities)
    if (length(members) > 0) membership[i, members] <- 1
  }
  # In the parish reading, the six ULS that share a municipality take the
  # share of it their parishes hold in the census.
  if (identical(mode, "parish") && any(areas %in% planning_parish_units())) {
    membership <- membership * planning_parish_weights(areas, municipalities, basis, mode, year)
  }
  membership
}

# Components of each area in each year. Portugal and Continente use INE's
# published row where the year's file has one; everything else is the sum of
# its municipalities. `members_found` counts members with a population row, so
# partial coverage is visible rather than read as a small area.
planning_components <- function(areas, years, lookup = get_nuts_lookup(), mode = PLANNING_DEFAULT_SPLIT_MODE) {
  membership <- planning_membership_matrix(areas, lookup, mode)
  areas <- rownames(membership)
  municipalities <- colnames(membership)
  columns <- planning_component_columns
  blocks <- vapply(columns, planning_column_block, character(1))

  frames <- lapply(as.integer(years), function(year) {
    table <- planning_year_components(year)
    present <- attr(table, "present")

    rows <- match(municipalities, table$area)
    found <- !is.na(rows)
    values <- matrix(0, length(municipalities), length(columns))
    if (any(found)) {
      values[found, ] <- as.matrix(table[rows[found], columns, drop = FALSE])
    }
    with_population <- found & !is.na(table$pop_total[ifelse(found, rows, 1L)])

    # A municipality that exists but was not reported in a complete block makes
    # every area containing it missing, instead of dropping out of the sum.
    complete <- blocks %in% names(PLANNING_COMPLETE_BLOCKS) & blocks %in% present
    unreported <- matrix(FALSE, length(municipalities), length(columns))
    unreported[, complete] <- with_population & is.na(values[, complete, drop = FALSE])
    split <- matrix(FALSE, length(areas), length(columns))
    joint_members <- intersect(PLANNING_JOINT_REPORTING$municipality, municipalities)
    for (block in unique(blocks[blocks %in% present])) {
      cols <- which(blocks == block)
      # "No value" for a joint pair: blank, or zero where the block counts
      # something no municipality has none of.
      own <- values[match(joint_members, municipalities), cols[[1]]]
      missing_own <- stats::setNames(with_population[match(joint_members, municipalities)] & (is.na(own) | own <= 0), joint_members)
      split[planning_joint_split(membership, year, block, missing_own), cols] <- TRUE
      # A joint municipality's blank is covered by its holder in the areas
      # that keep the pair together.
      joint <- PLANNING_JOINT_REPORTING[match(joint_members, PLANNING_JOINT_REPORTING$municipality), ]
      active <- block != "pop" & (year <= joint$until |
        (vapply(joint$blocks, function(b) block %in% unlist(b), logical(1)) & missing_own))
      unreported[match(joint_members[active], municipalities), cols] <- FALSE
    }
    values[is.na(values)] <- 0
    # Each column is summed with the weights of its own basis: population by
    # age group, births by women of childbearing age, deaths by the expected
    # deaths of the parishes (see R/planning_parish.R). Outside the parish
    # reading every basis is the same 0/1 membership.
    bases <- vapply(columns, planning_weight_basis, character(1))
    sums <- matrix(0, length(areas), length(columns), dimnames = list(areas, columns))
    for (basis in unique(bases)) {
      pick <- bases == basis
      sums[, pick] <- planning_membership_matrix(areas, lookup, mode, basis, year) %*% values[, pick, drop = FALSE]
    }
    sums[(membership %*% unreported) > 0 | split] <- NA
    colnames(sums) <- columns

    for (area in intersect(areas, planning_published_areas)) {
      if (!area %in% table$area) next
      own <- as.numeric(table[match(area, table$area), columns, drop = FALSE])
      sums[area, !is.na(own)] <- own[!is.na(own)]
    }
    sums[, !blocks %in% present] <- NA

    tibble::as_tibble(sums) %>%
      dplyr::mutate(
        area = areas,
        year = as.integer(year),
        members = as.integer(rowSums(membership)),
        members_found = as.integer(membership %*% as.numeric(with_population)),
        .before = 1
      )
  })
  dplyr::bind_rows(frames)
}

planning_components_year <- function(area, year, lookup = get_nuts_lookup()) {
  planning_components(area, year, lookup)
}

# ---------------------------------------------------------
# Indicators
# ---------------------------------------------------------
# Exact intervals in closed form, vectorised: the Poisson one is the interval
# stats::poisson.test() reports for the rounded count, the binomial one the
# Clopper-Pearson interval of stats::binom.test().
planning_poisson_ci <- function(events, confidence = 0.95) {
  alpha <- 1 - confidence
  n <- round(events)
  lower <- ifelse(is.na(n), NA_real_, ifelse(n == 0, 0, stats::qgamma(alpha / 2, pmax(n, 1e-12))))
  upper <- ifelse(is.na(n), NA_real_, stats::qgamma(1 - alpha / 2, n + 1))
  list(lower = lower, upper = upper)
}

planning_binomial_ci <- function(successes, trials, confidence = 0.95) {
  alpha <- 1 - confidence
  n <- round(trials)
  x <- pmin(round(successes), n)
  ok <- !is.na(x) & !is.na(n) & n > 0
  lower <- rep(NA_real_, length(x))
  upper <- rep(NA_real_, length(x))
  lower[ok] <- ifelse(x[ok] == 0, 0, stats::qbeta(alpha / 2, pmax(x[ok], 1e-12), n[ok] - x[ok] + 1))
  upper[ok] <- ifelse(x[ok] == n[ok], 1, stats::qbeta(1 - alpha / 2, x[ok] + 1, pmax(n[ok] - x[ok], 1e-12)))
  list(lower = lower * 100, upper = upper * 100)
}

planning_poisson_rate <- function(events, denominator, multiplier) {
  ok <- is.finite(events) & is.finite(denominator) & denominator > 0
  ci <- planning_poisson_ci(events)
  list(
    value = ifelse(ok, events / denominator * multiplier, NA_real_),
    lower = ifelse(ok, ci$lower / denominator * multiplier, NA_real_),
    upper = ifelse(ok, ci$upper / denominator * multiplier, NA_real_)
  )
}

# Every requested indicator for every area and year of a component table, in one
# pass. Pooled indicators sum their components over the window ending in the
# year, and are missing unless every year of the window is present.
#
# Population-based indices carry no interval: the population estimates are not
# a sample of events, and INE publishes no error for them. Event rates get an
# exact Poisson interval on the count, the denominator treated as fixed.
planning_compute_indicators <- function(components, ids, undercount_years = integer(0), education_min_age = 0L) {
  keys <- paste(components$area, components$year)
  lag_of <- function(column, k) components[[column]][match(paste(components$area, components$year - k), keys)]
  pooled <- function(column, window) {
    out <- components[[column]]
    if (window > 1) for (k in seq_len(window - 1)) out <- out + lag_of(column, k)
    out
  }
  # Events over a year are divided by the mean population of that year - the
  # mean of the end-of-year estimates before and after, as INE and the workbook
  # do for crude rates, RSI and waste (the end-of-year one alone reads low
  # while the population grows: 2024 birth rate -0.5% on average). The first
  # year of the series has no previous estimate and falls back to its own.
  # Stocks counted on 31 December (pensioners) keep the end-of-year estimate.
  mean_population <- function(column) {
    now <- components[[column]]
    before <- lag_of(column, 1)
    ifelse(is.finite(before), (now + before) / 2, now)
  }
  pooled_values <- function(values, window) {
    out <- values
    if (window > 1) for (k in seq_len(window - 1)) out <- out + values[match(paste(components$area, components$year - k), keys)]
    out
  }
  ratio <- function(num, den, multiplier) ifelse(is.finite(num) & is.finite(den) & den > 0, num / den * multiplier, NA_real_)
  unstable <- function(births) !is.na(births) & births < infant_stable_births_min

  one <- function(id) {
    window <- planning_indicator_spec(id)$window
    total <- function(column) pooled(column, window)
    value <- numerator <- denominator <- lower <- upper <- rep(NA_real_, nrow(components))
    flag <- rep("", nrow(components))

    set_ratio <- function(num_col, den_col, multiplier, mean_denominator = FALSE) {
      numerator <<- total(num_col)
      denominator <<- if (mean_denominator) pooled_values(mean_population(den_col), window) else total(den_col)
      value <<- ratio(numerator, denominator, multiplier)
    }
    set_rate <- function(num_col, den_col, multiplier, mean_denominator = FALSE) {
      numerator <<- total(num_col)
      denominator <<- if (mean_denominator) pooled_values(mean_population(den_col), window) else total(den_col)
      rate <- planning_poisson_rate(numerator, denominator, multiplier)
      value <<- rate$value
      lower <<- rate$lower
      upper <<- rate$upper
    }
    set_count <- function(column, interval = TRUE) {
      value <<- total(column)
      if (interval) {
        ci <- planning_poisson_ci(value)
        lower <<- ci$lower
        upper <<- ci$upper
      }
    }
    set_share <- function(num_col, den_col) {
      numerator <<- total(num_col)
      denominator <<- total(den_col)
      value <<- ratio(numerator, denominator, 100)
      ci <- planning_binomial_ci(numerator, denominator)
      ok <- is.finite(numerator) & is.finite(denominator) & denominator > 0
      lower <<- ifelse(ok, ci$lower, NA_real_)
      upper <<- ifelse(ok, ci$upper, NA_real_)
    }

    switch(
      id,
      pop_total = set_count("pop_total", interval = FALSE),
      pct_0_14 = set_ratio("pop_0_14", "pop_total", 100),
      pct_65_plus = set_ratio("pop_65_plus", "pop_total", 100),
      pct_75_plus = set_ratio("pop_75_plus", "pop_total", 100),
      ageing_index = set_ratio("pop_65_plus", "pop_0_14", 100),
      youth_dependency = set_ratio("pop_0_14", "pop_15_64", 100),
      old_dependency = set_ratio("pop_65_plus", "pop_15_64", 100),
      births = set_count("births"),
      birth_rate = set_rate("births", "pop_total", 1000, mean_denominator = TRUE),
      deaths = set_count("deaths"),
      death_rate = set_rate("deaths", "pop_total", 1000, mean_denominator = TRUE),
      infant_rate = {
        set_rate("infant_deaths", "births", 1000)
        flag <- ifelse(unstable(denominator), "*", "")
        # Municipal under-1 counts incomplete in the source (1995-2001).
        in_window <- Reduce(`|`, lapply(seq_len(window) - 1L, function(k) (components$year - k) %in% undercount_years))
        flag <- ifelse(components$area != "Portugal" & in_window, paste0(flag, "†"), flag)
      },
      neonatal_rate = , early_neonatal_rate = , postneonatal_rate = {
        column <- switch(id, neonatal_rate = "neonatal_deaths", early_neonatal_rate = "early_neonatal_deaths",
                         postneonatal_rate = "postneonatal_deaths")
        set_rate(column, "births", 1000)
        flag <- ifelse(unstable(denominator), "*", "")
      },
      fertility_index = {
        # Sum of age-specific fertility rates over the seven five-year groups
        # 15-49, times five: births to mothers of each group per woman of that
        # group. Women are counted at mid-year, as the mean of the estimates at
        # the end of the previous year and of this one - which reproduces INE's
        # published index (Portugal 2019 1.43, 2023 1.32, 2024 1.27); the
        # end-of-year estimate alone reads up to 0.03 low.
        rates <- lapply(seq(15, 45, by = 5), function(g) {
          column <- paste0("pop_f_", g)
          now <- components[[column]]
          before <- lag_of(column, 1)
          women <- ifelse(is.finite(before), (now + before) / 2, now)
          ratio(components[[paste0("births_mage_", g)]], women, 1)
        })
        value <- Reduce(`+`, rates) * 5
      },
      teen_births_pct = set_share("births_mother_lt20", "births_mother_total"),
      low_birth_weight_pct = set_share("births_low_weight", "births_weight_known"),
      earnings_mean = set_ratio("earnings_value", "employees_total", 1),
      employees = set_count("employees_total", interval = FALSE),
      pct_employees_primary = , pct_employees_secondary = , pct_employees_tertiary = {
        set_ratio(sub("^pct_", "", id), "employees_total", 100)
        # More than 1% of the area's employees in sectors INE suppressed; see
        # planning_estimate_suppressed_sectors().
        flag <- ifelse(is.finite(value) & total("employees_estimated") > PLANNING_ESTIMATED_SHARE * denominator,
                       PLANNING_ESTIMATED_FLAG, "")
      },
      census_population = set_count("census_pop", interval = FALSE),
      census_population_change = {
        # Censuses are ten years apart; the change is against the previous one.
        now <- total("census_pop")
        before <- lag_of("census_pop", 10)
        numerator <- now - before
        denominator <- before
        value <- ratio(numerator, denominator, 100)
      },
      pct_education_none = , pct_education_basic = , pct_education_secondary = , pct_education_higher = if (education_min_age > 0) {
        # Restricted to an age: the 2011 and 2021 censuses by age group, level
        # over everyone of that age. Children, who have not had the time to
        # complete a level, no longer count as "without education".
        level <- sub("^pct_education_", "", id)
        ages <- PLANNING_EDUCATION_AGES[PLANNING_EDUCATION_AGES >= education_min_age]
        sum_ages <- function(lvl) Reduce(`+`, lapply(ages, function(a) total(paste0("census_edu_", lvl, "_", a))))
        numerator <- sum_ages(level)
        denominator <- sum_ages("total")
        value <- ratio(numerator, denominator, 100)
      } else {
        # The series counts only people with a completed level, so those with
        # none are the census population less that total.
        column <- switch(id, pct_education_basic = "census_education_basic",
                         pct_education_secondary = "census_education_secondary",
                         pct_education_higher = "census_education_higher", NA_character_)
        numerator <- if (is.na(column)) total("census_pop") - total("census_education_total") else total(column)
        denominator <- total("census_pop")
        value <- ratio(numerator, denominator, 100)
      },
      illiteracy_rate = set_ratio("census_illiterate", "census_pop_10plus", 100),
      late_fetal_rate = , perinatal_rate = {
        # Stillbirths of 28 or more weeks are the perinatal deaths less the
        # deaths under 7 days; both denominators are live births plus those
        # stillbirths, as INE defines them.
        perinatal <- total("perinatal_deaths")
        stillbirths <- pmax(perinatal - total("early_neonatal_deaths"), 0)
        births_total <- total("births")
        events <- if (identical(id, "late_fetal_rate")) stillbirths else perinatal
        rate <- planning_poisson_rate(events, births_total + stillbirths, 1000)
        numerator <- events
        denominator <- births_total + stillbirths
        value <- rate$value
        lower <- rate$lower
        upper <- rate$upper
        flag <- ifelse(unstable(births_total), "*", "")
      },
      older_births_pct = set_share("births_mother_ge35", "births_mother_total"),
      preterm_pct = set_share("births_preterm", "births_gest_known"),
      rsi_beneficiaries = set_count("rsi", interval = FALSE),
      pensioners = set_count("pensioners", interval = FALSE),
      rsi_rate = set_ratio("rsi", "pop_15_plus", 1000, mean_denominator = TRUE),
      pensioners_rate = set_ratio("pensioners", "pop_15_plus", 1000),
      pension_mean = set_ratio("pension_value", "pensioners", 1),
      purchasing_power = set_ratio("pp_share", "pp_weight", 100),
      waste_per_capita = set_ratio("waste_total", "pop_total", 1000, mean_denominator = TRUE),
      waste_selective_per_capita = set_ratio("waste_selective", "pop_total", 1000, mean_denominator = TRUE),
      stop("No computation for planning indicator ", id, call. = FALSE)
    )

    tibble::tibble(
      area = components$area, year = components$year, indicator = id,
      value = value, lower = lower, upper = upper,
      numerator = numerator, denominator = denominator, flag = flag
    )
  }

  dplyr::bind_rows(lapply(ids, one))
}

planning_indicator_table <- function(areas, years, ids = PLANNING_INDICATORS$id, lookup = get_nuts_lookup(), education_min_age = 0L,
                                     benchmark = "Portugal", mode = PLANNING_DEFAULT_SPLIT_MODE) {
  years <- as.integer(years)
  areas <- unique(as.character(areas))
  max_window <- max(PLANNING_INDICATORS$window[PLANNING_INDICATORS$id %in% ids])
  # One year earlier than the widest window: the fertility index needs the
  # previous year's population for its mid-year denominator.
  component_years <- seq.int(min(years) - max_window, max(years))
  components <- if (length(setdiff(ids, c(life_expectancy_ids, standardised_ids))) > 0) planning_components(areas, component_years, lookup, mode) else NULL

  undercount <- if (any(c("infant_rate") %in% ids)) {
    infant_undercount_years(component_years, municipalities = lookup$municipality)
  } else {
    integer(0)
  }

  life_ids <- intersect(ids, life_expectancy_ids)
  standard_ids <- intersect(ids, standardised_ids)
  other_ids <- setdiff(ids, c(life_ids, standard_ids))
  results <- list()
  if (length(other_ids) > 0) {
    results$other <- planning_compute_indicators(components, other_ids, undercount_years = undercount, education_min_age = as.integer(education_min_age)) %>%
      dplyr::filter(.data$year %in% years)
  }
  if (length(life_ids) > 0) {
    results$life <- planning_life_expectancy_table(areas, years, ids = life_ids, lookup = lookup, mode = mode) %>%
      dplyr::select(-dplyr::any_of("reason"))
  }
  if (length(standard_ids) > 0) {
    # The SMR is against the benchmark's rates, so it follows the Portugal choice.
    results$standardised <- planning_standardised_table(areas, years, ids = standard_ids, lookup = lookup, benchmark = benchmark, mode = mode)
  }
  dplyr::bind_rows(results) %>%
    dplyr::arrange(match(.data$area, areas), .data$year, match(.data$indicator, ids))
}

# Notes attached to an indicator's label with an asterisk wherever it is shown.
PLANNING_INDICATOR_NOTES <- c(
  earnings_mean = "Quadros de Pessoal (MTSSS): trabalhadores por conta de outrem, contados no local de trabalho e não no de residência; não incluem a Administração Pública nem os trabalhadores por conta própria. Um concelho com muitos empregos mas poucos residentes (sede de distrito, zona industrial) aparece acima do que os seus residentes ganham; os concelhos-dormitório à volta, abaixo.",
  life_expectancy = "Tábua de mortalidade abreviada (Chiang II, método do PHE e do Eurostat). Reproduz o Eurostat para Portugal (2017-2019: 81,9 anos na aplicação; 82,0 no Eurostat), mas fica cerca de 0,8-0,9 anos acima dos valores publicados pelo INE, que usa outra metodologia (2007). A ordenação das regiões coincide com a do INE (correlação 0,97): compare valores da aplicação entre si, não com os do INE."
)
PLANNING_INDICATOR_NOTES[["smr_all"]] <- "Óbitos observados sobre os esperados se a área tivesse as taxas por idade de Portugal (a opção escolhida: total do INE ou soma dos municípios) no mesmo triénio, vezes 100. Portugal = 100. Os óbitos por idade e causa de cada município são primeiro completados até ao total de todas as idades publicado pelo INE (\u2021 quando mais de 2% foram redistribuídos)."
PLANNING_INDICATOR_NOTES[["dsr_all"]] <- "Taxa que a área teria com a estrutura etária da População Padrão Europeia de 2013, por 100.000 habitantes; intervalo de Dobson. Os óbitos por idade de cada município são completados até ao total publicado pelo INE (\u2021 quando mais de 2% foram redistribuídos)."
PLANNING_INDICATOR_NOTES[["dsr_premature"]] <- PLANNING_INDICATOR_NOTES[["dsr_all"]]
PLANNING_INDICATOR_NOTES[["dsr_preventable"]] <- "Listas Eurostat/OCDE (2019) adaptadas à lista sucinta europeia do INE, menos detalhada: só entram as causas cuja correspondência não exige juízo clínico. Seis causas (cerca de 18% dos óbitos antes dos 75) ficam de fora, pelo que é um limite inferior da mortalidade evitável e não reproduz os valores do Eurostat. Padronizada para a População Padrão Europeia de 2013 (0-74 anos)."
PLANNING_INDICATOR_NOTES[["dsr_treatable"]] <- PLANNING_INDICATOR_NOTES[["dsr_preventable"]]
PLANNING_INDICATOR_NOTES[["avoidable_deaths"]] <- PLANNING_INDICATOR_NOTES[["dsr_preventable"]]
PLANNING_INDICATOR_NOTES[["ypll_rate"]] <- "Soma, para cada óbito antes dos 70 anos, dos anos que faltavam até aos 70 (a partir do ponto médio do grupo etário; os óbitos com menos de 1 ano contam com 69,5 anos), por 100.000 residentes com menos de 70 anos. Taxa bruta: não padronizada, depende da estrutura etária abaixo dos 70."
for (id in c("employees", "pct_employees_primary", "pct_employees_secondary", "pct_employees_tertiary")) {
  PLANNING_INDICATOR_NOTES[[id]] <- PLANNING_INDICATOR_NOTES[["earnings_mean"]]
}
for (id in c("life_expectancy_men", "life_expectancy_women", "life_expectancy_65", "life_expectancy_65_men", "life_expectancy_65_women")) {
  PLANNING_INDICATOR_NOTES[[id]] <- PLANNING_INDICATOR_NOTES[["life_expectancy"]]
}
PLANNING_EDUCATION_NOTE <- "Com idade mínima, só existem os Censos de 2011 e 2021 (o INE não publica a escolaridade por idade e município em 1991 e 2001)."

# Label as displayed: the age restriction of the education indicators, and an
# asterisk when a note applies.
planning_indicator_label <- function(id, education_min_age = 0L, ref = TRUE, star = TRUE) {
  spec <- planning_indicator_spec(id)
  label <- spec$label
  if (id %in% PLANNING_EDUCATION_IDS && isTRUE(education_min_age > 0)) label <- paste0(label, " (", education_min_age, " e mais anos)")
  if (ref) label <- paste0(label, " [", spec$ref, "]")
  if (star && id %in% names(PLANNING_INDICATOR_NOTES)) label <- paste0(label, " *")
  label
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

planning_proportional <- function(area, end_year, window = 3L, sex = "HM", lookup = get_nuts_lookup(), mode = PLANNING_DEFAULT_SPLIT_MODE) {
  years <- seq.int(as.integer(end_year) - window + 1L, as.integer(end_year))
  members <- planning_area_members(area, lookup)
  causes <- c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)

  frames <- lapply(years, planning_year_causes, sex = sex)
  if (any(vapply(frames, is.null, logical(1)))) return(tibble::tibble())
  membership <- planning_membership_matrix(area, lookup)
  weights_of <- function(year) planning_membership_matrix(area, lookup, mode, "mortality", year)
  if (any(vapply(years, function(y) planning_joint_split(membership, y, "deaths")[[1]], logical(1)))) return(tibble::tibble())

  pooled <- dplyr::bind_rows(Map(function(frame, year) {
    own_row <- area %in% planning_published_areas && area %in% frame$area
    rows <- if (own_row) frame[frame$area == area, , drop = FALSE] else frame[frame$area %in% members, , drop = FALSE]
    # A ULS that shares a municipality takes its share of that year's deaths.
    rows$weight <- if (own_row) 1 else weights_of(year)[1, rows$area]
    rows
  }, frames, years))
  counts <- vapply(causes, function(cause) sum(pooled$deaths[pooled$cause == cause] * pooled$weight[pooled$cause == cause]), numeric(1))

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
  ci <- planning_binomial_ci(deaths, total)
  c(ci$lower, ci$upper)
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
# The ULS to show together. "groups" (the default) keeps the exact unions for
# the ULS that share a municipality and leaves the individual ones out, so a
# sum over the list is the country exactly. "units" lists those ULS
# individually instead - what a ranking or a funnel compares - which under the
# whole-municipality reading counts Lisboa, Loures and Porto once per ULS.
planning_uls_units <- function(which = c("groups", "units")) {
  which <- match.arg(which)
  lookup <- get_health_lookup()
  if (is.null(lookup) || nrow(lookup) == 0) return(character(0))
  keep <- if (identical(which, "groups")) c("ULS", "ULS (grupo)") else c("ULS", "ULS (partilhada)")
  sort(unique(as.character(lookup$unit[lookup$kind %in% keep])))
}

# ---------------------------------------------------------
# Levels and comparators
# ---------------------------------------------------------
PLANNING_LEVELS <- c("Local", "ULS", "ARS", "NUTS III", "NUTS II", "NUTS I", "Portugal")

# Fixed colour per level (the reference categorical palette, in its validated
# order), so a level keeps its colour whichever comparators are shown.
PLANNING_LEVEL_COLOURS <- c(
  "Local" = "#2a78d6", "ULS" = "#eb6834", "ARS" = "#1baf7a", "NUTS III" = "#eda100",
  "NUTS II" = "#e87ba4", "NUTS I" = "#008300", "Portugal" = "#4a3aa7"
)

# Every area the app can build, with its level, highest level first. A label
# naming several levels at once (the islands are NUTS I, II and III) keeps the
# highest.
planning_area_levels <- function(lookup = get_nuts_lookup(), health = get_health_lookup()) {
  nuts <- function(column, level) {
    if (!column %in% names(lookup)) return(tibble::tibble(area = character(0), level = character(0)))
    tibble::tibble(area = unique(stats::na.omit(as.character(lookup[[column]]))), level = level)
  }
  health_units <- function(kinds, level) {
    if (is.null(health) || nrow(health) == 0) return(tibble::tibble(area = character(0), level = character(0)))
    units <- health_unit_choices(health)
    tibble::tibble(area = units[units %in% health$unit[health$kind %in% kinds]], level = level)
  }
  dplyr::bind_rows(
    tibble::tibble(area = "Portugal", level = "Portugal"),
    nuts("nuts1", "NUTS I"),
    nuts("nuts2", "NUTS II"),
    nuts("nuts3", "NUTS III"),
    health_units("ARS", "ARS"),
    health_units(c("ULS", "ULS (grupo)"), "ULS"),
    tibble::tibble(area = sort(unique(as.character(lookup$municipality))), level = "Município")
  ) %>%
    dplyr::distinct(area, .keep_all = TRUE)
}

# A ULS can sit inside a NUTS III (ULS Matosinhos in the Área Metropolitana do
# Porto), and a NUTS III inside an ARS, so the ranks interleave the two trees.
planning_level_rank <- c("Município" = 1, "ULS" = 1.5, "NUTS III" = 2, "ARS" = 3, "NUTS II" = 3, "NUTS I" = 4, "Portugal" = 5)

# The areas that contain `area`, one per level above it, nearest first: for a
# municipality its ULS, ARS, NUTS III, NUTS II, NUTS I and Portugal. A containing
# unit with exactly the same municipalities is left out (ULS Algarve is also
# NUTS II and NUTS III Algarve), since it would repeat the area's own values.
planning_comparators <- function(area, lookup = get_nuts_lookup(), health = get_health_lookup()) {
  levels <- planning_area_levels(lookup, health)
  own_level <- levels$level[match(area, levels$area)]
  if (is.na(own_level) || identical(own_level, "Portugal")) {
    return(tibble::tibble(level = character(0), area = character(0)))
  }
  members <- planning_area_members(area, lookup)
  candidates <- levels[levels$level != "Município" &
                         planning_level_rank[levels$level] > planning_level_rank[[own_level]], , drop = FALSE]

  keep <- vapply(seq_len(nrow(candidates)), function(i) {
    if (identical(candidates$level[[i]], "Portugal")) return(TRUE)
    unit <- planning_area_members(candidates$area[[i]], lookup)
    all(members %in% unit) && !setequal(members, unit)
  }, logical(1))

  found <- candidates[keep, , drop = FALSE]
  found <- found[order(match(found$level, PLANNING_LEVELS)), , drop = FALSE]
  # One per level; drop a unit whose municipalities equal a nearer comparator's.
  found <- found[!duplicated(found$level), , drop = FALSE]
  sets <- lapply(found$area, planning_area_members, lookup = lookup)
  # Portugal always stays: it reads INE's national row, which also counts
  # residents of unknown municipality.
  repeated <- vapply(seq_along(sets), function(i) {
    i > 1 && !identical(found$level[[i]], "Portugal") &&
      any(vapply(seq_len(i - 1), function(j) setequal(sets[[i]], sets[[j]]), logical(1)))
  }, logical(1))
  tibble::tibble(level = found$level[!repeated], area = found$area[!repeated])
}

# ---------------------------------------------------------
# Proportional mortality for many areas
# ---------------------------------------------------------
# The same figures as planning_proportional(), for every area and triennium in
# one pass: the cause counts are summed with the membership matrix.
planning_proportional_table <- function(areas, end_years, window = 3L, sex = "HM", lookup = get_nuts_lookup(), mode = PLANNING_DEFAULT_SPLIT_MODE) {
  # 0/1 here: each year's deaths are weighted in the loop below.
  membership <- planning_membership_matrix(areas, lookup)
  areas <- rownames(membership)
  municipalities <- colnames(membership)
  causes <- c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)
  needed <- sort(unique(unlist(lapply(as.integer(end_years), function(y) seq.int(y - window + 1L, y)))))

  per_year <- lapply(needed, function(year) {
    frame <- planning_year_causes(year, sex = sex)
    if (is.null(frame)) return(NULL)
    wide <- matrix(0, length(municipalities), length(causes), dimnames = list(municipalities, causes))
    hit <- frame[frame$area %in% municipalities & frame$cause %in% causes, , drop = FALSE]
    wide[cbind(match(hit$area, municipalities), match(hit$cause, causes))] <- hit$deaths
    sums <- planning_membership_matrix(areas, lookup, mode, "mortality", year) %*% wide
    for (area in intersect(areas, planning_published_areas)) {
      own <- frame[frame$area == area & frame$cause %in% causes, , drop = FALSE]
      if (nrow(own) > 0) {
        sums[area, ] <- 0
        sums[area, match(own$cause, causes)] <- own$deaths
      }
    }
    sums
  })
  names(per_year) <- needed

  dplyr::bind_rows(lapply(as.integer(end_years), function(end_year) {
    years <- as.character(seq.int(end_year - window + 1L, end_year))
    if (!all(years %in% names(per_year)) || any(vapply(per_year[years], is.null, logical(1)))) return(NULL)
    pooled <- Reduce(`+`, per_year[years])
    total <- pooled[, planning_all_causes]
    groups <- pooled[, PLANNING_CAUSE_GROUPS$cause, drop = FALSE]
    other <- pmax(total - rowSums(groups), 0)
    counts <- cbind(total, groups, other)
    codes <- c("C00", PLANNING_CAUSE_GROUPS$code, "Outras")
    labels <- c("Todas as causas de morte", PLANNING_CAUSE_GROUPS$label, "Restantes causas")
    long <- tibble::tibble(
      area = rep(areas, times = length(codes)),
      code = rep(codes, each = length(areas)),
      group = rep(labels, each = length(areas)),
      deaths = as.numeric(counts),
      all_deaths = rep(unname(total), times = length(codes))
    )
    # Deaths up to 1998 sit with the parent of Odivelas, Trofa and Vizela.
    joint <- Reduce(`|`, lapply(as.integer(years), function(y) planning_joint_split(membership, y, "deaths")))
    long$deaths[rep(joint, times = length(codes))] <- NA_real_
    long$all_deaths[rep(joint, times = length(codes))] <- NA_real_
    ci <- planning_binomial_ci(long$deaths, long$all_deaths)
    long %>%
      dplyr::mutate(
        period = paste0(min(years), "-", max(years)),
        end_year = end_year,
        share = ifelse(!is.na(.data$all_deaths) & .data$all_deaths > 0, .data$deaths / .data$all_deaths * 100, NA_real_),
        lower = ifelse(is.na(.data$all_deaths), NA_real_, ifelse(.data$code == "C00", 100, ci$lower)),
        upper = ifelse(is.na(.data$all_deaths), NA_real_, ifelse(.data$code == "C00", 100, ci$upper))
      ) %>%
      dplyr::select(-all_deaths)
  })) %>%
    dplyr::arrange(match(.data$area, areas), .data$end_year, match(.data$code, c("C00", PLANNING_CAUSE_GROUPS$code, "Outras")))
}
