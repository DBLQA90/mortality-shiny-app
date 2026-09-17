# =========================================================
# Planning indicators (the "Indicadores de Planeamento" tab)
# =========================================================
# Demographic and mortality indicators that local health plans report for every
# ULS, ARS, NUTS region and municipality: the population structure, dependency
# and ageing indices, births, deaths, crude rates, infant mortality and
# proportional mortality by large cause groups. They mirror the indicators of
# the DRS/PNS2030 support workbook (I1-I8, I37-I39, I45), recomputed from the
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
  ~id,               ~label,                                        ~ref,  ~unit,               ~window, ~digits,
  "pop_total",       "População residente (estimativa)",            "I1",  "N.º",               1L,      0L,
  "pct_0_14",        "Proporção de jovens (0-14 anos)",             "I1",  "%",                 1L,      1L,
  "pct_65_plus",     "Proporção de idosos (65 e mais anos)",        "I1",  "%",                 1L,      1L,
  "pct_75_plus",     "Proporção de 75 e mais anos",                 "I1",  "%",                 1L,      1L,
  "ageing_index",    "Índice de envelhecimento",                    "I4",  "por 100 jovens",    1L,      1L,
  "youth_dependency","Índice de dependência de jovens",             "I5",  "por 100 em idade activa", 1L, 1L,
  "old_dependency",  "Índice de dependência de idosos",             "I6",  "por 100 em idade activa", 1L, 1L,
  "births",          "Nados-vivos",                                 "I7",  "N.º",               1L,      0L,
  "birth_rate",      "Taxa bruta de natalidade",                    "I8",  "‰",                 1L,      1L,
  "deaths",          "Óbitos",                                      "I37", "N.º",               1L,      0L,
  "death_rate",      "Taxa bruta de mortalidade",                   "I38", "‰",                 1L,      1L,
  "infant_rate",     "Taxa de mortalidade infantil (triénio)",      "I39", "‰ nados-vivos",     3L,      1L
)

planning_indicator_choices <- function() {
  stats::setNames(
    PLANNING_INDICATORS$id,
    paste0(PLANNING_INDICATORS$label, " (", PLANNING_INDICATORS$ref, ")")
  )
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
  "pop_total", "pop_0_14", "pop_15_64", "pop_65_plus", "pop_75_plus",
  "births", "deaths", "infant_deaths"
)

planning_year_components <- function(year) {
  key <- paste(infant_snapshot_root(), "components", year, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) {
    return(get(key, envir = planning_cache, inherits = FALSE))
  }

  parts <- list()

  pop <- read_year_file("population", year)
  if (!is.null(pop)) {
    pop <- pop[pop$sex == "HM", , drop = FALSE]
    lower <- planning_age_lower(pop$age_band)
    parts$pop <- tibble::tibble(area = pop$area, pop = pop$pop, lower = lower) %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(
        pop_total = sum(pop, na.rm = TRUE),
        pop_0_14 = sum(pop[lower < 15], na.rm = TRUE),
        pop_15_64 = sum(pop[lower >= 15 & lower < 65], na.rm = TRUE),
        pop_65_plus = sum(pop[lower >= 65], na.rm = TRUE),
        pop_75_plus = sum(pop[lower >= 75], na.rm = TRUE),
        .groups = "drop"
      )
  }

  births <- read_year_file("births", year)
  if (!is.null(births)) {
    parts$births <- births %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(births = sum(births, na.rm = TRUE), .groups = "drop")
  }

  totals <- read_death_totals_year(year)
  if (!is.null(totals)) {
    parts$deaths <- totals %>%
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
    parts$infant <- infant %>%
      dplyr::filter(.data$sex == "HM") %>%
      dplyr::group_by(area) %>%
      dplyr::summarise(infant_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  }

  table <- if (length(parts) == 0) {
    tibble::tibble(area = character(0))
  } else {
    Reduce(function(x, y) dplyr::full_join(x, y, by = "area"), parts)
  }
  for (column in planning_component_columns) {
    if (!column %in% names(table)) table[[column]] <- NA_real_
  }
  # Which datasets exist for the year at all: a dataset whose file is missing
  # stays NA for every area, while an area merely absent from a present file
  # counts as zero (INE omits empty cells).
  attr(table, "present") <- c(
    pop = !is.null(pop), births = !is.null(births), deaths = !is.null(totals), infant = !is.null(infant)
  )

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

  published <- area %in% planning_published_areas
  value_of <- function(column, dataset) {
    if (is.null(present) || !isTRUE(present[[dataset]])) return(NA_real_)
    values <- table[[column]]
    if (published && area %in% table$area && !is.na(values[table$area == area][[1]])) {
      return(values[table$area == area][[1]])
    }
    sum(values[table$area %in% members], na.rm = TRUE)
  }

  tibble::tibble(
    area = area,
    year = as.integer(year),
    members = length(members),
    members_found = sum(members %in% table$area[!is.na(table$pop_total)]),
    pop_total = value_of("pop_total", "pop"),
    pop_0_14 = value_of("pop_0_14", "pop"),
    pop_15_64 = value_of("pop_15_64", "pop"),
    pop_65_plus = value_of("pop_65_plus", "pop"),
    pop_75_plus = value_of("pop_75_plus", "pop"),
    births = value_of("births", "births"),
    deaths = value_of("deaths", "deaths"),
    infant_deaths = value_of("infant_deaths", "infant")
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
    }
  )

  out
}

planning_indicator_table <- function(areas, years, ids = PLANNING_INDICATORS$id, lookup = get_nuts_lookup()) {
  years <- as.integer(years)
  max_window <- max(PLANNING_INDICATORS$window[PLANNING_INDICATORS$id %in% ids])
  component_years <- seq.int(min(years) - max_window + 1L, max(years))
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
      Indicador = paste0(.data$label, " [", .data$ref, "]"),
      Unidade = .data$unit,
      order = match(.data$indicator, specs$id)
    ) %>%
    dplyr::select(order, Indicador, Unidade, Período, area, cell) %>%
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
