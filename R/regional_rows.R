# =========================================================
# Regional deaths: INE's own rows, or the sum of municipalities
# =========================================================
# Regions are expanded into municipalities so that one membership list can be
# applied to every year (R/regions.R). For deaths by cause that sum is short.
# INE publishes complete municipal totals but incomplete municipal age
# breakdowns, concentrated where counts are small, so a region rebuilt from its
# municipalities' age bands under-counts - lung cancer, 2013: Açores -18.4%,
# Madeira -11.8%, Alentejo -8.9%, Norte -1.2% - and in 2014 the loss is
# 31-84% in every region. Every year is affected, not only 2014.
#
# INE's regional rows do not have the defect: their age bands add up to their
# totals exactly. So two sources are offered, as a single app-wide choice:
#
#   "ine_rows"       INE's regional row wherever one exists for the selected
#                    territory in that year; the municipal sum only for the
#                    years without one. The default.
#   "municipal_sum"  Always the sum of municipalities. One consistent source and
#                    no seams, but biased low for cause-specific figures.
#
# Which INE rows build a region, by year:
#
#   Continente, Norte, Algarve, Açores and Madeira are the same territory under
#   NUTS 2013 and NUTS 2024, so their own rows serve every year under either
#   vintage (0008206 to 2021, 0013166 from 2022, the death archive's precedence).
#
#   The redrawn regions are composed from subregion rows in the years their own
#   row does not exist. Oeste, Médio Tejo and Lezíria do Tejo were NUTS III units
#   under both definitions, which is what makes this possible. Every composition
#   below was checked against the two municipality lookups, not assumed:
#
#     NUTS 2013, 2023-2024   Centro   = Centro + Oeste + Médio Tejo (2024)
#                            Alentejo = Alentejo + Lezíria do Tejo (2024)
#                            AML      = Grande Lisboa + Península de Setúbal
#     NUTS 2024, 1991-2021   Alentejo = its four subregions under NUTS 2013
#                            Centro   = six NUTS-2013 subregions
#                                       + Sertã + Vila de Rei
#                            OVT      = Oeste + Médio Tejo + Lezíria (2013)
#                                       - Sertã - Vila de Rei
#
#   Sertã and Vila de Rei moved from Médio Tejo to Beira Baixa in 2024 and are
#   the only reason two compositions are not pure subregion sums. Their own
#   municipal rows make the correction - about 1% of Centro - so those two
#   regions carry a trace of the municipal bias the rest avoid.
#
#   Grande Lisboa and Península de Setúbal before 2022 cannot be composed:
#   under NUTS 2013 the Lisbon metropolitan area was a single subregion. They
#   fall back to the municipal sum, and step where the source changes.
#   region_source_seam_warning() says so.
#
# Only deaths are substituted. Population is always the sum of municipalities,
# which is complete - the defect is in the age breakdown of deaths, not in the
# denominator.
#
# A municipality selected on its own cannot be repaired this way: there is no
# finer row to take. municipal_age_detail_warning() flags it instead.
#
# Rows are keyed by INE geography code, never by label: Algarve, the Lisbon
# metropolitan area and both autonomous regions carry the same name at two or
# three NUTS levels, and the municipal archive's label-keyed regional rows are
# multi-counted for exactly that reason.

region_source_choices <- c(
  "Linhas regionais do INE (recomendado)" = "ine_rows",
  "Soma dos municípios" = "municipal_sum"
)

normalize_region_source <- function(source) {
  source <- as.character(source)[1]
  if (isTRUE(source %in% unname(region_source_choices))) source else "ine_rows"
}

default_region_source <- function() {
  normalize_region_source(Sys.getenv("MORTALITY_REGION_SOURCE", unset = "ine_rows"))
}

# The INE rows that build each region, by the app's label and vintage.
# `vintage = NA` marks a territory identical under both. `codes` are summed;
# `plus` and `minus` are municipalities whose own rows are added or taken away,
# for the two compositions a municipality move keeps from being exact.
.moved_2024 <- "Sertã,Vila de Rei"

REGIONAL_ROW_TERRITORIES <- tibble::tribble(
  ~region,                        ~vintage, ~indicator, ~codes,                         ~plus,        ~minus,       ~from, ~to,
  # Same territory under both vintages.
  "Continente",                   NA,       "0008206",  "1",                            "",           "",           1991L, 2021L,
  "Continente",                   NA,       "0013166",  "1",                            "",           "",           2022L, 2024L,
  "Norte",                        NA,       "0008206",  "11",                           "",           "",           1991L, 2021L,
  "Norte",                        NA,       "0013166",  "11",                           "",           "",           2022L, 2024L,
  "Algarve",                      NA,       "0008206",  "15",                           "",           "",           1991L, 2021L,
  "Algarve",                      NA,       "0013166",  "15",                           "",           "",           2022L, 2024L,
  "Região Autónoma dos Açores",   NA,       "0008206",  "20",                           "",           "",           1991L, 2021L,
  "Região Autónoma dos Açores",   NA,       "0013166",  "20",                           "",           "",           2022L, 2024L,
  "Região Autónoma da Madeira",   NA,       "0008206",  "30",                           "",           "",           1991L, 2021L,
  "Região Autónoma da Madeira",   NA,       "0013166",  "30",                           "",           "",           2022L, 2024L,
  # NUTS 2013: own rows to 2022, composed from NUTS-2024 rows after.
  "Centro",                       "2013",   "0008206",  "16",                           "",           "",           1991L, 2022L,
  "Centro",                       "2013",   "0013166",  "19,1D1,1D2",                   "",           "",           2023L, 2024L,
  "Área Metropolitana de Lisboa", "2013",   "0008206",  "17",                           "",           "",           1991L, 2022L,
  "Área Metropolitana de Lisboa", "2013",   "0013166",  "1A,1B",                        "",           "",           2023L, 2024L,
  "Alentejo",                     "2013",   "0008206",  "18",                           "",           "",           1991L, 2022L,
  "Alentejo",                     "2013",   "0013166",  "1C,1D3",                       "",           "",           2023L, 2024L,
  # NUTS 2024: own rows from 2022, composed from NUTS-2013 subregions before.
  "Centro",                       "2024",   "0008206",  "16D,16E,16F,16G,16H,16J",      .moved_2024,  "",           1991L, 2021L,
  "Centro",                       "2024",   "0013166",  "19",                           "",           "",           2022L, 2024L,
  "Alentejo",                     "2024",   "0008206",  "181,184,186,187",              "",           "",           1991L, 2021L,
  "Alentejo",                     "2024",   "0013166",  "1C",                           "",           "",           2022L, 2024L,
  "Oeste e Vale do Tejo",         "2024",   "0008206",  "16B,16I,185",                  "",           .moved_2024,  1991L, 2021L,
  "Oeste e Vale do Tejo",         "2024",   "0013166",  "1D",                           "",           "",           2022L, 2024L,
  # No NUTS-2013 split of the Lisbon metropolitan area exists before 2022.
  "Grande Lisboa",                "2024",   "0013166",  "1A",                           "",           "",           2022L, 2024L,
  "Península de Setúbal",         "2024",   "0013166",  "1B",                           "",           "",           2022L, 2024L
)

split_list <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# Year-by-year plan for one region under one vintage: which indicator and code
# supply its row. Years absent from the result fall back to the municipal sum.
regional_row_plan <- function(region, vintage, years) {
  rows <- REGIONAL_ROW_TERRITORIES[
    REGIONAL_ROW_TERRITORIES$region == region &
      (is.na(REGIONAL_ROW_TERRITORIES$vintage) | REGIONAL_ROW_TERRITORIES$vintage == as.character(vintage)),
    ,
    drop = FALSE
  ]

  years <- sort(unique(as.integer(years)))
  empty <- tibble::tibble(
    year = integer(0), indicator = character(0),
    codes = character(0), plus = character(0), minus = character(0)
  )
  if (nrow(rows) == 0 || length(years) == 0) {
    return(empty)
  }

  plan <- lapply(years, function(year) {
    hit <- rows[rows$from <= year & rows$to >= year, , drop = FALSE]
    if (nrow(hit) == 0) return(NULL)
    tibble::tibble(
      year = year, indicator = hit$indicator[[1]],
      codes = hit$codes[[1]], plus = hit$plus[[1]], minus = hit$minus[[1]]
    )
  })

  dplyr::bind_rows(plan)
}

regional_rows_root <- function() {
  # Same resolution as the other parallel datasets, so tests can point it at the
  # repository copy without loading the snapshot machinery.
  infant_snapshot_root()
}

.regional_rows_cache <- new.env(parent = emptyenv())

read_regional_rows <- function(indicator, year) {
  key <- paste(regional_rows_root(), indicator, year)
  if (!is.null(.regional_rows_cache[[key]])) {
    return(.regional_rows_cache[[key]])
  }

  path <- file.path(regional_rows_root(), "regional_deaths", indicator, paste0("year_", year, ".rds"))
  rows <- if (file.exists(path)) readRDS(path) else NULL
  .regional_rows_cache[[key]] <- rows
  rows
}

# Replace each expanded region's municipal deaths with its INE row, for every
# year a row exists. The frame keeps its shape - one row per (year, area, sex,
# cause, band) - so every downstream calculation is unchanged: the region
# simply appears as one area in the substituted years and as its municipalities
# in the others.
#
# Nested selections are handled by covering each municipality once: regions are
# taken largest first, and a region whose municipalities are already covered in
# that year is skipped. Continente with Norte substitutes Continente only.
#
# Returns the frame with an attribute `regional_substitutions`, a table of the
# (region, year, indicator) actually applied, for source reporting and warnings.
# `load_municipal_deaths(areas, years, causes)` supplies rows for a `minus`
# municipality, which by construction is not a member of the region and so is
# not in `df`. Without it those compositions fall back to the municipal sum.
substitute_regional_deaths <- function(df,
                                       expanded_regions,
                                       vintage,
                                       lookup,
                                       source = "ine_rows",
                                       load_municipal_deaths = NULL) {
  applied <- tibble::tibble(region = character(0), year = integer(0), indicator = character(0))
  attr(df, "regional_substitutions") <- applied

  if (!identical(normalize_region_source(source), "ine_rows") ||
      length(expanded_regions) == 0 || nrow(df) == 0) {
    return(df)
  }

  has_pop <- "pop" %in% names(df)
  was_factor <- is.factor(df$age_band)
  df$age_band <- as.character(df$age_band)

  members <- lapply(expanded_regions, region_municipalities, lookup = lookup)
  names(members) <- expanded_regions
  regions_by_size <- expanded_regions[order(-lengths(members))]

  causes <- unique(as.character(df$cause))
  sexes <- unique(as.character(df$sex))

  for (year in sort(unique(as.integer(df$year)))) {
    covered <- character(0)

    for (region in regions_by_size) {
      region_members <- members[[region]]
      if (length(region_members) == 0 || all(region_members %in% covered)) next

      plan <- regional_row_plan(region, vintage, year)
      if (nrow(plan) == 0) next

      rows <- read_regional_rows(plan$indicator[[1]], year)
      if (is.null(rows)) next

      codes <- split_list(plan$codes[[1]])
      plus <- split_list(plan$plus[[1]])
      minus <- split_list(plan$minus[[1]])

      # A composition needs every one of its subregion rows; a partial set
      # would silently drop part of the territory.
      if (!all(codes %in% rows$region_code)) next

      key <- c("year", "sex", "cause", "age_band")
      regional <- rows[
        rows$region_code %in% codes & rows$cause %in% causes & rows$sex %in% sexes,
        c(key, "deaths"),
        drop = FALSE
      ] %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
        dplyr::summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop")
      if (nrow(regional) == 0) next

      municipal_part <- function(areas, frame) {
        frame[frame$year == year & frame$area %in% areas, c(key, "deaths"), drop = FALSE] %>%
          dplyr::mutate(age_band = as.character(.data$age_band)) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
          dplyr::summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop")
      }

      if (length(plus) > 0) {
        added <- municipal_part(plus, df)
        regional <- dplyr::bind_rows(regional, added) %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
          dplyr::summarise(deaths = sum(.data$deaths), .groups = "drop")
      }

      if (length(minus) > 0) {
        if (is.null(load_municipal_deaths)) next
        extra <- load_municipal_deaths(minus, year, causes)
        if (is.null(extra) || nrow(extra) == 0) next
        extra <- extra[extra$sex %in% sexes, , drop = FALSE]
        taken <- municipal_part(minus, extra) %>% dplyr::rename(minus_deaths = "deaths")
        regional <- dplyr::left_join(regional, taken, by = key) %>%
          dplyr::mutate(
            # Municipal rows under-count, so subtracting them can only leave the
            # composition slightly high, never negative - but guard anyway.
            deaths = pmax(.data$deaths - dplyr::coalesce(.data$minus_deaths, 0), 0)
          ) %>%
          dplyr::select(-"minus_deaths")
      }

      in_year <- df$year == year & df$area %in% region_members
      replacement <- tibble::as_tibble(regional) %>%
        dplyr::mutate(
          area = region,
          death_source = paste0(
            plan$indicator[[1]],
            if (length(codes) > 1 || length(plus) > 0 || length(minus) > 0) " (linhas regionais compostas)" else " (linha regional)"
          )
        )

      if (has_pop) {
        # Population repeats once per cause in a multi-cause frame; take it once
        # per (area, sex, band) before summing over the region.
        pop <- df[in_year, , drop = FALSE] %>%
          dplyr::distinct(.data$area, .data$sex, .data$age_band, .keep_all = TRUE) %>%
          dplyr::group_by(.data$sex, .data$age_band) %>%
          dplyr::summarise(
            pop = sum(.data$pop, na.rm = TRUE),
            population_source = paste(unique(.data$population_source), collapse = ", "),
            .groups = "drop"
          )
        replacement <- dplyr::inner_join(replacement, pop, by = c("sex", "age_band"))
      }

      keep_cols <- intersect(names(df), names(replacement))
      df <- dplyr::bind_rows(df[!in_year, , drop = FALSE], replacement[, keep_cols, drop = FALSE])

      covered <- c(covered, region_members)
      applied <- dplyr::bind_rows(applied, tibble::tibble(region = region, year = year, indicator = plan$indicator[[1]]))
    }
  }

  if (was_factor) {
    df$age_band <- factor(df$age_band, levels = age_levels, ordered = TRUE)
  }

  df <- df %>% dplyr::arrange(.data$year, .data$area, .data$sex, .data$age_band)
  attr(df, "regional_substitutions") <- applied
  df
}

# A region whose series mixes INE rows and municipal sums steps where the source
# changes, by roughly that region's under-count. Only the redrawn regions can do
# this; the five stable territories have rows for every year.
region_source_seam_warning <- function(expanded_regions, vintage, years, source = "ine_rows") {
  if (!identical(normalize_region_source(source), "ine_rows") || length(expanded_regions) == 0) {
    return(NULL)
  }

  years <- sort(unique(as.integer(years)))
  if (length(years) < 2) return(NULL)

  mixed <- vapply(expanded_regions, function(region) {
    covered <- regional_row_plan(region, vintage, years)$year
    length(covered) > 0 && length(covered) < length(years)
  }, logical(1))

  if (!any(mixed)) return(NULL)

  described <- vapply(expanded_regions[mixed], function(region) {
    covered <- regional_row_plan(region, vintage, years)$year
    as.character(glue::glue("{region} ({min(covered)}-{max(covered)})"))
  }, character(1))

  as.character(glue::glue(
    "Atenção: {paste(described, collapse = '; ')} usa as linhas regionais do INE ",
    "apenas nos anos indicados: antes de 2022 a Área Metropolitana de Lisboa era ",
    "uma única sub-região e não pode ser dividida. Nos restantes anos usa a soma ",
    "dos municípios, que subestima os óbitos por causa, pelo que a série pode dar ",
    "um salto onde a fonte muda, sem mudança real da mortalidade."
  ))
}

# A municipality on its own has no finer row to fall back on. Its cause-specific
# age bands can be incomplete in any year, and in 2014 they are badly so.
municipal_age_detail_warning <- function(areas, lookup, causes, years) {
  areas <- unique(as.character(areas))
  municipalities <- areas[!areas %in% c("Portugal", get_known_regions(lookup))]
  if (length(municipalities) == 0) return(NULL)

  cause_specific <- any(as.character(causes) != "Todas as causas de morte")
  includes_2014 <- 2014L %in% as.integer(years)
  if (!cause_specific && !includes_2014) return(NULL)

  parts <- character(0)
  if (cause_specific) {
    parts <- c(parts, paste0(
      "o INE publica os óbitos por causa de cada município com a repartição por ",
      "idade incompleta, sobretudo onde os números são pequenos, pelo que os ",
      "valores municipais podem estar subestimados"
    ))
  }
  if (includes_2014) {
    parts <- c(parts, paste0(
      "em 2014 essa perda é muito maior (cerca de 30% dos óbitos municipais por ",
      "cancro do pulmão sem idade), pelo que esse ano não é fiável ao nível municipal"
    ))
  }

  as.character(glue::glue(
    "Atenção: {paste(parts, collapse = '; e ')}. As regiões não têm este problema ",
    "quando usam as linhas regionais do INE; um município isolado não pode ser ",
    "corrigido desta forma."
  ))
}
