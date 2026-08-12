# =========================================================
# Regional aggregates and NUTS vintage handling
# =========================================================
# INE publishes mortality and population under three different NUTS vintages:
#
#   population 0003182  NUTS-2002   1991-2013
#   population 0008273  NUTS-2013   2011-2023
#   deaths     0008206  NUTS-2013   1980-2022
#   deaths     0013166  NUTS-2024   2022-2024
#
# Municipality boundaries are identical across all three, but the regional
# groupings are not: NUTS-2024 moved Lezíria do Tejo out of Alentejo into the
# new "Oeste e Vale do Tejo". Reading INE's own regional rows across the 2022
# seam therefore compares two different Alentejos - 2022 all-cause deaths are
# 11,327 under 0008206 and 7,898 under 0013166 - and dividing the later figure
# by NUTS-2013 population understates regional mortality by about 30%.
#
# Two honest responses, both offered to the user (see `region_mode_choices`):
#
#   "municipal"  Rebuild the region by summing its municipalities, using one
#                fixed NUTS-2024 membership list for every year. The series is
#                continuous and means "this region as defined today". It is an
#                approximation only where a municipality cannot be matched by
#                name (see below), not where boundaries moved.
#
#   "original"   Keep INE's own regional rows and leave the 2022 change of
#                definition visible, with a warning naming the affected region.
#
# Regions are unions of whole municipalities, so no parish-level data is needed.

nuts_lookup_path <- function() {
  file.path(app_dir_or_wd(), "data", "nuts_lookup.rds")
}

app_dir_or_wd <- function() {
  if (exists("app_dir", envir = globalenv(), inherits = TRUE)) {
    get("app_dir", envir = globalenv())
  } else if (exists("app_dir", inherits = TRUE)) {
    get("app_dir", inherits = TRUE)
  } else {
    getwd()
  }
}

# Cached because every area resolution consults it.
.nuts_lookup_cache <- new.env(parent = emptyenv())

get_nuts_lookup <- function(path = nuts_lookup_path()) {
  key <- normalizePath(path, mustWork = FALSE)

  if (!is.null(.nuts_lookup_cache[[key]])) {
    return(.nuts_lookup_cache[[key]])
  }

  lookup <- if (file.exists(path)) {
    readRDS(path)
  } else {
    tibble::tibble(
      municipality = character(0),
      municipality_code = character(0),
      nuts2 = character(0),
      nuts3 = character(0)
    )
  }

  .nuts_lookup_cache[[key]] <- lookup
  lookup
}

# Region labels the lookup can rebuild from municipalities.
get_known_regions <- function(lookup = get_nuts_lookup()) {
  sort(unique(stats::na.omit(c(lookup$nuts2, lookup$nuts3))))
}

is_region_label <- function(area, lookup = get_nuts_lookup()) {
  as.character(area) %in% get_known_regions(lookup)
}

# Municipalities making up a region, restricted to those the data actually
# knows about. `available_areas` is the area vocabulary of the active data
# source; passing it lets the caller distinguish "region not covered" from
# "region covered except for two islands".
region_municipalities <- function(region, lookup = get_nuts_lookup(), available_areas = NULL) {
  region <- as.character(region)

  members <- lookup %>%
    dplyr::filter(.data$nuts2 %in% region | .data$nuts3 %in% region) %>%
    dplyr::pull(.data$municipality) %>%
    unique()

  if (is.null(available_areas)) {
    return(sort(members))
  }

  sort(intersect(members, as.character(available_areas)))
}

# What a municipal rebuild of `region` would miss, given the available areas.
# Used to warn instead of silently returning an undercount: the Azores and
# Madeira both contain a "Calheta", and INE labels them plainly in some
# indicators and with a region suffix in others, so those municipalities cannot
# always be matched by name.
region_coverage <- function(region, lookup = get_nuts_lookup(), available_areas = NULL) {
  expected <- lookup %>%
    dplyr::filter(.data$nuts2 %in% region | .data$nuts3 %in% region) %>%
    dplyr::pull(.data$municipality) %>%
    unique()

  found <- if (is.null(available_areas)) expected else intersect(expected, as.character(available_areas))
  missing <- setdiff(expected, found)

  list(
    region = region,
    expected = length(expected),
    found = length(found),
    missing = sort(missing),
    complete = length(missing) == 0
  )
}

region_coverage_warning <- function(coverage) {
  if (isTRUE(coverage$complete) || coverage$expected == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "A reconstrução de {coverage$region} por municípios usa {coverage$found} de ",
    "{coverage$expected} municípios. Em falta: {paste(coverage$missing, collapse = ', ')}. ",
    "Os valores apresentados subestimam a região."
  ))
}

# Multiple selected areas are summed into one combined geography, which is only
# meaningful when they do not overlap. Selecting Portugal alongside anything, or
# a region alongside one of its own municipalities, double-counts those deaths
# and that population. Distinct NUTS II regions are disjoint, so combining them
# is legitimate and not flagged.
overlapping_selection_warning <- function(areas, lookup = get_nuts_lookup()) {
  areas <- unique(as.character(areas))

  if (length(areas) < 2) {
    return(NULL)
  }

  problems <- character(0)

  if ("Portugal" %in% areas) {
    problems <- c(problems, "Portugal já inclui todas as outras áreas seleccionadas")
  }

  regions <- setdiff(areas[is_region_label(areas, lookup)], "Portugal")
  for (region in regions) {
    inside <- intersect(setdiff(areas, region), region_municipalities(region, lookup))
    if (length(inside) > 0) {
      problems <- c(problems, as.character(glue::glue(
        "{region} já inclui {paste(inside, collapse = ', ')}"
      )))
    }
  }

  if (length(problems) == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "Selecção sobreposta: {paste(problems, collapse = '; ')}. As áreas ",
    "seleccionadas são somadas, pelo que estes óbitos e esta população seriam ",
    "contados duas vezes."
  ))
}

normalize_region_mode <- function(region_mode) {
  region_mode <- as.character(region_mode)[[1]]

  if (isTRUE(region_mode %in% unname(region_mode_choices))) {
    region_mode
  } else {
    "original"
  }
}

# Swap any selected region for its municipalities when the municipal mode is
# active. The rest of the app already sums multiple selected areas into one
# combined geography, so expanding the region here is all that is needed for
# every downstream metric to become vintage-consistent.
resolve_region_areas <- function(areas,
                                 region_mode = "original",
                                 lookup = get_nuts_lookup(),
                                 available_areas = NULL) {
  areas <- as.character(areas)
  region_mode <- normalize_region_mode(region_mode)

  if (identical(region_mode, "original") || length(areas) == 0) {
    return(list(areas = areas, expanded = character(0), warnings = character(0)))
  }

  regions <- areas[is_region_label(areas, lookup)]

  if (length(regions) == 0) {
    return(list(areas = areas, expanded = character(0), warnings = character(0)))
  }

  keep <- setdiff(areas, regions)
  members <- unlist(lapply(
    regions,
    region_municipalities,
    lookup = lookup,
    available_areas = available_areas
  ))

  warnings <- unlist(lapply(regions, function(region) {
    region_coverage_warning(region_coverage(region, lookup, available_areas))
  }))

  list(
    areas = sort(unique(c(keep, members))),
    expanded = regions,
    warnings = as.character(stats::na.omit(warnings))
  )
}

# Regions whose definition changed between the NUTS vintages the app mixes.
# Portugal and Norte are byte-identical across NUTS-2013 and NUTS-2024, so only
# the genuinely affected regions are flagged.
NUTS_BREAK_YEAR <- 2022L
NUTS_BREAK_REGIONS <- c("Alentejo", "Centro", "Oeste e Vale do Tejo", "Lezíria do Tejo", "Médio Tejo", "Oeste")

region_vintage_warning <- function(areas, years = NULL, region_mode = "original") {
  region_mode <- normalize_region_mode(region_mode)

  if (identical(region_mode, "municipal")) {
    return(NULL)
  }

  affected <- intersect(as.character(areas), NUTS_BREAK_REGIONS)
  if (length(affected) == 0) {
    return(NULL)
  }

  spans_break <- is.null(years) ||
    (any(as.integer(years) < NUTS_BREAK_YEAR) && any(as.integer(years) >= NUTS_BREAK_YEAR))

  if (!spans_break) {
    return(NULL)
  }

  as.character(glue::glue(
    "Atenção: {paste(affected, collapse = ', ')} mudou de definição em ",
    "{NUTS_BREAK_YEAR} (NUTS-2024). Os óbitos até {NUTS_BREAK_YEAR - 1} seguem a ",
    "definição NUTS-2013 e a partir de {NUTS_BREAK_YEAR} a NUTS-2024, pelo que a ",
    "série tem uma quebra que não é epidemiológica. Escolha ",
    "'Aproximação por municípios' para uma série contínua."
  ))
}
