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
#                fixed membership list for every year. The series is continuous
#                and means "this region as defined by the chosen vintage". It is
#                an approximation only where a municipality cannot be matched by
#                name (see below), not where boundaries moved.
#
#   "original"   Keep INE's own regional rows and leave the 2022 change of
#                definition visible, with a warning naming the affected region.
#
# Regions are unions of whole municipalities, so no parish-level data is needed.
#
# Which membership list is used is the user's choice: the app ships one lookup
# per NUTS vintage and the vintage is selectable. Both cover the same 308
# municipalities under the same labels, so switching vintage never changes what
# data is read - only how it is grouped:
#
#   NUTS 2013   7 regions.  Alentejo includes Lezíria do Tejo; Centro includes
#               Oeste and Médio Tejo; Lisboa is one region, Área Metropolitana
#               de Lisboa. Use it to line up with anything published before the
#               2024 reform.
#   NUTS 2024   9 regions.  Lezíria do Tejo, Oeste and Médio Tejo form Oeste e
#               Vale do Tejo; Área Metropolitana de Lisboa splits into Grande
#               Lisboa and Península de Setúbal. The current definition, and
#               the default.
#
# Six region names exist in both vintages and mean different things in each -
# NUTS-2013 Centro has 100 municipalities, NUTS-2024 Centro has 77 - so the
# active vintage is stated wherever regional figures are shown rather than left
# to be inferred from the numbers.

# Newest first, so the first entry is the default.
nuts_vintage_choices <- c(
  "NUTS 2024 (definição actual)" = "2024",
  "NUTS 2013" = "2013"
)

normalize_nuts_vintage <- function(vintage) {
  vintage <- as.character(vintage)[[1]]

  if (isTRUE(vintage %in% unname(nuts_vintage_choices))) {
    vintage
  } else {
    unname(nuts_vintage_choices)[[1]]
  }
}

default_nuts_vintage <- function() {
  normalize_nuts_vintage(
    Sys.getenv("MORTALITY_NUTS_VINTAGE", unset = unname(nuts_vintage_choices)[[1]])
  )
}

nuts_lookup_path <- function(vintage = default_nuts_vintage()) {
  file.path(
    app_dir_or_wd(),
    "data",
    paste0("nuts_lookup_", normalize_nuts_vintage(vintage), ".rds")
  )
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

get_nuts_lookup <- function(vintage = default_nuts_vintage()) {
  path <- if (file.exists(vintage)) vintage else nuts_lookup_path(vintage)
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

# Every column of the lookup that names a grouping of municipalities, coarsest
# first. NUTS I is `Continente` plus the two autonomous regions; the islands
# carry the same name at all three levels because they are the same territory.
NUTS_LEVEL_COLUMNS <- c("nuts1", "nuts2", "nuts3")

nuts_level_values <- function(lookup, columns = NUTS_LEVEL_COLUMNS) {
  columns <- intersect(columns, names(lookup))
  if (length(columns) == 0) {
    return(character(0))
  }
  as.character(unlist(lookup[columns], use.names = FALSE))
}

# Region labels the lookup can rebuild from municipalities.
get_known_regions <- function(lookup = get_nuts_lookup()) {
  sort(unique(stats::na.omit(nuts_level_values(lookup))))
}

# The regions of a vintage offered in the selectors: NUTS I first, then NUTS II,
# each in INE's own order rather than alphabetically. The geography code carries
# the hierarchy, so ordering by the lowest municipality code reproduces INE's
# ordering without a hand-maintained list to drift out of date.
#
# NUTS III is deliberately left out. It would add 21 more entries to a dropdown
# already 316 long, and every one of them sits inside a NUTS II region that is
# offered, so nothing is unreachable - `region_municipalities()` resolves a
# NUTS III name if one is supplied.
#
# The islands are named identically at NUTS I and NUTS II, so they appear once.
# `Continente` is the only label NUTS I contributes that NUTS II does not.
region_choices_for <- function(vintage = default_nuts_vintage(),
                               lookup = get_nuts_lookup(vintage),
                               columns = c("nuts1", "nuts2")) {
  if (nrow(lookup) == 0) {
    return(character(0))
  }

  columns <- intersect(columns, names(lookup))

  ordered <- lapply(columns, function(column) {
    lookup %>%
      dplyr::filter(!is.na(.data[[column]])) %>%
      dplyr::group_by(region = .data[[column]]) %>%
      dplyr::summarise(first_code = min(.data$municipality_code), .groups = "drop") %>%
      dplyr::arrange(.data$first_code) %>%
      dplyr::pull(.data$region)
  })

  unique(unlist(ordered, use.names = FALSE))
}

# Every region name any vintage knows. The area vocabulary has to admit all of
# them, because a selection made under one vintage must survive being read back
# under another long enough to be reported as unavailable.
all_region_choices <- function() {
  sort(unique(unlist(lapply(
    unname(nuts_vintage_choices),
    function(vintage) get_known_regions(get_nuts_lookup(vintage))
  ))))
}

# Areas that name a region under some other vintage but not under this one.
#
# Without this the area is passed through untouched and reaches the loader as an
# ordinary place name, which fails with "snapshot file not found" - an error
# about storage for what is really a mismatch between the selection and the
# active definition.
unknown_region_areas <- function(areas, lookup = get_nuts_lookup()) {
  areas <- unique(as.character(areas))
  setdiff(intersect(areas, all_region_choices()), get_known_regions(lookup))
}

# What an area selection becomes when the vintage changes.
#
# Kept out of the observer that applies it so the decision can be tested:
# `MockShinySession` does not record `updateSelectInput` messages, so an
# observer's effect on a dropdown is invisible to the test harness.
#
# Municipality choices are identical across vintages, so a municipal selection
# always survives. A selected region may not exist in the new vintage - there is
# no NUTS-2013 `Grande Lisboa` - and those are returned separately so the caller
# can name them rather than leave them in place to fail later as unknown areas.
reconcile_area_selection <- function(selected,
                                     vintage = default_nuts_vintage(),
                                     lookup = get_nuts_lookup(vintage),
                                     exclude = character(0)) {
  selected <- unique(as.character(selected))
  choices <- setdiff(area_choices_for(vintage, lookup), exclude)

  list(
    choices = choices,
    selected = intersect(selected, choices),
    dropped = setdiff(selected, choices)
  )
}

vintage_mismatch_message <- function(areas,
                                     vintage = default_nuts_vintage(),
                                     lookup = get_nuts_lookup(vintage),
                                     region_mode = default_region_mode()) {
  # Only rebuilt regions have to exist in the lookup. Under the "original"
  # escape hatch the label is passed straight to INE, which still publishes
  # rows for the older vintage's regions, so refusing them would narrow the
  # escape hatch for no gain.
  if (identical(normalize_region_mode(region_mode), "original")) {
    return(NULL)
  }

  offenders <- unknown_region_areas(areas, lookup)

  if (length(offenders) == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "{paste(offenders, collapse = ', ')} não existe na definição NUTS ",
    "{vintage}. Escolha outra definição das regiões no topo da página, ou ",
    "seleccione uma região desta definição."
  ))
}

# The area vocabulary offered by the selectors under a given vintage: the
# national total, that vintage's regions, then every municipality.
#
# Municipalities come from the lookup rather than from a separate list, so the
# two can never disagree about which places exist. Both vintages carry the same
# 308 under the same labels, which is what makes switching vintage a regrouping
# rather than a change of coverage.
area_choices_for <- function(vintage = default_nuts_vintage(),
                             lookup = get_nuts_lookup(vintage)) {
  if (nrow(lookup) == 0) {
    return(character(0))
  }

  c(
    "Portugal",
    region_choices_for(vintage, lookup),
    sort(unique(as.character(lookup$municipality)))
  )
}

is_region_label <- function(area, lookup = get_nuts_lookup()) {
  as.character(area) %in% get_known_regions(lookup)
}

# Municipalities making up a region, restricted to those the data actually
# knows about. `available_areas` is the area vocabulary of the active data
# source; passing it lets the caller distinguish "region not covered" from
# "region covered except for two islands".
# Municipalities of a region at whichever NUTS level names it. Written against
# the level columns generically so adding NUTS I needed no new special case.
region_members <- function(region, lookup = get_nuts_lookup()) {
  columns <- intersect(NUTS_LEVEL_COLUMNS, names(lookup))
  if (length(columns) == 0 || nrow(lookup) == 0) {
    return(character(0))
  }

  hit <- Reduce(`|`, lapply(columns, function(column) {
    as.character(lookup[[column]]) %in% as.character(region)
  }))

  unique(as.character(lookup$municipality[!is.na(hit) & hit]))
}

region_municipalities <- function(region, lookup = get_nuts_lookup(), available_areas = NULL) {
  region <- as.character(region)

  members <- region_members(region, lookup)

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
  expected <- region_members(region, lookup)

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

  # Two different consequences, so two different messages.
  #
  # `Portugal` is not a region label, so it is never expanded into
  # municipalities: its own row is loaded and summed with whatever else is
  # selected, and the overlap really is counted twice (Portugal 124,830 + Beja
  # 552 = 125,382 for 2021).
  #
  # A region *is* expanded, and the expansion is a sorted unique set, so
  # selecting a region together with something inside it de-duplicates instead:
  # `Alentejo + Beja` yields Alentejo, and `Continente + Centro` yields
  # Continente. Nothing is double-counted, but the user does not get the
  # combination they asked for, which is still worth saying.
  double_counted <- character(0)
  absorbed <- character(0)

  if ("Portugal" %in% areas) {
    double_counted <- setdiff(areas, "Portugal")
  }

  regions <- setdiff(areas[is_region_label(areas, lookup)], "Portugal")

  # Compare municipality *sets*, not region-against-municipality-name.
  #
  # The earlier version intersected a region's municipalities with the other
  # selected areas, which caught "Alentejo + Beja" but never "Continente +
  # Norte": Norte is a region label, so it is not in anyone's municipality list
  # and the intersection was empty. That was harmless only while every offered
  # region sat at the same NUTS level and so was disjoint from the others.
  # NUTS I regions contain NUTS II ones, so containment between two selected
  # regions has to be detected directly.
  for (region in regions) {
    members <- region_municipalities(region, lookup)

    contained_regions <- setdiff(regions, region)
    contained_regions <- contained_regions[vapply(
      contained_regions,
      function(other) {
        other_members <- region_municipalities(other, lookup)
        length(other_members) > 0 && all(other_members %in% members)
      },
      logical(1)
    )]

    inside <- c(
      contained_regions,
      intersect(setdiff(areas, regions), members)
    )

    if (length(inside) > 0) {
      absorbed <- c(absorbed, as.character(glue::glue(
        "{region} já inclui {paste(sort(inside), collapse = ', ')}"
      )))
    }
  }

  messages <- character(0)

  if (length(double_counted) > 0) {
    messages <- c(messages, as.character(glue::glue(
      "Portugal já inclui {paste(sort(double_counted), collapse = ', ')}, e o ",
      "total nacional é somado tal como está: estes óbitos e esta população ",
      "seriam contados duas vezes."
    )))
  }

  if (length(absorbed) > 0) {
    messages <- c(messages, as.character(glue::glue(
      "{paste(absorbed, collapse = '; ')}. As regiões são somadas a partir dos ",
      "seus municípios, pelo que a área contida é absorvida em vez de contada ",
      "duas vezes - o resultado é apenas a área maior."
    )))
  }

  if (length(messages) == 0) {
    return(NULL)
  }

  as.character(glue::glue("Selecção sobreposta. {paste(messages, collapse = ' ')}"))
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
