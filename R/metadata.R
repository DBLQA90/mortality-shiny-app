# =========================================================
# Available metadata
# =========================================================

# Years with a population denominator, tracked separately so the app can offer
# a death year that has no matching population estimate while still refusing to
# compute a rate for it.
population_years <- get_metadata_or_fallback(
  "population_years",
  fallback_population_years,
  sort(get_indicator_years(population_indicators))
)

# The union, so the latest death year stays selectable, but bounded below by the
# first year with a population estimate.
#
# The historical death indicator reaches back to 1980 while population starts in
# 1991, so an unbounded union offered 1980-1990: years with no denominator, no
# downloaded death chunks, and nothing the app can compute. Deaths-only years
# are useful going forward, where the denominator has yet to be published and
# counts are still meaningful; they are not useful going backward, where no
# denominator will ever exist. Rate metrics still check `population_years`
# before computing, which is what keeps 2024 honest.
# Now bounded at both ends, because the two sides no longer end together.
#
# Every metric the app computes needs deaths, so a year with population but no
# deaths - 2025, since the revised population series reaches further than the
# by-cause death indicator - is not selectable at all. Population running ahead
# is the reverse of the situation this used to guard against and would have
# offered a year that yields nothing.
# Years the by-cause death archive covers. Every metric except the infant pair
# reads it, so a year outside it can produce nothing at all.
death_archive_years <- get_metadata_or_fallback(
  "death_years",
  fallback_year_of_interest,
  sort(get_indicator_years(death_indicators))
)

year_of_interest <- get_metadata_or_fallback(
  "years",
  fallback_year_of_interest,
  {
    # Bounded above by what some metric can actually answer, not by population.
    # 2025 has population and live births but no deaths by cause, so it is
    # selectable for infant mortality and refused for everything else; 2026
    # would have nothing and is not offered.
    years <- sort(union(death_archive_years, infant_metric_years(infant_metric_id)))

    first_with_population <- suppressWarnings(min(population_years, na.rm = TRUE))
    if (is.finite(first_with_population)) years <- years[years >= first_with_population]

    years
  }
)

# Years that can produce a rate: both a numerator and a denominator. The
# observed and forecast tabs work only in rates, so their sliders use this
# rather than population_years, which now runs a year further than deaths.
# Years past the by-cause archive whose all-cause deaths by age were rebuilt from
# INE's weekly and parish deaths (tools/build_calendar_deaths_by_age.R): 2025
# until 0013166 reaches it. They serve all causes, both sexes, and areas made of
# whole NUTS II regions only; see calendar_selection_supported().
calendar_death_years <- local({
  files <- list.files(file.path(infant_snapshot_root(), "calendar_deaths_by_age"), pattern = "^year_\\d{4}\\.rds$")
  sort(setdiff(as.integer(sub("^year_(\\d{4})\\.rds$", "\\1", files)), as.integer(death_archive_years)))
})

rate_years <- sort(intersect(as.integer(population_years), union(as.integer(death_archive_years), calendar_death_years)))

# Years the user can select that the main death archive cannot serve.
# `calendar_ok`: the selection is one the rebuilt latest year can serve.
years_without_deaths <- function(years, calendar_ok = FALSE) {
  have <- as.integer(death_archive_years)
  if (isTRUE(calendar_ok)) have <- union(have, calendar_death_years)
  sort(setdiff(as.integer(years), have))
}

death_gap_message <- function(years, metric_label = "esta métrica") {
  missing <- years_without_deaths(years)

  if (length(missing) == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "O INE ainda não publica óbitos por causa de morte para ",
    "{paste(missing, collapse = ', ')}, pelo que {metric_label} não pode ser ",
    "calculada nesse(s) ano(s). A mortalidade infantil, que não depende da ",
    "lista de causas, continua disponível."
  ))
}

# Years the user can select that will not support a rate.
years_without_population <- function(years) {
  sort(setdiff(as.integer(years), as.integer(population_years)))
}

# Message for a selection that cannot produce a rate, naming the years and the
# reason, so this reads as a known limit of the source data rather than a fault.
population_gap_message <- function(years, metric_label = "esta métrica") {
  missing <- years_without_population(years)

  if (length(missing) == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "Não há população residente por idade e sexo disponível para ",
    "{paste(missing, collapse = ', ')}, pelo que {metric_label} não pode ser ",
    "calculada nesse(s) ano(s). Os óbitos, a mortalidade proporcional e os AVPP ",
    "continuam disponíveis."
  ))
}


# The area vocabulary of the default NUTS vintage. Built from the committed
# lookup so the regions offered and the regions the app can actually rebuild are
# the same list; the hand-written fallback only applies if no lookup shipped.
# Switching vintage at runtime updates the selectors in place, so this is the
# initial state rather than the only one.
local_area <- local({
  derived <- tryCatch(area_choices_for(), error = function(e) character(0))
  if (length(derived) == 0) fallback_local_area else derived
})

diseases <- get_metadata_or_fallback(
  "causes",
  fallback_diseases,
  get_available_cause_choices()
)

# The cause selectors: the short list, for every year, and the detailed malignant
# tumours of 0008461 (2002-2022) that the snapshot holds, in INE's order.
detailed_cause_choices <- local({
  held <- tryCatch(detailed_tumour_causes(), error = function(e) character(0))
  ordered <- tryCatch({
    cats <- get_dimension_categories(death_indicator_detailed, target_dim_num = 5)
    cats$categ_dsg[order(cats$categ_ord_num)]
  }, error = function(e) character(0))
  c(intersect(ordered, held), sort(setdiff(held, ordered)))
})

cause_choices <- if (length(detailed_cause_choices) > 0) {
  stats::setNames(
    list(diseases, detailed_cause_choices),
    c("Lista sucinta europeia (todos os anos)",
      paste0("Tumores malignos detalhados (INE 0008461, ", min(DETAILED_TUMOUR_YEARS), "-", max(DETAILED_TUMOUR_YEARS), ")"))
  )
} else {
  diseases
}
