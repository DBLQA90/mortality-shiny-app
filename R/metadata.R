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
year_of_interest <- get_metadata_or_fallback(
  "years",
  fallback_year_of_interest,
  {
    years <- sort(union(
      get_indicator_years(population_indicators),
      get_indicator_years(death_indicators)
    ))
    first_with_population <- suppressWarnings(min(population_years, na.rm = TRUE))
    if (is.finite(first_with_population)) years[years >= first_with_population] else years
  }
)

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
    "O INE ainda não publica população residente por idade e sexo para ",
    "{paste(missing, collapse = ', ')}, pelo que {metric_label} não pode ser ",
    "calculada nesse(s) ano(s). Os óbitos, a mortalidade proporcional e os AVPP ",
    "continuam disponíveis."
  ))
}

local_area <- fallback_local_area

diseases <- get_metadata_or_fallback(
  "causes",
  fallback_diseases,
  get_available_cause_choices()
)
