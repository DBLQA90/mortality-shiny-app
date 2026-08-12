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

# The union, so the latest death year stays selectable. Rate metrics check
# `population_years` before computing.
year_of_interest <- get_metadata_or_fallback(
  "years",
  fallback_year_of_interest,
  {
    sort(union(
      get_indicator_years(population_indicators),
      get_indicator_years(death_indicators)
    ))
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
