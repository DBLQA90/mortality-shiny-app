# =========================================================
# Available metadata
# =========================================================

year_of_interest <- get_metadata_or_fallback(
  "years",
  fallback_year_of_interest,
  {
    sort(intersect(
      get_indicator_years(population_indicators),
      get_indicator_years(death_indicators)
    ))
  }
)

local_area <- fallback_local_area

diseases <- get_metadata_or_fallback(
  "causes",
  fallback_diseases,
  get_available_cause_choices()
)
