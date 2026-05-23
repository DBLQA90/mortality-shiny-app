# =========================================================
# Lazy INE loader: get data only for selected area + cause
# =========================================================

prepare_population_data <- function(indicator, years, area, source_priority, year_order = "asc") {
  download_data_slices(
    indicator,
    years     = years,
    areas     = area,
    has_cause = FALSE,
    year_order = year_order
  ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      source_priority = source_priority,
      source_indicator = indicator
    )
}

prepare_death_data <- function(indicator, years, area, cause, source_priority, year_order = "asc") {
  download_data_slices(
    indicator,
    years     = years,
    areas     = area,
    cause     = cause,
    has_cause = TRUE,
    year_order = year_order
  ) %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(
      age_band = dplyr::case_when(
        age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ age_band
      )
    ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      source_priority = source_priority,
      source_indicator = indicator
    )
}


get_data_for_snapshot <- function(area, cause, years = year_of_interest) {
  df_pop <- get_snapshot_population_data(years, area)
  df_death <- get_snapshot_death_data(years, area, cause)

  df_full <- df_pop %>%
    dplyr::rename(population_source = source_indicator) %>%
    dplyr::right_join(
      df_death %>%
        dplyr::rename(death_source = source_indicator),
      by = c("year", "area", "sex", "age_band")
    ) %>%
    tidyr::replace_na(list(deaths = 0)) %>%
    dplyr::filter(!is.na(pop)) %>%
    dplyr::mutate(
      age_band = factor(age_band, levels = age_levels, ordered = TRUE)
    )

  if (nrow(df_full) == 0) {
    stop(
      "The snapshot population and death files do not have compatible rows for the selected filters.",
      call. = FALSE
    )
  }

  df_trunc <- df_full %>%
    dplyr::filter(!age_band %in% exclude_bands)

  list(full = df_full, trunc = df_trunc)
}

get_death_data_for <- function(area, cause, years = year_of_interest, year_order = "asc", data_source = "ine") {
  data_source <- normalize_data_source(data_source)
  if (identical(data_source, "snapshot")) {
    years <- order_years(years, year_order)
    return(get_snapshot_death_data(years, area, cause))
  }

  year_order <- normalize_year_order(year_order)
  years <- order_years(years, year_order)

  death_plan <- get_source_year_plan(
    indicators = c(death_indicator_legacy, death_indicator_current),
    priorities = c(1L, 2L),
    requested_years = years,
    year_order = year_order
  )

  df_death_sources <- purrr::pmap_dfr(
    death_plan,
    function(indicator, source_priority, years) {
      prepare_death_data(indicator, years, area, cause, source_priority, year_order = year_order)
    }
  )

  if (nrow(df_death_sources) == 0) {
    stop(
      "Não foi possível carregar dados de óbitos para a selecção actual.",
      call. = FALSE
    )
  }

  df_death_sources %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)
}

get_data_for <- function(area, cause, years = year_of_interest, year_order = "asc", data_source = "ine") {
  data_source <- normalize_data_source(data_source)
  if (identical(data_source, "snapshot")) {
    years <- order_years(years, year_order)
    return(get_data_for_snapshot(area, cause, years))
  }

  year_order <- normalize_year_order(year_order)
  years <- order_years(years, year_order)

  population_plan <- get_source_year_plan(
    indicators = c(population_indicator_current, population_indicator_legacy),
    priorities = c(1L, 2L),
    requested_years = years,
    year_order = year_order
  )
  death_plan <- get_source_year_plan(
    indicators = c(death_indicator_legacy, death_indicator_current),
    priorities = c(1L, 2L),
    requested_years = years,
    year_order = year_order
  )

  df_pop_sources <- purrr::pmap_dfr(
    population_plan,
    function(indicator, source_priority, years) {
      prepare_population_data(indicator, years, area, source_priority, year_order = year_order)
    }
  )

  if (nrow(df_pop_sources) == 0) {
    stop(
      "Não foi possível carregar dados de população para a selecção actual. Tente novamente ou escolha menos anos/locais.",
      call. = FALSE
    )
  }

  df_pop <- df_pop_sources %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)

  df_death_sources <- purrr::pmap_dfr(
    death_plan,
    function(indicator, source_priority, years) {
      prepare_death_data(indicator, years, area, cause, source_priority, year_order = year_order)
    }
  )

  if (nrow(df_death_sources) == 0) {
    stop(
      "Não foi possível carregar dados de óbitos para a selecção actual. Nada foi guardado na cache final; tente novamente ou escolha menos anos.",
      call. = FALSE
    )
  }

  df_death <- df_death_sources %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)

  # Combine
  df_full <- df_pop %>%
    dplyr::rename(population_source = source_indicator) %>%
    dplyr::right_join(
      df_death %>%
        dplyr::rename(death_source = source_indicator),
      by = c("year", "area", "sex", "age_band")
    ) %>%
    tidyr::replace_na(list(deaths = 0)) %>%
    dplyr::filter(!is.na(pop)) %>%
    dplyr::mutate(
      age_band = factor(age_band, levels = age_levels, ordered = TRUE)
    )

  if (nrow(df_full) == 0) {
    stop(
      "Os dados de população e óbitos carregados não têm anos/idades compatíveis para a selecção actual.",
      call. = FALSE
    )
  }

  df_trunc <- df_full %>%
    dplyr::filter(!age_band %in% exclude_bands)

  list(full = df_full, trunc = df_trunc)
}
