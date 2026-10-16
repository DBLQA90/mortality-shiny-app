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


# The population a rate divides by. INE's estimates refer to 31 December; with
# basis "mean" each year's estimate is replaced by the mean of it and the
# previous year's, as INE does for rates. `df_pop` must already hold the
# previous years where they exist (see population_years_to_load()); a year
# whose previous estimate is missing - the first of the series - keeps its own.
# The previous years are then dropped, so the frame covers exactly `years`.
apply_population_basis <- function(df_pop, years, basis = default_population_basis()) {
  years <- as.integer(years)
  if (identical(normalize_population_basis(basis), "end") || nrow(df_pop) == 0) {
    return(df_pop[as.integer(df_pop$year) %in% years, , drop = FALSE])
  }
  before <- df_pop %>%
    dplyr::transmute(year = as.integer(.data$year) + 1L, .data$area, .data$sex, .data$age_band, pop_before = .data$pop)
  df_pop %>%
    dplyr::mutate(year = as.integer(.data$year)) %>%
    dplyr::filter(.data$year %in% years) %>%
    dplyr::left_join(before, by = c("year", "area", "sex", "age_band")) %>%
    dplyr::mutate(pop = ifelse(is.finite(.data$pop_before), (.data$pop + .data$pop_before) / 2, .data$pop)) %>%
    dplyr::select(-"pop_before")
}

# The years whose population a selection needs: with the mean, each year's
# predecessor too. Years the sources do not have are simply absent when loaded.
population_years_to_load <- function(years, basis = default_population_basis()) {
  years <- as.integer(years)
  if (identical(normalize_population_basis(basis), "end")) return(years)
  sort(unique(c(years, years - 1L)))
}

# The previous years' estimates the mean needs, read one year at a time: a year
# the snapshot does not have - 1990, before the series starts - is skipped rather
# than failing the whole load (the chunked reader errors on a missing year), and
# the first year then keeps its own estimate.
load_previous_population <- function(years, basis, loader) {
  extra <- setdiff(population_years_to_load(years, basis), as.integer(years))
  dplyr::bind_rows(lapply(extra, function(year) tryCatch(loader(year), error = function(e) NULL)))
}

# Deaths without a published age, filled in.
#
# INE publishes each municipality's all-ages total, but its age bands do not
# always add up to it: some deaths have no recorded age, and cells are
# suppressed where numbers are small - in 2014 and 2024 above all (in 2024, 28
# municipalities lose more than 10% of their all-cause deaths). The planning
# module fills them in with planning_complete_by_age(); the other tabs used to
# drop them, so Portugal's 2024 standardised rate differed between tabs and
# Barrancos had 19 deaths in 2024 instead of 32. Every death frame the tabs load
# now gets the same fill, so an area and year read the same everywhere.
#
# The fill is computed once per year, cause and sex on the whole snapshot (it
# needs every municipality and the national row, its two margins) and returned
# as the deaths added to each area and band. Adding it, rather than replacing
# the frame's values, leaves frames loaded live from INE their own recorded
# deaths. NULL where the snapshot lacks the year or the cause.
death_totals_raw <- function(year) {
  key <- paste(infant_snapshot_root(), "totalsraw", year, sep = "|")
  if (!exists(key, envir = planning_cache, inherits = FALSE)) {
    value <- NULL
    for (indicator in c("0013166", "0008206")) {
      path <- file.path(death_totals_dir(), indicator, paste0("year_", year, ".rds"))
      if (file.exists(path)) {
        value <- readRDS(path)
        break
      }
    }
    assign(key, value, envir = planning_cache)
  }
  get(key, envir = planning_cache, inherits = FALSE)
}

death_age_fill <- function(year, cause) {
  key <- paste(infant_snapshot_root(), "agefill", year, cause, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  value <- NULL
  path <- planning_death_cause_file(year, cause)
  pop <- planning_read("population", year)
  # The totals as published, not repaired: repair_death_totals() reads the
  # age bands of every cause group to replace blank totals with their band sum,
  # which costs most of the time here and changes nothing - a blank total gives
  # no missing deaths in planning_complete_by_age() either way.
  totals <- death_totals_raw(year)
  if (!is.null(path) && !is.null(pop) && !is.null(totals)) {
    labels <- unique(c(sort(unique(get_nuts_lookup("2024")$municipality)), planning_published_areas))
    bands <- age_levels
    deaths <- readRDS(path)
    matrix_of <- function(frame, value, sex) {
      out <- matrix(0, length(labels), length(bands), dimnames = list(labels, bands))
      rows <- frame[frame$sex == sex & frame$area %in% labels & frame$age_band %in% bands, , drop = FALSE]
      if (nrow(rows) > 0) {
        summed <- rowsum(rows[[value]], paste(rows$area, rows$age_band, sep = "\r"), reorder = FALSE)
        parts <- do.call(rbind, strsplit(rownames(summed), "\r", fixed = TRUE))
        out[parts] <- summed[, 1]
      }
      out
    }
    value <- lapply(stats::setNames(c("HM", "H", "M"), c("HM", "H", "M")), function(sex) {
      recorded <- matrix_of(deaths, "deaths", sex)
      population <- matrix_of(pop, "pop", sex)
      if (sum(population["Portugal", ]) <= 0) return(NULL)
      total <- stats::setNames(rep(NA_real_, length(labels)), labels)
      rows <- totals[totals$sex == sex & totals$cause == cause & totals$area %in% labels, , drop = FALSE]
      total[rows$area] <- rows$deaths
      list(added = planning_complete_by_age(recorded, population, total, labels)$deaths - recorded,
           recorded = rowSums(recorded))
    })
  }
  assign(key, value, envir = planning_cache)
  value
}

# Adds the filled deaths to a death frame (year, area, sex, cause, age_band,
# deaths), and records them in `deaths_filled` so a table can say how much of a
# figure was estimated. Areas outside the snapshot's municipalities and
# national rows (INE's regional rows, which are complete) are left as they are.
#
# The fill is computed on the snapshot. A frame loaded live from INE is filled
# only where its recorded deaths are the snapshot's (`check_recorded`): if INE
# has revised a year since, adding the snapshot's gap could count deaths twice,
# so that area is left as published and a warning says so.
fill_death_ages <- function(df_death, check_recorded = FALSE) {
  df_death$deaths_filled <- 0
  if (nrow(df_death) == 0) return(df_death)
  keys <- unique(df_death[, c("year", "cause")])
  unfilled <- character(0)
  for (i in seq_len(nrow(keys))) {
    fill <- tryCatch(death_age_fill(keys$year[[i]], keys$cause[[i]]), error = function(e) {
      warning("Deaths without a published age not filled for ", keys$cause[[i]], " ", keys$year[[i]], ": ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(fill)) next
    for (sex in intersect(names(fill), unique(df_death$sex))) {
      part <- fill[[sex]]
      if (is.null(part)) next
      rows <- which(df_death$year == keys$year[[i]] & df_death$cause == keys$cause[[i]] & df_death$sex == sex &
                      df_death$area %in% rownames(part$added) & as.character(df_death$age_band) %in% colnames(part$added))
      if (check_recorded && length(rows) > 0) {
        loaded <- tapply(df_death$deaths[rows], df_death$area[rows], sum)
        differs <- names(loaded)[abs(loaded - part$recorded[names(loaded)]) > 0.5]
        unfilled <- c(unfilled, differs)
        rows <- rows[!df_death$area[rows] %in% differs]
      }
      df_death$deaths_filled[rows] <- part$added[cbind(df_death$area[rows], as.character(df_death$age_band[rows]))]
    }
  }
  if (length(unfilled) > 0) {
    warning("INE's live deaths differ from the snapshot for ", paste(utils::head(unique(unfilled), 5), collapse = ", "),
            "; their deaths without a published age were not filled.", call. = FALSE)
  }
  df_death$deaths <- df_death$deaths + df_death$deaths_filled
  df_death
}

# The latest year, before INE's deaths by cause reach it -----------------------
#
# calendar_death_years holds all-cause deaths by age for the NUTS II 2024
# regions, Portugal and the Continente only (see planning_calendar_deaths()).
# The tabs work on municipalities, so each region's deaths in each age band are
# shared among its municipalities in proportion to the deaths their population
# would have at the region's age-specific rates. Any area made of whole regions
# then sums back to the published figure exactly; a single municipality would
# only carry the region's age pattern, so calendar_selection_supported()
# refuses it rather than present an estimate as observed.
calendar_municipal_deaths <- function(year) {
  key <- paste(infant_snapshot_root(), "calendarmunicipal", year, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  calendar <- planning_calendar_deaths(year)
  pop <- read_year_file("population", year)
  out <- NULL
  if (!is.null(calendar) && !is.null(pop)) {
    lookup <- get_nuts_lookup("2024")
    pop <- pop[pop$sex == "HM" & pop$age_band %in% age_levels, , drop = FALSE]
    calendar <- calendar[calendar$sex == "HM" & calendar$age_band %in% age_levels, , drop = FALSE]
    parts <- lapply(intersect(unique(lookup$nuts2), calendar$area), function(region) {
      members <- lookup$municipality[lookup$nuts2 == region]
      p <- pop[pop$area %in% members, c("area", "age_band", "pop"), drop = FALSE]
      region_pop <- tapply(p$pop, p$age_band, sum)
      region_deaths <- tapply(calendar$deaths[calendar$area == region], calendar$age_band[calendar$area == region], sum)
      rate <- region_deaths[p$age_band] / region_pop[p$age_band]
      p$deaths <- ifelse(is.finite(rate), p$pop * rate, 0)
      p
    })
    published <- calendar[calendar$area %in% planning_published_areas, c("area", "age_band", "deaths"), drop = FALSE]
    out <- dplyr::bind_rows(lapply(parts, function(p) p[, c("area", "age_band", "deaths")]), published) %>%
      dplyr::mutate(year = as.integer(year), sex = "HM", cause = "Todas as causas de morte",
                    source_indicator = unique(calendar$source_indicator)[[1]], deaths_filled = 0)
  }
  assign(key, out, envir = planning_cache)
  out
}

calendar_gap_message <- function(years) {
  as.character(glue::glue(
    "Para {paste(years, collapse = ', ')} o INE ainda não publica óbitos por causa e idade. ",
    "A aplicação reconstrói os óbitos de todas as causas por idade a partir dos óbitos ",
    "semanais, mas só por região NUTS II (e Portugal, Continente e Lisboa e Vale do Tejo), ",
    "para ambos os sexos. Fora disso, {paste(years, collapse = ', ')} fica de fora."
  ))
}

# calendar_death_years is set by R/metadata.R when the app starts; modules loaded
# on their own (the tests, the tools) read the files instead.
active_calendar_years <- function() {
  if (exists("calendar_death_years", inherits = TRUE)) calendar_death_years else planning_calendar_years()
}

# Whether the calendar years can serve a selection: all causes, both sexes, and
# every area either Portugal, the Continente or a set of whole NUTS II regions.
calendar_selection_supported <- function(areas, causes = "Todas as causas de morte", sex = "HM") {
  if (!identical(as.character(sex), "HM") || !all(as.character(causes) == "Todas as causas de morte")) return(FALSE)
  areas <- setdiff(as.character(areas), planning_published_areas)
  if (length(areas) == 0) return(TRUE)
  lookup <- get_nuts_lookup("2024")
  # A region given by name (an indirect rate's reference) stands for its municipalities.
  named <- areas[areas %in% c(lookup$nuts1, lookup$nuts2) & !areas %in% lookup$municipality]
  areas <- unique(c(setdiff(areas, named), lookup$municipality[lookup$nuts1 %in% named | lookup$nuts2 %in% named]))
  regions <- unique(lookup$nuts2[lookup$municipality %in% areas])
  length(regions) > 0 && all(areas %in% lookup$municipality) &&
    setequal(areas, lookup$municipality[lookup$nuts2 %in% regions])
}

# Deaths of `years` for `area`, from the by-cause archive and, for the calendar
# years, from the rebuilt all-cause deaths (only for all causes).
load_snapshot_deaths <- function(years, area, cause) {
  calendar <- intersect(as.integer(years), active_calendar_years())
  archived <- setdiff(as.integer(years), calendar)
  parts <- list()
  if (length(archived) > 0) parts[[1]] <- fill_death_ages(get_snapshot_death_data(archived, area, cause))
  if (length(calendar) > 0) {
    if (!all(as.character(cause) == "Todas as causas de morte")) stop(calendar_gap_message(calendar), call. = FALSE)
    rows <- dplyr::bind_rows(lapply(calendar, calendar_municipal_deaths))
    parts[[length(parts) + 1L]] <- rows[rows$area %in% area, , drop = FALSE]
  }
  dplyr::bind_rows(parts)
}

# The same calendar rows for a frame loaded live from INE, which has no deaths
# by age for those years.
add_calendar_deaths <- function(df_death, years, area, cause) {
  calendar <- setdiff(intersect(as.integer(years), active_calendar_years()), as.integer(df_death$year))
  if (length(calendar) == 0 || !all(as.character(cause) == "Todas as causas de morte")) return(df_death)
  rows <- dplyr::bind_rows(lapply(calendar, calendar_municipal_deaths))
  dplyr::bind_rows(df_death, rows[rows$area %in% area, , drop = FALSE])
}

get_data_for_snapshot <- function(area, cause, years = year_of_interest, population_basis = default_population_basis()) {
  df_pop <- dplyr::bind_rows(
    get_snapshot_population_data(years, area),
    load_previous_population(years, population_basis, function(year) get_snapshot_population_data(year, area))
  ) %>%
    apply_population_basis(years, population_basis)
  df_death <- load_snapshot_deaths(years, area, cause)

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
    return(load_snapshot_deaths(years, area, cause))
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
    dplyr::select(-source_priority) %>%
    fill_death_ages(check_recorded = TRUE) %>%
    add_calendar_deaths(years, area, cause)
}

get_data_for <- function(area, cause, years = year_of_interest, year_order = "asc", data_source = "ine",
                         population_basis = default_population_basis()) {
  data_source <- normalize_data_source(data_source)
  if (identical(data_source, "snapshot")) {
    years <- order_years(years, year_order)
    return(get_data_for_snapshot(area, cause, years, population_basis = population_basis))
  }

  year_order <- normalize_year_order(year_order)
  years <- order_years(years, year_order)

  population_plan <- get_source_year_plan(
    indicators = names(population_source_priorities),
    priorities = unname(population_source_priorities),
    requested_years = population_years_to_load(years, population_basis),
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
    dplyr::select(-source_priority) %>%
    apply_population_basis(years, population_basis)

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
    dplyr::select(-source_priority) %>%
    fill_death_ages(check_recorded = TRUE) %>%
    add_calendar_deaths(years, area, cause)

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
