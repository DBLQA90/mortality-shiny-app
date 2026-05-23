# =========================================================
# INE client and live downloader
# =========================================================

ine_client <- ineptr2::INEClient$new(lang = "PT", timeout = 3600)

get_cat_code <- function(value, dimension_values, dim_name = "unknown") {
  requested <- as.character(value)
  requested_norm <- normalize_ine_label(requested)
  dim_number <- suppressWarnings(as.integer(sub("^dim", "", dim_name)))
  code_col <- if ("cat_id" %in% names(dimension_values)) "cat_id" else "categ_cod"

  dimension_values_for_dim <- dimension_values
  if (!is.na(dim_number) && "dim_num" %in% names(dimension_values_for_dim)) {
    dimension_values_for_dim <- dimension_values_for_dim %>%
      dplyr::filter(as.integer(.data$dim_num) == dim_number)
  }

  codes <- dimension_values_for_dim %>%
    dplyr::mutate(
      categ_dsg_chr = as.character(categ_dsg),
      categ_dsg_norm = normalize_ine_label(categ_dsg_chr),
      categ_cod_chr = if ("categ_cod" %in% names(.)) as.character(categ_cod) else NA_character_,
      cat_id_chr = if ("cat_id" %in% names(.)) as.character(cat_id) else NA_character_
    ) %>%
    dplyr::filter(
      categ_dsg_chr %in% requested |
        categ_dsg_norm %in% requested_norm |
        categ_cod_chr %in% requested |
        cat_id_chr %in% requested
    ) %>%
    dplyr::pull(dplyr::all_of(code_col)) %>%
    as.character() %>%
    unique()

  if (length(codes) == 0) {
    stop(
      glue::glue("Não foi possível mapear categorias para {dim_name}: {paste(requested, collapse = ', ')}"),
      call. = FALSE
    )
  }

  codes
}

# Bounded cache for INE dimension metadata (per indicator)
dim_values_cache <- cachem::cache_mem(
  max_size = 30 * 1024^2,
  max_age  = 24 * 60 * 60
)

get_dim_values_cached <- memoise::memoise(
  function(indicator) {
    with_persistent_cache(
      path = metadata_cache_file(indicator),
      max_age = persistent_metadata_cache_max_age,
      label = glue::glue("metadata for indicator {indicator}"),
      ine_client$get_dim_values(indicator)
    )
  },
  cache = dim_values_cache
)

get_dimension_categories <- function(indicators, target_dim_num) {
  purrr::map_dfr(indicators, function(indicator) {
    dv <- get_dim_values_cached(indicator)

    if (!"categ_cod" %in% names(dv)) {
      dv$categ_cod <- dv$cat_id
    }
    if (!"categ_ord" %in% names(dv)) {
      dv$categ_ord <- seq_len(nrow(dv))
    }
    if (!"categ_nivel" %in% names(dv)) {
      dv$categ_nivel <- NA_character_
    }

    dv %>%
      dplyr::filter(as.integer(.data$dim_num) == .env$target_dim_num) %>%
      dplyr::mutate(
        indicator = indicator,
        categ_cod = as.character(.data$categ_cod),
        cat_id = as.character(.data$cat_id),
        categ_dsg = as.character(.data$categ_dsg),
        categ_ord_num = suppressWarnings(as.numeric(.data$categ_ord)),
        categ_nivel_num = suppressWarnings(as.numeric(.data$categ_nivel))
      )
  })
}

get_indicator_years <- function(indicators) {
  get_dimension_categories(indicators, target_dim_num = 1) %>%
    dplyr::transmute(year = suppressWarnings(as.integer(.data$categ_dsg))) %>%
    dplyr::filter(!is.na(.data$year)) %>%
    dplyr::pull(.data$year) %>%
    unique()
}

normalize_year_order <- function(year_order = "asc") {
  if (is.null(year_order) || length(year_order) == 0 || is.na(year_order[[1]])) {
    return("asc")
  }

  year_order <- as.character(year_order[[1]])
  if (!year_order %in% c("asc", "desc")) {
    return("asc")
  }
  year_order
}

order_years <- function(years, year_order = "asc") {
  years <- unique(as.integer(years))
  years <- years[!is.na(years)]

  sort(years, decreasing = identical(normalize_year_order(year_order), "desc"))
}

get_source_year_plan <- function(indicators, priorities, requested_years, year_order = "asc") {
  requested_years <- order_years(requested_years, year_order)
  plan <- tibble(
    indicator = indicators,
    source_priority = priorities
  ) %>%
    dplyr::arrange(dplyr::desc(.data$source_priority))

  covered_years <- integer(0)
  plan$years <- vector("list", nrow(plan))

  for (i in seq_len(nrow(plan))) {
    available_years <- requested_years[requested_years %in% get_indicator_years(plan$indicator[[i]])]
    source_years <- available_years[!available_years %in% covered_years]
    plan$years[[i]] <- source_years
    covered_years <- unique(c(covered_years, source_years))
  }

  plan %>%
    dplyr::filter(lengths(.data$years) > 0)
}

get_available_cause_choices <- function() {
  cause_categories <- get_dimension_categories(death_indicators, target_dim_num = 5)
  cause_codes_by_indicator <- split(cause_categories$categ_cod, cause_categories$indicator)
  available_cause_codes <- Reduce(intersect, cause_codes_by_indicator)

  cause_categories %>%
    dplyr::filter(.data$categ_cod %in% available_cause_codes) %>%
    dplyr::mutate(
      source_priority = dplyr::case_when(
        .data$indicator == death_indicator_current ~ 1L,
        .data$indicator == death_indicator_legacy ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    dplyr::arrange(.data$source_priority, .data$categ_ord_num, .data$categ_dsg) %>%
    dplyr::group_by(.data$categ_cod) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$categ_ord_num, .data$categ_dsg) %>%
    dplyr::pull(.data$categ_dsg)
}

get_metadata_or_fallback <- function(label, fallback, expr) {
  tryCatch(
    {
      value <- force(expr)
      if (length(value) == 0) {
        stop(glue::glue("A metadata query returned no {label}."), call. = FALSE)
      }
      value
    },
    error = function(e) {
      warning(
        glue::glue("Using fallback {label} because INE metadata could not be read: {conditionMessage(e)}"),
        call. = FALSE
      )
      fallback
    }
  )
}

# General downloader: can include or exclude cause (dim5)
download_data <- function(indicator, dims, has_cause = FALSE) {
  dv   <- get_dim_values_cached(indicator)
  cats <- purrr::imap(dims, ~ get_cat_code(.x, dv, dim_name = .y))
  names(cats) <- names(dims)

  args <- c(list(indicator), cats)

  out <- with_persistent_cache(
    path = data_cache_file(indicator, cats, has_cause),
    max_age = persistent_data_cache_max_age,
    label = glue::glue("data for indicator {indicator}"),
    {
      raw <- do.call(ine_client$get_data, args)

      if (is.null(raw) || nrow(raw) == 0) {
        stop("API returned no data. Check that the indicator and dimension filters are correct.", call. = FALSE)
      }

      out <- raw %>%
        dplyr::transmute(
          year     = as.integer(dim_1),
          area     = geodsg,
          sex_raw  = dim_3_t,
          sex      = dplyr::recode(
            as.character(dim_3),
            "T" = "HM",
            "1" = "H",
            "2" = "M",
            .default = dplyr::recode(
              as.character(sex_raw),
              "MF"       = "HM",
              "Total"    = "HM",
              "Homens"   = "H",
              "Mulheres" = "M",
              .default   = as.character(sex_raw)
            )
          ),
          age_band = dim_4_t,
          value    = as.numeric(valor),
          source_indicator = indicator
        )

      if (has_cause) {
        out <- out %>%
          dplyr::mutate(cause = raw$dim_5_t)
      }

      if (nrow(out) == 0) {
        stop("API returned no rows after normalisation.", call. = FALSE)
      }

      out
    }
  )

  if (is.null(out) || nrow(out) == 0) {
    stop("Cached or live INE query returned no usable rows.", call. = FALSE)
  }

  out
}

expand_download_slices <- function(years, areas, cause = NULL, year_order = "asc") {
  slices <- tidyr::expand_grid(
    year = order_years(years, year_order),
    area = sort(unique(as.character(areas)))
  )

  if (!is.null(cause)) {
    slices$cause <- as.character(cause)
  }

  slices
}

empty_download_data <- function(has_cause = FALSE) {
  out <- tibble(
    year = integer(),
    area = character(),
    sex_raw = character(),
    sex = character(),
    age_band = character(),
    value = numeric(),
    source_indicator = character()
  )

  if (has_cause) {
    out$cause <- character()
  }

  out
}

data_load_control <- new.env(parent = emptyenv())
data_load_control$cancel_checker <- NULL
data_load_control$cancelled <- FALSE

service_pending_shiny_events <- function() {
  tryCatch(
    later::run_now(timeoutSecs = 0),
    error = function(e) NULL
  )
  invisible(NULL)
}

is_data_load_cancelled <- function() {
  service_pending_shiny_events()

  checker <- data_load_control$cancel_checker
  if (!is.function(checker)) {
    return(FALSE)
  }

  cancelled <- isTRUE(tryCatch(checker(), error = function(e) FALSE))
  if (cancelled) {
    data_load_control$cancelled <- TRUE
  }
  cancelled
}

with_data_load_cancel_checker <- function(cancel_checker, expr) {
  old_checker <- data_load_control$cancel_checker
  old_cancelled <- data_load_control$cancelled

  data_load_control$cancel_checker <- cancel_checker
  data_load_control$cancelled <- FALSE

  on.exit({
    data_load_control$cancel_checker <- old_checker
    data_load_control$cancelled <- old_cancelled
  }, add = TRUE)

  value <- force(expr)
  attr(value, "cancelled") <- isTRUE(data_load_control$cancelled)
  value
}

download_data_slices <- function(indicator, years, areas, cause = NULL, has_cause = FALSE, year_order = "asc") {
  slices <- expand_download_slices(years, areas, cause, year_order = year_order)
  if (nrow(slices) == 0) {
    return(empty_download_data(has_cause = has_cause))
  }

  failures <- list()
  out <- list()
  cancelled <- FALSE

  for (i in seq_len(nrow(slices))) {
    if (is_data_load_cancelled()) {
      cancelled <- TRUE
      break
    }

    year <- slices$year[[i]]
    area <- slices$area[[i]]
    cause_value <- if (has_cause) slices$cause[[i]] else NULL

    result <- tryCatch(
      {
        dims <- list(dim1 = year, dim2 = area)
        if (has_cause) {
          dims$dim5 <- cause_value
        }
        download_data(indicator, dims = dims, has_cause = has_cause)
      },
      error = function(e) {
        failures[[length(failures) + 1L]] <<- glue::glue(
          "{indicator} / {year} / {area}{if (has_cause) paste0(' / ', cause_value) else ''}: {conditionMessage(e)}"
        )
        NULL
      }
    )

    if (!is.null(result) && nrow(result) > 0) {
      out[[length(out) + 1L]] <- result
    }

    if (is_data_load_cancelled()) {
      cancelled <- TRUE
      break
    }
  }

  if (cancelled) {
    warning(
      glue::glue("INE loading was interrupted after {length(out)} completed slice{if (length(out) == 1) '' else 's'} for indicator {indicator}. Completed slices remain cached."),
      call. = FALSE
    )
  }

  if (length(failures) > 0) {
    warning(
      glue::glue(
        "Some INE slices could not be loaded and will be omitted from this run: {paste(unlist(failures), collapse = '; ')}"
      ),
      call. = FALSE
    )
  }
  if (length(out) == 0) {
    warning(
      glue::glue("No cached or live INE slices were available for indicator {indicator}; this source will be omitted from this run."),
      call. = FALSE
    )
    return(empty_download_data(has_cause = has_cause))
  }

  dplyr::bind_rows(out)
}
