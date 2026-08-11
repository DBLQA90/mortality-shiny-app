# =========================================================
# Indirect standardisation (SMR / ISR) and multi-year pooling
# =========================================================
# Direct standardisation (see R/metrics.R) needs stable age-specific rates in
# the area being standardised. For small municipalities and rare causes those
# rates are driven by one or two deaths per age band, so the DSR becomes
# unstable and its interval very wide - or `calculate_dsr()` cannot estimate an
# interval at all.
#
# Indirect standardisation avoids that: it applies the *reference* area's
# age-specific rates to the local age structure to obtain expected deaths, then
# compares the local observed total against it. Only the local total has to be
# large enough to be informative, not each age band. This is the conventional
# choice for small-area mortality comparison.
#
# Multi-year pooling addresses the same sparsity from the other direction, by
# summing deaths and person-years over a 3- or 5-year window.

# Odd windows only, so the pooled value can be centred on its target year.
POOLING_WINDOW_CHOICES <- c(
  "Ano único" = "1",
  "3 anos (média móvel)" = "3",
  "5 anos (média móvel)" = "5"
)

normalize_pooling_window <- function(window) {
  window <- suppressWarnings(as.integer(window)[[1]])

  if (!is.finite(window) || window < 1L) {
    return(1L)
  }
  if (window %% 2L == 0L) {
    stop(
      "A janela de agregação plurianual tem de ser ímpar (1, 3 ou 5 anos) para poder ser centrada.",
      call. = FALSE
    )
  }

  window
}

make_period_label <- function(period_start, period_end) {
  dplyr::if_else(
    period_start == period_end,
    as.character(period_start),
    paste0(period_start, "-", period_end)
  )
}

# Rolling multi-year pooling of an age-band level table.
#
# Deaths and population are both summed over the window, so the denominator
# becomes person-years and every downstream rate stays a valid rate per 100,000
# - no averaging of previously computed rates, which would weight small years
# equally with large ones.
#
# Windows are centred on each target year. At the ends of the series the window
# is truncated rather than dropped, so the most recent year stays visible;
# `n_years` records how many years actually contributed and `period_label`
# names the real span, so a truncated window is never presented as a full one.
pool_age_data <- function(df, window = 1L, keys = c("area", "sex", "cause", "age_band")) {
  window <- normalize_pooling_window(window)
  keys <- intersect(keys, names(df))

  if (!all(c("year", "deaths", "pop") %in% names(df))) {
    stop("A tabela a agregar tem de conter as colunas `year`, `deaths` e `pop`.", call. = FALSE)
  }

  if (window == 1L) {
    return(
      df %>%
        dplyr::mutate(
          period_start = as.integer(year),
          period_end = as.integer(year),
          n_years = 1L,
          period_label = as.character(year)
        )
    )
  }

  half <- (window - 1L) %/% 2L
  years <- sort(unique(as.integer(df$year)))

  windows <- dplyr::bind_rows(lapply(years, function(target) {
    tibble::tibble(
      target_year = target,
      source_year = years[years >= target - half & years <= target + half]
    )
  }))

  df %>%
    dplyr::mutate(year = as.integer(year)) %>%
    dplyr::inner_join(windows, by = c("year" = "source_year"), relationship = "many-to-many") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("target_year", keys)))) %>%
    dplyr::summarise(
      deaths = sum(deaths, na.rm = TRUE),
      pop = sum(pop, na.rm = TRUE),
      period_start = min(year),
      period_end = max(year),
      n_years = dplyr::n_distinct(year),
      .groups = "drop"
    ) %>%
    dplyr::rename(year = target_year) %>%
    dplyr::mutate(period_label = make_period_label(period_start, period_end))
}

describe_pooling_window <- function(window) {
  window <- normalize_pooling_window(window)

  if (window == 1L) {
    return("Valores anuais (sem agregação plurianual).")
  }

  as.character(glue::glue(
    "Agregação plurianual de {window} anos centrada em cada ano: óbitos e ",
    "população somados na janela (denominador em pessoas-ano). Nos extremos da ",
    "série a janela é truncada e o rótulo indica o período efectivo."
  ))
}

# Observed and expected deaths for one area/period against a reference schedule.
#
# `df_area` and `df_ref` must each hold one row per age band, with `deaths` and
# `pop`. Expected deaths are sum(local population x reference age-specific
# rate). Age bands present in the area but missing from the reference cannot
# contribute an expected count, so they are dropped and reported through
# `missing_bands` rather than silently treated as zero risk.
compute_expected_deaths <- function(df_area, df_ref) {
  ref <- df_ref %>%
    dplyr::group_by(age_band) %>%
    dplyr::summarise(
      ref_deaths = sum(deaths, na.rm = TRUE),
      ref_pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(ref_pop), ref_pop > 0)

  joined <- df_area %>%
    dplyr::group_by(age_band) %>%
    dplyr::summarise(
      deaths = sum(deaths, na.rm = TRUE),
      pop = sum(pop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(ref, by = "age_band")

  usable <- joined %>%
    dplyr::filter(is.finite(ref_pop), ref_pop > 0, is.finite(pop))

  list(
    observed = sum(usable$deaths, na.rm = TRUE),
    expected = sum(usable$pop * usable$ref_deaths / usable$ref_pop, na.rm = TRUE),
    ref_rate = if (sum(usable$ref_pop) > 0) sum(usable$ref_deaths) / sum(usable$ref_pop) else NA_real_,
    missing_bands = setdiff(as.character(joined$age_band), as.character(usable$age_band)),
    bands_used = nrow(usable)
  )
}

# Standardised mortality ratio for one area/period against one reference.
#
# `refvalue = 100` expresses the reference as 100, the usual presentation for
# SMRs. Intervals come from PHEindicatormethods, which applies Byar's method
# (exact Poisson below 10 observed deaths) - the same convention used across
# the PHE/OHID small-area indicator suite.
#
# `isr` is the matching indirectly standardised *rate* per 100,000: the SMR
# multiplied by the reference crude rate. It puts the same comparison on a rate
# scale, so it can be read next to the crude and directly standardised rates.
compute_smr <- function(df_area, df_ref, confidence = 0.95, refvalue = 100, multiplier = 1e5) {
  na_result <- function(observed = NA_real_, expected = NA_real_, note = NA_character_) {
    tibble::tibble(
      observed = observed,
      expected = expected,
      smr = NA_real_,
      smr_lower = NA_real_,
      smr_upper = NA_real_,
      isr = NA_real_,
      isr_lower = NA_real_,
      isr_upper = NA_real_,
      ref_rate = NA_real_,
      bands_used = NA_integer_,
      note = note
    )
  }

  if (nrow(df_area) == 0 || nrow(df_ref) == 0) {
    return(na_result(note = "Sem dados para a área ou para a referência."))
  }

  parts <- compute_expected_deaths(df_area, df_ref)

  if (!is.finite(parts$expected) || parts$expected <= 0) {
    return(na_result(
      observed = parts$observed,
      expected = parts$expected,
      note = "Óbitos esperados nulos: a referência não tem mortalidade nestas faixas etárias."
    ))
  }

  # calculate_ISRatio() works from age-band vectors: local deaths/population
  # alongside the reference deaths/population for the same bands.
  ratio <- tryCatch(
    {
      area_bands <- df_area %>%
        dplyr::group_by(age_band) %>%
        dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), pop = sum(pop, na.rm = TRUE), .groups = "drop")
      ref_bands <- df_ref %>%
        dplyr::group_by(age_band) %>%
        dplyr::summarise(ref_deaths = sum(deaths, na.rm = TRUE), ref_pop = sum(pop, na.rm = TRUE), .groups = "drop")

      aligned <- area_bands %>%
        dplyr::inner_join(ref_bands, by = "age_band") %>%
        dplyr::filter(is.finite(ref_pop), ref_pop > 0) %>%
        dplyr::arrange(age_band)

      aligned %>%
        PHEindicatormethods::calculate_ISRatio(
          x = deaths,
          n = pop,
          x_ref = aligned$ref_deaths,
          n_ref = aligned$ref_pop,
          refvalue = refvalue,
          confidence = confidence,
          type = "full"
        )
    },
    error = function(e) NULL
  )

  if (is.null(ratio) || nrow(ratio) == 0 || !is.finite(ratio$value[[1]])) {
    return(na_result(
      observed = parts$observed,
      expected = parts$expected,
      note = "O intervalo do SMR não pôde ser estimado para esta selecção."
    ))
  }

  scale_to_rate <- parts$ref_rate * multiplier / refvalue

  tibble::tibble(
    observed = parts$observed,
    expected = parts$expected,
    smr = ratio$value[[1]],
    smr_lower = ratio$lowercl[[1]],
    smr_upper = ratio$uppercl[[1]],
    isr = ratio$value[[1]] * scale_to_rate,
    isr_lower = ratio$lowercl[[1]] * scale_to_rate,
    isr_upper = ratio$uppercl[[1]] * scale_to_rate,
    ref_rate = parts$ref_rate * multiplier,
    bands_used = as.integer(parts$bands_used),
    note = if (length(parts$missing_bands) > 0) {
      paste0(
        "Faixas etárias sem referência utilizável, excluídas do cálculo: ",
        paste(parts$missing_bands, collapse = ", "), "."
      )
    } else {
      NA_character_
    }
  )
}

# Whether an SMR differs from the reference at the requested confidence level.
# Reported as a plain-language flag so the annual comparison table can mark the
# areas that are genuinely above or below the reference, rather than leaving
# every reader to compare intervals by eye.
classify_smr <- function(smr_lower, smr_upper, refvalue = 100) {
  dplyr::case_when(
    !is.finite(smr_lower) | !is.finite(smr_upper) ~ NA_character_,
    smr_lower > refvalue ~ "Acima da referência",
    smr_upper < refvalue ~ "Abaixo da referência",
    TRUE ~ "Sem diferença significativa"
  )
}

# SMR for many areas/periods at once, against a reference series matched on the
# same period, sex and cause. Used by the annual comparison tab so Portugal, the
# region and the selected local area are all standardised against one schedule.
compute_smr_series <- function(df, ref_df, confidence = 0.95, refvalue = 100, multiplier = 1e5) {
  match_keys <- intersect(c("year", "period_label", "sex", "cause"), names(df))
  match_keys <- intersect(match_keys, names(ref_df))
  group_keys <- unique(c(match_keys, intersect("area", names(df))))

  if (nrow(df) == 0) {
    return(tibble::tibble())
  }

  groups <- df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(group_keys)))

  dplyr::bind_rows(lapply(seq_len(nrow(groups)), function(i) {
    key <- groups[i, , drop = FALSE]

    area_rows <- df %>% dplyr::semi_join(key, by = group_keys)
    ref_rows <- ref_df %>% dplyr::semi_join(key[, match_keys, drop = FALSE], by = match_keys)

    result <- compute_smr(
      df_area = area_rows,
      df_ref = ref_rows,
      confidence = confidence,
      refvalue = refvalue,
      multiplier = multiplier
    )

    dplyr::bind_cols(key, result)
  })) %>%
    dplyr::mutate(smr_class = classify_smr(smr_lower, smr_upper, refvalue = refvalue))
}
