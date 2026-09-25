# =========================================================
# Planning indicators: primary care (SNS Transparency portal)
# =========================================================
# Monthly primary-care indicators per ULS from transparencia.sns.gov.pt
# (tools/fetch_sns.R), for the location's ULS, its ARS and the Continente.
#
# What the source is, and what that implies:
#
#   Unit      the ULS primary-care area, from January 2024 (ACES before, not
#             used here: ACES do not map one-to-one onto ULS). The portal
#             reports every ULS separately, including the six that share a
#             municipality - its units follow the parishes, so they do not
#             overlap and the 39 add up to the Continente. An area made of
#             whole ULS (an ARS, a group, the Continente) is their sum.
#             No data for the autonomous regions.
#   Base      people registered in primary care ("inscritos"), not residents.
#   Months    several indicators accumulate over a cycle and reset: blood
#             pressure and HbA1c by semester (≈10% in January, 50-65% in June,
#             ≈12% again in July), the foot exam and the screenings by calendar
#             year. A July value is not comparable with a December one, so each
#             indicator is compared only at the end of its cycle (December, or
#             June and December), and the current, incomplete cycle is shown as
#             year-to-date against the same month a year earlier.
#   Revision  the latest month is provisional and revised at the end of the
#             next one.
#
# Aggregation is exact: every proportion is rebuilt as a numerator and a
# denominator per ULS (the denominator from the published count and
# proportion when only those two are given) and summed.

SNS_INDICATORS <- tibble::tribble(
  ~id,                 ~dataset,                                            ~numerator,                                                                                             ~proportion,                                          ~denominator,            ~cycle,     ~label,
  "sns_no_gp",         "utentes-inscritos-em-cuidados-de-saude-primarios",  "total_utentes_sem_mdf_atribuido",                                                                  NA,                                                    "utentes_inscritos_csp", "month",    "Utentes inscritos sem médico de família atribuído",
  "sns_consultation",  "utentes-inscritos-em-cuidados-de-saude-primarios",  NA,                                                                                                   "taxa_de_utilizacao_consultas_medicas_1_ano_todos_os_utentes", "utentes_inscritos_csp", "month", "Utentes com pelo menos uma consulta médica no último ano",
  "sns_mammography",   "rastreios-oncologicos",                             "contagem_de_mulheres_com_registo_de_mamografia_nos_ultimos_dois_anos",                                "proporcao_mulheres_50_70_a_c_mamogr_2_anos",          NA,                      "year",     "Mulheres de 50-70 anos com mamografia nos últimos dois anos",
  "sns_cervical",      "rastreios-oncologicos",                             "contagem_de_mulheres_com_colpocitologia_atualizada",                                                  "proporcao_mulheres_25_60_a_c_colpoc_atuali",          NA,                      "year",     "Mulheres de 25-60 anos com colpocitologia atualizada",
  "sns_colorectal",    "rastreios-oncologicos",                             "contagem_de_utentes_inscritos_com_rastreio_do_cancro_do_colon_e_reto_efetuado",                       "proporcao_utentes_50_75_a_c_rastreio_cancro_cr",      NA,                      "year",     "Utentes de 50-75 anos com rastreio do cancro do cólon e reto",
  "sns_diabetes_feet", "diabetes",                                          "contagem_de_utentes_inscritos_com_diabetes_com_exame_dos_pes_realizado_no_ultimo_ano",                "proporcao_dm_com_exame_pes_ultimo_ano",               NA,                      "year",     "Diabéticos com exame dos pés no último ano",
  "sns_diabetes_hba1c", "diabetes",                                         "contagem_de_utentes_inscritos_com_diabetes_com_ultimo_resultado_de_hgba1c_inferior_ou_igual_a_8_0",   "proporcao_dm_c_ultima_hgba1c_8_0",                    NA,                      "semester", "Diabéticos com última HbA1c igual ou inferior a 8,0%",
  "sns_hypertension",  "hipertensao",                                       "contagem_de_utentes_inscritos_com_hipertensao_arterial_com_pressao_arterial_inferior_a_150_90_mmhg_n", "proporcao_hipertensos_65_a_com_pa_150_90",            NA,                      "semester", "Hipertensos com menos de 65 anos com pressão arterial inferior a 150/90 mmHg",
  "sns_newborn_visit", "saude-da-mulher-e-crianca",                         "contagem_de_recem_nascidos_com_pelo_menos_uma_consulta_medica_de_vigilancia_nos_primeiros_28_dias_de", "proporcao_rn_c_cons_med_vigil_ate_28_dias_vida",      NA,                      "month",    "Recém-nascidos com consulta médica de vigilância nos primeiros 28 dias",
  "sns_newborn_home",  "saude-da-mulher-e-crianca",                         "contagem_de_recem_nascidos_que_tiveram_pelo_menos_um_domicilio_de_enfermagem_durante_os_primeiros_15", "proporcao_rn_c_domicilio_enf_ate_15o_dia_de_vida",    NA,                      "month",    "Recém-nascidos com domicílio de enfermagem nos primeiros 15 dias"
)
sns_ids <- SNS_INDICATORS$id
SNS_FIRST_PERIOD <- "2024-01"

SNS_CYCLE_LABELS <- c(
  month = "comparável mês a mês",
  semester = "acumula ao longo de cada semestre: compara-se em Junho e Dezembro",
  year = "acumula ao longo do ano: compara-se em Dezembro"
)

planning_sns_choices <- function() stats::setNames(SNS_INDICATORS$id, SNS_INDICATORS$label)

# The portal's unit label -> the app's ULS name, where the two differ in
# spelling only.
SNS_UNIT_RENAMES <- c(
  "ULS Gaia/Espinho" = "ULS Vila Nova de Gaia/Espinho",
  "ULS Póvoa Varzim/Vila Conde" = "ULS Póvoa de Varzim/Vila do Conde",
  "ULS Póvoa Varzim/Vila do Conde" = "ULS Póvoa de Varzim/Vila do Conde",
  "ULS Trás-os-Montes Alto Douro" = "ULS Trás-os-Montes e Alto Douro"
)

planning_sns_unit <- function(label) {
  key <- sub("^(Área dos )?CSP da ", "", as.character(label))
  key <- gsub("\\s*/\\s*", "/", key)
  ifelse(key %in% names(SNS_UNIT_RENAMES), unname(SNS_UNIT_RENAMES[key]), key)
}

planning_sns_dir <- function() file.path(infant_snapshot_root(), "sns")

planning_sns_available <- function() {
  all(file.exists(file.path(planning_sns_dir(), paste0(unique(SNS_INDICATORS$dataset), ".rds"))))
}

# Numerators and denominators per app unit and month, for every indicator.
planning_sns_components <- function() {
  key <- paste(infant_snapshot_root(), "sns", sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  units <- planning_uls_units("units")
  out <- list()
  for (dataset in unique(SNS_INDICATORS$dataset)) {
    path <- file.path(planning_sns_dir(), paste0(dataset, ".rds"))
    if (!file.exists(path)) next
    raw <- readRDS(path)
    raw <- raw[raw$period >= SNS_FIRST_PERIOD, , drop = FALSE]
    wide <- tidyr::pivot_wider(raw[, c("period", "unit", "field", "value")], names_from = field, values_from = value, values_fn = sum)
    wide$app_unit <- planning_sns_unit(wide$unit)
    unmatched <- setdiff(unique(wide$app_unit), units)
    if (length(unmatched) > 0) warning("SNS units without an app ULS: ", paste(unmatched, collapse = "; "), call. = FALSE)
    wide <- wide[wide$app_unit %in% units, , drop = FALSE]
    for (i in which(SNS_INDICATORS$dataset == dataset)) {
      spec <- SNS_INDICATORS[i, ]
      # A field the portal renamed or dropped leaves its indicator out, with a
      # warning, rather than breaking the others.
      needed <- stats::na.omit(c(spec$numerator, spec$proportion, spec$denominator))
      if (!all(needed %in% names(wide))) {
        warning("SNS indicator ", spec$id, " skipped: missing ", paste(setdiff(needed, names(wide)), collapse = ", "), call. = FALSE)
        next
      }
      proportion <- if (!is.na(spec$proportion)) wide[[spec$proportion]] else NULL
      denominator <- if (!is.na(spec$denominator)) wide[[spec$denominator]] else ifelse(proportion > 0, wide[[spec$numerator]] / proportion * 100, NA_real_)
      numerator <- if (!is.na(spec$numerator)) wide[[spec$numerator]] else proportion * denominator / 100
      frame <- tibble::tibble(indicator = spec$id, period = wide$period, unit = wide$app_unit, numerator = numerator, denominator = denominator)
      out[[length(out) + 1L]] <- frame %>%
        dplyr::group_by(.data$indicator, .data$period, .data$unit) %>%
        dplyr::summarise(numerator = sum(.data$numerator), denominator = sum(.data$denominator), .groups = "drop")
    }
  }
  table <- dplyr::bind_rows(out)
  assign(key, table, envir = planning_cache)
  table
}

# The app ULS units an area is made of, or NULL when the area does not line up
# with whole ULS. A municipality reads the ULS that serves it - the six that
# share one read all the ULS serving it, since the portal cannot say which of
# a municipality's users belong to which. Portugal reads the Continente (the
# portal has no autonomous regions).
planning_sns_area_units <- function(area, lookup = get_nuts_lookup()) {
  units <- planning_uls_units("units")
  if (area %in% units) return(list(units = area, label = area, note = NULL))
  if (area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL, "Continente")) {
    return(list(units = units, label = "Continente", note = if (area != "Continente") "Os dados do SNS cobrem apenas o Continente." else NULL))
  }
  members <- planning_area_members(area, lookup)
  if (length(members) == 0) return(NULL)
  unit_members <- lapply(stats::setNames(units, units), planning_area_members, lookup = lookup)
  if (length(members) == 1) {
    holding <- Filter(function(u) members %in% unit_members[[u]], units)
    if (length(holding) == 0) return(NULL)
    if (length(holding) > 1) {
      # A municipality divided between ULS: all of them serve it.
      return(list(units = holding, label = paste(holding, collapse = " + "),
                  note = paste0(area, " está repartido por ", length(holding), " ULS ao nível da freguesia; os dados do SNS são por ULS, pelo que se somam todas.")))
    }
    return(list(units = holding, label = holding[[1]],
                note = paste0("Os dados do SNS são por ULS: ", area, " mostra ", holding[[1]], ".")))
  }
  inside <- units[vapply(units, function(u) all(unit_members[[u]] %in% members), logical(1))]
  covered <- unique(unlist(unit_members[inside]))
  if (length(inside) > 0 && setequal(covered, members)) return(list(units = inside, label = area, note = NULL))
  NULL
}

# Values for labelled sets of units: `sets` is a named list of unit vectors.
planning_sns_table <- function(sets, ids = sns_ids, lookup = get_nuts_lookup()) {
  components <- planning_sns_components()
  if (nrow(components) == 0) return(tibble::tibble())
  # A set may name an area rather than ULS (a group, an ARS): expand it, so a
  # caller never gets an empty answer for an area the portal does cover.
  known <- unique(components$unit)
  sets <- lapply(sets, function(units) {
    unknown <- setdiff(units, known)
    if (length(unknown) == 0) return(units)
    expanded <- unlist(lapply(unknown, function(a) {
      resolved <- planning_sns_area_units(a, lookup)
      if (is.null(resolved)) character(0) else resolved$units
    }))
    unique(c(intersect(units, known), expanded))
  })
  components <- components[components$indicator %in% ids, , drop = FALSE]
  latest <- tapply(components$period, components$indicator, max)
  dplyr::bind_rows(lapply(names(sets), function(label) {
    components[components$unit %in% sets[[label]], , drop = FALSE] %>%
      dplyr::group_by(.data$indicator, .data$period) %>%
      dplyr::summarise(numerator = sum(.data$numerator, na.rm = TRUE), denominator = sum(.data$denominator, na.rm = TRUE),
                       units = dplyr::n(), .groups = "drop") %>%
      # Every unit of the set must report that month.
      dplyr::filter(.data$units == length(sets[[label]])) %>%
      dplyr::mutate(area = label)
  })) %>%
    dplyr::mutate(
      value = ifelse(.data$denominator > 0, .data$numerator / .data$denominator * 100, NA_real_),
      month = as.integer(substr(.data$period, 6, 7)),
      cycle = SNS_INDICATORS$cycle[match(.data$indicator, SNS_INDICATORS$id)],
      complete = .data$cycle == "month" | (.data$cycle == "year" & .data$month == 12L) | (.data$cycle == "semester" & .data$month %in% c(6L, 12L)),
      provisional = .data$period == as.vector(latest[.data$indicator])
    ) %>%
    dplyr::mutate(
      lower = planning_binomial_ci(round(.data$numerator), round(.data$denominator))$lower,
      upper = planning_binomial_ci(round(.data$numerator), round(.data$denominator))$upper
    ) %>%
    dplyr::select(area, indicator, period, month, value, lower, upper, numerator, denominator, complete, provisional)
}

# The location and its comparators as SNS unit sets, dropping the ones that do
# not line up with ULS and repeats of the same set.
planning_sns_sets <- function(areas, lookup = get_nuts_lookup()) {
  sets <- list()
  notes <- character(0)
  levels <- character(0)
  for (i in seq_len(nrow(areas))) {
    resolved <- planning_sns_area_units(areas$area[[i]], lookup)
    if (is.null(resolved)) next
    if (any(vapply(sets, function(s) setequal(s, resolved$units), logical(1)))) next
    sets[[resolved$label]] <- resolved$units
    levels[[resolved$label]] <- if (identical(resolved$label, "Continente")) "Portugal" else if (resolved$label %in% planning_uls_units("units")) "ULS" else areas$level[[i]]
    if (!is.null(resolved$note)) notes <- c(notes, resolved$note)
  }
  if (!"Continente" %in% names(sets)) {
    sets[["Continente"]] <- planning_uls_units("units")
    levels[["Continente"]] <- "Portugal"
  }
  levels[[names(sets)[[1]]]] <- "Local"
  list(sets = sets, levels = levels, notes = unique(notes))
}
