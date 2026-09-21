# =========================================================
# Planning indicators: one-click location profile (Word)
# =========================================================
# An editable .docx portrait of one location against its comparators, for the
# teams writing a local health plan:
#
#   1. in short        the indicators significantly above or below Portugal
#   2. key indicators  every indicator in its latest period, the location beside
#                      its comparators and Portugal, with significance marks
#   3. population      pyramid against Portugal
#   4. evolution       small multiples of six indicators, location and comparators
#   5. among the ULS   the location's ULS ranked on the event rates
#   6. causes of death proportional mortality, all ages and under 75
#   7. notes           marks, significance method, the * notes, sources
#
# Word rather than PDF because the teams paste from it and edit it. Charts are
# static ggplot images with the tab's colours (one per level, fixed).

PLANNING_PROFILE_TRENDS <- c("ageing_index", "birth_rate", "death_rate", "infant_rate", "life_expectancy", "fertility_index")
PLANNING_PROFILE_RANKED <- c("birth_rate", "death_rate", "infant_rate", "perinatal_rate", "teen_births_pct",
                             "preterm_pct", "low_birth_weight_pct", "life_expectancy", "life_expectancy_65")

planning_profile_theme <- function() {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = PLANNING_INK$secondary),
      axis.text = ggplot2::element_text(colour = PLANNING_INK$muted),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = PLANNING_INK$grid, linewidth = 0.3),
      strip.text = ggplot2::element_text(colour = PLANNING_INK$primary, face = "bold", hjust = 0),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Write the profile of `areas$area[[1]]` (with `areas` its comparators, levels
# included, Portugal already swapped for the chosen benchmark) to `path`.
write_planning_profile <- function(path,
                                   areas,
                                   lookup = get_nuts_lookup(),
                                   vintage = default_nuts_vintage(),
                                   data_date = NA_character_,
                                   benchmark = "Portugal",
                                   education_min_age = 0L,
                                   last_year = NULL,
                                   progress = function(value, detail = NULL) invisible(NULL)) {
  for (pkg in c("officer", "flextable", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("The ", pkg, " package is required for the profile.", call. = FALSE)
  }
  local <- areas$area[[1]]
  if (!benchmark %in% areas$area) areas <- dplyr::bind_rows(areas, tibble::tibble(area = benchmark, level = "Portugal"))
  all_years <- sort(unique(unlist(lapply(PLANNING_INDICATORS$id, planning_indicator_years))))
  if (is.null(last_year)) last_year <- max(all_years)
  years <- all_years[all_years <= last_year & all_years >= last_year - 19L]
  census <- all_years[all_years <= last_year & all_years %% 10L == 1L]
  years <- sort(unique(c(years, utils::tail(census, 1))))

  progress(0.1, "indicadores")
  table <- planning_indicator_table(areas$area, years, lookup = lookup, education_min_age = education_min_age) %>%
    planning_add_significance(benchmark)
  summary <- planning_profile_summary(table, areas, education_min_age)

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, paste0("Perfil do local: ", local), style = "heading 1")
  doc <- officer::body_add_par(doc, paste0(
    "Indicadores de apoio aos Planos Locais de Saúde. Gerado em ", format(Sys.Date(), "%d/%m/%Y"),
    "; dados do INE importados até ", ifelse(is.na(data_date), "data desconhecida", data_date),
    "; regiões NUTS ", vintage, "; método ", PLANNING_METHOD_VERSION, ". Referência: ", benchmark, "."
  ), style = "Normal")

  # 1. In short --------------------------------------------------------------
  doc <- officer::body_add_par(doc, "Em resumo", style = "heading 2")
  latest <- summary$rows
  higher <- latest[latest$significance %in% "Superior", , drop = FALSE]
  lower <- latest[latest$significance %in% "Inferior", , drop = FALSE]
  describe <- function(rows) paste0(rows$Indicador, ": ", rows$local_text, " (", benchmark, " ", rows$benchmark_text, ")")
  doc <- officer::body_add_par(doc, paste0(
    "Dos ", sum(!is.na(latest$significance)), " indicadores com intervalo de confiança no último período disponível, ",
    nrow(higher), " estão significativamente acima e ", nrow(lower), " significativamente abaixo de ", benchmark,
    " (intervalo de 95% inteiramente acima ou abaixo). A diferença não diz se é boa ou má: depende do indicador."
  ), style = "Normal")
  if (nrow(higher) > 0) {
    doc <- officer::body_add_par(doc, paste0("Acima de ", benchmark, ":"), style = "Normal")
    for (line in describe(higher)) doc <- officer::body_add_par(doc, paste0("\u2022 ", line), style = "Normal")
  }
  if (nrow(lower) > 0) {
    doc <- officer::body_add_par(doc, paste0("Abaixo de ", benchmark, ":"), style = "Normal")
    for (line in describe(lower)) doc <- officer::body_add_par(doc, paste0("\u2022 ", line), style = "Normal")
  }

  # 2. Key indicators ---------------------------------------------------------
  progress(0.3, "quadro de indicadores")
  # The wide table sits on landscape pages of its own.
  doc <- officer::body_end_section_portrait(doc)
  doc <- officer::body_add_par(doc, "Indicadores no último período disponível", style = "heading 2")
  ft <- flextable::flextable(summary$display)
  ft <- flextable::fontsize(ft, size = 7.5, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::bg(ft, bg = "#e1e0d9", part = "header")
  ft <- flextable::valign(ft, valign = "top", part = "all")
  ft <- flextable::padding(ft, padding = 1.5, part = "all")
  ft <- flextable::merge_v(ft, j = "Tema")
  ft <- flextable::border_inner_h(ft, border = officer::fp_border(color = "#e1e0d9", width = 0.5))
  n_areas <- ncol(summary$display) - 4L
  ft <- flextable::width(ft, j = seq_len(ncol(summary$display)),
                         width = c(0.8, 2.7, 0.7, 0.75, rep(min(1.1, 4.7 / max(n_areas, 1)), n_areas)))
  ft <- flextable::set_table_properties(ft, layout = "fixed")
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_add_par(doc, paste0(
    "▲ / ▼ / = : intervalo de confiança de 95% acima, abaixo ou a incluir o valor de ", benchmark,
    ". * ver notas no fim. As contagens absolutas mostram-se só para o local."
  ), style = "Normal")
  doc <- officer::body_end_section_landscape(doc)

  # 3. Population ---------------------------------------------------------------
  progress(0.45, "pirâmide etária")
  pop_year <- max(snapshot_years_for("population")[snapshot_years_for("population") <= last_year])
  own <- planning_pyramid(local, pop_year, lookup)
  reference <- planning_pyramid(benchmark, pop_year, lookup)
  if (nrow(own) > 0) {
    doc <- officer::body_add_par(doc, paste0("Estrutura etária, ", pop_year), style = "heading 2")
    doc <- planning_profile_add_plot(doc, planning_profile_pyramid(own, reference, local, benchmark), width = 6, height = 3.6)
  }

  # 4. Evolution --------------------------------------------------------------------
  progress(0.6, "evolução")
  doc <- officer::body_add_par(doc, "Evolução", style = "heading 2")
  doc <- planning_profile_add_plot(doc, planning_profile_trends(table, areas), width = 6.5, height = 5.2)

  # 5. Among the ULS -----------------------------------------------------------------
  progress(0.7, "posição entre as ULS")
  ranked <- planning_profile_ranking(local, areas, lookup, last_year, benchmark)
  if (!is.null(ranked)) {
    doc <- officer::body_add_par(doc, paste0("Posição de ", ranked$unit, " entre as ULS do Continente"), style = "heading 2")
    ft <- flextable::flextable(ranked$table)
    ft <- flextable::fontsize(ft, size = 8, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::bg(ft, bg = "#e1e0d9", part = "header")
    ft <- flextable::width(ft, j = 1:5, width = c(2.9, 0.9, 0.8, 0.9, 1.1))
    ft <- flextable::set_table_properties(ft, layout = "fixed")
    doc <- flextable::body_add_flextable(doc, ft)
    doc <- officer::body_add_par(doc, "Posição 1 = valor mais alto. A posição não tem em conta o intervalo de confiança: veja a coluna de significância.", style = "Normal")
  }

  # 6. Causes of death --------------------------------------------------------------
  progress(0.8, "mortalidade proporcional")
  proportional_years <- planning_proportional_years()
  proportional_years <- proportional_years[proportional_years <= last_year]
  if (length(proportional_years) > 0) {
    end_year <- max(proportional_years)
    prop <- planning_proportional_table(c(local, benchmark), end_year, lookup = lookup)
    if (nrow(prop) > 0 && any(!is.na(prop$share[prop$area == local]))) {
      doc <- officer::body_add_par(doc, paste0("Mortalidade proporcional, todas as idades, ", end_year - 2L, "-", end_year, " [I45]"), style = "heading 2")
      doc <- planning_profile_add_plot(doc, planning_profile_proportional(prop, local, benchmark), width = 6.5, height = 4)
    }
  }

  # 7. Notes -------------------------------------------------------------------------
  doc <- officer::body_add_par(doc, "Notas", style = "heading 2")
  notes <- c(
    "Cada área é a soma dos seus municípios e cada indicador a razão dessas somas. As ULS e ARS usam a composição actual, aplicada a todos os anos.",
    if (identical(benchmark, PLANNING_PORTUGAL_MUNICIPAL)) "Portugal: soma dos 308 municípios, sem os acontecimentos de residência desconhecida (compara igual com igual)." else "Portugal: total publicado pelo INE, que inclui os acontecimentos de residência desconhecida (0,3-0,9% dos óbitos).",
    "Significância: um valor está acima (abaixo) de Portugal quando todo o seu intervalo de confiança de 95% fica acima (abaixo) do valor de Portugal no mesmo período, como no PHE Fingertips. Só para taxas, proporções e esperança de vida.",
    paste("* Ganho médio e trabalhadores por sector:", PLANNING_INDICATOR_NOTES[["earnings_mean"]]),
    paste("* Esperança de vida:", PLANNING_INDICATOR_NOTES[["life_expectancy"]]),
    if (education_min_age > 0) paste0("Escolaridade: população com ", education_min_age, " e mais anos. ", PLANNING_EDUCATION_NOTE) else "Escolaridade: toda a população, incluindo crianças, como no ficheiro de apoio.",
    "Marcas: * taxa sobre menos de 1.000 nados-vivos; † óbitos com menos de 1 ano incompletos (1995-2001); ‡ mais de 2% dos óbitos sem idade redistribuídos; ≈ sectores de actividade estimados.",
    "Fonte: INE; tratamento pela aplicação de monitorização da mortalidade (DGS/PNS2030). Método descrito na nota metodológica."
  )
  for (line in notes) doc <- officer::body_add_par(doc, line, style = "Normal")

  progress(0.95, "a gravar")
  print(doc, target = path)
  invisible(path)
}

# The summary as rows (for the text) and as a display table (for Word).
planning_profile_summary <- function(table, areas, education_min_age = 0L) {
  local <- areas$area[[1]]
  benchmark_area <- areas$area[areas$level == "Portugal"][1]
  rows <- list()
  display <- list()
  for (id in PLANNING_INDICATORS$id) {
    spec <- planning_indicator_spec(id)
    own <- table[table$indicator == id & table$area == local & !is.na(table$value), , drop = FALSE]
    if (nrow(own) == 0) next
    year <- max(own$year)
    values <- table[table$indicator == id & table$year == year, , drop = FALSE]
    own <- values[values$area == local, , drop = FALSE]
    bench <- values[values$area == benchmark_area, , drop = FALSE]
    label <- planning_indicator_label(id, education_min_age)
    cell <- function(a) {
      if (!identical(a, local) && !isTRUE(spec$comparable)) return("")
      v <- values[values$area == a, , drop = FALSE]
      if (nrow(v) == 0 || is.na(v$value[[1]])) return("—")
      paste0(planning_format_value(v$value[[1]], spec$digits), v$flag[[1]], planning_significance_mark(v$significance[[1]]))
    }
    rows[[length(rows) + 1L]] <- tibble::tibble(
      Indicador = label, significance = own$significance[[1]],
      local_text = paste0(planning_format_value(own$value[[1]], spec$digits), " ", spec$unit, ", ", planning_period_label(id, year)),
      benchmark_text = if (nrow(bench) > 0) planning_format_value(bench$value[[1]], spec$digits) else "—"
    )
    display[[length(display) + 1L]] <- tibble::as_tibble(c(
      list(Tema = spec$theme, Indicador = label, Unidade = spec$unit, `Período` = planning_period_label(id, year)),
      stats::setNames(lapply(areas$area, cell), planning_series_name(areas$area, areas$level))
    ))
  }
  list(rows = dplyr::bind_rows(rows), display = dplyr::bind_rows(display))
}

planning_profile_add_plot <- function(doc, plot, width, height) {
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file))
  ggplot2::ggsave(file, plot, width = width, height = height, dpi = 200, bg = "white")
  officer::body_add_img(doc, file, width = width, height = height)
}

planning_profile_pyramid <- function(own, reference, local, benchmark) {
  own$signed <- ifelse(own$sex == "H", -own$share, own$share)
  own$band <- factor(own$age_band, levels = unique(own$age_band[order(own$lower)]))
  p <- ggplot2::ggplot(own, ggplot2::aes(x = .data$signed, y = .data$band)) +
    ggplot2::geom_col(ggplot2::aes(fill = ifelse(.data$sex == "H", "Homens", "Mulheres")), width = 0.85)
  if (nrow(reference) > 0) {
    reference$signed <- ifelse(reference$sex == "H", -reference$share, reference$share)
    reference$band <- factor(reference$age_band, levels = levels(own$band))
    p <- p + ggplot2::geom_step(data = reference[reference$sex == "H", ], ggplot2::aes(group = 1, colour = benchmark),
                                direction = "mid", linewidth = 0.6, orientation = "y") +
      ggplot2::geom_step(data = reference[reference$sex == "M", ], ggplot2::aes(group = 1, colour = benchmark),
                         direction = "mid", linewidth = 0.6, orientation = "y")
  }
  p + ggplot2::scale_fill_manual(values = c(Homens = PLANNING_LEVEL_COLOURS[["Local"]], Mulheres = "#8fb8ec")) +
    ggplot2::scale_colour_manual(values = stats::setNames(PLANNING_INK$primary, benchmark)) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(format(abs(x), decimal.mark = ","), "%")) +
    ggplot2::labs(x = paste0("% da população de ", local, " (contorno: ", benchmark, ")"), y = NULL) +
    planning_profile_theme()
}

planning_profile_trends <- function(table, areas) {
  rows <- table[table$indicator %in% PLANNING_PROFILE_TRENDS & table$area %in% areas$area & !is.na(table$value), , drop = FALSE]
  rows$level <- areas$level[match(rows$area, areas$area)]
  rows$series <- planning_series_name(rows$area, rows$level)
  rows$panel <- vapply(rows$indicator, function(id) {
    spec <- planning_indicator_spec(id)
    paste0(spec$label, " (", spec$unit, ")")
  }, character(1))
  rows$panel <- factor(rows$panel, levels = unique(rows$panel[order(match(rows$indicator, PLANNING_PROFILE_TRENDS))]))
  colours <- stats::setNames(PLANNING_LEVEL_COLOURS[areas$level], planning_series_name(areas$area, areas$level))
  rows$series <- factor(rows$series, levels = names(colours))
  local_rows <- rows[rows$level == "Local" & !is.na(rows$lower), , drop = FALSE]
  ggplot2::ggplot(rows, ggplot2::aes(x = .data$year, y = .data$value, colour = .data$series)) +
    ggplot2::geom_ribbon(data = local_rows, ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$series),
                         colour = NA, alpha = 0.15, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~panel, scales = "free_y", ncol = 2) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::labs(x = NULL, y = NULL, caption = "Faixa: intervalo de confiança de 95% do local. Indicadores de triénio identificados pelo último ano.") +
    planning_profile_theme()
}

# The location's ULS (the location itself when it is one) ranked among all ULS
# on the event rates. NULL when the location has no single ULS.
planning_profile_ranking <- function(local, areas, lookup, last_year, benchmark) {
  units <- planning_uls_units()
  # The smallest ULS holding every municipality of the location (a ULS with
  # the same municipalities is not offered as a comparator, but ranks here).
  members <- planning_area_members(local, lookup)
  holding <- Filter(function(u) all(members %in% planning_area_members(u, lookup)), units)
  if (length(holding) == 0) return(NULL)
  unit <- holding[[which.min(vapply(holding, function(u) length(planning_area_members(u, lookup)), integer(1)))]]
  rows <- lapply(PLANNING_PROFILE_RANKED, function(id) {
    years <- planning_indicator_years(id)
    years <- years[years <= last_year]
    if (length(years) == 0) return(NULL)
    year <- max(years)
    t <- planning_indicator_table(c(benchmark, units), year, ids = id, lookup = lookup) %>% planning_add_significance(benchmark)
    ranked <- t[t$area %in% units & !is.na(t$value), , drop = FALSE]
    ranked <- ranked[order(-ranked$value), , drop = FALSE]
    position <- match(unit, ranked$area)
    if (is.na(position)) return(NULL)
    spec <- planning_indicator_spec(id)
    tibble::tibble(
      Indicador = planning_indicator_label(id), `Período` = planning_period_label(id, year),
      Valor = planning_format_value(ranked$value[[position]], spec$digits),
      `Posição` = paste0(position, " de ", nrow(ranked)),
      `Face a Portugal` = ifelse(is.na(ranked$significance[[position]]), "—", ranked$significance[[position]])
    )
  })
  table <- dplyr::bind_rows(rows)
  if (nrow(table) == 0) return(NULL)
  list(unit = unit, table = table)
}

planning_profile_proportional <- function(prop, local, benchmark) {
  rows <- prop[prop$code != "C00" & !is.na(prop$share), , drop = FALSE]
  order <- rows$group[rows$area == local][order(rows$share[rows$area == local])]
  rows$group <- factor(rows$group, levels = unique(order))
  own <- rows[rows$area == local, , drop = FALSE]
  ref <- rows[rows$area == benchmark, , drop = FALSE]
  ggplot2::ggplot(own, ggplot2::aes(y = .data$group, x = .data$share)) +
    ggplot2::geom_col(fill = PLANNING_LEVEL_COLOURS[["Local"]], width = 0.7) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), orientation = "y", width = 0.25, colour = PLANNING_INK$muted, linewidth = 0.3) +
    ggplot2::geom_point(data = ref, ggplot2::aes(colour = benchmark), size = 2) +
    ggplot2::scale_colour_manual(values = stats::setNames(PLANNING_LEVEL_COLOURS[["Portugal"]], benchmark)) +
    ggplot2::labs(x = "% dos óbitos", y = NULL, caption = paste0("Barras: ", local, ", com o intervalo de confiança de 95%. Pontos: ", benchmark, ".")) +
    planning_profile_theme()
}
