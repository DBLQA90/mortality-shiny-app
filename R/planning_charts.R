# =========================================================
# Planning indicators: charts and display tables
# =========================================================
# Chart choices, by what each view has to show:
#
#   rates, shares, indices   lines over time, the location against its comparators
#                            (one colour per level, fixed), with the location's
#                            95% interval as a band; direct labels at the line ends
#   absolute counts          bars over time for the location alone, with the
#                            interval as error bars; no comparators, since a
#                            count depends on the size of the area
#   ranking of ULS           sorted horizontal bars, the location's unit in its
#                            colour and Portugal as a reference line
#   population pyramid       the location's bars against the nearest comparator's
#                            outline, as shares so different sizes compare
#   cause groups             the location's share as bars, comparators as dots
#
# Text and axes stay in neutral ink; colour carries only the identity of the
# series. Every chart has a table view in the tab.

PLANNING_INK <- list(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781", grid = "#e1e0d9", surface = "#fcfcfb")

# Shared look; `xaxis`, `yaxis` and the rest are merged into the defaults rather
# than passed twice, which plotly would resolve in favour of the defaults.
planning_plotly_layout <- function(p, xaxis = list(), yaxis = list(), margin = list(), ...) {
  axis <- list(gridcolor = PLANNING_INK$grid, zeroline = FALSE, linecolor = PLANNING_INK$grid,
               tickfont = list(color = PLANNING_INK$muted), titlefont = list(color = PLANNING_INK$secondary))
  p <- plotly::layout(
    p,
    paper_bgcolor = PLANNING_INK$surface,
    plot_bgcolor = PLANNING_INK$surface,
    font = list(family = "system-ui, -apple-system, 'Segoe UI', sans-serif", color = PLANNING_INK$secondary, size = 12),
    xaxis = utils::modifyList(axis, xaxis),
    yaxis = utils::modifyList(axis, yaxis),
    legend = list(orientation = "h", y = -0.2, font = list(color = PLANNING_INK$secondary)),
    margin = utils::modifyList(list(t = 30, r = 30, l = 60, b = 60), margin),
    separators = ",.",
    ...
  )
  plotly::config(p, displaylogo = FALSE, modeBarButtonsToRemove = list("lasso2d", "select2d", "autoScale2d"))
}

planning_series_name <- function(area, level) {
  ifelse(level %in% c("Local", "Portugal"), area, paste0(area, " (", level, ")"))
}

planning_hover_value <- function(value, lower, upper, digits, flag) {
  paste0(
    planning_format_value(value, digits), flag,
    ifelse(is.na(lower), "", paste0(" (IC 95%: ", planning_format_value(lower, digits), " - ", planning_format_value(upper, digits), ")"))
  )
}

# Lines for a location and its comparators. `areas` is a tibble of area and level
# in display order, the location first.
planning_trend_chart <- function(series, spec, areas) {
  series <- series %>%
    dplyr::filter(!is.na(.data$value)) %>%
    dplyr::mutate(
      period = vapply(.data$year, function(y) planning_period_label(spec$id, y), character(1)),
      hover = paste0(
        "<b>", .data$area, "</b><br>", .data$period, ": ",
        planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag)
      )
    )
  p <- plotly::plot_ly()
  local <- areas$area[[1]]
  own <- series[series$area == local, , drop = FALSE]

  if (nrow(own) > 0 && any(!is.na(own$lower))) {
    p <- plotly::add_ribbons(
      p, data = own, x = ~year, ymin = ~lower, ymax = ~upper,
      fillcolor = paste0(PLANNING_LEVEL_COLOURS[["Local"]], "26"), line = list(width = 0),
      name = paste0(local, " (IC 95%)"), hoverinfo = "skip", showlegend = TRUE
    )
  }

  shown <- areas[areas$area %in% series$area, , drop = FALSE]
  # Direct labels only while they cannot crowd each other; the legend always
  # names every line.
  label_ends <- nrow(shown) <= 4
  for (i in seq_len(nrow(areas))) {
    rows <- series[series$area == areas$area[[i]], , drop = FALSE]
    if (nrow(rows) == 0) next
    level <- areas$level[[i]]
    colour <- PLANNING_LEVEL_COLOURS[[level]]
    flagged <- nzchar(rows$flag)
    p <- plotly::add_trace(
      p, data = rows, x = ~year, y = ~value, type = "scatter", mode = "lines+markers",
      name = planning_series_name(areas$area[[i]], level),
      line = list(color = colour, width = if (identical(level, "Local")) 3 else 2),
      marker = list(color = ifelse(flagged, PLANNING_INK$surface, colour), size = 8,
                    line = list(color = colour, width = 2)),
      text = ~hover, hoverinfo = "text"
    )
    if (label_ends) {
      # Direct label at the end of the line, in ink rather than the series colour.
      last <- rows[which.max(rows$year), , drop = FALSE]
      p <- plotly::add_annotations(
        p, x = last$year, y = last$value, text = areas$area[[i]], xanchor = "left", xshift = 8,
        showarrow = FALSE, font = list(color = PLANNING_INK$secondary, size = 11)
      )
    }
  }

  breaks <- PLANNING_SERIES_BREAKS[PLANNING_SERIES_BREAKS$indicator == spec$id, , drop = FALSE]
  shapes <- lapply(breaks$year, function(y) list(
    type = "line", x0 = y - 0.5, x1 = y - 0.5, yref = "paper", y0 = 0, y1 = 1,
    line = list(color = PLANNING_INK$muted, dash = "dot", width = 1)
  ))

  span <- range(series$year)
  planning_plotly_layout(
    p,
    hovermode = "closest",
    shapes = shapes,
    margin = list(r = if (label_ends) 170 else 30),
    xaxis = list(title = if (spec$window > 1) "Último ano do triénio" else "Ano",
                 range = c(span[[1]] - 0.5, span[[2]] + 0.5), dtick = if (diff(span) > 12) 2 else 1),
    yaxis = list(title = spec$unit, rangemode = "tozero", tickformat = if (spec$digits == 0) ",.0f" else "")
  )
}

# Bars over time for an absolute count, the location alone.
planning_count_chart <- function(series, spec, local) {
  rows <- series %>%
    dplyr::filter(.data$area == local, !is.na(.data$value)) %>%
    dplyr::mutate(hover = paste0("<b>", local, "</b><br>", .data$year, ": ",
                                 planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag)))
  error <- if (any(!is.na(rows$lower))) {
    list(type = "data", symmetric = FALSE, array = rows$upper - rows$value, arrayminus = rows$value - rows$lower,
         color = PLANNING_INK$secondary, thickness = 1, width = 3)
  } else {
    NULL
  }
  p <- plotly::plot_ly(
    rows, x = ~year, y = ~value, type = "bar", name = local,
    marker = list(color = PLANNING_LEVEL_COLOURS[["Local"]]),
    error_y = error, text = ~hover, hoverinfo = "text", textposition = "none"
  )
  planning_plotly_layout(
    p, bargap = 0.35, showlegend = FALSE,
    margin = list(l = 80, b = 50),
    xaxis = list(title = "Ano", dtick = if (diff(range(rows$year)) > 12) 2 else 1),
    yaxis = list(title = paste0(spec$label, " (", spec$unit, ")"), tickformat = ",.0f", rangemode = "tozero")
  )
}

# All ULS for one indicator and year. `highlight` is the unit to mark (the
# location itself, or the ULS containing it).
# Significance against Portugal is a polarity: a warm and a cool pole with a
# grey midpoint (validated for colour-vision deficiency: worst adjacent OKLab
# distance 16.8). Orange rather than red, because which direction is desirable
# depends on the indicator (life expectancy against the death rate).
PLANNING_SIGNIFICANCE_COLOURS <- c(Superior = "#eb6834", Inferior = "#2a78d6", Semelhante = "#c3c2b7")

planning_ranking_chart <- function(ranking, spec, year, highlight = character(0), benchmark_area = "Portugal") {
  reference <- ranking$value[ranking$area == benchmark_area]
  ranking <- planning_add_significance(ranking, benchmark_area)
  units <- ranking %>%
    dplyr::filter(!.data$area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL), !is.na(.data$value)) %>%
    dplyr::arrange(.data$value) %>%
    dplyr::mutate(
      label = factor(.data$area, levels = .data$area),
      colour = ifelse(is.na(.data$significance), "#c3c2b7", PLANNING_SIGNIFICANCE_COLOURS[.data$significance]),
      outline = ifelse(.data$area %in% highlight, PLANNING_INK$primary, "rgba(0,0,0,0)"),
      hover = paste0("<b>", .data$area, "</b><br>", planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag),
                     ifelse(is.na(.data$significance), "", paste0("<br>Face a ", benchmark_area, ": ", tolower(.data$significance))))
    )
  error <- if (any(!is.na(units$lower))) {
    list(type = "data", symmetric = FALSE, array = units$upper - units$value, arrayminus = units$value - units$lower,
         color = PLANNING_INK$muted, thickness = 1, width = 2)
  } else {
    NULL
  }
  p <- plotly::plot_ly(
    units, y = ~label, x = ~value, type = "bar", orientation = "h",
    marker = list(color = ~colour, line = list(color = ~outline, width = 2)), error_x = error, text = ~hover, hoverinfo = "text", textposition = "none"
  )
  shapes <- if (length(reference) == 1 && is.finite(reference)) {
    list(list(type = "line", x0 = reference, x1 = reference, yref = "paper", y0 = 0, y1 = 1,
              line = list(color = PLANNING_LEVEL_COLOURS[["Portugal"]], dash = "dash", width = 2)))
  } else {
    list()
  }
  annotations <- if (length(shapes) > 0) {
    list(list(x = reference, y = 1.01, yref = "paper", text = paste0(benchmark_area, ": ", planning_format_value(reference, spec$digits)),
              showarrow = FALSE, xanchor = "left", font = list(color = PLANNING_INK$secondary, size = 11)))
  } else {
    list()
  }
  planning_plotly_layout(
    p, showlegend = FALSE, shapes = shapes, annotations = annotations, bargap = 0.25,
    margin = list(t = 40, l = 260, b = 50),
    xaxis = list(title = paste0(spec$label, " (", spec$unit, "), ", planning_period_label(spec$id, year))),
    yaxis = list(title = "", tickfont = list(color = PLANNING_INK$secondary, size = 11))
  )
}

# The location's population as bars (men left, women right), against a
# comparator's outline, both as shares of their own population.
planning_pyramid_chart <- function(local_pyramid, comparator_pyramid = NULL, local, comparator = NULL) {
  shape <- function(pyr) {
    pyr %>%
      dplyr::mutate(
        band = factor(.data$age_band, levels = unique(.data$age_band[order(.data$lower)])),
        signed = ifelse(.data$sex == "H", -.data$share, .data$share),
        sex_label = ifelse(.data$sex == "H", "Homens", "Mulheres")
      )
  }
  own <- shape(local_pyramid)
  p <- plotly::plot_ly()
  for (sex in c("H", "M")) {
    rows <- own[own$sex == sex, , drop = FALSE]
    p <- plotly::add_bars(
      p, data = rows, y = ~band, x = ~signed, orientation = "h",
      name = paste0(local, " - ", rows$sex_label[[1]]),
      marker = list(color = if (sex == "H") PLANNING_LEVEL_COLOURS[["Local"]] else "#8fb4e6"),
      text = ~paste0("<b>", local, "</b><br>", sex_label, ", ", age_band, ": ",
                     planning_format_value(share, 1), "% (", planning_format_value(pop, 0), ")"),
      hoverinfo = "text", textposition = "none"
    )
  }
  if (!is.null(comparator_pyramid) && nrow(comparator_pyramid) > 0) {
    other <- shape(comparator_pyramid)
    for (sex in c("H", "M")) {
      rows <- other[other$sex == sex, , drop = FALSE]
      rows <- rows[order(rows$lower), , drop = FALSE]
      p <- plotly::add_trace(
        p, data = rows, y = ~band, x = ~signed, type = "scatter", mode = "lines",
        line = list(color = PLANNING_INK$primary, width = 2, shape = "vh"),
        name = comparator, legendgroup = "comparator", showlegend = identical(sex, "H"),
        text = ~paste0("<b>", comparator, "</b><br>", sex_label, ", ", age_band, ": ", planning_format_value(share, 1), "%"),
        hoverinfo = "text"
      )
    }
  }
  limit <- max(abs(c(own$signed, if (!is.null(comparator_pyramid)) shape(comparator_pyramid)$signed)), na.rm = TRUE)
  ticks <- pretty(c(0, limit), n = 4)
  ticks <- sort(unique(c(-ticks, ticks)))
  planning_plotly_layout(
    p, barmode = "overlay", bargap = 0.08,
    margin = list(l = 110, b = 70),
    xaxis = list(title = "% da população: homens à esquerda, mulheres à direita", range = c(-limit, limit) * 1.08,
                 tickvals = ticks, ticktext = paste0(planning_format_value(abs(ticks), 1), "%")),
    yaxis = list(title = "", tickfont = list(color = PLANNING_INK$secondary))
  )
}

# Cause-group shares: the location as bars, comparators as dots on the same row.
planning_proportional_chart <- function(proportional, areas) {
  rows <- proportional %>%
    dplyr::filter(.data$code != "C00") %>%
    dplyr::left_join(areas, by = "area")
  local <- areas$area[[1]]
  order <- rows %>%
    dplyr::filter(.data$area == local) %>%
    dplyr::arrange(.data$code == "Outras", .data$share) %>%
    dplyr::pull(.data$group)
  rows$group <- factor(rows$group, levels = order)
  rows$hover <- paste0("<b>", rows$area, "</b><br>", rows$group, ": ", planning_format_value(rows$share, 1),
                       "% (", planning_format_value(rows$deaths, 0), " óbitos)")

  own <- rows[rows$area == local, , drop = FALSE]
  p <- plotly::plot_ly(
    own, y = ~group, x = ~share, type = "bar", orientation = "h", name = local,
    marker = list(color = PLANNING_LEVEL_COLOURS[["Local"]]), text = ~hover, hoverinfo = "text", textposition = "none"
  )
  for (i in seq_len(nrow(areas))[-1]) {
    comp <- rows[rows$area == areas$area[[i]], , drop = FALSE]
    if (nrow(comp) == 0) next
    p <- plotly::add_trace(
      p, data = comp, y = ~group, x = ~share, type = "scatter", mode = "markers",
      name = planning_series_name(areas$area[[i]], areas$level[[i]]),
      marker = list(color = PLANNING_LEVEL_COLOURS[[areas$level[[i]]]], size = 10, line = list(color = PLANNING_INK$surface, width = 2)),
      text = ~hover, hoverinfo = "text"
    )
  }
  planning_plotly_layout(
    p, bargap = 0.3,
    margin = list(l = 330, b = 60),
    xaxis = list(title = "% dos óbitos", rangemode = "tozero"),
    yaxis = list(title = "", tickfont = list(color = PLANNING_INK$secondary, size = 11))
  )
}

# The series as a table, most recent period first, one column per area.
planning_series_display <- function(series, spec, areas) {
  if (nrow(series) == 0) return(tibble::tibble())
  shown <- series %>%
    dplyr::filter(.data$area %in% areas$area) %>%
    dplyr::mutate(
      cell = ifelse(
        is.na(.data$value), "—",
        paste0(planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag),
               if ("significance" %in% names(series)) planning_significance_mark(.data$significance) else "")
      ),
      `Período` = vapply(.data$year, function(y) planning_period_label(spec$id, y), character(1))
    )
  columns <- planning_series_name(areas$area, areas$level)
  shown$column <- columns[match(shown$area, areas$area)]
  shown %>%
    dplyr::select(year, `Período`, column, cell) %>%
    tidyr::pivot_wider(names_from = column, values_from = cell) %>%
    dplyr::arrange(dplyr::desc(.data$year)) %>%
    dplyr::select(-year) %>%
    dplyr::select(dplyr::any_of(c("Período", columns)))
}

# Every indicator in its latest period for the location, with comparators beside
# the comparable ones.
planning_summary_display <- function(table, areas, education_min_age = 0L) {
  local <- areas$area[[1]]
  has_significance <- "significance" %in% names(table)
  columns <- planning_series_name(areas$area, areas$level)
  rows <- lapply(PLANNING_INDICATORS$id, function(id) {
    spec <- planning_indicator_spec(id)
    own <- table[table$indicator == id & table$area == local & !is.na(table$value), , drop = FALSE]
    if (nrow(own) == 0) return(NULL)
    year <- max(own$year)
    values <- table[table$indicator == id & table$year == year, , drop = FALSE]
    cells <- vapply(areas$area, function(a) {
      if (!identical(a, local) && !isTRUE(spec$comparable)) return("")
      v <- values[values$area == a, , drop = FALSE]
      if (nrow(v) == 0 || is.na(v$value[[1]])) return("—")
      paste0(planning_format_value(v$value[[1]], spec$digits), v$flag[[1]],
             if (has_significance) planning_significance_mark(v$significance[[1]]) else "")
    }, character(1))
    tibble::as_tibble(c(
      list(Tema = spec$theme, Indicador = planning_indicator_label(id, education_min_age),
           Unidade = spec$unit, `Período` = planning_period_label(id, year)),
      stats::setNames(as.list(cells), columns)
    ))
  })
  dplyr::bind_rows(rows)
}

# ---------------------------------------------------------
# Funnel plot
# ---------------------------------------------------------
# Each unit's value against the size of its denominator, with control limits
# around the benchmark: what a unit of that size would show by chance alone,
# 95% and 99.8% (two and three standard deviations). Points outside the limits
# differ by more than chance; small units scatter widely inside them. Only for
# indicators with a count model: Poisson for the event rates, binomial for
# the birth proportions.
PLANNING_FUNNEL_MODELS <- tibble::tribble(
  ~id,                    ~model,     ~multiplier, ~denominator_label,
  "birth_rate",           "poisson",  1000, "População média",
  "death_rate",           "poisson",  1000, "População média",
  "infant_rate",          "poisson",  1000, "Nados-vivos no triénio",
  "neonatal_rate",        "poisson",  1000, "Nados-vivos no triénio",
  "early_neonatal_rate",  "poisson",  1000, "Nados-vivos no triénio",
  "postneonatal_rate",    "poisson",  1000, "Nados-vivos no triénio",
  "late_fetal_rate",      "poisson",  1000, "Nados-vivos e fetos mortos no triénio",
  "perinatal_rate",       "poisson",  1000, "Nados-vivos e fetos mortos no triénio",
  "teen_births_pct",      "binomial", 100,  "Nados-vivos no triénio",
  "older_births_pct",     "binomial", 100,  "Nados-vivos no triénio",
  "preterm_pct",          "binomial", 100,  "Nascimentos com duração da gestação conhecida, no triénio",
  "low_birth_weight_pct", "binomial", 100,  "Nascimentos com peso conhecido, no triénio"
)

# Limits at `probability` (two-sided) for denominators `n` around the rate
# `reference` (in the indicator's unit): the exact quantiles of the count a unit
# of that size would show under the benchmark rate - Poisson for event rates,
# binomial for proportions - interpolated between integers so the lines are
# smooth (Spiegelhalter, Stat Med 2005). A count of zero in a small unit is then
# never "significantly low", which an interval around the expected count would
# wrongly make it.
planning_funnel_limits <- function(n, reference, model, multiplier, probability) {
  rate <- reference / multiplier
  quantile_of <- function(p) {
    if (identical(model, "poisson")) {
      lambda <- rate * n
      r <- stats::qpois(p, lambda)
      below <- stats::ppois(r - 1, lambda); at <- stats::ppois(r, lambda)
    } else {
      size <- round(n)
      r <- stats::qbinom(p, size, rate)
      below <- stats::pbinom(r - 1, size, rate); at <- stats::pbinom(r, size, rate)
    }
    alpha <- ifelse(at > below, (at - p) / (at - below), 0)
    pmax(r - alpha, 0) / n * multiplier
  }
  list(lower = quantile_of((1 - probability) / 2), upper = quantile_of(1 - (1 - probability) / 2))
}

# Units with their position against the limits.
planning_funnel_data <- function(table, spec, benchmark_area) {
  model <- PLANNING_FUNNEL_MODELS[PLANNING_FUNNEL_MODELS$id == spec$id, , drop = FALSE]
  if (nrow(model) == 0) return(NULL)
  reference <- table$value[table$area == benchmark_area]
  if (length(reference) != 1 || !is.finite(reference)) return(NULL)
  units <- table[!table$area %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL) & is.finite(table$value) &
                   is.finite(table$denominator) & table$denominator > 0, , drop = FALSE]
  l95 <- planning_funnel_limits(units$denominator, reference, model$model, model$multiplier, 0.95)
  l998 <- planning_funnel_limits(units$denominator, reference, model$model, model$multiplier, 0.998)
  units$position <- dplyr::case_when(
    units$value > l998$upper ~ "Acima do limite de 99,8%",
    units$value > l95$upper ~ "Acima do limite de 95%",
    units$value < l998$lower ~ "Abaixo do limite de 99,8%",
    units$value < l95$lower ~ "Abaixo do limite de 95%",
    TRUE ~ "Dentro dos limites"
  )
  attr(units, "reference") <- reference
  attr(units, "model") <- model
  units
}

# The same three colours; beyond 99.8% is told apart from beyond 95% by the
# marker's shape, not by a lighter step (a five-step ramp fails the
# normal-vision distance check against the grey).
PLANNING_FUNNEL_COLOURS <- c(
  "Acima do limite de 99,8%" = "#eb6834", "Acima do limite de 95%" = "#eb6834",
  "Dentro dos limites" = "#c3c2b7",
  "Abaixo do limite de 95%" = "#2a78d6", "Abaixo do limite de 99,8%" = "#2a78d6"
)

planning_funnel_chart <- function(units, spec, year, benchmark_area, highlight = character(0), unit_label = "Unidades", denominator_label = "Denominador") {
  reference <- attr(units, "reference")
  model <- attr(units, "model")
  grid <- exp(seq(log(max(min(units$denominator) * 0.8, 1)), log(max(units$denominator) * 1.1), length.out = 200))
  l95 <- planning_funnel_limits(grid, reference, model$model, model$multiplier, 0.95)
  l998 <- planning_funnel_limits(grid, reference, model$model, model$multiplier, 0.998)
  units$hover <- paste0("<b>", units$area, "</b><br>", planning_hover_value(units$value, units$lower, units$upper, spec$digits, units$flag),
                        "<br>", denominator_label, ": ", planning_format_value(units$denominator, 0), "<br>", units$position)
  units$highlighted <- units$area %in% highlight

  p <- plotly::plot_ly()
  for (bound in list(list(l998, "Limites de 99,8%", "dot"), list(l95, "Limites de 95%", "dash"))) {
    p <- plotly::add_lines(p, x = grid, y = bound[[1]]$upper, name = bound[[2]], legendgroup = bound[[2]],
                           line = list(color = PLANNING_INK$muted, dash = bound[[3]], width = 1), hoverinfo = "skip")
    p <- plotly::add_lines(p, x = grid, y = bound[[1]]$lower, name = bound[[2]], legendgroup = bound[[2]], showlegend = FALSE,
                           line = list(color = PLANNING_INK$muted, dash = bound[[3]], width = 1), hoverinfo = "skip")
  }
  p <- plotly::add_lines(p, x = range(grid), y = c(reference, reference), name = paste0(benchmark_area, ": ", planning_format_value(reference, spec$digits)),
                         line = list(color = PLANNING_LEVEL_COLOURS[["Portugal"]], width = 2), hoverinfo = "skip")
  for (position in names(PLANNING_FUNNEL_COLOURS)) {
    rows <- units[units$position == position, , drop = FALSE]
    if (nrow(rows) == 0) next
    p <- plotly::add_markers(p, data = rows, x = ~denominator, y = ~value, name = position, text = ~hover, hoverinfo = "text",
                             marker = list(color = PLANNING_FUNNEL_COLOURS[[position]], size = if (grepl("99,8", position)) 12 else 9,
                                           symbol = if (grepl("99,8", position)) "diamond" else "circle",
                                           line = list(color = ifelse(rows$highlighted, PLANNING_INK$primary, "#ffffff"),
                                                       width = ifelse(rows$highlighted, 2.5, 0.5))))
  }
  marked <- units[units$highlighted, , drop = FALSE]
  annotations <- lapply(seq_len(nrow(marked)), function(i) list(
    x = log10(marked$denominator[[i]]), y = marked$value[[i]], text = marked$area[[i]], showarrow = TRUE, arrowhead = 0,
    ax = 30, ay = -25, font = list(color = PLANNING_INK$primary, size = 11)
  ))
  planning_plotly_layout(
    p, annotations = annotations,
    xaxis = list(title = paste0(denominator_label, " (escala logarítmica)"), type = "log"),
    # The limits of the smallest units run far above the data; the axis
    # follows the points, and the lines leave through the top.
    yaxis = list(title = paste0(spec$label, " (", spec$unit, "), ", planning_period_label(spec$id, year)),
                 range = c(min(0, min(units$value)), max(units$value, reference) * 1.15)),
    margin = list(b = 90)
  )
}
