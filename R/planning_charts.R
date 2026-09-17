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
planning_ranking_chart <- function(ranking, spec, year, highlight = character(0)) {
  reference <- ranking$value[ranking$area == "Portugal"]
  units <- ranking %>%
    dplyr::filter(.data$area != "Portugal", !is.na(.data$value)) %>%
    dplyr::arrange(.data$value) %>%
    dplyr::mutate(
      label = factor(.data$area, levels = .data$area),
      colour = ifelse(.data$area %in% highlight, PLANNING_LEVEL_COLOURS[["Local"]], "#c3c2b7"),
      hover = paste0("<b>", .data$area, "</b><br>", planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag))
    )
  error <- if (any(!is.na(units$lower))) {
    list(type = "data", symmetric = FALSE, array = units$upper - units$value, arrayminus = units$value - units$lower,
         color = PLANNING_INK$muted, thickness = 1, width = 2)
  } else {
    NULL
  }
  p <- plotly::plot_ly(
    units, y = ~label, x = ~value, type = "bar", orientation = "h",
    marker = list(color = ~colour), error_x = error, text = ~hover, hoverinfo = "text", textposition = "none"
  )
  shapes <- if (length(reference) == 1 && is.finite(reference)) {
    list(list(type = "line", x0 = reference, x1 = reference, yref = "paper", y0 = 0, y1 = 1,
              line = list(color = PLANNING_LEVEL_COLOURS[["Portugal"]], dash = "dash", width = 2)))
  } else {
    list()
  }
  annotations <- if (length(shapes) > 0) {
    list(list(x = reference, y = 1.01, yref = "paper", text = paste0("Portugal: ", planning_format_value(reference, spec$digits)),
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
        planning_hover_value(.data$value, .data$lower, .data$upper, spec$digits, .data$flag)
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
planning_summary_display <- function(table, areas) {
  local <- areas$area[[1]]
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
      paste0(planning_format_value(v$value[[1]], spec$digits), v$flag[[1]])
    }, character(1))
    tibble::as_tibble(c(
      list(Tema = spec$theme, Indicador = paste0(spec$label, " [", spec$ref, "]"),
           Unidade = spec$unit, `Período` = planning_period_label(id, year)),
      stats::setNames(as.list(cells), columns)
    ))
  })
  dplyr::bind_rows(rows)
}
