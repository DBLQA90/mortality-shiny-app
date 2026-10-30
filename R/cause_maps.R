# =========================================================
# Maps: any cause, sex and metric, drawn interactively
# =========================================================
#
# The Mapas tab maps one cause of death (the short list or INE 0008461's
# detailed tumours) for every area of a level, over a triennium: the SMR against
# Portugal, the standardised rate at all ages or under 75. The values come from
# the planning tab's standardised module (planning_cause_standardised(), with
# the cause as an extra measure), so a map and the profile of an area agree.
#
# interactive_choropleth() draws any such table as a plotly choropleth on the
# app's own boundaries (data/boundaries, from geopt), one trace per map: hover
# for the values, zoom and pan, and a click returns the area's name. The
# planning tab's maps use it too; their PNG stays the ggplot drawing
# (planning_map_draw()).

CAUSE_MAP_METRICS <- c(
  "SMR face a Portugal (= 100)" = "smr",
  "Taxa padronizada (por 100 000)" = "dsr",
  "Taxa padronizada, menos de 75 anos (por 100 000)" = "dsr75"
)

# End years of the triennia a cause can be mapped for: all three years need
# deaths by cause, and the detailed tumours exist only in 2002-2022.
cause_map_end_years <- function(cause) {
  years <- planning_cause_years()
  if (is_detailed_tumour_cause(cause)) years <- intersect(years, DETAILED_TUMOUR_YEARS)
  years <- sort(unique(years))
  years[(years - 1L) %in% years & (years - 2L) %in% years]
}

# One row per area: the value, its interval, observed and expected deaths, the
# index against the benchmark (100 = benchmark), significance and quintile.
cause_map_values <- function(cause, sex, end_year, metric, level, lookup, benchmark = "Portugal",
                             basis = default_population_basis(), smoothed = FALSE) {
  units <- planning_map_units(level, lookup)
  # The ULS that share a municipality are drawn with their own parishes.
  mode <- if (identical(level, "uls")) "parish" else PLANNING_DEFAULT_SPLIT_MODE
  table <- planning_cause_standardised(c(units, benchmark), end_year, lookup = lookup, benchmark = benchmark, sex = sex,
                                       mode = mode, basis = basis, causes = cause)
  if (nrow(table) == 0) return(tibble::tibble())
  pick <- function(suffix = "") table[[paste0(metric, suffix)]]
  value <- pick(); lower <- pick("_lower"); upper <- pick("_upper")
  reference <- value[table$area == benchmark][1]
  method <- "observado"
  if (isTRUE(smoothed) && identical(metric, "smr")) {
    units_tab <- data.frame(value = table$smr, denominator = table$expected)
    ok <- is.finite(units_tab$value) & is.finite(units_tab$denominator) & units_tab$denominator > 0 & table$area %in% units
    if (sum(ok) >= 3) {
      value[ok] <- planning_eb_smooth(units_tab[ok, , drop = FALSE], list(model = "poisson", multiplier = 100))$smoothed
      method <- "suavizado"
    }
  }
  significance <- planning_significance(lower, upper, rep(reference, length(lower)))
  out <- tibble::tibble(
    area = table$area, value = value, observed_value = pick(), lower = lower, upper = upper,
    observed = table$observed, expected = table$expected, flag = table$flag,
    benchmark_value = reference,
    index = ifelse(is.finite(value) & is.finite(reference) & reference > 0, value / reference * 100, NA_real_),
    significance = significance, method = method
  )
  out <- out[out$area %in% units, , drop = FALSE]
  out$quintile <- map_quintile(out$value)
  out$class <- map_significance_class(out$significance)
  out
}

cause_map_metric_label <- function(metric) names(CAUSE_MAP_METRICS)[match(metric, CAUSE_MAP_METRICS)]

# ---------------------------------------------------------
# Significance, trend and change between two triennia
# ---------------------------------------------------------

# Colours of three classes, lowest first, by what a higher value means:
# mortality ("worse"), an indicator where higher is better, or neither.
map_class_colours <- function(labels, direction = "worse") {
  pal <- switch(direction,
                better = c(PLANNING_MAP_DIVERGING[["high"]], "#e6e4dd", PLANNING_MAP_DIVERGING[["low"]]),
                neutral = c(PLANNING_MAP_NEUTRAL[["low"]], "#e6e4dd", PLANNING_MAP_NEUTRAL[["high"]]),
                c(PLANNING_MAP_DIVERGING[["low"]], "#e6e4dd", PLANNING_MAP_DIVERGING[["high"]]))
  stats::setNames(pal, labels)
}

# The class of each area against the benchmark, from planning_significance().
MAP_SIGNIFICANCE_LABELS <- c(Inferior = "\u25bc abaixo de Portugal", Semelhante = "= semelhante", Superior = "\u25b2 acima de Portugal")
map_significance_class <- function(significance) unname(MAP_SIGNIFICANCE_LABELS[significance])

# Non-overlapping triennia ending at `end_year`, `n` of them, oldest first;
# only those the cause can be mapped for.
cause_map_triennia <- function(cause, end_year, n) {
  ends <- end_year - 3L * (seq_len(n) - 1L)
  sort(intersect(ends, cause_map_end_years(cause)))
}

MAP_TREND_LABELS <- c("\u25bc a descer", "= sem tendência significativa", "\u25b2 a subir")
MAP_CHANGE_LABELS <- c("\u25bc desceu", "= sem mudança significativa", "\u25b2 subiu")
MAP_TREND_LIMIT <- 5        # % per year at either end of the colour scale
MAP_PROJECTION_LIMIT <- 20  # % change to the target year, likewise

# Annual percentage change of each area over `n` non-overlapping triennia ending
# at `end_year` (trend_apc(): log-linear, each triennium weighted by its
# precision, as the planning tab's trend), with its interval and class.
cause_map_trend <- function(cause, sex, end_year, n, metric, level, lookup, benchmark = "Portugal",
                            basis = default_population_basis(), progress = NULL) {
  ends <- cause_map_triennia(cause, end_year, n)
  if (length(ends) < 4) return(tibble::tibble())
  periods <- dplyr::bind_rows(lapply(seq_along(ends), function(i) {
    if (is.function(progress)) progress(i, length(ends))
    v <- cause_map_values(cause, sex, ends[[i]], metric, level, lookup, benchmark, basis)
    if (nrow(v) == 0) return(NULL)
    dplyr::mutate(v, end_year = ends[[i]])
  }))
  if (nrow(periods) == 0) return(tibble::tibble())
  rows <- lapply(split(periods, periods$area), function(a) {
    a <- a[order(a$end_year), , drop = FALSE]
    trend <- trend_apc(a$end_year - 1L, a$value, a$lower, a$upper, min_points = 4L)
    first <- a$value[[1]]; last <- a$value[[nrow(a)]]
    if (is.null(trend)) {
      return(tibble::tibble(area = a$area[[1]], z = NA_real_, apc_lower = NA_real_, apc_upper = NA_real_, n = sum(is.finite(a$value) & a$value > 0),
                            first = first, last = last))
    }
    tibble::tibble(area = a$area[[1]], z = trend$aapc$apc, apc_lower = trend$aapc$lower, apc_upper = trend$aapc$upper, n = trend$n,
                   first = first, last = last)
  })
  out <- dplyr::bind_rows(rows)
  out$class <- ifelse(!is.finite(out$z), NA_character_,
                      ifelse(out$apc_lower > 0, MAP_TREND_LABELS[[3]], ifelse(out$apc_upper < 0, MAP_TREND_LABELS[[1]], MAP_TREND_LABELS[[2]])))
  out$quintile <- map_quintile(out$z)
  attr(out, "periods") <- ends
  out
}

# The metric in triennium `end_year` against triennium `base_year` (= 100): the
# ratio, with an interval from the two intervals on the log scale, and its class.
cause_map_change <- function(cause, sex, end_year, base_year, metric, level, lookup, benchmark = "Portugal",
                             basis = default_population_basis()) {
  now <- cause_map_values(cause, sex, end_year, metric, level, lookup, benchmark, basis)
  before <- cause_map_values(cause, sex, base_year, metric, level, lookup, benchmark, basis)
  if (nrow(now) == 0 || nrow(before) == 0) return(tibble::tibble())
  both <- dplyr::inner_join(now[, c("area", "value", "lower", "upper")], before[, c("area", "value", "lower", "upper")],
                            by = "area", suffix = c("", "_base"))
  z <- stats::qnorm(0.975)
  se <- function(lo, hi) ifelse(lo > 0 & hi > 0, (log(hi) - log(lo)) / (2 * z), NA_real_)
  ok <- is.finite(both$value) & is.finite(both$value_base) & both$value > 0 & both$value_base > 0
  log_ratio <- ifelse(ok, log(both$value / both$value_base), NA_real_)
  spread <- sqrt(se(both$lower, both$upper)^2 + se(both$lower_base, both$upper_base)^2)
  out <- tibble::tibble(
    area = both$area, value = both$value, value_base = both$value_base,
    index = exp(log_ratio) * 100, change_lower = exp(log_ratio - z * spread) * 100, change_upper = exp(log_ratio + z * spread) * 100
  )
  out$class <- ifelse(!is.finite(out$index) | !is.finite(out$change_lower), NA_character_,
                      ifelse(out$change_lower > 100, MAP_CHANGE_LABELS[[3]], ifelse(out$change_upper < 100, MAP_CHANGE_LABELS[[1]], MAP_CHANGE_LABELS[[2]])))
  out$quintile <- map_quintile(out$index)
  out
}

map_quintile <- function(x) {
  out <- rep(NA_integer_, length(x))
  finite <- is.finite(x)
  if (sum(finite) >= 5) {
    breaks <- unique(stats::quantile(x[finite], seq(0, 1, 0.2), names = FALSE))
    out[finite] <- if (length(breaks) > 2) as.integer(cut(x[finite], breaks, include.lowest = TRUE)) else 3L
  }
  out
}

cause_map_trend_hover <- function(values, metric) {
  digits <- if (identical(metric, "smr")) 0 else 1
  pct <- function(x) formatC(x, format = "f", digits = 1, decimal.mark = ",", flag = "+")
  paste0("<b>", values$area, "</b>",
         ifelse(is.finite(values$z),
                paste0("<br>", pct(values$z), "% por ano (IC 95% ", pct(values$apc_lower), " a ", pct(values$apc_upper), ")",
                       "<br>", values$class,
                       "<br>Primeiro triénio: ", planning_format_value(values$first, digits),
                       " · último: ", planning_format_value(values$last, digits), " (", values$n, " triénios)"),
                "<br>Sem tendência (menos de 4 triénios com valor)"))
}

cause_map_change_hover <- function(values, metric, base_label, end_label) {
  digits <- if (identical(metric, "smr")) 0 else 1
  pct <- function(x) formatC(x - 100, format = "f", digits = 0, decimal.mark = ",", flag = "+")
  paste0("<b>", values$area, "</b>",
         ifelse(is.finite(values$index),
                paste0("<br>", base_label, ": ", planning_format_value(values$value_base, digits),
                       "<br>", end_label, ": ", planning_format_value(values$value, digits),
                       "<br>Variação: ", pct(values$index), "% (IC 95% ", pct(values$change_lower), "% a ", pct(values$change_upper), "%)",
                       ifelse(is.na(values$class), "", paste0("<br>", values$class))),
                "<br>Sem valor num dos triénios"))
}

cause_map_title <- function(cause, sex, end_year, metric, smoothed = FALSE) {
  sex_label <- c(HM = "ambos os sexos", H = "homens", M = "mulheres")[[sex]]
  paste0(cause, " · ", sex_label, " · ", end_year - 2L, "-", end_year, "<br><sup>",
         cause_map_metric_label(metric), if (isTRUE(smoothed)) ", suavizada" else "", "</sup>")
}

# ---------------------------------------------------------
# The interactive drawing
# ---------------------------------------------------------

# The boundaries as a GeoJSON FeatureCollection keyed by area name, one feature
# per area (its pieces a MultiPolygon). The tables are already in degrees, with
# the islands moved beside the mainland, so they serve as longitude/latitude.
map_geojson <- function(boundaries) {
  key <- paste("geojson", digest_or_names(boundaries))
  if (exists(key, envir = .map_geojson_cache, inherits = FALSE)) return(get(key, envir = .map_geojson_cache))
  p <- boundaries$polygons
  # Plotly's geo drawing (d3-geo, on the sphere) takes a ring wound
  # anticlockwise as everything outside it: exterior rings go clockwise, holes
  # anticlockwise. Each piece is one polygon, its first ring the exterior.
  ring <- function(r, exterior) {
    xy <- unname(as.matrix(r[, c("x", "y")]))
    if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
    n <- nrow(xy)
    anticlockwise <- sum(xy[-n, 1] * xy[-1, 2] - xy[-1, 1] * xy[-n, 2]) > 0
    if (anticlockwise == exterior) xy <- xy[n:1, , drop = FALSE]
    xy
  }
  features <- lapply(split(p, factor(p$id, levels = unique(p$id))), function(area) {
    polygons <- lapply(split(area, factor(area$piece, levels = unique(area$piece))), function(piece) {
      rings <- split(piece, factor(piece$hole, levels = unique(piece$hole)))
      unname(lapply(seq_along(rings), function(k) ring(rings[[k]], exterior = k == 1L)))
    })
    list(type = "Feature", id = area$id[[1]], properties = list(name = area$name[[1]]),
         geometry = list(type = "MultiPolygon", coordinates = unname(polygons)))
  })
  value <- list(type = "FeatureCollection", features = unname(features))
  assign(key, value, envir = .map_geojson_cache)
  value
}
.map_geojson_cache <- new.env(parent = emptyenv())
digest_or_names <- function(boundaries) paste(nrow(boundaries$polygons), utils::head(boundaries$polygons$id, 1), utils::tail(boundaries$polygons$id, 1))

# Plotly colour scale for three colours, from `low` through `mid` to `high`.
map_colorscale <- function(colours) list(list(0, colours[["low"]]), list(0.5, colours[["mid"]]), list(1, colours[["high"]]))

# Grey for areas without a value: darker than the scales' near-white midpoint,
# which it was too close to (#d9d8d3 against #f0efec).
MAP_NA_COLOUR <- "#b4b3ad"

# A choropleth of `values` (columns area, hover, and the scale's column) on
# `boundaries`. `scale`:
#   "index":    log(index / 100) on the planning maps' range, 67 to 150, with
#               `colours` (low, mid, high) - red worse by default
#   "quintile": the five quintile colours
#   "classes":  column `class`, one colour per class (`classes`, named colours
#               in legend order)
#   "linear":   column `z` on `limits`, with `colours`, ticks `ticks`
# Areas without a value are grey and say so on hover. `highlight` outlines
# areas; `source` names the plotly source whose clicks the server reads.
# `legend` "bottom" lays the colour bar under the map, for maps side by side.
interactive_choropleth <- function(values, boundaries, scale = "index", colours = PLANNING_MAP_DIVERGING,
                                   title = NULL, legend_title = "", source = "map", highlight = character(0),
                                   height = NULL, classes = NULL, limits = NULL, ticks = NULL, tick_labels = NULL,
                                   legend = "right") {
  # Levels whose island units have no id (ULS, ACES) are drawn as the mainland
  # alone: their inset boxes would stay empty.
  if (any(is.na(boundaries$polygons$id))) {
    boundaries$polygons <- boundaries$polygons[!is.na(boundaries$polygons$id), , drop = FALSE]
    boundaries$frames <- boundaries$frames[0, , drop = FALSE]
  }
  geo <- map_geojson(boundaries)
  units <- unique(boundaries$polygons$id)
  v <- values[values$area %in% units, , drop = FALSE]
  has_value <- switch(scale, quintile = is.finite(v$quintile), classes = !is.na(v$class) & v$class %in% names(classes),
                      linear = is.finite(v$z), is.finite(v$index))
  missing <- setdiff(units, v$area[has_value])
  shown <- v[!v$area %in% missing, , drop = FALSE]
  hover_missing <- values$hover[match(missing, values$area)]
  hover_missing[is.na(hover_missing)] <- paste0("<b>", missing[is.na(hover_missing)], "</b><br>Sem valor")

  bar <- function(title_text, tickvals, ticktext) {
    common <- list(title = list(text = title_text), tickvals = tickvals, ticktext = ticktext, thickness = 12)
    if (identical(legend, "bottom")) {
      common$title$side <- "top"
      c(common, list(orientation = "h", x = 0.5, xanchor = "center", y = 0.16, yanchor = "top", len = 0.9,
                     tickangle = 0, tickfont = list(size = 10)))
    } else {
      c(common, list(len = 0.6))
    }
  }
  stepped_scale <- function(cols) {
    n <- length(cols)
    unlist(lapply(seq_len(n), function(i) list(list((i - 1) / n, cols[[i]]), list(i / n, cols[[i]]))), recursive = FALSE)
  }
  outline <- list(line = list(color = "#ffffff", width = 0.4))

  p <- plotly::plot_ly(source = source, height = height)
  if (length(missing) > 0) {
    p <- plotly::add_trace(
      p, type = "choropleth", geojson = geo, featureidkey = "id", locations = missing, z = rep(0, length(missing)),
      colorscale = list(list(0, MAP_NA_COLOUR), list(1, MAP_NA_COLOUR)), showscale = FALSE,
      text = hover_missing, hoverinfo = "text", marker = outline, customdata = missing, name = "Sem valor"
    )
  }
  if (nrow(shown) > 0) {
    trace <- switch(
      scale,
      quintile = list(z = shown$quintile, zmin = 0.5, zmax = 5.5, colorscale = stepped_scale(PLANNING_MAP_QUINTILES),
                      colorbar = bar("Quintil", 1:5, c("1 (mais baixo)", "2", "3", "4", "5 (mais alto)"))),
      classes = list(z = match(shown$class, names(classes)), zmin = 0.5, zmax = length(classes) + 0.5,
                     colorscale = stepped_scale(unname(classes)),
                     colorbar = bar(legend_title, seq_along(classes), names(classes))),
      linear = list(z = pmin(pmax(shown$z, limits[1]), limits[2]), zmin = limits[1], zmax = limits[2],
                    colorscale = map_colorscale(colours), colorbar = bar(legend_title, ticks, tick_labels %||% ticks)),
      {
        lim <- log(PLANNING_MAP_INDEX_LIMITS)
        index_ticks <- c(67, 80, 100, 125, 150)
        list(z = pmin(pmax(log(shown$index / 100), lim[1]), lim[2]), zmin = lim[1], zmax = lim[2],
             colorscale = map_colorscale(colours), colorbar = bar(legend_title, log(index_ticks / 100), index_ticks))
      }
    )
    p <- do.call(plotly::add_trace, c(list(
      p, type = "choropleth", geojson = geo, featureidkey = "id", locations = shown$area,
      text = shown$hover, hoverinfo = "text", customdata = shown$area, marker = outline
    ), trace))
  }
  # Frames around the island insets, and the outline of the chosen areas.
  frames <- boundaries$frames
  if (!is.null(frames) && nrow(frames) > 0) {
    path <- do.call(rbind, lapply(split(frames, frames$group), function(f) rbind(f[, c("x", "y")], data.frame(x = NA, y = NA))))
    p <- plotly::add_trace(p, type = "scattergeo", mode = "lines", lon = path$x, lat = path$y,
                           line = list(color = "#898781", width = 1), hoverinfo = "skip", showlegend = FALSE)
  }
  marked <- boundaries$polygons[boundaries$polygons$id %in% highlight, , drop = FALSE]
  if (nrow(marked) > 0) {
    path <- do.call(rbind, lapply(split(marked, marked$hole), function(f) rbind(f[, c("x", "y")], f[1, c("x", "y")], data.frame(x = NA, y = NA))))
    p <- plotly::add_trace(p, type = "scattergeo", mode = "lines", lon = path$x, lat = path$y,
                           line = list(color = "#0b0b0b", width = 2), hoverinfo = "skip", showlegend = FALSE)
  }
  p <- plotly::layout(
    p,
    title = if (is.null(title)) NULL else list(text = title, font = list(size = 13), x = 0.5),
    geo = list(fitbounds = "locations", visible = FALSE, projection = list(type = "mercator"), bgcolor = "rgba(0,0,0,0)",
               domain = if (identical(legend, "bottom")) list(x = c(0, 1), y = c(0.24, 1)) else NULL),
    margin = list(l = 0, r = 0, t = if (is.null(title)) 10 else 64, b = 0),
    # Vertical, at the right edge: across the top it covered long titles.
    modebar = list(orientation = "v"),
    paper_bgcolor = "rgba(0,0,0,0)", showlegend = FALSE,
    annotations = list(list(text = boundaries$meta$attribution %||% "", showarrow = FALSE, xref = "paper", yref = "paper",
                            x = 1, y = 0, xanchor = "right", yanchor = "bottom",
                            font = list(size = 9, color = "#52514e")))
  )
  p <- plotly::event_register(p, "plotly_click")
  # Plotly fetches a world topojson for every geo map, drawn or not; an empty
  # one served by the app (map_topojson_url()) keeps the maps off the internet.
  plotly::config(p, displaylogo = FALSE, scrollZoom = TRUE, topojsonURL = map_topojson_url(),
                 modeBarButtonsToRemove = list("lasso2d", "select2d", "toggleHover", "hoverClosestGeo"))
}

# Where the empty world topojson is served: data/boundaries/topojson, registered
# as a resource path by the app (shiny::addResourcePath()); outside Shiny the
# folder itself.
map_topojson_url <- function() {
  if ("map-topojson" %in% names(shiny::resourcePaths())) return("map-topojson/")
  paste0("file://", normalizePath(file.path(app_dir_or_wd(), "data", "boundaries", "topojson")), "/")
}

# Hover text for a cause map.
cause_map_hover <- function(values, metric, benchmark = "Portugal") {
  digits <- if (identical(metric, "smr")) 0 else 1
  fmt <- function(x) planning_format_value(x, digits)
  sig <- c(Superior = paste("acima de", benchmark), Inferior = paste("abaixo de", benchmark), Semelhante = paste("semelhante a", benchmark))
  paste0(
    "<b>", values$area, "</b>",
    ifelse(is.finite(values$value),
           paste0("<br>", ifelse(values$method == "suavizado", "Suavizado: ", ""), fmt(values$value), values$flag,
                  ifelse(values$method == "suavizado", paste0(" (observado ", fmt(values$observed_value), ")"), ""),
                  "<br>IC 95%: ", fmt(values$lower), " a ", fmt(values$upper),
                  ifelse(is.na(values$significance), "", paste0("<br>", planning_significance_mark(values$significance), " ", sig[values$significance])),
                  "<br>Óbitos: ", planning_format_value(values$observed, 0),
                  ifelse(is.finite(values$expected), paste0(" (esperados ", planning_format_value(values$expected, 1), ")"), ""),
                  if (!identical(metric, "smr")) paste0("<br>", benchmark, ": ", fmt(values$benchmark_value)) else ""),
           "<br>Sem valor")
  )
}

# Hover text for a planning map (planning_map_values()).
planning_map_hover <- function(values, benchmark = "Portugal") {
  paste0(
    "<b>", values$area, "</b>",
    ifelse(is.finite(values$value),
           paste0("<br>", planning_format_value(values$value, values$digits), ifelse(is.na(values$flag), "", values$flag),
                  "<br>", benchmark, ": ", planning_format_value(values$benchmark_value, values$digits),
                  "<br>Índice (", benchmark, " = 100): ", planning_format_value(values$index, 0),
                  ifelse(is.na(values$quintile), "", paste0("<br>Quintil: ", values$quintile))),
           "<br>Sem valor")
  )
}

cause_map_projection_hover <- function(values, target) {
  pct <- function(x) formatC(x - 100, format = "f", digits = 0, decimal.mark = ",", flag = "+")
  paste0("<b>", values$area, "</b>",
         ifelse(is.finite(values$value),
                paste0("<br>", values$last_year, ": ", planning_format_value(values$last_value, 1),
                       "<br>", target, ": ", planning_format_value(values$value, 1),
                       " (", planning_format_value(values$lower, 1), " a ", planning_format_value(values$upper, 1), ")",
                       "<br>Variação: ", pct(values$index), "%",
                       ifelse(is.na(values$class), "", paste0("<br>", values$class)),
                       "<br>Portugal em ", target, ": ", planning_format_value(values$benchmark_value, 1)),
                "<br>Sem projecção"))
}
