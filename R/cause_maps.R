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
  finite <- is.finite(out$value)
  out$quintile <- NA_integer_
  if (sum(finite) >= 5) {
    breaks <- unique(stats::quantile(out$value[finite], seq(0, 1, 0.2), names = FALSE))
    out$quintile[finite] <- if (length(breaks) > 2) as.integer(cut(out$value[finite], breaks, include.lowest = TRUE)) else 3L
  }
  out
}

cause_map_metric_label <- function(metric) names(CAUSE_MAP_METRICS)[match(metric, CAUSE_MAP_METRICS)]

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
  features <- lapply(split(p, factor(p$id, levels = unique(p$id))), function(area) {
    rings <- lapply(split(area, factor(area$piece, levels = unique(area$piece))), function(piece) {
      xy <- unname(as.matrix(piece[, c("x", "y")]))
      if (!all(xy[1, ] == xy[nrow(xy), ])) xy <- rbind(xy, xy[1, ])
      # Plotly's geo drawing (d3-geo, on the sphere) takes a ring wound
      # anticlockwise as everything outside it: exterior rings go clockwise.
      n <- nrow(xy)
      if (sum(xy[-n, 1] * xy[-1, 2] - xy[-1, 1] * xy[-n, 2]) > 0) xy <- xy[n:1, , drop = FALSE]
      list(xy)
    })
    list(type = "Feature", id = area$id[[1]], properties = list(name = area$name[[1]]),
         geometry = list(type = "MultiPolygon", coordinates = unname(rings)))
  })
  value <- list(type = "FeatureCollection", features = unname(features))
  assign(key, value, envir = .map_geojson_cache)
  value
}
.map_geojson_cache <- new.env(parent = emptyenv())
digest_or_names <- function(boundaries) paste(nrow(boundaries$polygons), utils::head(boundaries$polygons$id, 1), utils::tail(boundaries$polygons$id, 1))

# Plotly colour scale for three colours, from `low` through `mid` to `high`.
map_colorscale <- function(colours) list(list(0, colours[["low"]]), list(0.5, colours[["mid"]]), list(1, colours[["high"]]))

# A choropleth of `values` (columns area, index, quintile, hover) on `boundaries`.
#   scale "index":    log(index / 100) on the planning maps' range, 67 to 150,
#                     with `colours` (low, mid, high) - red worse by default
#   scale "quintile": the five quintile colours
# Areas without a value are grey and say so on hover. `highlight` outlines
# areas; `source` names the plotly source whose clicks the server reads.
interactive_choropleth <- function(values, boundaries, scale = "index", colours = PLANNING_MAP_DIVERGING,
                                   title = NULL, legend_title = "", source = "map", highlight = character(0),
                                   height = NULL) {
  geo <- map_geojson(boundaries)
  units <- unique(boundaries$polygons$id)
  v <- values[values$area %in% units, , drop = FALSE]
  missing <- setdiff(units, v$area[is.finite(if (identical(scale, "quintile")) v$quintile else v$index)])
  shown <- v[!v$area %in% missing, , drop = FALSE]
  hover_missing <- values$hover[match(missing, values$area)]
  hover_missing[is.na(hover_missing)] <- paste0("<b>", missing[is.na(hover_missing)], "</b><br>Sem valor")

  p <- plotly::plot_ly(source = source, height = height)
  if (length(missing) > 0) {
    p <- plotly::add_trace(
      p, type = "choropleth", geojson = geo, featureidkey = "id", locations = missing, z = rep(0, length(missing)),
      colorscale = list(list(0, "#d9d8d3"), list(1, "#d9d8d3")), showscale = FALSE,
      text = hover_missing, hoverinfo = "text", marker = list(line = list(color = "#ffffff", width = 0.4)),
      customdata = missing, name = "Sem valor"
    )
  }
  if (nrow(shown) > 0) {
    if (identical(scale, "quintile")) {
      steps <- PLANNING_MAP_QUINTILES
      stepped <- unlist(lapply(seq_along(steps), function(i) list(list((i - 1) / 5, steps[[i]]), list(i / 5, steps[[i]]))), recursive = FALSE)
      p <- plotly::add_trace(
        p, type = "choropleth", geojson = geo, featureidkey = "id", locations = shown$area, z = shown$quintile,
        zmin = 0.5, zmax = 5.5, colorscale = stepped, text = shown$hover, hoverinfo = "text", customdata = shown$area,
        marker = list(line = list(color = "#ffffff", width = 0.4)),
        colorbar = list(title = list(text = "Quintil"), tickvals = 1:5, ticktext = c("1 (mais baixo)", "2", "3", "4", "5 (mais alto)"),
                        len = 0.6, thickness = 14)
      )
    } else {
      limits <- log(PLANNING_MAP_INDEX_LIMITS)
      z <- pmin(pmax(log(shown$index / 100), limits[1]), limits[2])
      ticks <- c(67, 80, 100, 125, 150)
      p <- plotly::add_trace(
        p, type = "choropleth", geojson = geo, featureidkey = "id", locations = shown$area, z = z,
        zmin = limits[1], zmax = limits[2], colorscale = map_colorscale(colours), text = shown$hover, hoverinfo = "text",
        customdata = shown$area, marker = list(line = list(color = "#ffffff", width = 0.4)),
        colorbar = list(title = list(text = legend_title), tickvals = log(ticks / 100), ticktext = ticks, len = 0.6, thickness = 14)
      )
    }
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
    path <- do.call(rbind, lapply(split(marked, marked$piece), function(f) rbind(f[, c("x", "y")], f[1, c("x", "y")], data.frame(x = NA, y = NA))))
    p <- plotly::add_trace(p, type = "scattergeo", mode = "lines", lon = path$x, lat = path$y,
                           line = list(color = "#0b0b0b", width = 2), hoverinfo = "skip", showlegend = FALSE)
  }
  p <- plotly::layout(
    p,
    title = if (is.null(title)) NULL else list(text = title, font = list(size = 13), x = 0.5),
    geo = list(fitbounds = "locations", visible = FALSE, projection = list(type = "mercator"), bgcolor = "rgba(0,0,0,0)"),
    margin = list(l = 0, r = 0, t = if (is.null(title)) 10 else 60, b = 0),
    paper_bgcolor = "rgba(0,0,0,0)", showlegend = FALSE,
    annotations = list(list(text = boundaries$meta$attribution %||% "", showarrow = FALSE, xref = "paper", yref = "paper",
                            x = 1, y = 0, xanchor = "right", yanchor = "bottom", font = list(size = 9, color = "#52514e")))
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
