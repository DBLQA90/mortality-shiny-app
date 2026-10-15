# =========================================================
# Maps of Portugal: several indicators side by side
# =========================================================
#
# Small multiples of the country by municipality, NUTS III, NUTS II or ULS, one
# map per indicator. By default every indicator is shown against the benchmark
# (Portugal = 100) on one common diverging scale - blue below, red above, grey
# at Portugal - so maps of different indicators, in different units, can be read
# together as a heatmap: where several are high at once. The alternative is
# each indicator's own quintiles, which shows its internal ranking whatever its
# spread.
#
# Above Portugal is not worse: a high life expectancy and a high mortality rate
# are both red. The note under the maps says so.
#
# The boundaries are plain coordinate tables built once by
# tools/build_boundaries.R from the geopt package (CAOP, DGT); drawing them
# needs only ggplot2.

PLANNING_MAP_LEVELS <- c("Municípios" = "municipio", "NUTS III" = "nuts3", "ULS" = "uls", "NUTS II" = "nuts2")
PLANNING_MAP_MAX_INDICATORS <- 4L
# The diverging pair and neutral midpoint of the chart palette; the index is
# drawn on a log scale so that 150 and 67 sit at the same distance from 100.
PLANNING_MAP_DIVERGING <- c(low = "#256abf", mid = "#f0efec", high = "#e34948")
PLANNING_MAP_INDEX_LIMITS <- c(2 / 3, 3 / 2)
PLANNING_MAP_QUINTILES <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")

planning_map_boundaries <- function(level, vintage = default_nuts_vintage()) {
  file <- if (level %in% c("nuts2", "nuts3")) paste0(level, "_", normalize_nuts_vintage(vintage)) else level
  key <- paste("boundaries", file, sep = "|")
  if (!exists(key, envir = planning_cache, inherits = FALSE)) {
    path <- file.path(app_dir_or_wd(), "data", "boundaries", paste0(file, ".rds"))
    assign(key, if (file.exists(path)) readRDS(path) else NULL, envir = planning_cache)
  }
  get(key, envir = planning_cache, inherits = FALSE)
}

# The areas each level is drawn with, as the app names them.
planning_map_units <- function(level, lookup) {
  switch(
    level,
    municipio = sort(unique(as.character(lookup$municipality))),
    nuts3 = sort(unique(as.character(lookup$nuts3))),
    nuts2 = sort(unique(as.character(lookup$nuts2))),
    uls = planning_uls_units("units")
  )
}

# The comparable indicators: a count only measures size, so it is not mapped.
planning_map_indicator_choices <- function() {
  choices <- planning_indicator_choices()
  comparable <- PLANNING_INDICATORS$id[PLANNING_INDICATORS$comparable]
  choices <- lapply(choices, function(group) group[unlist(group) %in% comparable])
  choices[lengths(choices) > 0]
}

# The latest year of each indicator up to `max_year`. The rebuilt latest year
# exists only for whole regions, so for life expectancy and standardised rates
# the smaller levels keep the last year INE published by cause.
planning_map_year <- function(id, max_year, level) {
  years <- planning_indicator_years(id)
  if (id %in% c(standardised_ids, life_expectancy_ids) && level != "nuts2") years <- setdiff(years, planning_calendar_years())
  years <- years[years <= max_year]
  if (length(years) == 0) NA_integer_ else max(years)
}

# One row per indicator and area: the value shown, its index against the
# benchmark (100 = benchmark) and its quintile among the areas.
planning_map_values <- function(ids, level, lookup, max_year, benchmark = "Portugal", smoothed = FALSE,
                                mode = PLANNING_DEFAULT_SPLIT_MODE, basis = default_population_basis()) {
  units <- planning_map_units(level, lookup)
  # Each of the six ULS that share a municipality is drawn with its own
  # parishes, so its values must be its share of the municipality too.
  if (identical(level, "uls")) mode <- "parish"
  dplyr::bind_rows(lapply(ids, function(id) {
    spec <- planning_indicator_spec(id)
    year <- planning_map_year(id, max_year, level)
    if (is.na(year)) return(NULL)
    table <- planning_indicator_table(c(benchmark, units), year, ids = id, lookup = lookup, benchmark = benchmark,
                                      mode = mode, basis = basis)
    reference <- table$value[table$area == benchmark][1]
    rows <- table[table$area %in% units, , drop = FALSE]
    method <- "observado"
    model <- PLANNING_FUNNEL_MODELS[PLANNING_FUNNEL_MODELS$id == id, , drop = FALSE]
    rows$observed <- rows$value
    if (isTRUE(smoothed) && nrow(model) == 1) {
      ok <- is.finite(rows$value) & is.finite(rows$denominator) & rows$denominator > 0
      if (sum(ok) >= 3) {
        rows$value[ok] <- planning_eb_smooth(rows[ok, , drop = FALSE], model)$smoothed
        method <- "suavizado"
      }
    }
    finite <- is.finite(rows$value)
    quintile <- rep(NA_integer_, nrow(rows))
    if (sum(finite) >= 5) {
      breaks <- unique(stats::quantile(rows$value[finite], seq(0, 1, 0.2), names = FALSE))
      quintile[finite] <- if (length(breaks) > 2) as.integer(cut(rows$value[finite], breaks, include.lowest = TRUE)) else 3L
    }
    tibble::tibble(
      indicator = id,
      panel = paste0(paste(strwrap(spec$label, 46), collapse = "\n"), "\n", spec$unit, ", ", planning_period_label(id, year),
                     if (identical(method, "suavizado")) ", suavizado" else ""),
      area = rows$area,
      value = rows$value,
      observed = rows$observed,
      flag = rows$flag,
      digits = spec$digits,
      benchmark_value = reference,
      index = ifelse(is.finite(rows$value) & is.finite(reference) & reference > 0, rows$value / reference * 100, NA_real_),
      quintile = quintile
    )
  }))
}

# Small multiples, two per row. `scale`: "index" (Portugal = 100, common to
# every map) or "quintile" (each indicator's own fifths).
planning_map_chart <- function(values, boundaries, scale = "index", highlight = character(0), benchmark = "Portugal") {
  polygons <- boundaries$polygons
  panels <- unique(values$panel)
  drawn <- dplyr::bind_rows(lapply(panels, function(p) {
    v <- values[values$panel == p, c("area", "index", "quintile"), drop = FALSE]
    out <- dplyr::left_join(polygons, v, by = c("id" = "area"))
    out$panel <- p
    out
  }))
  drawn$panel <- factor(drawn$panel, levels = panels)
  frames <- boundaries$frames
  if (identical(scale, "quintile")) {
    drawn$quintile <- factor(drawn$quintile, levels = 1:5)
    plot <- ggplot2::ggplot(drawn, ggplot2::aes(x = .data$x, y = .data$y, group = .data$piece, subgroup = .data$hole)) +
      ggplot2::geom_polygon(ggplot2::aes(fill = .data$quintile), colour = "#ffffff", linewidth = 0.08) +
      ggplot2::scale_fill_manual(values = stats::setNames(PLANNING_MAP_QUINTILES, 1:5), na.value = "#d9d8d3", drop = FALSE,
                                 breaks = as.character(1:5), labels = c("1 (mais baixo)", "2", "3", "4", "5 (mais alto)"),
                                 name = "Quintil (cinzento: sem valor)")
  } else {
    drawn$log_index <- log(pmin(pmax(drawn$index / 100, PLANNING_MAP_INDEX_LIMITS[1]), PLANNING_MAP_INDEX_LIMITS[2]))
    breaks <- c(67, 80, 100, 125, 150)
    plot <- ggplot2::ggplot(drawn, ggplot2::aes(x = .data$x, y = .data$y, group = .data$piece, subgroup = .data$hole)) +
      ggplot2::geom_polygon(ggplot2::aes(fill = .data$log_index), colour = "#ffffff", linewidth = 0.08) +
      ggplot2::scale_fill_gradient2(
        low = PLANNING_MAP_DIVERGING[["low"]], mid = PLANNING_MAP_DIVERGING[["mid"]], high = PLANNING_MAP_DIVERGING[["high"]],
        midpoint = 0, limits = log(PLANNING_MAP_INDEX_LIMITS), breaks = log(breaks / 100), labels = breaks,
        na.value = "#d9d8d3", name = paste0(benchmark, " = 100 (cinzento: sem valor)")
      )
  }
  if (length(highlight) > 0) {
    marked <- drawn[drawn$id %in% highlight, , drop = FALSE]
    plot <- plot + ggplot2::geom_polygon(data = marked, fill = NA, colour = "#0b0b0b", linewidth = 0.6)
  }
  if (!is.null(frames) && nrow(frames) > 0) {
    plot <- plot + ggplot2::geom_path(data = frames, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                                      inherit.aes = FALSE, colour = "#898781", linewidth = 0.3)
  }
  plot +
    ggplot2::facet_wrap(~panel, ncol = min(2L, length(panels))) +
    ggplot2::coord_quickmap() +
    ggplot2::labs(caption = boundaries$meta$attribution) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 11, colour = "#0b0b0b", margin = ggplot2::margin(4, 0, 4, 0), lineheight = 1.1),
      legend.position = "bottom", legend.key.width = grid::unit(2.2, "cm"),
      plot.caption = ggplot2::element_text(colour = "#52514e", size = 9),
      panel.spacing = grid::unit(1, "lines")
    )
}

# The area under a click, by ray casting on the polygons (outline in, holes out).
planning_map_hit <- function(boundaries, x, y) {
  if (is.null(x) || is.null(y)) return(NULL)
  p <- boundaries$polygons
  near <- p[abs(p$x - x) < 1.5 & abs(p$y - y) < 1.5, "id", drop = TRUE]
  for (id in unique(stats::na.omit(near))) {
    rings <- split(p[p$id == id, c("x", "y")], p$hole[p$id == id])
    crossings <- sum(vapply(rings, function(r) {
      xi <- r$x; yi <- r$y; xj <- c(xi[length(xi)], xi[-length(xi)]); yj <- c(yi[length(yi)], yi[-length(yi)])
      sum(((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
    }, numeric(1)))
    if (crossings %% 2 == 1) return(id)
  }
  NULL
}
