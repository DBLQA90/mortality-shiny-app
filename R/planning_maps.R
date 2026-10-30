# =========================================================
# Maps of Portugal: several indicators side by side
# =========================================================
#
# Small multiples of the country by municipality, NUTS III, NUTS II or ULS, one
# map per indicator. By default every indicator is shown against the benchmark
# (Portugal = 100) on the same log-symmetric range, 67 to 150, coloured by what
# the difference means (PLANNING_INDICATOR_DIRECTION): red where the area is
# worse than Portugal and blue where it is better - so a high mortality and a
# low life expectancy are both red, and maps of different indicators, in
# different units, read together as a heatmap of where things are worse. An
# indicator with no better or worse direction (age structure, births,
# employment sectors) is drawn violet below and yellow above, a pair that
# carries no judgement. The alternative is each indicator's own quintiles,
# which shows its internal ranking whatever its spread.
#
# Each map has its own legend, so the panels are separate plots laid out with
# grid (planning_map_draw()).
#
# The boundaries are plain coordinate tables built once by
# tools/build_boundaries.R from the geopt package (CAOP, DGT); drawing them
# needs only ggplot2.

PLANNING_MAP_LEVELS <- c("Municípios" = "municipio", "NUTS III" = "nuts3", "ULS" = "uls", "NUTS II" = "nuts2")
PLANNING_MAP_MAX_INDICATORS <- 4L
# The diverging pair and neutral midpoint of the chart palette; the index is
# drawn on a log scale so that 150 and 67 sit at the same distance from 100.
PLANNING_MAP_DIVERGING <- c(low = "#256abf", mid = "#f0efec", high = "#e34948")
# Neutral indicators: violet and yellow, distinct from the blue/red of better
# and worse (OKLab distance x100 between them: 42 in normal vision, 33 or more
# under each colour-vision deficiency).
PLANNING_MAP_NEUTRAL <- c(low = "#6250d6", mid = "#f0efec", high = "#eda100")

# Which way is better, for the indicators where that is not a matter of
# opinion. Everything else - age structure and dependency, births and
# fertility, employment sectors, pensioners, waste per inhabitant, the crude
# death rate (which mostly measures how old the population is) - is neutral.
PLANNING_INDICATOR_DIRECTION <- c(
  smr_all = "worse", dsr_all = "worse", dsr_premature = "worse", dsr_preventable = "worse", dsr_treatable = "worse",
  ypll_rate = "worse", infant_rate = "worse", neonatal_rate = "worse", early_neonatal_rate = "worse",
  postneonatal_rate = "worse", late_fetal_rate = "worse", perinatal_rate = "worse", preterm_pct = "worse",
  low_birth_weight_pct = "worse", teen_births_pct = "worse", pct_education_none = "worse", illiteracy_rate = "worse",
  rsi_rate = "worse",
  life_expectancy = "better", life_expectancy_men = "better", life_expectancy_women = "better",
  life_expectancy_65 = "better", life_expectancy_65_men = "better", life_expectancy_65_women = "better",
  purchasing_power = "better", earnings_mean = "better", pension_mean = "better", pct_education_higher = "better",
  waste_selective_per_capita = "better"
)
planning_indicator_direction <- function(id) {
  direction <- unname(PLANNING_INDICATOR_DIRECTION[id])
  ifelse(is.na(direction), "neutral", direction)
}
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
planning_map_year <- function(id, max_year, level, vintage = "2024") {
  # Only the NUTS II regions of 2024 are made of the rebuilt year's regions;
  # under NUTS 2013, Centro and Alentejo are not.
  whole_regions <- identical(level, "nuts2") && identical(normalize_nuts_vintage(vintage), "2024")
  years <- if (whole_regions) planning_indicator_years(id) else planning_small_area_years(id)
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
    year <- planning_map_year(id, max_year, level, planning_lookup_vintage(lookup))
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
      direction = planning_indicator_direction(id),
      panel = paste0(paste(strwrap(spec$label, 46), collapse = "\n"), "\n", spec$unit, ", ", planning_period_label(id, year),
                     if (identical(method, "suavizado")) ", suavizado" else ""),
      area = rows$area,
      value = rows$value,
      observed = rows$observed,
      flag = rows$flag,
      digits = spec$digits,
      benchmark_value = reference,
      index = ifelse(is.finite(rows$value) & is.finite(reference) & reference > 0, rows$value / reference * 100, NA_real_),
      quintile = quintile,
      # Of the observed value's interval, smoothed or not.
      significance = if (all(c("lower", "upper") %in% names(rows))) planning_significance(rows$lower, rows$upper, rep(reference, nrow(rows))) else NA_character_,
      class = map_significance_class(significance)
    )
  }))
}

# One map per indicator, each with its own legend, as a list that
# planning_map_draw() lays out two per row. `scale`: "index" (against the
# benchmark, coloured by meaning) or "quintile" (each indicator's own fifths).
#
# For the Mapas tab and the significance colours, `scale` may also be "classes"
# (column `class`, coloured by `classes`, named colours lowest first) or
# "linear" (column `z` on `limits`, labelled `ticks`/`tick_labels`), with
# `legend_title`. "significance" draws each panel's significance classes in the
# colours of its direction (map_class_colours()).
planning_map_chart <- function(values, boundaries, scale = "index", highlight = character(0), benchmark = "Portugal",
                               classes = NULL, limits = NULL, ticks = NULL, tick_labels = NULL, legend_title = NULL) {
  polygons <- boundaries$polygons
  frames <- boundaries$frames
  panels <- unique(values$panel)
  plots <- lapply(panels, function(p) {
    v <- values[values$panel == p, , drop = FALSE]
    keep <- intersect(c("area", "index", "quintile", "class", "z"), names(v))
    drawn <- dplyr::left_join(polygons, v[, keep], by = c("id" = "area"))
    direction <- v$direction[[1]]
    base <- function(data) ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y, group = .data$piece, subgroup = .data$hole))
    if (identical(scale, "significance") && all(is.na(v$significance))) {
      scale <- "index"
      legend_title <- paste0(benchmark, " = 100 (sem intervalo)\n", switch(direction, neutral = "sem melhor nem pior", "vermelho: pior"))
    }
    if (identical(scale, "significance")) {
      scale <- "classes"
      classes <- map_class_colours(unname(MAP_SIGNIFICANCE_LABELS), direction)
      legend_title <- paste0("Face a ", benchmark)
    }
    if (identical(scale, "classes")) {
      drawn$class <- factor(drawn$class, levels = names(classes))
      plot <- base(drawn) +
        ggplot2::geom_polygon(ggplot2::aes(fill = .data$class), colour = "#ffffff", linewidth = 0.08) +
        ggplot2::scale_fill_manual(values = classes, na.value = MAP_NA_COLOUR, drop = FALSE, breaks = names(classes),
                                   name = legend_title %||% "")
    } else if (identical(scale, "linear")) {
      plot <- base(drawn) +
        ggplot2::geom_polygon(ggplot2::aes(fill = pmin(pmax(.data$z, limits[1]), limits[2])), colour = "#ffffff", linewidth = 0.08) +
        ggplot2::scale_fill_gradient2(low = PLANNING_MAP_DIVERGING[["low"]], mid = PLANNING_MAP_DIVERGING[["mid"]],
                                      high = PLANNING_MAP_DIVERGING[["high"]], midpoint = 0, limits = limits,
                                      breaks = ticks, labels = tick_labels %||% ticks, na.value = MAP_NA_COLOUR,
                                      name = legend_title %||% "")
    } else if (identical(scale, "quintile")) {
      drawn$quintile <- factor(drawn$quintile, levels = 1:5)
      plot <- base(drawn) +
        ggplot2::geom_polygon(ggplot2::aes(fill = .data$quintile), colour = "#ffffff", linewidth = 0.08) +
        ggplot2::scale_fill_manual(values = stats::setNames(PLANNING_MAP_QUINTILES, 1:5), na.value = MAP_NA_COLOUR, drop = FALSE,
                                   breaks = as.character(1:5), labels = c("1 (mais baixo)", "2", "3", "4", "5 (mais alto)"),
                                   name = "Quintil")
    } else {
      drawn$log_index <- log(pmin(pmax(drawn$index / 100, PLANNING_MAP_INDEX_LIMITS[1]), PLANNING_MAP_INDEX_LIMITS[2]))
      colours <- switch(direction,
                        worse = PLANNING_MAP_DIVERGING,
                        better = c(low = PLANNING_MAP_DIVERGING[["high"]], mid = PLANNING_MAP_DIVERGING[["mid"]], high = PLANNING_MAP_DIVERGING[["low"]]),
                        PLANNING_MAP_NEUTRAL)
      meaning <- switch(direction, worse = "vermelho: pior", better = "vermelho: pior", "sem melhor nem pior")
      breaks <- c(67, 80, 100, 125, 150)
      plot <- base(drawn) +
        ggplot2::geom_polygon(ggplot2::aes(fill = .data$log_index), colour = "#ffffff", linewidth = 0.08) +
        ggplot2::scale_fill_gradient2(
          low = colours[["low"]], mid = colours[["mid"]], high = colours[["high"]],
          midpoint = 0, limits = log(PLANNING_MAP_INDEX_LIMITS), breaks = log(breaks / 100), labels = breaks,
          na.value = MAP_NA_COLOUR, name = legend_title %||% paste0(benchmark, " = 100\n", meaning)
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
      ggplot2::coord_quickmap() +
      ggplot2::labs(title = p) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 11, colour = "#0b0b0b", hjust = 0.5, lineheight = 1.1),
        legend.position = "bottom", legend.key.width = grid::unit(1.4, "cm"), legend.title.position = "top",
        plot.margin = ggplot2::margin(4, 4, 14, 4),
        legend.title = ggplot2::element_text(size = 9, colour = "#52514e", lineheight = 1),
        legend.text = ggplot2::element_text(size = 8, colour = "#52514e")
      )
  })
  structure(list(plots = plots, caption = boundaries$meta$attribution), class = "planning_map")
}

# Draws the maps two per row, with the boundaries' attribution underneath.
planning_map_draw <- function(chart) {
  n <- length(chart$plots)
  columns <- min(2L, n)
  rows <- ceiling(n / columns)
  grid::grid.newpage()
  layout <- grid::grid.layout(rows + 1L, columns, heights = grid::unit(c(rep(1, rows), 1.5), c(rep("null", rows), "lines")))
  grid::pushViewport(grid::viewport(layout = layout))
  for (i in seq_len(n)) {
    print(chart$plots[[i]], vp = grid::viewport(layout.pos.row = (i - 1L) %/% columns + 1L, layout.pos.col = (i - 1L) %% columns + 1L))
  }
  grid::pushViewport(grid::viewport(layout.pos.row = rows + 1L, layout.pos.col = seq_len(columns)))
  grid::grid.text(chart$caption, x = 0.98, hjust = 1, gp = grid::gpar(fontsize = 9, col = "#52514e"))
  grid::popViewport(2)
  invisible(chart)
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
