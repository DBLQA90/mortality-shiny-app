# =========================================================
# General helpers
# =========================================================

normalize_ine_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:punct:]]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

get_default_area_selection <- function() {
  "Portugal"
}

get_selection_label <- function(selected_areas, custom_label = NULL) {
  label <- if (is.null(custom_label)) "" else trimws(custom_label)
  if (nzchar(label)) {
    return(label)
  }
  if (length(selected_areas) == 1) {
    return(selected_areas[[1]])
  }
  "Soma de locais"
}

safe_filename_token <- function(x) {
  gsub("[^[:alnum:]_]+", "_", iconv(x, to = "ASCII//TRANSLIT"))
}

get_model_labels <- function(model_ids) {
  labels <- names(forecast_model_choices)[match(model_ids, forecast_model_choices)]
  labels[is.na(labels)] <- model_ids[is.na(labels)]
  labels
}

format_source_indicators <- function(x) {
  x <- sort(unique(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0) {
    return("N/D")
  }

  paste(x, collapse = ", ")
}

get_loaded_source_summary <- function(data) {
  population_source <- if ("population_source" %in% names(data)) {
    format_source_indicators(data$population_source)
  } else {
    "N/D"
  }

  death_source <- if ("death_source" %in% names(data)) {
    format_source_indicators(data$death_source)
  } else if ("source_indicator" %in% names(data)) {
    format_source_indicators(data$source_indicator)
  } else {
    "N/D"
  }

  list(
    population_source = population_source,
    death_source = death_source
  )
}

# =========================================================
# Provenance of exported files
# =========================================================
# INE revises published data, so the same analysis can give different figures
# months apart. Every export says which version of the data produced it: the
# import date, the region definition in force, and when the file was made. A CSV
# carries it as a comment line after the data, a chart as a caption.

export_stamp <- function(data_date = NA_character_, vintage = NULL) {
  parts <- c(
    if (is.na(data_date) || !nzchar(data_date)) "Dados do INE: data de importação desconhecida" else paste0("Dados do INE importados até ", data_date),
    if (!is.null(vintage) && nzchar(vintage)) paste0("regiões NUTS ", vintage) else NULL,
    paste0("exportado em ", format(Sys.Date(), "%Y-%m-%d"))
  )
  paste(parts, collapse = "; ")
}

helpers_write_csv_utf8 <- function(x, file, stamp = NULL) {
  utils::write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
  if (!is.null(stamp) && nzchar(stamp)) {
    # After the data, so a reader that ignores the last line still parses.
    connection <- file(file, open = "a", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    writeLines(paste0("# ", stamp), connection, useBytes = FALSE)
  }
  invisible(file)
}

# Add the stamp under a ggplot, keeping any caption the plot already has.
stamp_ggplot <- function(plot_obj, stamp = NULL) {
  if (is.null(stamp) || !nzchar(stamp) || !inherits(plot_obj, "ggplot")) return(plot_obj)
  existing <- plot_obj$labels$caption
  caption <- if (is.null(existing) || !nzchar(existing)) stamp else paste0(existing, "\n", stamp)
  plot_obj +
    ggplot2::labs(caption = caption) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 7, colour = "grey35", hjust = 0))
}

helpers_save_ggplot_png <- function(file, plot_obj, width = 1200, height = 800, res = 150, stamp = NULL) {
  grDevices::png(file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(stamp_ggplot(plot_obj, stamp))
}

helpers_save_base_plot_png <- function(file, plot_expr, width = 1200, height = 800, res = 150, stamp = NULL) {
  grDevices::png(file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(plot_expr)
  if (!is.null(stamp) && nzchar(stamp)) {
    graphics::mtext(stamp, side = 1, line = -0.2, outer = TRUE, adj = 0, cex = 0.55, col = "grey35")
  }
}
