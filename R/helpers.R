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

# Models the guided forecast builds from the others, not fitted on their own.
forecast_derived_model_labels <- c("Combinação (média dos 3 melhores)" = "combo",
                                   "Média de 6 modelos" = "mean6")

get_model_labels <- function(model_ids) {
  choices <- c(forecast_model_choices, forecast_derived_model_labels)
  labels <- names(choices)[match(model_ids, choices)]
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

# Corrections to values INE publishes wrongly (data/source_corrections.rds,
# built and documented by tools/build_source_corrections.R; reasons in
# ESCOLHAS_METODOLOGICAS.md). They are applied when a death file is read, so
# the snapshots keep INE's data as published. options(mortality.source_corrections
# = FALSE) reads the data exactly as INE publishes them.
source_corrections_enabled <- function() !identical(getOption("mortality.source_corrections"), FALSE)

source_corrections_cache <- new.env(parent = emptyenv())

# The registry beside the snapshot directory that `path` lies in, or NULL.
source_corrections_for <- function(path) {
  root <- sub("/snapshots/.*$", "", path)
  if (identical(root, path)) return(NULL)
  file <- file.path(root, "source_corrections.rds")
  if (!file.exists(file)) return(NULL)
  stamp <- paste(file, file.mtime(file))
  if (!identical(source_corrections_cache$stamp, stamp)) {
    source_corrections_cache$value <- tryCatch(readRDS(file), error = function(e) NULL)
    source_corrections_cache$stamp <- stamp
  }
  source_corrections_cache$value
}

# Replace, in a death file read from `path`, the rows the registry corrects:
# by (area, sex, cause, age band) in the files by age, by (area, sex, cause) in
# the all-ages totals. Other files are returned unchanged.
apply_source_corrections <- function(x, path) {
  if (!source_corrections_enabled() || !is.data.frame(x) || nrow(x) == 0 || !is.character(path) || grepl("^https?://", path)) return(x)
  m <- regmatches(path, regexec("/snapshots/(deaths|death_totals)/([0-9]+)/year_([0-9]{4})", path))[[1]]
  if (length(m) == 0) return(x)
  registry <- source_corrections_for(path)
  if (is.null(registry)) return(x)
  part <- if (identical(m[[2]], "deaths")) registry$deaths else registry$totals
  part <- part[part$indicator == m[[3]] & part$year == as.integer(m[[4]]) & part$cause %in% unique(x$cause), , drop = FALSE]
  if (nrow(part) == 0) return(x)
  keys <- intersect(c("area", "sex", "cause", "age_band"), names(x))
  replaced <- do.call(paste, c(x[keys], sep = "\r")) %in% do.call(paste, c(part[keys], sep = "\r"))
  if (!any(replaced)) return(x)
  template <- x[replaced, , drop = FALSE]
  at <- match(do.call(paste, c(part[keys], sep = "\r")), do.call(paste, c(template[keys], sep = "\r")))
  part <- part[!is.na(at), , drop = FALSE]
  new_rows <- template[at[!is.na(at)], , drop = FALSE]
  new_rows$deaths <- part$deaths
  out <- rbind(x[!replaced, , drop = FALSE], new_rows)
  attr(out, "source_corrected") <- TRUE
  out
}

read_death_file <- function(path) apply_source_corrections(readRDS(path), path)
