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
