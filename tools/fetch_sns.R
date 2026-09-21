#!/usr/bin/env Rscript
# Fetch primary-care indicators from the SNS Transparency portal
# (transparencia.sns.gov.pt, Opendatasoft) into data/snapshots/sns.
#
#   Rscript tools/fetch_sns.R [datasets=all]
#
# The portal's Explore API is public and needs no key; each dataset is exported
# whole as JSON (a few thousand rows) and stored in long form, one file per
# dataset, through R/data_versions.R - so a month the portal revises (the most
# recent one is provisional and revised at the end of the following month) is
# archived, not overwritten silently.
#
# Output: data/snapshots/sns/<dataset>.rds
#   columns: period ("YYYY-MM"), unit (the portal's label: "CSP da ULS Guarda",
#   "Área dos CSP da ULS Algarve", or an ACES before 2024), region, field,
#   value, dataset

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "fetch_sns.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

# Dataset -> its period and unit fields.
SNS_DATASETS <- list(
  "utentes-inscritos-em-cuidados-de-saude-primarios" = c(period = "periodo", unit = "aces", region = "ars"),
  "rastreios-oncologicos" = c(period = "tempo", unit = "area_csp", region = "regiao"),
  "diabetes" = c(period = "tempo", unit = "area_csp", region = "regiao"),
  "hipertensao" = c(period = "tempo", unit = "area_csp", region = "regiao"),
  "saude-da-mulher-e-crianca" = c(period = "tempo", unit = "area_csp", region = "regiao")
)
wanted <- get_arg("datasets", "all")
if (!identical(wanted, "all")) SNS_DATASETS <- SNS_DATASETS[strsplit(wanted, ",", fixed = TRUE)[[1]]]

base <- "https://transparencia.sns.gov.pt/api/explore/v2.1/catalog/datasets"
failed <- character(0)
for (dataset in names(SNS_DATASETS)) {
  spec <- SNS_DATASETS[[dataset]]
  url <- paste0(base, "/", dataset, "/exports/json")
  raw <- NULL
  for (attempt in 1:4) {
    raw <- tryCatch(jsonlite::fromJSON(url, flatten = TRUE), error = function(e) {
      message("  ", dataset, ": ", conditionMessage(e))
      NULL
    })
    if (is.data.frame(raw) && nrow(raw) > 0) break
    Sys.sleep(20 * attempt)
  }
  if (!is.data.frame(raw) || nrow(raw) == 0) {
    failed <- c(failed, dataset)
    next
  }
  numeric_fields <- names(raw)[vapply(raw, is.numeric, logical(1))]
  numeric_fields <- numeric_fields[!grepl("^(localizacao|ponto)", numeric_fields)]
  tidy <- dplyr::bind_rows(lapply(numeric_fields, function(name) {
    values <- as.numeric(raw[[name]])
    tibble::tibble(
      period = substr(as.character(raw[[spec[["period"]]]]), 1, 7),
      unit = trimws(as.character(raw[[spec[["unit"]]]])),
      region = trimws(as.character(raw[[spec[["region"]]]])),
      field = name,
      value = values
    )
  })) %>%
    dplyr::filter(!is.na(.data$period), nzchar(.data$unit)) %>%
    dplyr::mutate(dataset = dataset) %>%
    dplyr::arrange(.data$period, .data$unit, .data$field)
  message(sprintf("%s: %d rows, %s to %s, %d units", dataset, nrow(raw), min(tidy$period), max(tidy$period), dplyr::n_distinct(tidy$unit)))
  save_rds_atomic(tidy, file.path("data/snapshots/sns", paste0(dataset, ".rds")))
  Sys.sleep(2)
}
message("Done. Failed: ", length(failed), if (length(failed)) paste0(" (", paste(failed, collapse = ", "), ")") else "")
if (length(failed) > 0) quit(status = 1)
