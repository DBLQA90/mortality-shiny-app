#!/usr/bin/env Rscript
# Fetch the INE municipal indicators the app is validated against.
#
#   Rscript tools/fetch_validation_refs.R
#
# These are INE's own published values for each municipality - ageing index,
# dependency ratios, crude birth and death rates, five-year infant and teenage
# birth rates, waste per inhabitant, RSI beneficiaries and pensioners per 1.000
# residents aged 15-64. The app computes all of them from its own snapshots, so
# comparing the two is the strongest end-to-end check there is: it exercises the
# fetchers, the municipality matching and the arithmetic at once.
#
# tools/regression_check.R reads the result and reports any municipality whose
# value drifts away from INE's.
#
# Output: data/snapshots/validation/ine_municipal_refs.rds
#   columns: indicator, code, label, value, period, source_indicator

suppressMessages(library(dplyr))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())

# The period is part of the check: it pins the year the baseline in
# data/regression_baseline.csv was measured on. Moving to a newer year means
# refetching here and re-recording the baseline.
REFS <- list(
  ageing_index     = list(id = "0012909", period = "S7A2024"),
  old_dependency   = list(id = "0012910", period = "S7A2024"),
  youth_dependency = list(id = "0012906", period = "S7A2024"),
  birth_rate       = list(id = "0013044", period = "S7A2024"),
  death_rate       = list(id = "0013046", period = "S7A2024"),
  infant_5y        = list(id = "0013578", period = "S11A20202024"),
  teen_5y          = list(id = "0013336", period = "S11A20202024"),
  waste            = list(id = "0012765", period = "S7A2024"),
  rsi_1564         = list(id = "0013420", period = "S7A2024"),
  pens_1564        = list(id = "0014599", period = "S7A2024")
)

client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)
out <- list()
failed <- character(0)
for (name in names(REFS)) {
  raw <- NULL
  for (attempt in 1:4) {
    raw <- tryCatch(client$get_data(REFS[[name]]$id, dim1 = REFS[[name]]$period), error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) break
    Sys.sleep(30 * attempt)
  }
  if (is.null(raw) || nrow(raw) == 0) {
    message(name, ": FAILED")
    failed <- c(failed, name)
    next
  }
  # Keep the total of every dimension other than the period (sex, typology).
  label_cols <- grep("^dim_[0-9]+_t$", names(raw), value = TRUE)
  keep <- rep(TRUE, nrow(raw))
  for (col in setdiff(label_cols, "dim_1_t")) {
    values <- trimws(as.character(raw[[col]]))
    if (any(values %in% c("Total", "HM"))) keep <- keep & values %in% c("Total", "HM")
  }
  out[[name]] <- tibble::tibble(
    indicator = name, code = as.character(raw$geocod[keep]), label = as.character(raw$geodsg[keep]),
    value = suppressWarnings(as.numeric(raw$valor[keep])),
    period = REFS[[name]]$period, source_indicator = REFS[[name]]$id
  )
  message(sprintf("  %-17s %s %s: %d rows", name, REFS[[name]]$id, REFS[[name]]$period, sum(keep)))
  Sys.sleep(5)
}

if (length(out) == 0) {
  message("Nothing fetched.")
  quit(status = 1)
}
refs <- dplyr::bind_rows(out)
versioned_save_rds(refs, "data/snapshots/validation/ine_municipal_refs.rds",
                   tool = "fetch_validation_refs.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))
message("Done: ", nrow(refs), " rows, ", dplyr::n_distinct(refs$indicator), " indicators. Failed: ", length(failed))
if (length(failed) > 0) quit(status = 1)
