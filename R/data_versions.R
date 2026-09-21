# =========================================================
# Data versions: import dates, archived previous versions, reconstruction
# =========================================================
# INE revises published figures: provisional years become final, population
# series are re-estimated, indicators are replaced by new editions. An analysis
# run today and repeated next year can therefore give different numbers for the
# same selection, with nothing in the app to say why. This module keeps the
# history needed to explain it:
#
#   data/import_log.csv   one row per file written: when, by which tool, whether
#                         it was new or replaced a different version, and the
#                         row count and value total before and after
#   data/archive/<run>/   the previous version of every replaced file, under its
#                         original relative path
#
# Writing identical content is a no-op, so re-fetching a year that INE has not
# revised leaves no trace and no git change. `data_as_of()` rebuilds the data
# directory as it stood on any date, which the app can then be pointed at.
#
# Versions from before this log existed are recorded from git history by
# tools/backfill_import_log.R; their previous versions stay in git
# (`archived_as` = "git:<commit>:<path>") rather than being copied.
#
# Only files under a directory named `data` are versioned. A tool writing to a
# scratch directory elsewhere writes plainly.

import_log_columns <- c(
  "run_id", "imported_at", "relpath", "dataset", "year", "action",
  "rows_old", "rows_new", "value_column", "total_old", "total_new",
  "archived_as", "tool", "note", "rows_changed"
)

# The `data` directory a path lives under, or NULL.
data_root_for <- function(path) {
  full <- normalizePath(path, mustWork = FALSE)
  parts <- strsplit(full, "/", fixed = TRUE)[[1]]
  hits <- which(parts == "data")
  if (length(hits) == 0) return(NULL)
  paste(parts[seq_len(max(hits))], collapse = "/")
}

import_log_path <- function(data_root) file.path(data_root, "import_log.csv")

# One identifier per refresh, shared by every tool the driver launches through
# the environment, so a run's changes can be read together.
current_run_id <- function() {
  configured <- Sys.getenv("DATA_RUN_ID", unset = "")
  if (nzchar(configured)) return(configured)
  if (!exists(".data_run_id", envir = .data_versions_state, inherits = FALSE)) {
    assign(".data_run_id", format(Sys.time(), "%Y-%m-%dT%H%M%S"), envir = .data_versions_state)
  }
  get(".data_run_id", envir = .data_versions_state)
}
.data_versions_state <- new.env(parent = emptyenv())

# Dataset key and year of a relative path: "snapshots/deaths/0013166/year_2024/
# cause_x.rds" -> dataset "deaths", year 2024.
relpath_dataset <- function(relpath) {
  parts <- strsplit(relpath, "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 2 && identical(p[[1]], "snapshots")) {
      if (p[[2]] %in% c("planning_extra", "sns", "weekly_deaths") && length(p) >= 3) return(paste0(p[[2]], "/", sub("\\.rds$", "", p[[3]])))
      return(sub("\\.rds$", "", p[[2]]))
    }
    sub("\\.rds$", "", p[[length(p)]])
  }, character(1))
}

relpath_year <- function(relpath) {
  out <- rep(NA_integer_, length(relpath))
  has <- grepl("year_\\d{4}", relpath)
  out[has] <- as.integer(sub("year_", "", regmatches(relpath, regexpr("year_\\d{4}", relpath))))
  out
}

# The column that carries the measurement, for the before/after totals.
value_column_of <- function(x) {
  if (!is.data.frame(x)) return(NA_character_)
  hit <- intersect(c("deaths", "pop", "births", "value"), names(x))
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

# A total comparable across versions of a file: Portugal, both sexes, all
# causes, the "Total" category - whichever of those the file has - summed over
# the rest (age bands, municipalities). Summing every row would count regional
# rows and sex breakdowns, which change between INE editions even when no value
# does: the population revision would read +110% instead of +3.8%.
value_total <- function(x, column) {
  if (!is.data.frame(x) || is.na(column) || !column %in% names(x)) return(NA_real_)
  keep <- rep(TRUE, nrow(x))
  narrow <- function(col, value) {
    if (col %in% names(x) && any(as.character(x[[col]]) == value, na.rm = TRUE)) {
      keep <<- keep & as.character(x[[col]]) %in% value
    }
  }
  narrow("area", "Portugal")
  narrow("sex", "HM")
  narrow("cause", "Todas as causas de morte")
  narrow("category", "Total")
  sum(suppressWarnings(as.numeric(x[[column]][keep])), na.rm = TRUE)
}

# Rows whose value differs between two versions, matched on every other column
# except the source indicator (an edition change alone is not a revision). Rows
# present in only one version count too. This is what shows a correction that
# leaves the national total unchanged, such as births moved between
# municipalities.
rows_changed_count <- function(old, new, column) {
  if (!is.data.frame(old) || !is.data.frame(new) || is.na(column) ||
      !column %in% names(old) || !column %in% names(new)) {
    return(NA_real_)
  }
  keys <- setdiff(intersect(names(old), names(new)), c(column, "source_indicator"))
  if (length(keys) == 0) return(NA_real_)
  key_of <- function(df) do.call(paste, c(lapply(df[keys], as.character), sep = "\r"))
  a <- stats::setNames(suppressWarnings(as.numeric(old[[column]])), key_of(old))
  b <- stats::setNames(suppressWarnings(as.numeric(new[[column]])), key_of(new))
  shared <- intersect(names(a), names(b))
  differs <- sum(abs(a[shared] - b[shared]) > 1e-9 | xor(is.na(a[shared]), is.na(b[shared])), na.rm = TRUE)
  differs + length(setdiff(names(a), names(b))) + length(setdiff(names(b), names(a)))
}

# Same content regardless of row order and attributes.
snapshot_content_equal <- function(old, new) {
  if (!is.data.frame(old) || !is.data.frame(new)) {
    return(isTRUE(all.equal(old, new, check.attributes = FALSE)))
  }
  if (!setequal(names(old), names(new)) || nrow(old) != nrow(new)) return(FALSE)
  columns <- sort(names(new))
  normalise <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)[, columns, drop = FALSE]
    for (col in columns) if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    df <- df[do.call(order, unname(as.list(df))), , drop = FALSE]
    rownames(df) <- NULL
    df
  }
  isTRUE(all.equal(normalise(old), normalise(new), check.attributes = FALSE, tolerance = 1e-9))
}

read_import_log <- function(data_root) {
  path <- import_log_path(data_root)
  empty <- as.data.frame(stats::setNames(replicate(length(import_log_columns), character(0), simplify = FALSE), import_log_columns))
  if (!grepl("^https?://", path) && !file.exists(path)) return(tibble::as_tibble(empty))
  log <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character", encoding = "UTF-8"),
    error = function(e) empty
  )
  for (col in setdiff(import_log_columns, names(log))) log[[col]] <- NA_character_
  log <- tibble::as_tibble(log[, import_log_columns])
  log$year <- suppressWarnings(as.integer(log$year))
  for (col in c("rows_old", "rows_new", "total_old", "total_new", "rows_changed")) log[[col]] <- suppressWarnings(as.numeric(log[[col]]))
  log
}

append_import_log <- function(data_root, entries) {
  if (nrow(entries) == 0) return(invisible(NULL))
  path <- import_log_path(data_root)
  entries <- as.data.frame(entries, stringsAsFactors = FALSE)
  for (col in setdiff(import_log_columns, names(entries))) entries[[col]] <- NA
  entries <- entries[, import_log_columns]
  exists_already <- file.exists(path)
  if (exists_already) {
    header <- strsplit(gsub('"', "", readLines(path, n = 1, warn = FALSE)), ",", fixed = TRUE)[[1]]
    if (!identical(header, import_log_columns)) {
      # A log written before a column was added: rewrite it with the current
      # columns so appended rows line up.
      previous <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character", encoding = "UTF-8")
      for (col in setdiff(import_log_columns, names(previous))) previous[[col]] <- NA
      utils::write.table(previous[, import_log_columns], path, sep = ",", row.names = FALSE,
                         qmethod = "double", na = "", fileEncoding = "UTF-8")
    }
  }
  utils::write.table(
    entries, path, sep = ",", row.names = FALSE, col.names = !exists_already,
    append = exists_already, qmethod = "double", na = "", fileEncoding = "UTF-8"
  )
  invisible(NULL)
}

write_rds_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Could not move temporary file into ", path, call. = FALSE)
}

# Write `x` to `path`, keeping the history. Returns "added", "replaced" or
# "unchanged" invisibly.
versioned_save_rds <- function(x, path, tool = NA_character_, note = NA_character_) {
  data_root <- data_root_for(path)
  if (is.null(data_root)) {
    write_rds_atomic(x, path)
    return(invisible("unversioned"))
  }

  full <- normalizePath(path, mustWork = FALSE)
  relpath <- substring(full, nchar(data_root) + 2)
  column <- value_column_of(x)
  run_id <- current_run_id()

  entry <- list(
    run_id = run_id,
    imported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    relpath = relpath,
    dataset = relpath_dataset(relpath),
    year = relpath_year(relpath),
    rows_new = if (is.data.frame(x)) nrow(x) else NA,
    value_column = column,
    total_new = value_total(x, column),
    tool = tool,
    note = note
  )

  if (!file.exists(full)) {
    write_rds_atomic(x, full)
    append_import_log(data_root, as.data.frame(c(entry, action = "added"), stringsAsFactors = FALSE))
    return(invisible("added"))
  }

  old <- tryCatch(readRDS(full), error = function(e) NULL)
  if (!is.null(old) && snapshot_content_equal(old, x)) {
    return(invisible("unchanged"))
  }

  archived <- file.path("archive", run_id, relpath)
  dir.create(dirname(file.path(data_root, archived)), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(full, file.path(data_root, archived), overwrite = FALSE, copy.date = TRUE)) {
    # Already archived in this run (the same file written twice): keep the
    # first, which is the version from before the run.
    if (!file.exists(file.path(data_root, archived))) stop("Could not archive ", relpath, call. = FALSE)
  }
  write_rds_atomic(x, full)
  append_import_log(data_root, as.data.frame(c(
    entry,
    action = "replaced",
    rows_old = if (is.data.frame(old)) nrow(old) else NA,
    total_old = value_total(old, value_column_of(old)),
    rows_changed = rows_changed_count(old, x, column),
    archived_as = archived
  ), stringsAsFactors = FALSE))
  invisible("replaced")
}

# ---------------------------------------------------------
# Reading the history
# ---------------------------------------------------------
import_timestamp <- function(x) {
  as.POSIXct(sub("([+-]\\d{2})(\\d{2})$", "\\1\\2", x), format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
}

# The entry that produced each file's current version.
current_file_versions <- function(log) {
  if (nrow(log) == 0) return(log)
  log$ts <- import_timestamp(log$imported_at)
  log <- log[order(log$relpath, log$ts), , drop = FALSE]
  log[!duplicated(log$relpath, fromLast = TRUE), , drop = FALSE]
}

# Per dataset: files, years, when the current data was imported and when a
# published value last changed.
dataset_import_summary <- function(log) {
  if (nrow(log) == 0) return(tibble::tibble())
  current <- current_file_versions(log)
  changes <- log[log$action == "replaced", , drop = FALSE]
  last_change <- if (nrow(changes) > 0) {
    tapply(substr(changes$imported_at, 1, 10), changes$dataset, max)
  } else {
    character(0)
  }
  current %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(
      files = dplyr::n(),
      years = if (all(is.na(year))) "" else paste0(min(year, na.rm = TRUE), "-", max(year, na.rm = TRUE)),
      first_import = min(substr(imported_at, 1, 10)),
      last_import = max(substr(imported_at, 1, 10)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(last_change = unname(last_change[dataset])) %>%
    dplyr::arrange(dataset)
}

# One row per run that added or replaced something.
import_runs_summary <- function(log) {
  if (nrow(log) == 0) return(tibble::tibble())
  log %>%
    dplyr::group_by(run_id) %>%
    dplyr::summarise(
      date = min(substr(imported_at, 1, 10)),
      tools = paste(sort(unique(stats::na.omit(tool[nzchar(tool)]))), collapse = ", "),
      datasets = paste(sort(unique(dataset)), collapse = ", "),
      added = sum(action == "added"),
      replaced = sum(action == "replaced"),
      note = paste(unique(stats::na.omit(note[nzchar(note)])), collapse = " | "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(date), dplyr::desc(run_id))
}

# What a run changed, by dataset and year, with the value totals before and
# after, so a revision can be sized without opening files.
run_changes <- function(log, run_id) {
  # Files without a measured value (the inventory, reports) are not revisions.
  rows <- log[log$run_id == run_id & log$action == "replaced" & !is.na(log$value_column) & nzchar(log$value_column), , drop = FALSE]
  if (nrow(rows) == 0) return(tibble::tibble())
  # Death chunks are one file per cause, and causes nest (chapters contain their
  # sub-causes), so their totals cannot be added up: where the all-cause file is
  # among a year's files, the total is that file's alone.
  rows$all_cause <- grepl("cause_todas_as_causas", rows$relpath)
  rows %>%
    dplyr::group_by(dataset, year, value_column) %>%
    dplyr::summarise(
      files = dplyr::n(),
      rows_changed = sum(rows_changed, na.rm = TRUE),
      total_old = if (any(all_cause)) sum(total_old[all_cause], na.rm = TRUE) else sum(total_old, na.rm = TRUE),
      total_new = if (any(all_cause)) sum(total_new[all_cause], na.rm = TRUE) else sum(total_new, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(change_pct = ifelse(total_old > 0, (total_new - total_old) / total_old * 100, NA_real_)) %>%
    dplyr::arrange(dataset, year)
}

# The latest import of the files currently in use, for display.
latest_import_date <- function(log) {
  if (nrow(log) == 0) return(NA_character_)
  max(substr(current_file_versions(log)$imported_at, 1, 10))
}

# ---------------------------------------------------------
# Reconstruction
# ---------------------------------------------------------
# Rebuild `data_root` as it stood at the end of `as_of` (a date or timestamp)
# into `out_dir`. Files unchanged since are hard-linked where the filesystem
# allows, so a reconstruction costs little space. A file first added after the
# date is left out; a file replaced after the date is restored from the version
# archived by its first later replacement - from data/archive, or from git for
# versions recorded by the backfill.
data_as_of <- function(as_of, data_root, out_dir, repo_dir = dirname(data_root)) {
  cutoff <- if (inherits(as_of, "POSIXct")) as_of else as.POSIXct(paste0(as.character(as_of), " 23:59:59"), tz = "UTC")
  log <- read_import_log(data_root)
  log$ts <- import_timestamp(log$imported_at)

  files <- list.files(data_root, recursive = TRUE, all.files = FALSE)
  files <- files[!startsWith(files, "archive/") & files != "import_log.csv"]

  later <- log[!is.na(log$ts) & log$ts > cutoff, , drop = FALSE]
  later <- later[order(later$relpath, later$ts), , drop = FALSE]
  first_later <- later[!duplicated(later$relpath), , drop = FALSE]

  restored <- 0L
  dropped <- 0L
  linked <- 0L
  missing <- character(0)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  place <- function(src, dest) {
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(dest)) unlink(dest)
    ok <- suppressWarnings(file.link(src, dest))
    if (!ok) ok <- file.copy(src, dest, copy.date = TRUE)
    ok
  }

  for (relpath in files) {
    change <- first_later[first_later$relpath == relpath, , drop = FALSE]
    dest <- file.path(out_dir, relpath)
    if (nrow(change) == 0) {
      if (place(file.path(data_root, relpath), dest)) linked <- linked + 1L
      next
    }
    if (identical(change$action, "added")) {
      # Not present yet at the cutoff.
      dropped <- dropped + 1L
      next
    }
    source <- change$archived_as
    if (startsWith(source, "git:")) {
      spec <- sub("^git:", "", source)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      status <- suppressWarnings(system2("git", c("-C", shQuote(repo_dir), "show", shQuote(spec)), stdout = dest, stderr = FALSE))
      if (!identical(as.integer(status), 0L)) {
        unlink(dest)
        missing <- c(missing, relpath)
        next
      }
    } else if (!place(file.path(data_root, source), dest)) {
      missing <- c(missing, relpath)
      next
    }
    restored <- restored + 1L
  }

  list(out_dir = out_dir, unchanged = linked, restored = restored, not_yet_imported = dropped, missing = missing)
}

# ---------------------------------------------------------
# For the app
# ---------------------------------------------------------
# The data directory holding the snapshots in use: the parent of the snapshot
# directory, which is how data_as_of() lays a reconstruction out as well.
app_data_root <- function(snapshot_dir) {
  snapshot_dir <- sub("/+$", "", snapshot_dir)
  if (identical(basename(snapshot_dir), "snapshots")) dirname(snapshot_dir) else snapshot_dir
}

DATASET_LABELS <- c(
  "deaths" = "Óbitos por causa e idade",
  "death_totals" = "Óbitos por causa, todas as idades",
  "regional_deaths" = "Óbitos por causa e idade, linhas regionais do INE",
  "population" = "População residente",
  "births" = "Nados-vivos",
  "infant_deaths" = "Óbitos com menos de 1 ano, por causa",
  "infant_totals" = "Óbitos com menos de 1 ano (contagens completas)",
  "planning_extra/rsi_beneficiaries" = "Beneficiários do RSI",
  "planning_extra/pensioners" = "Pensionistas",
  "planning_extra/pension_mean" = "Valor médio das pensões",
  "planning_extra/purchasing_power_per_capita" = "Poder de compra per capita",
  "planning_extra/purchasing_power_share" = "Proporção do poder de compra",
  "planning_extra/waste_collected" = "Resíduos urbanos recolhidos",
  "planning_extra/births_by_mother_age" = "Nados-vivos por idade da mãe",
  "planning_extra/births_by_gestation" = "Nados-vivos por duração da gestação",
  "planning_extra/infant_deaths_by_age" = "Óbitos com menos de 1 ano, por idade",
  "sns/utentes-inscritos-em-cuidados-de-saude-primarios" = "SNS: utentes inscritos e médico de família",
  "sns/rastreios-oncologicos" = "SNS: rastreios oncológicos",
  "sns/diabetes" = "SNS: programa de diabetes",
  "sns/hipertensao" = "SNS: programa de hipertensão",
  "sns/saude-da-mulher-e-crianca" = "SNS: saúde da mulher e da criança",
  "weekly_deaths/0012100" = "Óbitos semanais por NUTS III e idade (NUTS 2024)",
  "weekly_deaths/0010112" = "Óbitos semanais por NUTS III e idade (NUTS 2013)",
  "nuts_lookup_2013" = "Municípios por região (NUTS 2013)",
  "nuts_lookup_2024" = "Municípios por região (NUTS 2024)",
  "uls_lookup" = "Municípios por ULS e ARS",
  "snapshot_inventory" = "Inventário dos ficheiros",
  "ambiguous_areas" = "Relatório de nomes ambíguos"
)

dataset_label <- function(dataset) {
  label <- unname(DATASET_LABELS[dataset])
  ifelse(is.na(label), dataset, label)
}
