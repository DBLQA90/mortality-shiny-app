#!/usr/bin/env Rscript
# Build data/import_log.csv from git history, for data written before the log
# existed.
#
#   Rscript tools/backfill_import_log.R [since=2026-08-01] [force=false]
#
# For every file under data/ it records:
#   - "added"    at the commit that created it - or, for files that were built
#                up over many commits before `since` (the 0008206 death chunks
#                were assembled area by area in May), at the last of those
#                commits, which is when the file reached the form it had when
#                the data was first complete;
#   - "replaced" at every later commit from `since` on that changed it, with
#                archived_as = "git:<parent commit>:<path>" so the previous
#                version can be restored from git by data_as_of() without being
#                copied into data/archive, and with the value totals of both
#                versions.
#
# The import date of a backfilled version is its commit date, the closest record
# of when it was fetched. Run once; later writes are logged as they happen.

suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
repo <- normalizePath(file.path(script_dir, ".."))
setwd(repo)
sys.source(file.path(repo, "R", "data_versions.R"), envir = environment())

since <- as.Date(get_arg("since", "2026-08-01"))
force <- tolower(get_arg("force", "false")) %in% c("true", "1", "yes")
data_root <- file.path(repo, "data")

if (file.exists(import_log_path(data_root)) && !force) {
  stop("data/import_log.csv already exists; pass force=true to rebuild it from git.", call. = FALSE)
}

git <- function(...) system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = FALSE)

raw <- git("log", "--reverse", "--no-renames", shQuote("--format=__%H|%P|%cI|%s"), "--name-status", "--", "data")
commits <- list()
current <- NULL
for (line in raw) {
  if (startsWith(line, "__")) {
    parts <- strsplit(sub("^__", "", line), "|", fixed = TRUE)[[1]]
    current <- list(hash = parts[[1]], parent = strsplit(parts[[2]], " ")[[1]][[1]], date = parts[[3]],
                    subject = paste(parts[-(1:3)], collapse = "|"))
  } else if (nzchar(line) && !is.null(current)) {
    status <- substr(line, 1, 1)
    path <- sub("^[A-Z]\t", "", line)
    commits[[length(commits) + 1L]] <- c(current, status = status, path = path)
  }
}
history <- bind_rows(lapply(commits, as.data.frame, stringsAsFactors = FALSE)) %>%
  filter(startsWith(path, "data/"), !startsWith(path, "data/archive/"), path != "data/import_log.csv") %>%
  mutate(relpath = sub("^data/", "", path), day = as.Date(substr(date, 1, 10)))

present <- list.files(data_root, recursive = TRUE)
# Data files only: the refresh status note is a report about the data, not data.
present <- present[!startsWith(present, "archive/") & present != "import_log.csv" & grepl("\\.rds$", present)]
history <- history %>% filter(relpath %in% present)
message(length(present), " files; ", nrow(history), " file-commit records in git.")

read_at <- function(commit, path) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  status <- suppressWarnings(system2("git", c("-C", shQuote(repo), "show", shQuote(paste0(commit, ":", path))), stdout = tmp, stderr = FALSE))
  if (!identical(as.integer(status), 0L)) return(NULL)
  tryCatch(readRDS(tmp), error = function(e) NULL)
}

totals_of <- function(x) {
  column <- value_column_of(x)
  c(rows = if (is.data.frame(x)) nrow(x) else NA, total = value_total(x, column), column = column)
}

totals_at <- function(commit, path) totals_of(read_at(commit, path))

git_date <- function(iso) format(as.POSIXct(sub("([+-]\\d{2}):(\\d{2})$", "\\1\\2", iso), format = "%Y-%m-%dT%H:%M:%S%z", tz = "UTC"), "%Y-%m-%dT%H:%M:%S+0000")

entries <- list()
for (relpath in present) {
  rows <- history[history$relpath == relpath, , drop = FALSE]
  if (nrow(rows) == 0) {
    # Not in git: date it by the file itself.
    info <- file.info(file.path(data_root, relpath))
    x <- tryCatch(readRDS(file.path(data_root, relpath)), error = function(e) NULL)
    column <- value_column_of(x)
    entries[[length(entries) + 1L]] <- data.frame(
      run_id = paste0("file-", format(info$mtime, "%Y-%m-%d")), imported_at = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
      relpath = relpath, action = "added", rows_new = if (is.data.frame(x)) nrow(x) else NA,
      value_column = column, total_new = value_total(x, column), tool = "not in git",
      note = "sem histórico git: data do ficheiro", stringsAsFactors = FALSE
    )
    next
  }

  before <- rows[rows$day < since, , drop = FALSE]
  after <- rows[rows$day >= since, , drop = FALSE]
  after_mods <- after[after$status == "M", , drop = FALSE]

  # The version that stood at `since`, or the file's creation after it.
  base <- if (nrow(before) > 0) before[nrow(before), ] else after[after$status == "A", , drop = FALSE][1, ]
  if (nrow(base) == 0 || is.na(base$hash)) base <- rows[1, ]
  base_totals <- totals_at(base$hash, base$path)
  entries[[length(entries) + 1L]] <- data.frame(
    run_id = paste0("git-", substr(base$hash, 1, 7)), imported_at = git_date(base$date), relpath = relpath,
    action = "added", rows_new = base_totals[["rows"]], value_column = base_totals[["column"]],
    total_new = base_totals[["total"]], tool = "git",
    note = if (nrow(before) > 1) paste0(base$subject, " (construído em ", nrow(before), " commits)") else base$subject,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(after_mods))) {
    mod <- after_mods[i, ]
    if (identical(mod$hash, base$hash)) next
    old_x <- read_at(mod$parent, mod$path)
    new_x <- read_at(mod$hash, mod$path)
    old <- totals_of(old_x)
    new <- totals_of(new_x)
    entries[[length(entries) + 1L]] <- data.frame(
      run_id = paste0("git-", substr(mod$hash, 1, 7)), imported_at = git_date(mod$date), relpath = relpath,
      action = "replaced", rows_old = old[["rows"]], rows_new = new[["rows"]], value_column = new[["column"]],
      total_old = old[["total"]], total_new = new[["total"]],
      archived_as = paste0("git:", mod$parent, ":", mod$path), tool = "git", note = mod$subject,
      rows_changed = rows_changed_count(old_x, new_x, new[["column"]]),
      stringsAsFactors = FALSE
    )
  }
}

log <- bind_rows(entries) %>%
  mutate(dataset = relpath_dataset(relpath), year = relpath_year(relpath)) %>%
  arrange(imported_at, relpath)
for (col in setdiff(import_log_columns, names(log))) log[[col]] <- NA
log <- log[, import_log_columns]

unlink(import_log_path(data_root))
append_import_log(data_root, log)
message("Wrote ", nrow(log), " entries (", sum(log$action == "replaced"), " replacements) to data/import_log.csv")
print(import_runs_summary(read_import_log(data_root)), n = 40, width = 200)
