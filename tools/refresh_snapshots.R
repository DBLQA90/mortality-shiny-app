#!/usr/bin/env Rscript
# Resumable snapshot refresh driver, intended for unattended CI runs.
#
#   Rscript tools/refresh_snapshots.R [task=all] [minutes=300] [recent=2] [note=...]
#
# Tasks (each fetches what is missing, and re-checks the last `recent` calendar
# years for INE revisions):
#   deaths         deaths by cause and age (0013166)
#   population     resident population by age and sex (0012918)
#   deathtotals    municipal death totals by cause, all ages
#   regional       INE's regional death rows (NUTS I/II/III)
#   infant         live births, under-1 deaths by cause, complete under-1 counts
#   planning       socio-economic, birth and neonatal components of the planning tab
#   weekly         INE weekly deaths by NUTS III and age (0012100, 0010112)
#   parish         births and deaths by parish, for the ULS that share a
#                  municipality (0012450/0012542 and earlier editions)
#   sns            primary-care indicators from the SNS Transparency portal
#   validation     INE's own municipal indicators, which tools/regression_check.R
#                  compares the app's values against
#   current        weekly + sns + the latest deaths and population: what the
#                  scheduled refresh runs (tools/scheduled_refresh.sh)
#   ambiguous      report municipalities INE labels ambiguously (Calheta, Lagoa)
#   inventory      rebuild data/snapshots/snapshot_inventory.rds
#   all            every task above, in that order
# Run explicitly only:
#   fixareas       re-run the Lisboa/Calheta/Lagoa repair of old death chunks
#   nuts2          backfill regional rows into area-by-area death chunks
#
# Every file is written through R/data_versions.R. All tools launched by one run
# share its run id (DATA_RUN_ID), so the import log groups the run's changes: a
# re-fetched file identical to the stored one is left untouched, a revised one is
# archived under data/archive/<run id>/ before being replaced. REFRESH_STATUS.md
# ends with what the run changed. `recent=0` fetches missing years only.
#
# The driver is deadline-aware: it stops starting new work once the time budget
# is spent and leaves the archive in a consistent state, because the underlying
# builders write each chunk atomically and skip areas already present. A run
# that is cut short simply makes progress; the next run continues where it
# stopped. That is what makes a scheduled workflow viable against an API as
# slow as INE's. Tasks run one at a time: INE answers parallel requests with
# "429 Too Many Requests" and then blocks the address for hours.

suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  value <- if (length(hit) == 0) Sys.getenv(toupper(name), unset = "") else sub(paste0("^", name, "="), "", hit[[1]])
  if (!nzchar(value)) default else value
}

task <- tolower(get_arg("task", "all"))
budget_minutes <- suppressWarnings(as.numeric(get_arg("minutes", "300")))
if (!is.finite(budget_minutes) || budget_minutes <= 0) budget_minutes <- 300
out_dir <- get_arg("out", "data/snapshots")
recent <- suppressWarnings(as.integer(get_arg("recent", "2")))
if (is.na(recent) || recent < 0) recent <- 2L
run_note <- get_arg("note", paste0("refresh task=", task))

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
repo_root <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

started_at <- Sys.time()

# One run id for every tool this run launches; child processes inherit it.
run_id <- format(started_at, "%Y-%m-%dT%H%M%S")
Sys.setenv(DATA_RUN_ID = run_id, DATA_RUN_NOTE = run_note)
sys.source(file.path(repo_root, "R", "data_versions.R"), envir = environment())
recheck_from <- as.integer(format(Sys.Date(), "%Y")) - recent
deadline <- started_at + budget_minutes * 60
log_lines <- character(0)

say <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  message(line)
  log_lines <<- c(log_lines, line)
}

minutes_left <- function() as.numeric(difftime(deadline, Sys.time(), units = "mins"))

# Every task checks this before starting, and passes the remaining budget down
# so a long builder run cannot overshoot the workflow's own timeout.
have_time <- function(need_minutes = 5) {
  left <- minutes_left()
  if (left < need_minutes) {
    say("Time budget exhausted (", round(left, 1), " min left); stopping cleanly.")
    return(FALSE)
  }
  TRUE
}

run_builder <- function(script, env = character(0), label = script, args = character(0)) {
  say("-> ", label)

  # system2(env=) prepends NAME=value to the command line, where the shell word
  # splits it. Region labels contain spaces ("Oeste e Vale do Tejo"), which the
  # shell then tried to execute: "sh: 1: Oeste: not found". Quote the value of
  # each assignment, leaving the name and the "=" alone.
  quoted_env <- vapply(as.character(env), function(assignment) {
    at <- regexpr("=", assignment, fixed = TRUE)
    if (at < 1) return(assignment)
    paste0(substr(assignment, 1, at), shQuote(substring(assignment, at + 1)))
  }, character(1), USE.NAMES = FALSE)

  status <- tryCatch(
    system2(
      "Rscript",
      args = c(shQuote(file.path("tools", script)), shQuote(args)),
      env = c(quoted_env, "R_PROGRESSR_ENABLE=FALSE"),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )

  exit_code <- attr(status, "status")
  tail_lines <- utils::tail(as.character(status), 15)
  log_lines <<- c(log_lines, paste0("   ", tail_lines))

  if (!is.null(exit_code) && exit_code != 0) {
    say("   FAILED (exit ", exit_code, ")")
    return(FALSE)
  }

  say("   done")
  TRUE
}

# ---- Which death years does the archive still lack? -----------------------
# Asks INE what exists rather than hard-coding a year, so the same workflow
# keeps working when 2025 is published.
missing_death_years <- function(indicator = "0013166") {
  available <- tryCatch(
    {
      client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)
      years <- client$get_dim_values(indicator) %>%
        dplyr::filter(as.integer(dim_num) == 1) %>%
        dplyr::pull(categ_dsg)
      sort(unique(suppressWarnings(as.integer(as.character(years)))))
    },
    error = function(e) {
      say("Could not read INE years for ", indicator, ": ", conditionMessage(e))
      integer(0)
    }
  )

  present <- list.dirs(file.path(out_dir, "deaths", indicator), recursive = FALSE, full.names = FALSE)
  present <- suppressWarnings(as.integer(sub("^year_", "", present)))
  present <- present[!is.na(present)]

  setdiff(available[is.finite(available)], present)
}

# ---- Regional labels, per indicator vintage -------------------------------
# NUTS II labels differ between vintages (NUTS-2013 "Area Metropolitana de
# Lisboa" vs NUTS-2024 "Grande Lisboa"), so each indicator is asked for its own
# level-3/4 geographies instead of reusing one hard-coded list.
regional_labels <- function(indicator) {
  tryCatch(
    {
      client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)
      client$get_dim_values(indicator) %>%
        dplyr::filter(as.integer(dim_num) == 2) %>%
        dplyr::mutate(level = suppressWarnings(as.integer(categ_nivel))) %>%
        dplyr::filter(level %in% c(3L, 4L), nchar(as.character(categ_cod)) <= 2) %>%
        dplyr::pull(categ_dsg) %>%
        as.character() %>%
        setdiff(c("Total", "Ignorado", "Estrangeiro")) %>%
        unique()
    },
    error = function(e) {
      say("Could not read INE geographies for ", indicator, ": ", conditionMessage(e))
      character(0)
    }
  )
}

# ---- Tasks ----------------------------------------------------------------

# Years INE publishes for an indicator.
ine_years <- function(indicator) {
  tryCatch(
    {
      client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)
      years <- client$get_dim_values(indicator) %>%
        dplyr::filter(as.integer(dim_num) == 1) %>%
        dplyr::pull(categ_dsg)
      sort(unique(suppressWarnings(as.integer(as.character(years)))))
    },
    error = function(e) {
      say("Could not read INE years for ", indicator, ": ", conditionMessage(e))
      integer(0)
    }
  )
}

task_deaths_latest <- function() {
  say("== Task: deaths by cause and age ==")
  missing <- missing_death_years("0013166")
  published <- ine_years("0013166")
  recheck <- setdiff(published[published >= recheck_from], missing)

  if (length(missing) == 0 && length(recheck) == 0) {
    say("Nothing to fetch or re-check for 0013166.")
    return(invisible(TRUE))
  }

  say("Missing years: ", if (length(missing)) paste(missing, collapse = ", ") else "none",
      "; re-checking: ", if (length(recheck)) paste(recheck, collapse = ", ") else "none")
  for (year in sort(c(missing, recheck))) {
    if (!have_time(10)) break
    # fetch_death_year.R issues one request per cause across all areas, which
    # measured no slower than a single-area request: 66 calls for a year rather
    # than the 20,394 that per-(area, cause) slicing would need.
    run_builder(
      "fetch_death_year.R",
      env = c("INDICATOR=0013166", paste0("YEAR=", year), paste0("MINUTES=", max(5, floor(minutes_left() - 3))),
              if (year %in% recheck) "OVERWRITE=true"),
      label = paste0("deaths 0013166 ", year, if (year %in% recheck) " (re-check)" else "")
    )
  }
  invisible(TRUE)
}

task_population <- function() {
  say("== Task: population ==")
  published <- ine_years("0012918")
  present <- suppressWarnings(as.integer(sub("^year_(\\d+)\\.rds$", "\\1",
    list.files(file.path(out_dir, "population"), pattern = "^year_\\d+\\.rds$"))))
  years <- sort(union(setdiff(published, present), published[published >= recheck_from]))
  if (length(years) == 0) {
    say("Nothing to fetch or re-check for 0012918.")
    return(invisible(TRUE))
  }
  invisible(run_builder(
    "fetch_population_year.R",
    args = c("indicator=0012918", paste0("years=", paste(years, collapse = ",")), "overwrite=true"),
    label = paste0("population 0012918 ", paste(years, collapse = ", "))
  ))
}

task_simple <- function(script, label, extra = character(0)) {
  invisible(run_builder(script, args = c(paste0("recent=", recent), extra), label = label))
}

task_death_totals <- function() {
  say("== Task: municipal death totals ==")
  task_simple("fetch_death_totals.R", "death totals")
}

task_regional <- function() {
  say("== Task: regional death rows ==")
  task_simple("fetch_regional_deaths.R", "regional deaths")
}

task_planning <- function() {
  say("== Task: planning-tab components ==")
  task_simple("fetch_planning_extra.R", "planning extra (RSI, pensions, purchasing power, waste, births, neonatal)")
}

task_weekly <- function() {
  say("== Task: weekly deaths ==")
  invisible(run_builder("fetch_weekly_deaths.R", label = "weekly deaths 0012100 / 0010112"))
}

task_parish <- function() {
  say("== Task: births and deaths by parish ==")
  invisible(run_builder("fetch_parish_vitals.R", label = "parish births and deaths"))
}

task_sns <- function() {
  say("== Task: SNS Transparency portal ==")
  invisible(run_builder("fetch_sns.R", label = "SNS primary care"))
}

task_nuts2 <- function() {
  say("== Task: regional (NUTS II) rows ==")
  say("Note: rarely needed. Years fetched by fetch_death_year.R already include ",
      "INE's regional rows, because the all-areas response carries every NUTS ",
      "level; and the municipal region mode derives regions by summing ",
      "municipalities. This task only backfills regions into older chunks that ",
      "were built area-by-area.")

  for (indicator in c("0013166", "0008206")) {
    if (!have_time(25)) break
    labels <- regional_labels(indicator)
    if (length(labels) == 0) next

    say(indicator, " regions: ", paste(labels, collapse = " | "))
    run_builder(
      "build_death_snapshot_chunks.R",
      env = c(
        paste0("INDICATOR=", indicator),
        "YEARS=ALL",
        paste0("AREAS=", paste(labels, collapse = "|")),
        "CAUSES=ALL",
        "AREA_BATCH_SIZE=12",
        paste0("MAX_BATCHES=", max(1L, floor(minutes_left() / 4)))
      ),
      label = paste0("regional deaths ", indicator)
    )
  }

  if (have_time(25)) {
    labels <- regional_labels("0008273")
    if (length(labels) > 0) {
      run_builder(
        "build_population_snapshot_chunks.R",
        env = c(
          "YEARS=ALL",
          paste0("AREAS=", paste(labels, collapse = "|")),
          "AREA_BATCH_SIZE=12",
          paste0("MAX_BATCHES=", max(1L, floor(minutes_left() / 4)))
        ),
        label = "regional population"
      )
    }
  }
  invisible(TRUE)
}

# Calheta (Azores / Madeira) and Lagoa (Algarve / Azores) share a label in the
# older indicators, so they cannot be resolved by name. Rather than guess -
# which would silently attribute deaths to the wrong island - this reports the
# state so the fix can be made deliberately.
task_ambiguous <- function() {
  say("== Task: ambiguous municipality labels ==")

  report <- lapply(c("0008206", "0013166", "0008273", "0003182"), function(indicator) {
    labels <- tryCatch(
      {
        client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)
        client$get_dim_values(indicator) %>%
          dplyr::filter(as.integer(dim_num) == 2) %>%
          dplyr::mutate(label = as.character(categ_dsg), code = as.character(categ_cod)) %>%
          dplyr::filter(grepl("^(Calheta|Lagoa)", label)) %>%
          dplyr::select(code, label)
      },
      error = function(e) NULL
    )

    if (is.null(labels) || nrow(labels) == 0) {
      return(NULL)
    }

    ambiguous <- labels %>% dplyr::count(label) %>% dplyr::filter(n > 1) %>% dplyr::pull(label)
    say(indicator, ": ", paste(paste0(labels$label, " (", labels$code, ")"), collapse = "; "),
        if (length(ambiguous) > 0) paste0("  <-- AMBIGUOUS: ", paste(ambiguous, collapse = ", ")) else "  <-- unambiguous")

    labels %>% dplyr::mutate(indicator = indicator, ambiguous = label %in% ambiguous)
  })

  report <- dplyr::bind_rows(report)
  if (nrow(report) > 0) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(report, file.path(out_dir, "ambiguous_areas.rds"))
    say("Wrote ", file.path(out_dir, "ambiguous_areas.rds"))
  }
  invisible(TRUE)
}

# Repairs geographies that were resolved through an ambiguous INE label and
# stored as the sum of two different places. See tools/fix_ambiguous_areas.R for
# the full diagnosis; the headline case is "Lisboa", where 1991-2013 population
# is region + municipio and mortality is understated roughly six-fold.
task_fixareas <- function() {
  say("== Task: repair ambiguous geographies ==")
  run_builder(
    "fix_ambiguous_areas.R",
    env = c(paste0("MINUTES=", max(5, floor(minutes_left() - 5)))),
    label = "fix ambiguous areas (Lisboa, Calheta, Lagoa)"
  )
}

# Live births and under-1 deaths back the infant mortality rate. Both are small
# - one request per year each - so they are cheap to keep current.
task_infant <- function() {
  say("== Task: births and under-1 deaths ==")
  run_builder(
    "fetch_births.R",
    env = c("YEARS=ALL", paste0("MINUTES=", max(5, floor(minutes_left() / 3)))),
    args = paste0("recent=", recent),
    label = "live births"
  )
  if (have_time(5)) {
    run_builder(
      "fetch_infant_deaths.R",
      env = c("YEARS=ALL", paste0("MINUTES=", max(5, floor(minutes_left() - 5)))),
      args = paste0("recent=", recent),
      label = "under-1 deaths by cause"
    )
  }
  if (have_time(5)) task_simple("fetch_infant_totals.R", "complete under-1 death counts")
  invisible(TRUE)
}

task_validation <- function() {
  say("== Task: INE reference indicators for validation ==")
  invisible(run_builder("fetch_validation_refs.R", label = "municipal reference indicators"))
}

task_inventory <- function() {
  say("== Task: rebuild inventory ==")
  invisible(run_builder("update_snapshot_inventory.R", label = "snapshot inventory"))
}

# ---- Run ------------------------------------------------------------------

say("Refresh started; task=", task, ", budget=", budget_minutes, " min")

# The Lisboa/Calheta/Lagoa repair is complete (2026-08-11); it no longer runs
# as part of "all", where it re-read 8,580 chunks every time.
if (identical(task, "fixareas")) if (have_time(5)) task_fixareas()
if (task %in% c("all", "current", "weekly")) if (have_time(5)) task_weekly()
if (task %in% c("all", "current", "parish")) if (have_time(5)) task_parish()
if (task %in% c("all", "current", "sns")) if (have_time(5)) task_sns()
if (task %in% c("all", "current", "deaths2024", "deaths")) if (have_time(10)) task_deaths_latest()
if (task %in% c("all", "current", "population")) if (have_time(5)) task_population()
if (task %in% c("all", "current", "deathtotals")) if (have_time(10)) task_death_totals()
if (task %in% c("all", "regional")) if (have_time(5)) task_regional()
if (task %in% c("all", "current", "infant")) if (have_time(10)) task_infant()
if (task %in% c("all", "planning")) if (have_time(10)) task_planning()
if (task %in% c("all", "validation")) if (have_time(5)) task_validation()
if (task %in% c("all", "ambiguous")) if (have_time(5)) task_ambiguous()
# Not part of "all": regions are built by summing municipalities, so INE's own
# regional rows are needed only to reproduce a published regional figure via
# MORTALITY_REGION_MODE=original. Request it explicitly when that is the aim.
if (identical(task, "nuts2")) if (have_time(20)) task_nuts2()
if (task %in% c("all", "inventory")) if (have_time(3)) task_inventory()

elapsed <- round(as.numeric(difftime(Sys.time(), started_at, units = "mins")), 1)
say("Finished in ", elapsed, " min.")

# What this run changed, from the import log.
import_log <- read_import_log(file.path(repo_root, "data"))
run_rows <- import_log[import_log$run_id == run_id, , drop = FALSE]
changes <- run_changes(import_log, run_id)
change_lines <- c(
  "",
  "## Changes in this run",
  "",
  paste0("Run id: `", run_id, "`  |  Re-checked years from ", recheck_from, "  |  ",
         sum(run_rows$action == "added"), " file(s) added, ", sum(run_rows$action == "replaced"), " replaced"),
  ""
)
if (nrow(changes) > 0) {
  change_lines <- c(
    change_lines,
    paste0("Revised by INE since the previous import (previous versions in `data/archive/", run_id, "/`):"),
    "",
    "Totals are Portugal, both sexes, all causes where the files carry them; a",
    "correction between municipalities shows in the changed rows, not the total.",
    "",
    "| Dataset | Year | Files | Rows changed | Total before | Total after | Change |",
    "|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %d | %s | %s | %s | %s |", changes$dataset, changes$year, changes$files, format(changes$rows_changed, big.mark = " "),
            format(changes$total_old, big.mark = " "), format(changes$total_new, big.mark = " "),
            ifelse(is.na(changes$change_pct), "", sprintf("%+.2f%%", changes$change_pct)))
  )
} else {
  change_lines <- c(change_lines, "No stored value was revised.")
}

status_path <- file.path(out_dir, "REFRESH_STATUS.md")
dir.create(dirname(status_path), recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    "# Snapshot refresh status",
    "",
    paste0("Last run: ", format(started_at, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Task: `", task, "`  |  Duration: ", elapsed, " min"),
    "",
    "```",
    log_lines,
    "```",
    change_lines
  ),
  status_path
)
message("Wrote ", status_path)
