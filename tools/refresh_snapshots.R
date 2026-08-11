#!/usr/bin/env Rscript
# Resumable snapshot refresh driver, intended for unattended CI runs.
#
#   Rscript tools/refresh_snapshots.R [task=all] [minutes=300] [out=data/snapshots]
#
# Tasks:
#   deaths2024  fetch the latest death years missing from the archive
#   nuts2       backfill regional (NUTS II) rows into existing chunks
#   ambiguous   report municipalities INE labels ambiguously (Calheta, Lagoa)
#   inventory   rebuild data/snapshots/snapshot_inventory.rds
#   all         every task above, in that order
#
# The driver is deadline-aware: it stops starting new work once the time budget
# is spent and leaves the archive in a consistent state, because the underlying
# builders write each chunk atomically and skip areas already present. A run
# that is cut short simply makes progress; the next run continues where it
# stopped. That is what makes a scheduled workflow viable against an API as
# slow as INE's.

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

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
repo_root <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

started_at <- Sys.time()
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

run_builder <- function(script, env = character(0), label = script) {
  say("-> ", label)
  status <- tryCatch(
    system2(
      "Rscript",
      args = shQuote(file.path("tools", script)),
      env = c(env, paste0("R_PROGRESSR_ENABLE=FALSE")),
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

task_deaths_latest <- function() {
  say("== Task: latest death years ==")
  missing <- missing_death_years("0013166")

  if (length(missing) == 0) {
    say("No missing death years for 0013166.")
    return(invisible(TRUE))
  }

  say("Missing years: ", paste(missing, collapse = ", "))
  for (year in missing) {
    if (!have_time(20)) break
    run_builder(
      "build_death_snapshot_chunks.R",
      env = c(
        "INDICATOR=0013166",
        paste0("YEARS=", year),
        "AREAS=ALL",
        "CAUSES=ALL",
        "AREA_BATCH_SIZE=40"
      ),
      label = paste0("deaths 0013166 ", year)
    )
  }
  invisible(TRUE)
}

task_nuts2 <- function() {
  say("== Task: regional (NUTS II) rows ==")
  say("Note: the municipal region mode derives regions by summing municipalities ",
      "and needs none of this. These rows are only for 'Dados originais INE'.")

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

task_inventory <- function() {
  say("== Task: rebuild inventory ==")
  run_builder("update_snapshot_inventory.R", label = "snapshot inventory")
}

# ---- Run ------------------------------------------------------------------

say("Refresh started; task=", task, ", budget=", budget_minutes, " min")

# Order matters: the geography repair corrects existing rows and is the highest
# value work, so it runs before any new bulk download can consume the budget.
if (task %in% c("all", "fixareas")) if (have_time(5)) task_fixareas()
if (task %in% c("all", "deaths2024", "deaths")) if (have_time(10)) task_deaths_latest()
if (task %in% c("all", "ambiguous")) if (have_time(5)) task_ambiguous()
if (task %in% c("all", "nuts2")) if (have_time(20)) task_nuts2()
if (task %in% c("all", "inventory")) if (have_time(3)) task_inventory()

elapsed <- round(as.numeric(difftime(Sys.time(), started_at, units = "mins")), 1)
say("Finished in ", elapsed, " min.")

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
    "```"
  ),
  status_path
)
message("Wrote ", status_path)
