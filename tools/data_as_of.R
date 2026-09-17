#!/usr/bin/env Rscript
# Rebuild the app's data as it stood on a given date, to reproduce an earlier
# analysis or to explain why a figure changed.
#
#   Rscript tools/data_as_of.R date=2026-09-01 [out=.mortality-shiny-cache/data_as_of/2026-09-01]
#
# Files unchanged since the date are hard-linked (no extra space where the
# filesystem allows it); files revised after the date are restored from
# data/archive, or from git for revisions recorded before the import log
# existed; files first imported after the date are left out. The snapshot
# inventory is rebuilt for the reconstructed files.
#
# Then run the app on it:
#
#   MORTALITY_SNAPSHOT_DIR=<out>/snapshots Rscript -e 'shiny::runApp()'
#
# The NUTS and ULS lookups are restored into <out> as well, but the app reads
# its lookups from data/ in the repository, so they apply only if copied there.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = "") {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
repo <- normalizePath(file.path(script_dir, ".."))
setwd(repo)

date <- get_arg("date")
if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) stop("Pass date=YYYY-MM-DD.", call. = FALSE)
out <- get_arg("out", file.path(".mortality-shiny-cache", "data_as_of", date))

suppressMessages(library(dplyr))
sys.source(file.path(repo, "R", "data_versions.R"), envir = environment())

if (dir.exists(out)) unlink(out, recursive = TRUE)
result <- data_as_of(date, file.path(repo, "data"), out, repo_dir = repo)

message(sprintf(
  "Data as of %s written to %s\n  unchanged since: %d files\n  restored older version: %d\n  not yet imported then: %d",
  date, normalizePath(out), result$unchanged, result$restored, result$not_yet_imported
))
if (length(result$missing) > 0) {
  message("  could not restore ", length(result$missing), " file(s): ", paste(utils::head(result$missing, 5), collapse = ", "))
}

# The inventory lists which chunks exist; rebuild it for this version.
snapshot_dir <- file.path(out, "snapshots")
if (dir.exists(snapshot_dir)) {
  status <- tryCatch({
    env <- new.env()
    assign("get_app_dir", function() repo, envir = env)
    for (f in c("R/config.R", "R/helpers.R", "R/cache.R", "R/snapshots.R")) sys.source(file.path(repo, f), envir = env)
    inventory_path <- file.path(snapshot_dir, "snapshot_inventory.rds")
    if (file.exists(inventory_path)) unlink(inventory_path)
    inventory <- env$build_snapshot_inventory(snapshot_dir = snapshot_dir)
    saveRDS(inventory, inventory_path, version = 2)
    "rebuilt"
  }, error = function(e) paste("not rebuilt:", conditionMessage(e)))
  message("  snapshot inventory ", status)
}

message("\nRun the app on it with:\n  MORTALITY_SNAPSHOT_DIR=", normalizePath(snapshot_dir), " Rscript -e 'shiny::runApp()'")
