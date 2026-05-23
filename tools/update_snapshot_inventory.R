#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(memoise)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("tools/update_snapshot_inventory.R", mustWork = FALSE)
}
repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

get_app_dir <- function() {
  repo_dir
}

sys.source(file.path(repo_dir, "R/config.R"), envir = globalenv())
sys.source(file.path(repo_dir, "R/helpers.R"), envir = globalenv())
sys.source(file.path(repo_dir, "R/cache.R"), envir = globalenv())
sys.source(file.path(repo_dir, "R/snapshots.R"), envir = globalenv())

inventory <- write_snapshot_inventory()

cat(glue(
  "Wrote {nrow(inventory)} snapshot inventory rows to {get_snapshot_inventory_file()}\n"
))
cat(glue(
  "Population chunks: {sum(inventory$dataset == 'population')}; death chunks: {sum(inventory$dataset == 'deaths')}\n"
))
