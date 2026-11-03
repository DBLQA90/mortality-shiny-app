#!/usr/bin/env Rscript
# The fingerprint of the data the INLA forecasts are built from (for
# tools/scheduled_inla.sh): the municipal deaths by cause and age, the
# population and the death totals snapshots, file by file, and the source
# corrections without their build date, which changes at every refresh.
#
#   Rscript tools/inla_inputs_fingerprint.R
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..", "data")))
files <- sort(list.files(c("snapshots/deaths", "snapshots/population", "snapshots/death_totals"),
                         recursive = TRUE, full.names = TRUE))
sums <- unname(tools::md5sum(files))
corrections <- readRDS("source_corrections.rds")
corrections$notes$built <- NULL
path <- tempfile(fileext = ".rds")
writeBin(serialize(corrections, NULL, version = 3), path)
key <- paste(c(paste(files, sums), unname(tools::md5sum(path))), collapse = "\n")
writeLines(key, path)
cat(unname(tools::md5sum(path)), "\n", sep = "")
