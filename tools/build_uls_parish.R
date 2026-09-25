#!/usr/bin/env Rscript
# Build the parish -> ULS lookup for the municipalities the ULS reform split.
#
#   Rscript tools/build_uls_parish.R [out=data/uls_parish.rds]
#
# Reads data-raw/uls_parish.csv (the assignment, with its sources) and matches
# it to INE's parish list in data/snapshots/census_parish, by normalised name
# within the municipality: INE writes "União das freguesias de Aldoar, Foz do
# Douro e Nevogilde" where the decree and the ULS write the parish names alone.
# The build fails unless every parish of each municipality is matched exactly
# once, so a boundary or naming change is caught here and not in the app.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())
out_path <- get_arg("out", "data/uls_parish.rds")

normalise <- function(x) {
  x <- tolower(iconv(as.character(x), to = "ASCII//TRANSLIT"))
  x <- gsub("^uniao das freguesias d[eao]s? ", "", x)
  gsub("[^a-z0-9]+", "", x)
}

assignment <- utils::read.csv("data-raw/uls_parish.csv", comment.char = "#", encoding = "UTF-8", stringsAsFactors = FALSE)
census <- readRDS("data/snapshots/census_parish/year_2021.rds")
lookup <- readRDS("data/nuts_lookup_2024.rds") %>%
  transmute(dico = substr(as.character(municipality_code), 4, 7), municipality)

parishes <- census %>%
  distinct(code, parish, dico) %>%
  left_join(lookup, by = "dico") %>%
  filter(.data$municipality %in% unique(assignment$municipality)) %>%
  mutate(key = normalise(.data$parish))

wanted <- assignment %>% mutate(key = normalise(.data$parish))
joined <- parishes %>% left_join(wanted[, c("municipality", "key", "uls")], by = c("municipality", "key"))

missing <- joined[is.na(joined$uls), , drop = FALSE]
extra <- wanted[!paste(wanted$municipality, wanted$key) %in% paste(parishes$municipality, parishes$key), , drop = FALSE]
if (nrow(missing) > 0 || nrow(extra) > 0) {
  if (nrow(missing) > 0) message("Parishes in INE with no ULS: ", paste(missing$parish, collapse = "; "))
  if (nrow(extra) > 0) message("Rows in the csv with no INE parish: ", paste(extra$parish, collapse = "; "))
  stop("The parish assignment does not cover the municipalities exactly.", call. = FALSE)
}

out <- joined %>%
  transmute(code = .data$code, parish = .data$parish, dico = .data$dico, municipality = .data$municipality, unit = .data$uls) %>%
  arrange(.data$municipality, .data$unit, .data$parish)

message(sprintf("%d parishes across %d municipalities and %d ULS", nrow(out), n_distinct(out$municipality), n_distinct(out$unit)))
print(out %>% count(municipality, unit), n = 20)
versioned_save_rds(out, out_path, tool = "build_uls_parish.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))
message("Wrote ", out_path)
