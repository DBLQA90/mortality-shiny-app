#!/usr/bin/env Rscript
# Build the map boundaries the planning tab draws, into data/boundaries/.
#
#   Rscript tools/build_boundaries.R
#
# The boundaries come from the geopt package (CRAN, version 0.1.0 of 2026-10-21;
# github.com/c-matos/geopt, MIT): install.packages("geopt"). It
# bundles DGT's official CAOP limits pre-simplified, from parish to NUTS I
# and for the health geography, with the Azores and Madeira composed as insets
# next to the mainland. geopt and sf are needed here only: the output is plain
# coordinate tables, so the app draws them with ggplot2 and needs neither
# package, nor a network connection. Rebuild when CAOP, the NUTS or the ULS
# change.
#
# Output, one file per level: data/boundaries/<level>.rds, a list of
#   polygons  id, name, piece, hole, x, y - id is the app's own area name;
#   frames    group, x, y - the boxes around the island insets;
#   meta      source, geopt version, CAOP attribution, build date.
#
# Levels: municipio (308), nuts2 and nuts3 for both NUTS vintages, and uls (the
# 39 ULS of the mainland, drawn from their parishes, so the ULS that share
# Lisboa, Loures and Porto are exact; the islands' health units have no ULS
# in the app and are kept without an id).
#
# Checked on 2026-10-15 against the app's lookups: NUTS 2013 and 2024 agree
# for every municipality, the parish split of Lisboa, Loures and Porto for all
# 41 parishes. geopt puts Lomba (Gondomar) in ULS Vila Nova de Gaia/Espinho; the
# app keeps Gondomar whole in ULS Santo António (unverified; see METHODOLOGY).

suppressMessages({
  library(sf)
  library(geopt)
})

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
out_dir <- "data/boundaries"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

normalize <- function(x) {
  x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
  x <- gsub("[/-]", " ", x)
  gsub("\\s+", " ", trimws(x))
}

# Polygon rings as a data frame. Every polygon of a feature is a `piece`; its
# first ring is the outline and the others are holes (geom_polygon's subgroup).
polygon_table <- function(map, ids, names) {
  rows <- lapply(seq_len(nrow(map)), function(i) {
    geometry <- sf::st_cast(sf::st_geometry(map)[i], "MULTIPOLYGON")
    coords <- sf::st_coordinates(geometry)
    if (nrow(coords) == 0) return(NULL)
    piece <- paste(i, coords[, "L3"], sep = ".")
    data.frame(id = ids[[i]], name = names[[i]], piece = piece,
               hole = paste(piece, coords[, "L1"], sep = "."),
               x = round(coords[, "X"], 5), y = round(coords[, "Y"], 5), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

frame_table <- function(map) {
  frames <- map[map$geopt_type == "frame", ]
  if (nrow(frames) == 0) return(data.frame(group = character(0), x = numeric(0), y = numeric(0)))
  coords <- sf::st_coordinates(sf::st_cast(sf::st_geometry(frames), "MULTILINESTRING"))
  data.frame(group = paste(coords[, "L2"], coords[, "L1"], sep = "."), x = coords[, "X"], y = coords[, "Y"])
}

write_level <- function(level, map, ids, names) {
  boundary <- map$geopt_type == "boundary"
  value <- list(
    polygons = polygon_table(map[boundary, ], ids[boundary], names[boundary]),
    frames = frame_table(map),
    meta = list(
      source = "CAOP, Direção-Geral do Território, via o pacote geopt",
      attribution = "Limites administrativos: CAOP (DGT)",
      geopt_version = as.character(utils::packageVersion("geopt")),
      built_at = format(Sys.time(), "%Y-%m-%d")
    )
  )
  path <- file.path(out_dir, paste0(level, ".rds"))
  saveRDS(value, path, compress = "xz")
  message(sprintf("%-14s %3d areas, %6d vertices, %4.0f KB", level, length(unique(value$polygons$id)),
                  nrow(value$polygons), file.size(path) / 1024))
  invisible(value)
}

# Municipalities: the app names them as its NUTS lookup does, keyed by DICO.
lookup <- readRDS("data/nuts_lookup_2024.rds")
lookup$dtmn <- substr(lookup$municipality_code, nchar(lookup$municipality_code) - 3L, nchar(lookup$municipality_code))
m <- get_map("municipio")
ids <- lookup$municipality[match(m$dtmn, lookup$dtmn)]
stopifnot(sum(!is.na(ids[m$geopt_type == "boundary"])) == 308)
write_level("municipio", m, ids, ifelse(is.na(ids), m$municipio, ids))

# NUTS II and III, per vintage, named as in the app's lookups.
for (vintage in c("2013", "2024")) {
  lk <- readRDS(sprintf("data/nuts_lookup_%s.rds", vintage))
  for (level in c("nuts2", "nuts3")) {
    map <- get_map(level, nuts_vintage = vintage)
    app_names <- unique(lk[[level]])
    ids <- app_names[match(normalize(map[[level]]), normalize(app_names))]
    unmatched <- setdiff(normalize(app_names), normalize(ids))
    if (length(unmatched) > 0) stop(level, " ", vintage, ": no boundary for ", paste(unmatched, collapse = ", "))
    write_level(paste0(level, "_", vintage), map, ids, ifelse(is.na(ids), map[[level]], ids))
  }
}

# ULS: geopt's names, matched to the app's after normalising punctuation and
# the one name that differs ("ULS Dão-Lafões" is the app's "ULS Viseu Dão-Lafões").
uls_lookup <- readRDS("data/uls_lookup.rds")
app_uls <- unique(uls_lookup$unit[uls_lookup$kind %in% c("ULS", "ULS (partilhada)")])
aliases <- c("uls dao lafoes" = "uls viseu dao lafoes")
u <- get_map("uls")
key <- normalize(u$uls)
key <- ifelse(key %in% names(aliases), aliases[key], key)
ids <- app_uls[match(key, normalize(app_uls))]
missing <- setdiff(app_uls, ids)
if (length(missing) > 0) stop("ULS without a boundary: ", paste(missing, collapse = ", "))
write_level("uls", u, ids, ifelse(is.na(ids), u$uls, ids))
