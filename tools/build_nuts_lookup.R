#!/usr/bin/env Rscript
# Build the municipality -> NUTS II lookup used by the "municipal
# approximation" region mode.
#
# INE's geography codes are hierarchical: a municipality carries a 7-character
# code whose leading characters are its NUTS III and NUTS II parents (for
# example 1111601 Arcos de Valdevez -> 111 Alto Minho -> 11 Norte). Regions can
# therefore be rebuilt from municipalities by prefix, without needing a
# hand-maintained list.
#
# The lookup is written once and committed, so the app never needs a network
# call to aggregate a region:
#
#   Rscript tools/build_nuts_lookup.R [indicator=0013166] [out=data/nuts_lookup.rds]
#
# Default indicator is 0013166, the current NUTS-2024 deaths indicator. Using
# one vintage for every year is the point: the same municipality set is summed
# across the whole series, so the regional series stays continuous even though
# INE's own regional rows change definition in 2022.

suppressMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

indicator <- get_arg("indicator", "0013166")
out_path <- get_arg("out", "data/nuts_lookup.rds")

message("Reading geography metadata for indicator ", indicator, " ...")
client <- ineptr2::INEClient$new(lang = "PT", timeout = 600)

geo <- client$get_dim_values(indicator) %>%
  dplyr::filter(as.integer(dim_num) == 2) %>%
  dplyr::mutate(
    code = as.character(categ_cod),
    label = as.character(categ_dsg),
    level = suppressWarnings(as.integer(categ_nivel))
  ) %>%
  dplyr::filter(!is.na(level))

# INE lists "Ignorado" and "Estrangeiro" alongside the real geographies at
# every level. They are not places and must not join a regional aggregate.
pseudo_areas <- c("Total", "Ignorado", "Estrangeiro")

municipalities <- geo %>%
  dplyr::filter(level == max(level, na.rm = TRUE), !label %in% pseudo_areas)
regions <- geo %>%
  dplyr::filter(
    level < max(geo$level, na.rm = TRUE),
    nchar(code) <= 3,
    !label %in% pseudo_areas
  )

message("Found ", nrow(municipalities), " municipalities and ", nrow(regions), " candidate regions.")

# For each municipality, the NUTS II parent is the longest region code that
# prefixes it while still being shorter than a NUTS III code. Islands sit
# directly under a single-character code ("2" Açores, "3" Madeira), mainland
# municipalities under a two-character one, so the prefix search handles both
# without special-casing.
parent_of <- function(muni_code, max_nchar) {
  # startsWith() recycles over the whole candidate vector; substr() would not,
  # because it returns one element per `x` and ignores the remaining widths.
  candidates <- regions %>%
    dplyr::filter(nchar(code) <= max_nchar, startsWith(muni_code, code))

  if (nrow(candidates) == 0) {
    return(NA_character_)
  }

  candidates %>%
    dplyr::arrange(dplyr::desc(nchar(code))) %>%
    dplyr::slice(1) %>%
    dplyr::pull(label)
}

lookup <- municipalities %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    nuts2 = parent_of(code, 2L),
    nuts3 = parent_of(code, 3L)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    municipality = label,
    municipality_code = code,
    nuts2 = nuts2,
    nuts3 = nuts3,
    indicator = indicator
  )

unresolved <- lookup %>% dplyr::filter(is.na(nuts2))
if (nrow(unresolved) > 0) {
  warning(
    "Municipalities without a NUTS II parent: ",
    paste(unresolved$municipality, collapse = ", "),
    call. = FALSE
  )
}

# Municipality labels that INE repeats (Calheta, Lagoa) cannot be resolved by
# name alone. Record them so the app can warn instead of silently aggregating
# the wrong island.
duplicated_labels <- lookup %>%
  dplyr::count(municipality) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(municipality)

attr(lookup, "duplicated_labels") <- duplicated_labels
attr(lookup, "indicator") <- indicator
attr(lookup, "built_at") <- Sys.time()

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(lookup, out_path)

message("Wrote ", out_path, " (", nrow(lookup), " municipalities).")
message("NUTS II regions: ", paste(sort(unique(stats::na.omit(lookup$nuts2))), collapse = ", "))
if (length(duplicated_labels) > 0) {
  message("Ambiguous municipality labels: ", paste(duplicated_labels, collapse = ", "))
}
