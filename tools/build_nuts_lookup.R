#!/usr/bin/env Rscript
# Build the municipality -> NUTS II/III lookup used to aggregate regions.
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
#   Rscript tools/build_nuts_lookup.R [indicator=] [out=] [canonical=]
#
# One file per NUTS vintage. The app offers both as a user selection:
#
#   Rscript tools/build_nuts_lookup.R indicator=0013166 out=data/nuts_lookup_2024.rds
#   Rscript tools/build_nuts_lookup.R indicator=0008206 out=data/nuts_lookup_2013.rds \
#     canonical=data/nuts_lookup_2024.rds
#
# Whichever vintage is selected, the same municipality set is summed across the
# whole series, so a regional series stays continuous even though INE's own
# regional rows change definition in 2022. The vintage chooses which grouping
# the municipalities are summed by, not which years are read.
#
# `canonical=` reconciles municipality labels against another lookup, so both
# vintages speak the archive's vocabulary. INE labels the same place
# differently between indicators: 0008206 publishes two municipalities called
# "Calheta", while 0013166 distinguishes "Calheta (R.A.A.)" from
# "Calheta (R.A.M.)", and the archive stores the distinguished form.
#
# Matching is by code first, then by name. The code carries the hierarchy, so
# it is stable only for municipalities whose parents did not change: 1111601 is
# Arcos de Valdevez in both vintages, but every municipality NUTS-2024 moved
# between regions was renumbered with it. 132 of the 308 match by code and the
# remaining 176 by name; the ambiguous island labels are all in the first
# group, so no name is ever matched against two candidates. Anything that
# resolves by neither route is an error rather than a guess.

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
out_path <- get_arg("out", "data/nuts_lookup_2024.rds")
canonical_path <- get_arg("canonical", "")

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

# For each municipality, its parent at a given level is the longest region code
# that prefixes it while staying within that level's code width: 1 character for
# NUTS I ("1" Continente, "2" Açores, "3" Madeira), 2 for NUTS II, 3 for NUTS
# III. The islands carry the same name at all three levels, which is correct -
# they are the same territory - and the mainland is the only place NUTS I says
# something the other levels do not.
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
    nuts1 = parent_of(code, 1L),
    nuts2 = parent_of(code, 2L),
    nuts3 = parent_of(code, 3L)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    municipality = label,
    municipality_code = code,
    nuts1 = nuts1,
    nuts2 = nuts2,
    nuts3 = nuts3,
    indicator = indicator
  )

# Relabel municipalities to match another vintage's vocabulary, by code.
if (nzchar(canonical_path)) {
  if (!file.exists(canonical_path)) {
    stop("canonical lookup not found: ", canonical_path, call. = FALSE)
  }

  canonical <- readRDS(canonical_path) %>%
    dplyr::mutate(
      municipality_code = as.character(municipality_code),
      municipality = as.character(municipality)
    )

  by_code <- canonical %>%
    dplyr::transmute(municipality_code, label_by_code = municipality)

  # A name is only usable when it identifies one place on both sides.
  unique_names <- function(df) {
    df %>%
      dplyr::count(municipality) %>%
      dplyr::filter(n == 1) %>%
      dplyr::pull(municipality)
  }
  usable <- intersect(unique_names(canonical), unique_names(lookup))

  by_name <- canonical %>%
    dplyr::filter(municipality %in% usable) %>%
    dplyr::transmute(municipality, label_by_name = municipality)

  lookup <- lookup %>%
    dplyr::left_join(by_code, by = "municipality_code") %>%
    dplyr::left_join(by_name, by = "municipality") %>%
    dplyr::mutate(
      canonical_label = dplyr::coalesce(label_by_code, label_by_name)
    )

  unmatched <- lookup %>% dplyr::filter(is.na(canonical_label))
  if (nrow(unmatched) > 0) {
    stop(
      "Could not match to ", canonical_path, " by code or by unambiguous name: ",
      paste(unmatched$municipality, collapse = ", "),
      call. = FALSE
    )
  }

  message(
    "Matched ", sum(!is.na(lookup$label_by_code)), " by code and ",
    sum(is.na(lookup$label_by_code)), " by name; ",
    sum(lookup$municipality != lookup$canonical_label), " relabelled."
  )

  lookup <- lookup %>%
    dplyr::mutate(municipality = canonical_label) %>%
    dplyr::select(-label_by_code, -label_by_name, -canonical_label)
}

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
message("NUTS I regions: ", paste(sort(unique(stats::na.omit(lookup$nuts1))), collapse = ", "))
message("NUTS II regions: ", paste(sort(unique(stats::na.omit(lookup$nuts2))), collapse = ", "))
if (length(duplicated_labels) > 0) {
  message("Ambiguous municipality labels: ", paste(duplicated_labels, collapse = ", "))
}
