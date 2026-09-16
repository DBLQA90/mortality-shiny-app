#!/usr/bin/env Rscript
# Build the municipality -> ULS / ARS lookup from the PNS2030 planning workbook.
#
#   Rscript tools/build_uls_lookup.R workbook=/path/to/Indicadores_Apoio_PLS_vNN.xlsx \
#     [out=data/uls_lookup.rds]
#
# Source
# ------
# The `Var` sheet of "Indicadores de Apoio ao Planeamento Local em Saúde"
# (Equipa DRS / Equipa PNS2030) lists every ULS with its municipalities, keyed
# by INE's 2024 municipality codes ("refgeo2024"). The workbook itself is not
# committed; only the derived lookup is.
#
# The ARS regions are not listed per municipality in the workbook. They follow
# from the first digit of the ULS code - 1 Norte, 2 Centro, 3 Lisboa e Vale do
# Tejo, 4 Alentejo, 5 Algarve - which is how the workbook's own regional rows
# group the ULS.
#
# Split municipalities
# --------------------
# Three municipalities are divided between two ULS at parish level: Lisboa
# (Santa Maria, São José), Loures (Loures/Odivelas, São José) and Porto (Santo
# António, São João). Nothing below municipality level exists in the app's data,
# so the five ULS that touch them cannot be built individually. The smallest
# unions of those ULS that contain only whole municipalities are exact, and are
# offered instead:
#
#   ULS Santo António + São João                  Gondomar, Maia, Porto, Valongo
#   ULS Loures/Odivelas + São José + Santa Maria  Lisboa, Loures, Mafra, Odivelas
#
# Those groups are derived here from the mapping, not hard-coded, so a future
# workbook that moves a parish boundary produces the right groups automatically.

suppressMessages({
  library(dplyr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = "") {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

workbook <- get_arg("workbook")
out_path <- get_arg("out", "data/uls_lookup.rds")
if (!nzchar(workbook) || !file.exists(workbook)) {
  stop("workbook= must point to the planning indicators .xlsx", call. = FALSE)
}

var <- suppressMessages(read_excel(workbook, sheet = "Var", col_names = FALSE, col_types = "text"))
names(var) <- paste0("c", seq_len(ncol(var)))

header <- which(var$c1 == "cod_ULS" & var$c3 == "cod1_conc")[1]
if (is.na(header)) stop("Could not find the ULS-municipality table in the Var sheet", call. = FALSE)

ars_names <- c(
  "1" = "ARS Norte",
  "2" = "ARS Centro",
  "3" = "ARS Lisboa e Vale do Tejo",
  "4" = "ARS Alentejo",
  "5" = "ARS Algarve"
)

# The workbook spells one ULS two ways ("Loures-Odivelas" and
# "Loures/Odivelas"); settle on the slash form it uses in its data sheets. Only
# that name: a blanket hyphen-to-slash rule would also rewrite genuinely
# hyphenated names such as "Viseu Dão-Lafões".
normalise_uls <- function(x) sub("^ULS Loures\\s*-\\s*Odivelas$", "ULS Loures/Odivelas", x)

mapping <- var[(header + 1):nrow(var), ] %>%
  filter(!is.na(c3), !is.na(c5)) %>%
  transmute(
    uls_code = as.character(c1),
    uls = normalise_uls(trimws(c2)),
    municipality_code = as.character(c3),
    municipality = trimws(c5)
  ) %>%
  mutate(ars_code = substr(uls_code, 1, 1), ars = unname(ars_names[ars_code]))

if (any(is.na(mapping$ars))) stop("Unexpected ULS code prefix", call. = FALSE)

# Municipality names must be the app's; the codes decide.
nuts <- readRDS("data/nuts_lookup_2024.rds")
mapping <- mapping %>%
  select(-municipality) %>%
  inner_join(nuts %>% select(municipality_code, municipality), by = "municipality_code")

split <- mapping %>% count(municipality) %>% filter(n > 1) %>% pull(municipality)

# Smallest unions of ULS containing only whole municipalities.
touching <- unique(mapping$uls[mapping$municipality %in% split])
groups <- list()
left <- touching
while (length(left) > 0) {
  group <- left[[1]]
  repeat {
    members <- mapping$municipality[mapping$uls %in% group]
    add <- setdiff(unique(mapping$uls[mapping$municipality %in% intersect(members, split)]), group)
    if (length(add) == 0) break
    group <- c(group, add)
  }
  groups[[length(groups) + 1]] <- sort(group)
  left <- setdiff(left, group)
}

group_label <- function(group) {
  paste0("ULS ", paste(sub("^ULS ", "", group), collapse = " + "))
}

units <- bind_rows(
  mapping %>%
    filter(!uls %in% touching) %>%
    distinct(unit = uls, municipality, ars) %>%
    mutate(kind = "ULS"),
  bind_rows(lapply(groups, function(group) {
    mapping %>%
      filter(uls %in% group) %>%
      distinct(municipality, ars) %>%
      mutate(unit = group_label(group), kind = "ULS (grupo)")
  })),
  mapping %>%
    distinct(unit = ars, municipality, ars) %>%
    mutate(kind = "ARS")
) %>%
  distinct(unit, kind, municipality, .keep_all = TRUE) %>%
  arrange(kind, unit, municipality)

attr(units, "split_municipalities") <- split
attr(units, "split_groups") <- vapply(groups, group_label, character(1))
attr(units, "source") <- basename(workbook)
attr(units, "built_at") <- Sys.time()

saveRDS(units, out_path)

message("Wrote ", out_path)
message("  ULS (exact): ", n_distinct(units$unit[units$kind == "ULS"]))
message("  ULS groups:  ", paste(attr(units, "split_groups"), collapse = " | "))
message("  ARS:         ", n_distinct(units$unit[units$kind == "ARS"]))
message("  split municipalities: ", paste(split, collapse = ", "))
message("  municipalities covered: ", n_distinct(units$municipality))
