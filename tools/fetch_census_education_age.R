#!/usr/bin/env Rscript
# Fetch the census population by age group and highest completed level of
# education, per municipality, for the education indicator restricted to an age
# (I24 "com N e mais anos").
#
#   Rscript tools/fetch_census_education_age.R [overwrite=false]
#
# Only the 2011 and 2021 censuses publish this by age at municipal level
# (0006350, 0012364); the 1991-2021 series used for the whole population
# (0014380, census_education) has no age. Levels are harmonised with that
# series: post-secondary counts as secondary.
#
# Output: data/snapshots/planning_extra/census_education_by_age/year_<year>.rds
#   columns: year, area, age, category, value, source_indicator
#   `age` is the lower bound of the group (0 for "Menos de 15 anos", 75 for
#   "75 ou mais anos"); `category` is Total, Nenhum, Básico, Secundário or
#   Superior.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
overwrite <- any(tolower(args) %in% c("overwrite=true", "overwrite=1"))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))

dico_lookup <- readRDS("data/nuts_lookup_2024.rds") %>%
  transmute(dico = substr(as.character(municipality_code), 4, 7), municipality)
area_for_code <- function(code) {
  code <- as.character(code)
  out <- rep(NA_character_, length(code))
  out[code == "PT"] <- "Portugal"
  out[code == "1"] <- "Continente"
  # Municipal codes: NUTS III + DICO (7), one character + DICO (5), or the
  # bare DICO (4, the 2021 census by age); parishes (6) are left out.
  municipal <- nchar(code) %in% c(4L, 5L, 7L)
  out[municipal] <- dico_lookup$municipality[match(substr(code[municipal], nchar(code[municipal]) - 3L, nchar(code[municipal])), dico_lookup$dico)]
  out
}

sys.source(file.path("R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "fetch_census_education_age.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

# Level codes per edition, mapped to the harmonised categories.
EDITIONS <- list(
  `2011` = list(id = "0006350", levels = c(T = "Total", A = "Nenhum", `2` = "Básico", `3` = "Secundário", `4` = "Secundário", `5` = "Superior")),
  `2021` = list(id = "0012364", levels = c(T = "Total", `1` = "Nenhum", `2` = "Básico", `3` = "Secundário", `4` = "Secundário", `5` = "Superior"))
)

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)
failed <- character(0)
for (year in names(EDITIONS)) {
  spec <- EDITIONS[[year]]
  path <- file.path("data/snapshots/planning_extra/census_education_by_age", paste0("year_", year, ".rds"))
  if (file.exists(path) && !overwrite) next
  # One request per level: the API takes a single category per dimension.
  parts <- list()
  for (level in names(spec$levels)) {
    got <- NULL
    for (attempt in 1:4) {
      got <- tryCatch(client$get_data(spec$id, dim1 = paste0("S7A", year), dim3 = "T", dim5 = level), error = function(e) NULL)
      if (!is.null(got) && nrow(got) > 0) break
      Sys.sleep(60 * attempt)
    }
    if (is.null(got) || nrow(got) == 0) { parts <- NULL; break }
    message(sprintf("  %s level %s: %d rows", year, level, nrow(got)))
    parts[[level]] <- got
    Sys.sleep(5)
  }
  if (is.null(parts)) { message(year, ": FAILED"); failed <- c(failed, year); next }
  raw <- dplyr::bind_rows(lapply(parts, function(x) dplyr::mutate(x, dplyr::across(dplyr::everything(), as.character))))

  age_label <- trimws(as.character(raw$dim_4_t))
  age <- ifelse(age_label == "Total", NA_integer_,
                ifelse(grepl("^Menos de 15", age_label), 0L, suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", age_label)))))
  chunk <- tibble::tibble(
    code = as.character(raw$geocod), age = age, age_label = age_label,
    level = as.character(raw$dim_5), value = suppressWarnings(as.numeric(raw$valor))
  ) %>%
    filter(age_label != "Total", level %in% names(spec$levels)) %>%
    mutate(area = area_for_code(code), category = unname(spec$levels[level])) %>%
    filter(!is.na(area)) %>%
    group_by(area, age, category) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(year = as.integer(year), source_indicator = spec$id) %>%
    select(year, area, age, category, value, source_indicator) %>%
    arrange(area, age, category)

  n_mun <- n_distinct(setdiff(chunk$area, c("Portugal", "Continente")))
  pt <- chunk[chunk$area == "Portugal", ]
  parts <- sum(pt$value[pt$category != "Total"]); total <- sum(pt$value[pt$category == "Total"])
  message(sprintf("%s %s: %d municipalities, %d age groups; Portugal total %s, levels sum %s",
                  year, spec$id, n_mun, n_distinct(chunk$age), format(total, big.mark = " "), format(parts, big.mark = " ")))
  if (n_mun < 300) { message("  REJECTED"); failed <- c(failed, year); next }
  save_rds_atomic(chunk, path)
}
message("Done. Failed: ", length(failed))
