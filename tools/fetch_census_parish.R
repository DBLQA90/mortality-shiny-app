#!/usr/bin/env Rscript
# Fetch the 2021 census population by parish, age group and sex, for the
# weights that split a municipality between ULS (R/planning_parish.R).
#
#   Rscript tools/fetch_census_parish.R [overwrite=false]
#
# Source: 0012364 (Censos 2021, resident population by place of residence, sex,
# age group and highest completed level of education), read at its parish rows
# (six-character codes) with the education dimension pinned to the total. It is
# the only parish-level population INE publishes by age; there is no annual
# parish estimate, which is why the weights are census shares applied to the
# municipal estimates of each year.
#
# Output: data/snapshots/census_parish/year_2021.rds
#   columns: code (6), parish, dico (4), age_group, sex, pop, source_indicator

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
overwrite <- any(tolower(args) %in% c("overwrite=true", "overwrite=1"))
script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))
sys.source(file.path("R", "data_versions.R"), envir = environment())

path <- "data/snapshots/census_parish/year_2021.rds"
if (file.exists(path) && !overwrite) {
  message("Already present: ", path)
  quit(status = 0)
}

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)
parts <- list()
for (sex in c(T = "T", M = "2")) {
  got <- NULL
  for (attempt in 1:4) {
    got <- tryCatch(client$get_data("0012364", dim1 = "S7A2021", dim3 = sex, dim5 = "T"), error = function(e) NULL)
    if (!is.null(got) && nrow(got) > 0) break
    Sys.sleep(60 * attempt)
  }
  if (is.null(got) || nrow(got) == 0) stop("Could not fetch 0012364 for sex ", sex, call. = FALSE)
  message("sex ", sex, ": ", nrow(got), " rows")
  parts[[sex]] <- dplyr::mutate(got, sex = ifelse(sex == "T", "HM", "M"))
  Sys.sleep(5)
}

raw <- dplyr::bind_rows(parts)
tidy <- tibble::tibble(
  code = as.character(raw$geocod),
  parish = trimws(as.character(raw$geodsg)),
  age_group = trimws(as.character(raw$dim_4_t)),
  sex = raw$sex,
  pop = suppressWarnings(as.numeric(raw$valor)),
  source_indicator = "0012364"
) %>%
  # Parishes are the six-character codes; municipalities and regions are
  # shorter and are rebuilt by summing.
  dplyr::filter(nchar(.data$code) == 6L, .data$age_group != "Total") %>%
  dplyr::mutate(dico = substr(.data$code, 1, 4)) %>%
  dplyr::arrange(.data$code, .data$sex, .data$age_group)

message(sprintf("%d rows, %d parishes, %d municipalities, %d age groups; total population %s",
                nrow(tidy), dplyr::n_distinct(tidy$code), dplyr::n_distinct(tidy$dico), dplyr::n_distinct(tidy$age_group),
                format(sum(tidy$pop[tidy$sex == "HM"], na.rm = TRUE), big.mark = " ")))
if (dplyr::n_distinct(tidy$code) < 2000) stop("Too few parishes: ", dplyr::n_distinct(tidy$code), call. = FALSE)
versioned_save_rds(tidy, path, tool = "fetch_census_parish.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))
message("Wrote ", path)
