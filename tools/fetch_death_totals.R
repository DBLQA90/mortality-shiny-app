#!/usr/bin/env Rscript
# Fetch municipal death totals (all ages) by cause into data/snapshots/death_totals.
#
#   Rscript tools/fetch_death_totals.R [indicator=all] [years=ALL] [overwrite=false] [recent=0]
#
# Why: INE publishes each municipality's cause-specific deaths with a complete
# total but an incomplete breakdown by age, and the main death archive keeps only
# the age bands (it discards the "Total" row). Any figure that does not need ages
# - a count of deaths, a crude rate, proportional mortality at all ages - should
# come from the total instead. For lung cancer in 2014 the municipal age bands sum
# to 2,921 deaths; the totals sum to 4,193, against a national 4,288. Guarda had 14
# lung-cancer deaths that year; its age bands show none.
#
# One request per year returns every municipality, cause and sex.
#
# Municipalities are identified by geography code, mapped to the app's canonical
# names through the lookup of the matching vintage: 0008206 publishes two
# municipalities called "Calheta" and two called "Lagoa", and its codes are the
# NUTS-2013 ones. Portugal and Continente are kept as published.

suppressMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^", name, "="), "", hit[[1]])
}

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))

indicator_arg <- get_arg("indicator", "all")
years_arg <- get_arg("years", "ALL")
overwrite <- tolower(get_arg("overwrite", "false")) %in% c("true", "1", "yes")
# recent=N re-fetches the last N calendar years even when present, so a refresh
# picks up INE's revisions of provisional years; unchanged files stay untouched.
recent_years <- suppressWarnings(as.integer(get_arg("recent", "0")))
if (is.na(recent_years) || recent_years < 0) recent_years <- 0L
recheck_from <- as.integer(format(Sys.Date(), "%Y")) - recent_years

# Same precedence as the main death archive: 0013166 from 2022.
sources <- list(
  "0008206" = list(years = 1991:2021, lookup = "data/nuts_lookup_2013.rds"),
  "0013166" = list(years = 2022:2100, lookup = "data/nuts_lookup_2024.rds")
)
if (!identical(indicator_arg, "all")) sources <- sources[intersect(names(sources), strsplit(indicator_arg, ",")[[1]])]

parse_years <- function(value, default) {
  if (identical(toupper(value), "ALL")) return(default)
  if (grepl(":", value, fixed = TRUE)) { b <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]]); return(seq.int(min(b), max(b))) }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

client <- ineptr2::INEClient$new(lang = "PT", timeout = 900)
fetch <- function(indicator, year) {
  # INE answers bursts with "429 Too Many Requests" and then refuses connections
  # for a while, so back off between attempts and pause between years.
  for (attempt in 1:5) {
    r <- tryCatch(client$get_data(indicator, dim1 = paste0("S7A", year), dim4 = "T"), error = function(e) NULL)
    if (!is.null(r) && nrow(r) > 0) {
      Sys.sleep(5)
      return(r)
    }
    Sys.sleep(60 * attempt)
  }
  NULL
}

find_dim <- function(raw, pattern) {
  cols <- grep("^dim_[0-9]+_t$", names(raw), value = TRUE)
  hit <- cols[vapply(cols, function(col) any(grepl(pattern, raw[[col]])), logical(1))]
  if (length(hit) != 1) stop("Could not identify dimension ", pattern, call. = FALSE)
  hit
}

# Writes go through R/data_versions.R: identical content is left alone, a
# revised file is archived under data/archive before being replaced, and
# every write is recorded in data/import_log.csv with its date.
sys.source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]))), "..", "R", "data_versions.R"), envir = environment())
save_rds_atomic <- function(x, path) versioned_save_rds(x, path, tool = "fetch_death_totals.R", note = Sys.getenv("DATA_RUN_NOTE", unset = NA))

written <- 0L; failed <- character(0)
for (indicator in names(sources)) {
  spec <- sources[[indicator]]
  lookup <- readRDS(spec$lookup) %>% transmute(code = as.character(municipality_code), municipality)
  # Only years INE actually publishes: the ranges above are open-ended so a new
  # year is picked up without editing this file.
  published <- tryCatch({
    dv <- client$get_dim_values(indicator)
    sort(unique(suppressWarnings(as.integer(as.character(dv$categ_dsg[as.integer(dv$dim_num) == 1])))))
  }, error = function(e) spec$years)
  spec$years <- intersect(spec$years, published)
  for (year in intersect(parse_years(years_arg, spec$years), spec$years)) {
    path <- file.path("data/snapshots/death_totals", indicator, paste0("year_", year, ".rds"))
    if (file.exists(path) && !overwrite && year < recheck_from) next

    raw <- fetch(indicator, year)
    if (is.null(raw)) { message("  ", year, ": FAILED"); failed <- c(failed, paste(indicator, year)); next }

    sex_col <- find_dim(raw, "^HM$")
    cause_col <- find_dim(raw, "^Todas as causas de morte$")

    tidy <- raw %>%
      transmute(
        code = as.character(geocod),
        label = as.character(geodsg),
        sex = as.character(.data[[sex_col]]),
        cause = as.character(.data[[cause_col]]),
        deaths = suppressWarnings(as.numeric(valor))
      )

    municipal <- tidy %>%
      inner_join(lookup, by = "code") %>%
      transmute(area = municipality, sex, cause, deaths)
    national <- tidy %>%
      filter(code %in% c("PT", "1")) %>%
      transmute(area = if_else(code == "PT", "Portugal", "Continente"), sex, cause, deaths)

    chunk <- bind_rows(municipal, national) %>%
      mutate(year = as.integer(year), deaths = coalesce(deaths, 0), source_indicator = indicator) %>%
      select(year, area, sex, cause, deaths, source_indicator) %>%
      arrange(area, cause, sex)

    n_mun <- n_distinct(municipal$area)
    ac <- chunk %>% filter(sex == "HM", cause == "Todas as causas de morte")
    save_rds_atomic(chunk, path)
    written <- written + 1L
    message(sprintf("  %s %d: %d municipalities | all-cause HM: municipal sum %s, Portugal %s",
      indicator, year, n_mun,
      format(sum(ac$deaths[!ac$area %in% c("Portugal", "Continente")]), big.mark = " "),
      format(sum(ac$deaths[ac$area == "Portugal"]), big.mark = " ")))
    if (n_mun != 308) message("    ! expected 308 municipalities")
  }
}
message("Done: ", written, " written, ", length(failed), " failed.")
