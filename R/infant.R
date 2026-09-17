# =========================================================
# Infant mortality
# =========================================================
# The infant mortality rate is deaths under one year of age per 1,000 live
# births. Neither part comes from the app's main pipeline:
#
#   Numerator   Both death indicators publish `Menos de 1 ano` as its own age
#               band, but prepare_death_data() folds it into `0 - 4 anos` on
#               ingest, so the main death archive cannot separate an infant
#               death from a death at age four. tools/fetch_infant_deaths.R
#               writes a parallel dataset holding only the under-1 band.
#
#   Denominator Live births, not population. Neither population indicator has
#               an under-1 age band, so infant deaths per under-1 population is
#               not computable at all; the conventional denominator is live
#               births, which INE publishes separately.
#               tools/fetch_births.R writes those.
#
# Keeping both outside the main pipeline means no existing rate changes, and the
# `0 - 4 anos` band keeps behaving exactly as before everywhere else.
#
# Coverage of the rate is bounded by births, currently 1995-2025 with no gap,
# assembled from three INE vintages (see tools/fetch_births.R). The count of
# infant deaths needs only the numerator and reaches 1980-2025. Years outside
# those report as unavailable rather than being filled in, and both spans are
# read from the files present rather than hard-coded, so they follow whatever
# has been fetched.
#
# 2025 is a thinner year than the rest: the by-cause death indicators stop at
# 2024, so it comes from 0012540, which publishes under-1 deaths by
# municipality and month only. That year therefore holds a single all-cause,
# both-sexes total per area, and a cause- or sex-specific request has to be
# refused rather than answered with zero - see infant_year_has_detail().
#
# The under-1 death counts are also used outside this module, to correct the
# years-of-life-lost weight of the `0 - 4 anos` band - see
# split_infant_age_band() in R/metrics.R.

# Resolve the snapshot root without depending on R/snapshots.R.
#
# The infant datasets are read by a calculation module, which the test runner
# loads on its own precisely so the calculations can be tested without the
# INE/snapshot machinery. Calling get_snapshot_dir() directly would drag that
# machinery in, so honour the same environment variable and fall back to it only
# when it happens to be loaded.
infant_snapshot_root <- function() {
  configured <- Sys.getenv("MORTALITY_SNAPSHOT_DIR", unset = "")
  if (nzchar(configured)) {
    return(configured)
  }

  if (exists("get_snapshot_dir", mode = "function", inherits = TRUE)) {
    return(get_snapshot_dir())
  }

  file.path(getwd(), "data", "snapshots")
}

infant_snapshot_dir <- function(dataset) {
  file.path(infant_snapshot_root(), dataset)
}

# Years each dataset actually holds, read from the files rather than assumed, so
# coverage follows whatever has been fetched.
snapshot_years_for <- function(dataset) {
  dir_path <- tryCatch(infant_snapshot_dir(dataset), error = function(e) NULL)

  if (is.null(dir_path) || !dir.exists(dir_path)) {
    return(integer(0))
  }

  files <- list.files(dir_path, pattern = "^year_\\d+\\.rds$")
  sort(as.integer(sub("^year_(\\d+)\\.rds$", "\\1", files)))
}

births_years <- function() snapshot_years_for("births")
infant_death_years <- function() snapshot_years_for("infant_deaths")

# Years where both halves of the ratio exist. The count of infant deaths needs
# only the numerator, so it is available over a wider span - see
# infant_metric_years().
infant_mortality_years <- function() {
  sort(intersect(births_years(), infant_death_years()))
}

# Coverage of whichever infant metric is being asked for: the rate needs both
# datasets, the count only the deaths.
infant_metric_years <- function(metric_id = infant_metric_id) {
  if (identical(metric_id, infant_count_metric_id)) {
    infant_death_years()
  } else {
    infant_mortality_years()
  }
}

read_year_file <- function(dataset, year) {
  path <- file.path(infant_snapshot_dir(dataset), paste0("year_", year, ".rds"))
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

get_births_data <- function(years, areas) {
  rows <- lapply(years, function(year) {
    chunk <- read_year_file("births", year)
    if (is.null(chunk)) return(NULL)
    chunk %>% dplyr::filter(area %in% areas)
  })

  dplyr::bind_rows(rows)
}

# All-cause requests read data/snapshots/infant_totals where a year has it: INE's
# own under-1 death counts (0008181/0012541), complete at municipal level. The
# by-cause dataset takes the "Menos de 1 ano" band of the cause-of-death
# indicators, whose municipal age breakdown is incomplete - in 2014 its
# municipalities add up to 112 infant deaths against a national 236 - so it is
# used only for a specific cause, or where no complete count exists.
#
# `complete = FALSE` keeps the band-derived count. The AVPP correction needs it:
# it splits the `0 - 4 anos` band of the same age breakdown, and a complete count
# larger than that band would be capped rather than split.
infant_totals_years <- function() snapshot_years_for("infant_totals")

get_infant_death_data <- function(years, areas, cause = "Todas as causas de morte", sex = "HM", complete = TRUE) {
  use_totals <- isTRUE(complete) && identical(cause, "Todas as causas de morte")

  rows <- lapply(years, function(year) {
    if (use_totals) {
      totals <- read_year_file("infant_totals", year)
      if (!is.null(totals)) {
        return(
          totals %>%
            dplyr::filter(area %in% areas, sex == .env$sex) %>%
            dplyr::mutate(cause = .env$cause) %>%
            dplyr::select(year, area, sex, cause, deaths, source_indicator)
        )
      }
    }
    chunk <- read_year_file("infant_deaths", year)
    if (is.null(chunk)) return(NULL)
    chunk %>% dplyr::filter(area %in% areas, cause == .env$cause, sex == .env$sex)
  })

  dplyr::bind_rows(rows)
}

# Share of the national under-1 count that the municipalities account for, in a
# year read from the band-derived dataset. INE's complete municipal counts start
# in 2011; before that the "Menos de 1 ano" band of the cause-of-death
# indicators is all there is, and in 1995-2001 its municipalities add up to about
# 85% of Portugal. A year with complete counts returns 1.
infant_municipal_coverage <- function(year,
                                      cause = "Todas as causas de morte",
                                      sex = "HM",
                                      municipalities = get_nuts_lookup()$municipality) {
  if (identical(cause, "Todas as causas de morte") && !is.null(read_year_file("infant_totals", year))) {
    return(1)
  }
  chunk <- read_year_file("infant_deaths", year)
  if (is.null(chunk)) return(NA_real_)
  rows <- chunk[chunk$cause == cause & chunk$sex == sex, , drop = FALSE]
  national <- sum(rows$deaths[rows$area == "Portugal"])
  if (!is.finite(national) || national <= 0) return(NA_real_)
  sum(rows$deaths[rows$area %in% municipalities]) / national
}

infant_undercount_threshold <- 0.97

infant_undercount_years <- function(years, cause = "Todas as causas de morte", sex = "HM",
                                    municipalities = get_nuts_lookup()$municipality) {
  years <- intersect(as.integer(years), infant_death_years())
  coverage <- vapply(years, infant_municipal_coverage, numeric(1),
                     cause = cause, sex = sex, municipalities = municipalities)
  sort(years[!is.na(coverage) & coverage < infant_undercount_threshold])
}

# Warning for a selection that reads municipal under-1 counts in such a year.
# Portugal is read from its own published row and is unaffected.
infant_undercount_message <- function(years, areas, cause = "Todas as causas de morte", sex = "HM",
                                      municipalities = get_nuts_lookup()$municipality) {
  if (length(setdiff(as.character(areas), "Portugal")) == 0) return(NULL)
  short <- infant_undercount_years(years, cause, sex, municipalities)
  if (length(short) == 0) return(NULL)
  as.character(glue::glue(
    "Em {paste(short, collapse = ', ')}, os óbitos com menos de 1 ano por ",
    "município estão incompletos no INE (a soma dos municípios fica abaixo de ",
    "{round(infant_undercount_threshold * 100)}% do total nacional). As taxas e ",
    "contagens de regiões, ULS e municípios nesses anos estão subestimadas; ",
    "Portugal não é afectado. Contagens completas por município existem desde 2011."
  ))
}

# Total under-1 deaths for a selection, or NA when the window is not fully
# covered by the archive.
#
# Partial coverage has to report NA rather than a partial sum. A caller that
# subtracts this from a wider band - as the AVPP correction does - would
# otherwise silently attribute the uncovered years' infant deaths to ages one to
# four, which is a worse error than declining to make the correction at all.
infant_deaths_total <- function(years,
                                areas,
                                cause = "Todas as causas de morte",
                                sex = "HM",
                                complete = TRUE) {
  years <- as.integer(years)

  if (length(years) == 0 || !all(years %in% infant_death_years())) {
    return(NA_real_)
  }

  rows <- get_infant_death_data(years, areas, cause = cause, sex = sex, complete = complete)

  if (nrow(rows) == 0) {
    return(NA_real_)
  }

  sum(rows$deaths, na.rm = TRUE)
}

get_births_total <- function(years, areas) {
  rows <- get_births_data(years, areas)

  if (nrow(rows) == 0) {
    return(NA_real_)
  }

  sum(rows$births, na.rm = TRUE)
}

# Infant mortality rate for one area/period, with an exact Poisson interval on
# the death count scaled by births. Births are treated as a fixed denominator,
# which is the usual convention: the sampling variation that matters is in the
# small number of deaths, not in the birth count.
compute_infant_mortality <- function(deaths, births, confidence = 0.95, multiplier = 1000) {
  na_result <- function(note) {
    tibble::tibble(
      deaths = deaths,
      births = births,
      value = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      note = note
    )
  }

  if (!is.finite(deaths) || !is.finite(births)) {
    return(na_result("Sem dados de óbitos infantis ou de nados-vivos."))
  }
  if (births <= 0) {
    return(na_result("Sem nados-vivos registados no período, pelo que a taxa não é calculável."))
  }

  ci <- compute_count_interval(deaths, confidence = confidence)

  tibble::tibble(
    deaths = deaths,
    births = births,
    value = deaths / births * multiplier,
    lower = unname(ci[[1]]) / births * multiplier,
    upper = unname(ci[[2]]) / births * multiplier,
    note = NA_character_
  )
}

# Below this many live births the rate should not be read as a rate.
#
# The threshold is not a significance test but a statement about resolution: with
# fewer than 1,000 births in the period, one extra death moves the rate by more
# than one whole unit of the scale it is reported in - larger than the entire
# national rate of about 3 per 1,000. Ranking such places, or reading a change
# between years, is reading noise. Most Portuguese municipalities are below it,
# which is the point: the flag marks the ordinary case for what it is rather
# than singling out a few outliers.
#
# Nothing is suppressed. The value is still shown, still exact, and the interval
# already states the uncertainty; the mark only stops a reader skimming the
# table from treating the number as comparable.
infant_stable_births_min <- 1000L

infant_rate_is_unstable <- function(births, threshold = infant_stable_births_min) {
  isTRUE(is.finite(births) && births < threshold)
}

infant_instability_mark <- "*"

infant_instability_footnote <- function(threshold = infant_stable_births_min) {
  as.character(glue::glue(
    "* Menos de {format(threshold, big.mark = ' ')} nados-vivos no período: um ",
    "único óbito altera a taxa em mais de 1 por 1.000, pelo que o valor não é ",
    "comparável entre locais nem entre anos. Prefira a métrica «Óbitos ",
    "infantis (< 1 ano)», ou agregue vários anos."
  ))
}

# Whether a year holds the requested cause/sex breakdown.
#
# The by-cause indicators end at 2024. 2025 comes from 0012540, which publishes
# under-1 deaths by municipality and month only, so that year holds a single
# all-cause, both-sexes total per area. Asking it for a specific cause or sex
# returns no rows, which without this check reads as "no deaths" rather than
# "not published at that detail".
infant_year_has_detail <- function(year,
                                   cause = "Todas as causas de morte",
                                   sex = "HM") {
  chunk <- read_year_file("infant_deaths", year)

  if (is.null(chunk) || nrow(chunk) == 0) {
    return(FALSE)
  }

  any(chunk$cause == cause & chunk$sex == sex)
}

# Years in the selection that can answer at the requested detail.
infant_detail_years <- function(years,
                                cause = "Todas as causas de morte",
                                sex = "HM") {
  years <- intersect(as.integer(years), infant_death_years())

  if (length(years) == 0) {
    return(integer(0))
  }

  years[vapply(years, infant_year_has_detail, logical(1), cause = cause, sex = sex)]
}

# Years in the selection that cannot answer at the requested detail.
infant_detail_message <- function(years,
                                  cause = "Todas as causas de morte",
                                  sex = "HM") {
  years <- intersect(as.integer(years), infant_death_years())

  if (length(years) == 0) {
    return(NULL)
  }

  lacking <- years[!vapply(years, infant_year_has_detail, logical(1), cause = cause, sex = sex)]

  if (length(lacking) == 0) {
    return(NULL)
  }

  as.character(glue::glue(
    "Para {paste(sort(lacking), collapse = ', ')} o INE publica óbitos com ",
    "menos de 1 ano apenas no total (todas as causas, ambos os sexos). ",
    "Seleccione «Todas as causas de morte» e o sexo «HM» para incluir ",
    "esse(s) ano(s)."
  ))
}

# Message for a selection whose years fall outside the coverage of the dataset
# the chosen metric needs, naming the gap and its reason so it reads as a known
# limit of the source rather than a fault.
infant_gap_message <- function(years, metric_id = infant_metric_id) {
  available <- infant_metric_years(metric_id)
  missing <- setdiff(as.integer(years), available)

  if (length(missing) == 0) {
    return(NULL)
  }

  missing_text <- paste(sort(missing), collapse = ", ")

  if (identical(metric_id, infant_count_metric_id)) {
    return(as.character(glue::glue(
      "O arquivo de óbitos com menos de 1 ano não cobre {missing_text}, pelo ",
      "que a contagem não pode ser apresentada nesse(s) ano(s)."
    )))
  }

  as.character(glue::glue(
    "O INE não publica nados-vivos por município para {missing_text}, pelo que ",
    "a mortalidade infantil não pode ser calculada nesse(s) ano(s). A contagem ",
    "de óbitos infantis e as restantes métricas mantêm-se disponíveis."
  ))
}
