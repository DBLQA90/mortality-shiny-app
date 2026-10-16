# =========================================================
# ULS that share a municipality: whole municipality or parish weights
# =========================================================
# Four municipalities are divided between ULS at parish level - Lisboa (ULS
# Santa Maria, ULS São José, ULS Lisboa Ocidental), Loures (ULS Loures/Odivelas,
# ULS São José), Porto (ULS Santo António, ULS São João) and Gondomar (ULS Santo
# António, and the parish of Lomba in ULS Vila Nova de Gaia/Espinho) - and INE
# publishes nothing below the municipality except the census and annual counts
# of births and deaths. The tab offers two readings of those ULS, and says which
# one is in use:
#
#   "whole"   the ULS takes each municipality it serves whole. Nothing is
#             estimated, and it answers the question a ULS actually faces (it
#             cannot turn away a resident of the municipality it serves), but
#             the six overlap: summing them counts Lisboa, Loures and Porto
#             twice, 11,613 deaths in 2024.
#   "parish"  each shared municipality is split by the parishes' share of it in
#             the 2021 census, by age group. The parts add up to the
#             municipality and to the country, at the cost of assuming the
#             census shares still hold and that, within a municipality and age
#             group, the parishes of one ULS behave like the others.
#
# The exact groups (unions of ULS that contain only whole municipalities) stay
# available and assume nothing at all.
#
# Lomba is the exception to "whole": 1,284 of Gondomar's 164,257 residents
# (0.8%). Giving all of Gondomar to ULS Vila Nova de Gaia/Espinho for it would
# add half that ULS's population, so a part flagged `whole = FALSE` in
# data-raw/uls_parish.csv is left out of the whole reading (Gondomar counts
# whole in ULS Santo António only, as in every other tab) and enters only
# through the parish weights (planning_parish_extend()). The exact group of
# Santo António and São João therefore includes Lomba.
#
# Weights. Births and deaths are not estimated at all where INE publishes them
# by parish - annually since 2014, counts only (0012450/0012542 and the earlier
# editions). Each year's own parish split is used, and it is exact: the
# parishes of a municipality sum to its total. So:
#   births and everything per birth   the parishes' share of that year's births
#   deaths, infant and perinatal      their share of that year's deaths
#   population and everything else    the census share (the only parish
#                                     population there is)
# Age is used wherever it is known: a population or death band is weighted by
# the parishes' census share of that age group, and for deaths those band
# weights are then scaled so the bands add up to the year's real parish share -
# the census supplies the shape over ages, the registers the level.
# Before 2014 the parishes are the pre-reform ones and do not match, so those
# years fall back to the census shares: women 15-49 for births, and the census
# age structure weighted by national death rates for deaths.

PLANNING_SPLIT_MODES <- c(
  "Município inteiro" = "whole",
  "Ponderação por freguesias (Censos 2021)" = "parish"
)
PLANNING_DEFAULT_SPLIT_MODE <- "whole"

planning_parish_path <- function() file.path(app_dir_or_wd(), "data", "uls_parish.rds")
planning_parish_available <- function() file.exists(planning_parish_path())

planning_parish_lookup <- function() {
  if (!planning_parish_available()) return(NULL)
  key <- "parish|lookup"
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  value <- readRDS(planning_parish_path())
  assign(key, value, envir = planning_cache)
  value
}

# The ULS that share a municipality, and the municipalities each of them serves
# (whole, in the "whole" reading).
planning_parish_units <- function() {
  lookup <- planning_parish_lookup()
  if (is.null(lookup)) return(character(0))
  sort(unique(lookup$unit))
}

# The municipalities a ULS takes whole in the "whole" reading: not those where
# it holds only a part flagged `whole = FALSE`.
planning_parish_members <- function(unit) {
  lookup <- planning_parish_lookup()
  if (is.null(lookup)) return(character(0))
  whole <- if ("whole" %in% names(lookup)) lookup$whole else rep(TRUE, nrow(lookup))
  sort(unique(lookup$municipality[lookup$unit == unit & whole]))
}

# The (unit, municipality) pairs a ULS holds only in part, too small to take
# whole: they enter the parish reading only.
planning_parish_partial_pairs <- function() {
  lookup <- planning_parish_lookup()
  if (is.null(lookup) || !"whole" %in% names(lookup)) return(data.frame(unit = character(0), municipality = character(0)))
  pairs <- unique(lookup[, c("unit", "municipality")])
  whole <- vapply(seq_len(nrow(pairs)), function(i) {
    any(lookup$whole[lookup$unit == pairs$unit[[i]] & lookup$municipality == pairs$municipality[[i]]])
  }, logical(1))
  as.data.frame(pairs[!whole, , drop = FALSE])
}

# The 0/1 membership with those parts added, for the parish reading: the
# weights then give each its share.
planning_parish_extend <- function(membership) {
  pairs <- planning_parish_partial_pairs()
  for (i in seq_len(nrow(pairs))) {
    a <- match(pairs$unit[[i]], rownames(membership)); m <- match(pairs$municipality[[i]], colnames(membership))
    if (!is.na(a) && !is.na(m)) membership[a, m] <- 1
  }
  membership
}

# Census population by parish, age group and sex.
planning_parish_census <- function() {
  path <- file.path(infant_snapshot_root(), "census_parish", "year_2021.rds")
  if (!file.exists(path)) return(NULL)
  key <- paste(infant_snapshot_root(), "parishcensus", sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  value <- readRDS(path)
  assign(key, value, envir = planning_cache)
  value
}

planning_census_age_lower <- function(group) {
  ifelse(grepl("^Menos de 15", group), 0L, suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", group))))
}

# National age-specific death rates by census age group, for the mortality
# weight. Taken from the latest year with deaths by age and population.
planning_mortality_weights <- function() {
  key <- paste(infant_snapshot_root(), "mortweights", sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  census <- planning_parish_census()
  years <- intersect(life_death_years(), snapshot_years_for("population"))
  value <- NULL
  if (length(years) > 0 && !is.null(census)) {
    year <- max(years)
    deaths <- readRDS(life_death_file(year))
    pop <- read_year_file("population", year)
    groups <- sort(unique(planning_census_age_lower(census$age_group)))
    band_group <- function(bands) groups[findInterval(planning_age_lower(bands), groups)]
    d <- stats::aggregate(list(deaths = deaths$deaths[deaths$sex == "HM" & deaths$area == "Portugal"]),
                          by = list(group = band_group(deaths$age_band[deaths$sex == "HM" & deaths$area == "Portugal"])), FUN = sum)
    p <- stats::aggregate(list(pop = pop$pop[pop$sex == "HM" & pop$area == "Portugal"]),
                          by = list(group = band_group(pop$age_band[pop$sex == "HM" & pop$area == "Portugal"])), FUN = sum)
    merged <- merge(d, p, by = "group")
    value <- stats::setNames(merged$deaths / merged$pop, merged$group)
  }
  assign(key, value, envir = planning_cache)
  value
}

# Share of each shared municipality that belongs to each ULS, per basis:
# "total", "female_15_49", "mortality", or the lower bound of a census age
# group as a character ("0", "15", ... "75") for age-specific weights.
planning_parish_shares <- function() {
  key <- paste(infant_snapshot_root(), "parishshares", sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  lookup <- planning_parish_lookup()
  census <- planning_parish_census()
  if (is.null(lookup) || is.null(census)) {
    assign(key, NULL, envir = planning_cache)
    return(NULL)
  }
  rows <- census %>%
    dplyr::inner_join(lookup[, c("code", "municipality", "unit")], by = "code") %>%
    dplyr::mutate(lower = planning_census_age_lower(.data$age_group))

  share_of <- function(frame, basis) {
    if (nrow(frame) == 0) return(NULL)
    frame %>%
      dplyr::group_by(.data$municipality, .data$unit) %>%
      dplyr::summarise(part = sum(.data$pop, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::mutate(weight = .data$part / sum(.data$part), basis = basis) %>%
      dplyr::ungroup() %>%
      dplyr::select(municipality, unit, basis, weight)
  }

  both <- rows[rows$sex == "HM", , drop = FALSE]
  mortality <- planning_mortality_weights()
  weighted <- both
  weighted$pop <- weighted$pop * unname(mortality[as.character(weighted$lower)])
  ages <- sort(unique(both$lower))

  value <- dplyr::bind_rows(
    share_of(both, "total"),
    share_of(both[both$lower >= 15, , drop = FALSE], "15_plus"),
    share_of(both[both$lower >= 15 & both$lower < 65, , drop = FALSE], "15_64"),
    share_of(both[both$lower >= 65, , drop = FALSE], "65_plus"),
    share_of(rows[rows$sex == "M" & rows$lower >= 15 & rows$lower < 50, , drop = FALSE], "female_15_49"),
    # Women by age group, for the fertility index's denominators.
    dplyr::bind_rows(lapply(seq(15, 45, by = 5), function(a) {
      share_of(rows[rows$sex == "M" & rows$lower == a, , drop = FALSE], paste0("f", a))
    })),
    share_of(weighted, "mortality"),
    # One basis per census age group, for values that carry an age band.
    dplyr::bind_rows(lapply(ages, function(a) share_of(both[both$lower == a, , drop = FALSE], as.character(a))))
  )
  assign(key, value, envir = planning_cache)
  value
}

# Births and deaths by parish, as published (counts, both sexes).
planning_parish_vitals <- function(kind) {
  path <- file.path(infant_snapshot_root(), "parish_vitals", paste0(kind, ".rds"))
  if (!file.exists(path)) return(NULL)
  key <- paste(infant_snapshot_root(), "parishvitals", kind, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  value <- readRDS(path)
  assign(key, value, envir = planning_cache)
  value
}

# Each ULS's share of a shared municipality's births or deaths, per year, from
# the parish registers. NULL where the parishes of the year do not match the
# current ones (before the 2013 reform).
planning_parish_actual <- function(kind) {
  key <- paste(infant_snapshot_root(), "parishactual", kind, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
  lookup <- planning_parish_lookup()
  vitals <- planning_parish_vitals(kind)
  value <- NULL
  if (!is.null(lookup) && !is.null(vitals)) {
    rows <- vitals[vitals$dico %in% unique(lookup$dico), , drop = FALSE]
    totals <- rows %>%
      dplyr::group_by(.data$year, .data$dico) %>%
      dplyr::summarise(total = sum(.data$value, na.rm = TRUE), covered = sum(.data$value[.data$code %in% lookup$code], na.rm = TRUE), .groups = "drop")
    # A year counts only if its parishes are the current ones, so that the
    # parts of the municipality are complete.
    usable <- totals[totals$total > 0 & abs(totals$covered / totals$total - 1) < 1e-9, c("year", "dico", "total")]
    value <- rows %>%
      dplyr::inner_join(lookup[, c("code", "municipality", "unit")], by = "code") %>%
      dplyr::inner_join(usable, by = c("year", "dico")) %>%
      dplyr::group_by(.data$year, .data$municipality, .data$unit) %>%
      dplyr::summarise(weight = sum(.data$value, na.rm = TRUE) / dplyr::first(.data$total), .groups = "drop")
  }
  assign(key, value, envir = planning_cache)
  value
}

# The years the registers cover for both births and deaths.
planning_parish_actual_years <- function() {
  years <- lapply(c("births", "deaths"), function(kind) {
    actual <- planning_parish_actual(kind)
    if (is.null(actual)) integer(0) else sort(unique(actual$year))
  })
  Reduce(intersect, years)
}

# The weight basis for each component column of the planning tab.
planning_weight_basis <- function(column) {
  if (identical(column, "pop_total")) return("total")
  if (identical(column, "pop_15_plus")) return("15_plus")
  if (identical(column, "pop_0_14")) return("0")
  if (identical(column, "pop_15_64")) return("15_64")
  if (identical(column, "pop_65_plus")) return("65_plus")
  if (identical(column, "pop_75_plus")) return("75")
  if (startsWith(column, "pop_f_")) return(paste0("f", sub("^pop_f_", "", column)))
  if (startsWith(column, "births") || column %in% c("infant_deaths", "neonatal_deaths", "early_neonatal_deaths",
                                                    "postneonatal_deaths", "perinatal_deaths")) return("female_15_49")
  if (identical(column, "deaths")) return("mortality")
  "total"
}

# Weights of `areas` over municipalities for one basis, in the chosen mode.
# Every area that is not one of the six split ULS keeps whole municipalities.
planning_parish_weights <- function(areas, municipalities, basis = "total", mode = PLANNING_DEFAULT_SPLIT_MODE, year = NULL) {
  weights <- matrix(1, length(areas), length(municipalities), dimnames = list(areas, municipalities))
  if (!identical(mode, "parish")) return(weights)
  shares <- planning_parish_shares()
  if (is.null(shares)) return(weights)
  # Births and deaths of a year INE publishes by parish are not estimated.
  registered <- if (is.null(year)) NULL else switch(basis, female_15_49 = "births", mortality = "deaths", NULL)
  rows <- NULL
  if (!is.null(registered)) {
    actual <- planning_parish_actual(registered)
    if (!is.null(actual)) rows <- actual[actual$year == as.integer(year), c("municipality", "unit", "weight"), drop = FALSE]
  }
  if (is.null(rows) || nrow(rows) == 0) rows <- shares[shares$basis == basis, , drop = FALSE]
  if (nrow(rows) == 0) rows <- shares[shares$basis == "total", , drop = FALSE]
  for (i in seq_len(nrow(rows))) {
    a <- match(rows$unit[[i]], areas); m <- match(rows$municipality[[i]], municipalities)
    if (!is.na(a) && !is.na(m)) weights[a, m] <- rows$weight[[i]]
  }
  weights
}

# The weight of one age band, for the modules that carry deaths and population
# by band (life expectancy, standardised mortality).
planning_parish_band_weights <- function(areas, municipalities, band, mode = PLANNING_DEFAULT_SPLIT_MODE) {
  shares <- planning_parish_shares()
  if (!identical(mode, "parish") || is.null(shares)) {
    return(matrix(1, length(areas), length(municipalities), dimnames = list(areas, municipalities)))
  }
  groups <- sort(unique(suppressWarnings(as.integer(shares$basis[grepl("^[0-9]+$", shares$basis)]))))
  lower <- planning_age_lower(band)
  basis <- as.character(groups[findInterval(lower, groups)])
  planning_parish_weights(areas, municipalities, basis, mode)
}

# Sum a municipalities x bands matrix into areas x bands, weighting each band
# by the parishes' share of that age group. Outside the parish reading it is
# the plain membership product.
planning_band_product <- function(membership, values, mode = PLANNING_DEFAULT_SPLIT_MODE, bands = colnames(values), target = NULL) {
  municipalities <- colnames(membership)
  if (!identical(mode, "parish") || !any(rownames(membership) %in% planning_parish_units())) {
    return(membership %*% values[municipalities, , drop = FALSE])
  }
  weights <- lapply(bands, function(band) planning_parish_band_weights(rownames(membership), municipalities, band, mode))
  if (!is.null(target)) {
    # The census gives the shape over ages; the registers give the level. Scale
    # each municipality's band weights so they add up to the real share of its
    # births or deaths that year.
    totals <- rowSums(values[municipalities, bands, drop = FALSE])
    implied <- Reduce(`+`, lapply(seq_along(bands), function(j) weights[[j]] * matrix(values[municipalities, bands[[j]]], nrow(membership), length(municipalities), byrow = TRUE)))
    scale <- ifelse(implied > 0, target * matrix(totals, nrow(membership), length(municipalities), byrow = TRUE) / implied, 1)
    weights <- lapply(weights, function(w) w * scale)
  }
  membership <- planning_parish_extend(membership)
  out <- matrix(0, nrow(membership), length(bands), dimnames = list(rownames(membership), bands))
  for (j in seq_along(bands)) {
    out[, j] <- (membership * weights[[j]]) %*% values[municipalities, bands[[j]], drop = FALSE]
  }
  out
}

# Sum a vector over municipalities with the weights of one basis.
planning_weighted_sum <- function(membership, values, basis, mode = PLANNING_DEFAULT_SPLIT_MODE, year = NULL) {
  municipalities <- colnames(membership)
  weights <- planning_parish_weights(rownames(membership), municipalities, basis, mode, year)
  if (identical(mode, "parish")) membership <- planning_parish_extend(membership)
  as.numeric((membership * weights) %*% values[municipalities])
}

# A note for each selected ULS that shares a municipality: what it shares,
# with which ULS, and how the reading in use treats it. `mode` "whole" is also
# how every tab but the planning tab's parish weights reads a ULS.
uls_split_note <- function(units, mode = PLANNING_DEFAULT_SPLIT_MODE) {
  lookup <- planning_parish_lookup()
  if (is.null(lookup)) return(character(0))
  units <- intersect(as.character(units), unique(lookup$unit))
  if (length(units) == 0) return(character(0))
  whole <- if ("whole" %in% names(lookup)) lookup$whole else rep(TRUE, nrow(lookup))
  census <- planning_parish_census()
  resident <- function(codes) {
    if (is.null(census)) return(NA_real_)
    sum(census$pop[census$code %in% codes & census$sex == "HM"], na.rm = TRUE)
  }
  fmt <- function(x) formatC(x, format = "f", digits = 0, big.mark = ".", decimal.mark = ",")
  pct <- function(x) formatC(100 * x, format = "f", digits = 1, decimal.mark = ",")
  parish_reading <- identical(mode, "parish")
  unlist(lapply(units, function(unit) {
    municipalities <- sort(unique(lookup$municipality[lookup$unit == unit]))
    vapply(municipalities, function(m) {
      here <- lookup$municipality == m
      others <- sort(setdiff(unique(lookup$unit[here]), unit))
      own <- here & lookup$unit == unit
      share <- resident(lookup$code[own]) / resident(lookup$code[here])
      partial_here <- !any(whole[own])
      partial_other <- unique(lookup$unit[here & !whole])
      if (partial_here) {
        main <- setdiff(unique(lookup$unit[here & whole]), unit)
        paste0(unit, ": de ", m, " pertence-lhe só a freguesia ", paste(lookup$parish[own], collapse = ", "), " (",
               pct(share), "% dos residentes nos Censos de 2021). ",
               if (parish_reading) {
                 "Com a ponderação por freguesias, essa parte entra nesta ULS. "
               } else {
                 paste0("Com «Município inteiro», como nos outros separadores, ", m, " conta todo na ", paste(main, collapse = " e "),
                        " e não nesta: dar-lhe o município inteiro por uma parte tão pequena somaria ", fmt(resident(lookup$code[here])),
                        " residentes a esta ULS. Só a ponderação por freguesias a inclui. ")
               })
      } else if (length(partial_other) > 0) {
        paste0(unit, ": ", m, " é seu, excepto a freguesia ", paste(lookup$parish[here & lookup$unit %in% partial_other], collapse = ", "),
               " (", pct(1 - share), "%), que pertence à ", paste(partial_other, collapse = " e "), ". ",
               if (parish_reading) {
                 "Com a ponderação por freguesias, essa freguesia fica de fora. "
               } else {
                 paste0("Com «Município inteiro», como nos outros separadores, ", m, " conta todo nesta ULS, incluindo essa freguesia. ")
               })
      } else {
        paste0(unit, ": ", m, " está repartido por freguesias com ", if (length(others) > 1) "as " else "a ",
               paste(others, collapse = " e "), " (esta ULS tem ",
               pct(share), "% dos residentes). ",
               if (parish_reading) {
                 "Com a ponderação por freguesias, cada ULS leva a sua parte, pelas freguesias e pela idade. "
               } else {
                 paste0("Com «Município inteiro», como nos outros separadores, ", m, " conta inteiro em cada uma: nada é estimado, ",
                        "mas as ULS sobrepõem-se e não devem ser somadas. ")
               })
      }
    }, character(1), USE.NAMES = FALSE)
  }))
}
