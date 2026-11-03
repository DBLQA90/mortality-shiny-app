# =========================================================
# Parish mortality: the SMR of each parish (Mapas, "Freguesias")
# =========================================================
# INE publishes the deaths of each parish of residence every year, all causes
# and both sexes together, with no age (0008235 for 2014-2020, 0012542 from
# 2021: data/snapshots/parish_vitals), and the population of each parish by age
# only in the 2021 census (data/snapshots/census_parish). So the one comparable
# measure is the SMR by indirect standardisation: the deaths each parish would
# have if every age group died at the national rate - each age group's national
# deaths in the triennium, shared out by the parish's share of that group in
# the census - against the deaths it had. Portugal is 100; the expected deaths
# add up to the observed ones over the country (the national deaths by age are
# scaled to the parishes' total, the deaths of no known parish left out).
#
# Two assumptions follow from the data and are said on the map: the census
# shares hold for the triennium (hence only triennia within three years of
# 2021), and the parishes are those of 2013 - INE's 2025 figures use the map
# with the restored parishes, which the census does not have. Half the parishes
# expect fewer than 46 deaths in a triennium (2022-2024; 112 fewer than 10), so
# much of the observed SMR is chance and the map smooths it by default
# (empirical Bayes, as the municipal maps can).

PARISH_MAP_LEVEL <- "freguesia"
PARISH_MAP_END_YEARS <- 2020:2024
PARISH_CENSUS_YEAR <- 2021L

# The triennia the parish map can show: within the years above, with national
# deaths by age for all three years.
parish_map_end_years <- function() intersect(PARISH_MAP_END_YEARS, cause_map_end_years(planning_all_causes))

# Why the parish map cannot show a selection, or NULL.
parish_map_unavailable_reason <- function(cause, sex, metric, view) {
  if (!identical(cause, planning_all_causes) || !identical(sex, "HM")) {
    return("Por freguesia, o INE publica só os óbitos de todas as causas, sem separar os sexos: escolha «Todas as causas de morte» e «Ambos».")
  }
  if (!identical(metric, "smr")) {
    return("Por freguesia não há óbitos por idade, por isso não há taxa padronizada: escolha a SMR.")
  }
  if (!identical(view, "period")) {
    return("Por freguesia só há população no Censos 2021: o mapa mostra um triénio próximo, sem tendência, variação ou projecção.")
  }
  NULL
}

# The parishes of the census: code, name, municipality (by DICO, as the NUTS
# lookup names it) and the label shown on the map.
parish_map_units <- function(lookup) {
  census <- planning_parish_census()
  if (is.null(census)) return(NULL)
  units <- unique(census[census$sex == "HM", c("code", "parish", "dico")])
  dico <- substr(lookup$municipality_code, nchar(lookup$municipality_code) - 3L, nchar(lookup$municipality_code))
  units$municipality <- as.character(lookup$municipality)[match(units$dico, dico)]
  units$label <- paste0(units$parish, " (", units$municipality, ")")
  units[order(units$code), , drop = FALSE]
}

# One row per parish, with the columns of cause_map_values().
parish_map_values <- function(end_year, lookup, smoothed = TRUE) {
  units <- parish_map_units(lookup)
  deaths <- planning_parish_vitals("deaths")
  if (is.null(units) || is.null(deaths) || !end_year %in% parish_map_end_years()) return(tibble::tibble())
  window <- seq.int(end_year - 2L, end_year)

  # Observed: the parish deaths of the three years (missing if any year is).
  deaths <- deaths[deaths$year %in% window & deaths$code %in% units$code, , drop = FALSE]
  by_year <- tapply(deaths$value, list(deaths$code, deaths$year), sum)
  observed <- rowSums(by_year[units$code, as.character(window), drop = FALSE])

  # Expected: the national deaths of each census age group, shared by the
  # parish's census share of the group.
  census <- planning_parish_census()
  census <- census[census$sex == "HM" & census$code %in% units$code, , drop = FALSE]
  census$lower <- planning_census_age_lower(census$age_group)
  pop <- tapply(census$pop, list(census$code, census$lower), sum)
  pop[is.na(pop)] <- 0
  pop <- pop[units$code, , drop = FALSE]
  groups <- as.integer(colnames(pop))
  pool <- planning_standardised_pool(sort(unique(as.character(lookup$municipality))), end_year, "HM", lookup, span = 3L)
  if (is.null(pool)) return(tibble::tibble())
  by_band <- colSums(pool$deaths[, , "all"])
  group_of <- groups[findInterval(planning_age_lower(names(by_band)), groups)]
  national <- tapply(by_band, group_of, sum)[as.character(groups)]
  national <- national * sum(observed, na.rm = TRUE) / sum(national)
  expected <- as.vector(sweep(pop, 2, colSums(pop), "/") %*% national)

  ci <- planning_poisson_ci(observed)
  ok <- is.finite(observed) & expected > 0
  smr <- ifelse(ok, observed / expected * 100, NA_real_)
  lower <- ifelse(ok, ci$lower / expected * 100, NA_real_)
  upper <- ifelse(ok, ci$upper / expected * 100, NA_real_)
  value <- smr
  method <- "observado"
  if (isTRUE(smoothed) && sum(ok) >= 3) {
    value[ok] <- planning_eb_smooth(data.frame(value = smr[ok], denominator = expected[ok]),
                                    list(model = "poisson", multiplier = 100))$smoothed
    method <- "suavizado"
  }
  significance <- planning_significance(lower, upper, rep(100, length(lower)))
  out <- tibble::tibble(
    area = units$code, label = units$label, municipality = units$municipality,
    value = value, observed_value = smr, lower = lower, upper = upper,
    observed = unname(observed), expected = expected, flag = "", benchmark_value = 100,
    index = value, significance = significance, method = method
  )
  out$quintile <- map_quintile(out$value)
  out$class <- map_significance_class(out$significance)
  out
}

# The parishes of the areas chosen elsewhere (a municipality, region or ULS),
# for the outline.
parish_map_highlight <- function(focus, lookup) {
  units <- parish_map_units(lookup)
  if (is.null(units) || length(focus) == 0) return(character(0))
  municipalities <- unique(c(intersect(focus, units$municipality),
                             unlist(lapply(setdiff(focus, units$municipality), planning_area_members, lookup = lookup))))
  units$code[units$municipality %in% municipalities]
}
