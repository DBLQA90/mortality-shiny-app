# =========================================================
# Proportional mortality under 75 (planning indicator I46)
# =========================================================
# The same 13 cause groups as I45, restricted to deaths before 75. Unlike I45 it
# cannot read the all-ages totals, so it needs INE's breakdown of deaths by
# cause and age, which is complete for the rows INE publishes per region and
# very nearly complete when summed from municipalities.
#
# Each area therefore takes the better of the two sources, the same rule the
# app's regional deaths follow:
#
#   INE's own row      Portugal, Continente, every NUTS region, and the ULS that
#                      coincide with a NUTS III unit - exact
#   municipal sum      everything else
#
# How good the municipal sum is, measured against INE's own rows: for 2020-2022
# it reproduces the under-75 deaths of Alto Minho and Algarve exactly, and the
# shares to within 0.08 and 0.29 percentage points; for Portugal the shares are
# within 0.21 points. The exception is 2014, where INE published only 79.7% of
# municipal deaths with an age (52.7% for the worst cause group). Areas not on
# INE rows are therefore marked in the three triennia that contain 2014, where
# Alto Minho's shares are off by up to 2.2 points.
#
# The deaths are not rescaled to the complete totals. Tested against INE's rows,
# rescaling makes the shares worse, not better: deaths with no published age sit
# mostly at older ages, so spreading them proportionally moves too many of them
# below 75 (Portugal 2020-2022, malignant tumours: 38.69% exact, 38.48% from
# municipal sums, 39.20% rescaled).

UNDER75_INCOMPLETE_YEARS <- 2014L
UNDER75_FLAG <- "§"

# The death archive names each file after its cause; the same rule as
# snapshot_file_token(), kept here so this module needs no snapshot machinery.
planning_cause_file_token <- function(cause) {
  token <- tolower(iconv(as.character(cause), to = "ASCII//TRANSLIT"))
  token <- gsub("[^a-z0-9]+", "_", token)
  gsub("^_+|_+$", "", token)
}

planning_under75_bands <- function() {
  age_levels[seq_len(which(age_levels == "75 - 79 anos") - 1L)]
}

# Deaths before 75 for every cause group, per area label and per INE region
# code, in one year. Read once per year and cached.
planning_under75_year <- function(year, sex = "HM") {
  key <- paste(infant_snapshot_root(), "under75", year, sex, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) {
    return(get(key, envir = planning_cache, inherits = FALSE))
  }

  causes <- c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)
  bands <- planning_under75_bands()
  indicator <- NULL
  municipal <- list()
  regional <- list()

  for (cause in causes) {
    path <- NULL
    for (candidate in c("0013166", "0008206")) {
      file <- file.path(infant_snapshot_root(), "deaths", candidate, paste0("year_", year),
                        paste0("cause_", planning_cause_file_token(cause), ".rds"))
      if (file.exists(file)) {
        path <- file
        indicator <- candidate
        break
      }
    }
    if (is.null(path)) next

    rows <- readRDS(path)
    rows <- rows[rows$sex == sex & rows$age_band %in% bands, , drop = FALSE]
    municipal[[cause]] <- tapply(rows$deaths, rows$area, sum)

    regional_rows <- read_regional_rows(indicator, year)
    if (!is.null(regional_rows)) {
      hit <- regional_rows[regional_rows$sex == sex & regional_rows$cause == cause &
                             regional_rows$age_band %in% bands, , drop = FALSE]
      if (nrow(hit) > 0) regional[[cause]] <- tapply(hit$deaths, hit$region_code, sum)
    }
  }

  value <- list(municipal = municipal, regional = regional, indicator = indicator, causes = causes)
  assign(key, value, envir = planning_cache)
  value
}

# An area INE publishes a row for, whose municipalities are exactly this area's.
# Several territories carry two names - NUTS III Alto Minho is also ULS Alto
# Minho - and the table of rows lists only one of them.
planning_row_alias <- function(area, members, vintage, lookup) {
  key <- paste(vintage, "alias", area, sep = "|")
  if (exists(key, envir = planning_cache, inherits = FALSE)) {
    return(get(key, envir = planning_cache, inherits = FALSE))
  }
  candidates <- unique(REGIONAL_ROW_TERRITORIES$region)
  alias <- NA_character_
  # Asked for as a municipal sum on purpose: never swap it for INE's row.
  if (identical(area, PLANNING_PORTUGAL_MUNICIPAL)) candidates <- character(0)
  for (candidate in candidates) {
    if (identical(candidate, area)) next
    if (setequal(members, planning_area_members(candidate, lookup))) {
      alias <- candidate
      break
    }
  }
  assign(key, alias, envir = planning_cache)
  alias
}

# Deaths before 75 of one area in one year, per cause: INE's row where the area
# has one, otherwise the sum of its municipalities. Returns the counts and
# whether they came from INE's rows.
planning_under75_area_year <- function(area, year, members, vintage, lookup) {
  block <- planning_under75_year(year)
  causes <- block$causes
  from_row <- FALSE

  counts <- stats::setNames(rep(0, length(causes)), causes)
  plan <- if (area %in% planning_published_areas) {
    tibble::tibble()
  } else {
    own <- regional_row_plan(area, vintage, year)
    if (nrow(own) > 0) {
      own
    } else {
      alias <- planning_row_alias(area, members, vintage, lookup)
      if (is.na(alias)) own else regional_row_plan(alias, vintage, year)
    }
  }

  # Only the causes the year's archive actually holds; a cause with no file is
  # zero on either path. Every one of them must be in INE's rows, or the area
  # would mix sources between causes and its shares would not add up.
  present <- intersect(causes, names(block$municipal))

  if (nrow(plan) == 1 && identical(plan$indicator[[1]], block$indicator) && length(present) > 0) {
    codes <- split_list(plan$codes[[1]])
    plus <- split_list(plan$plus[[1]])
    minus <- split_list(plan$minus[[1]])
    available <- vapply(present, function(cause) {
      values <- block$regional[[cause]]
      !is.null(values) && all(codes %in% names(values))
    }, logical(1))
    if (all(available)) {
      for (cause in present) {
        values <- block$regional[[cause]]
        municipal <- block$municipal[[cause]]
        total <- sum(values[codes], na.rm = TRUE)
        if (length(plus) > 0) total <- total + sum(municipal[intersect(plus, names(municipal))], na.rm = TRUE)
        if (length(minus) > 0) total <- total - sum(municipal[intersect(minus, names(municipal))], na.rm = TRUE)
        counts[[cause]] <- total
      }
      return(list(counts = counts, from_row = TRUE))
    }
  }

  # Portugal and Continente publish their own rows in the death files.
  published <- area %in% planning_published_areas
  for (cause in causes) {
    values <- block$municipal[[cause]]
    if (is.null(values)) next
    if (published && area %in% names(values)) {
      counts[[cause]] <- values[[area]]
      from_row <- TRUE
    } else {
      counts[[cause]] <- sum(values[intersect(members, names(values))], na.rm = TRUE)
    }
  }
  list(counts = counts, from_row = from_row)
}

# Proportional mortality under 75 for each area and triennium, in the shape
# planning_proportional_table() returns, plus a flag for areas whose municipal
# sums cover a year INE published incompletely.
planning_under75_table <- function(areas, end_years, window = 3L, sex = "HM",
                                   lookup = get_nuts_lookup(), vintage = default_nuts_vintage()) {
  areas <- unique(as.character(areas))
  members <- stats::setNames(lapply(areas, planning_area_members, lookup = lookup), areas)
  codes <- c("C00", PLANNING_CAUSE_GROUPS$code, "Outras")
  labels <- c("Todas as causas de morte", PLANNING_CAUSE_GROUPS$label, "Restantes causas")

  # Municipal sums for every area at once, one matrix product per year; the
  # areas INE publishes a row for are then overwritten one by one (about thirty
  # of them).
  membership <- planning_membership_matrix(areas, lookup)
  municipalities <- colnames(membership)
  causes <- c(planning_all_causes, PLANNING_CAUSE_GROUPS$cause)

  municipal_year <- function(year) {
    key <- paste(infant_snapshot_root(), "under75sum", year, paste(areas, collapse = "|"), sep = "|")
    if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
    block <- planning_under75_year(year)
    values <- matrix(0, length(municipalities), length(causes), dimnames = list(municipalities, causes))
    for (cause in intersect(causes, names(block$municipal))) {
      counts <- block$municipal[[cause]]
      hit <- intersect(municipalities, names(counts))
      values[hit, cause] <- counts[hit]
    }
    summed <- membership %*% values
    colnames(summed) <- causes
    # Portugal and Continente publish their own rows in the death files - but
    # not in every year: Continente has none before 2022, and then falls back to
    # its municipalities like any other area.
    own_row <- stats::setNames(rep(FALSE, length(areas)), areas)
    for (area in intersect(areas, planning_published_areas)) {
      present <- intersect(causes, names(block$municipal))
      has_own <- length(present) > 0 && all(vapply(present, function(cause) area %in% names(block$municipal[[cause]]), logical(1)))
      if (!has_own) next
      for (cause in present) summed[area, cause] <- block$municipal[[cause]][[area]]
      own_row[[area]] <- TRUE
    }
    attr(summed, "own_row") <- own_row
    assign(key, summed, envir = planning_cache)
    summed
  }

  # Only areas that can have an INE row are looked up, and once per year rather
  # than once per triennium: about two dozen areas instead of every municipality.
  candidates <- areas[!areas %in% planning_published_areas]
  candidates <- candidates[vapply(candidates, function(area) {
    nrow(REGIONAL_ROW_TERRITORIES[REGIONAL_ROW_TERRITORIES$region == area, , drop = FALSE]) > 0 ||
      !is.na(planning_row_alias(area, members[[area]], vintage, lookup))
  }, logical(1))]

  row_values_year <- function(year) {
    key <- paste(infant_snapshot_root(), "under75rows", year, vintage, paste(candidates, collapse = "|"), sep = "|")
    if (exists(key, envir = planning_cache, inherits = FALSE)) return(get(key, envir = planning_cache, inherits = FALSE))
    value <- stats::setNames(lapply(candidates, function(area) {
      got <- planning_under75_area_year(area, year, members[[area]], vintage, lookup)
      if (got$from_row) got$counts else NULL
    }), candidates)
    assign(key, value, envir = planning_cache)
    value
  }

  rows <- lapply(as.integer(end_years), function(end_year) {
    years <- seq.int(end_year - window + 1L, end_year)
    if (!all(years %in% death_totals_years())) return(NULL)
    summed <- lapply(years, municipal_year)
    from_rows <- lapply(years, row_values_year)
    # Deaths up to 1998 sit with the parent of Odivelas, Trofa and Vizela.
    joint <- stats::setNames(Reduce(`|`, lapply(years, function(y) planning_joint_split(membership, y, "deaths"))), rownames(membership))

    dplyr::bind_rows(lapply(areas, function(area) {
      pooled <- stats::setNames(rep(0, length(causes)), causes)
      published <- area %in% planning_published_areas
      on_rows <- TRUE
      for (j in seq_along(years)) {
        got <- if (published) NULL else from_rows[[j]][[area]]
        if (!is.null(got)) {
          pooled <- pooled + got[causes]
        } else {
          # Portugal and Continente read their own rows in `summed` where the
          # year has one; otherwise this is the municipal sum, like any area.
          pooled <- pooled + summed[[j]][area, causes]
          if (!published || !isTRUE(attr(summed[[j]], "own_row")[[area]])) on_rows <- FALSE
        }
      }
      complete <- on_rows
      if (isTRUE(joint[[area]])) pooled[] <- NA_real_
      total <- pooled[[planning_all_causes]]
      groups <- pooled[PLANNING_CAUSE_GROUPS$cause]
      other <- max(total - sum(groups), 0)
      deaths <- c(total, unname(groups), other)
      ci <- planning_binomial_ci(deaths, rep(total, length(deaths)))

      tibble::tibble(
        area = area, code = codes, group = labels, deaths = deaths,
        period = paste0(min(years), "-", max(years)), end_year = end_year,
        share = if (!is.na(total) && total > 0) deaths / total * 100 else NA_real_,
        lower = if (is.na(total)) NA_real_ else ifelse(codes == "C00", 100, ci$lower),
        upper = if (is.na(total)) NA_real_ else ifelse(codes == "C00", 100, ci$upper),
        flag = if (!complete && any(years %in% UNDER75_INCOMPLETE_YEARS)) UNDER75_FLAG else ""
      )
    }))
  })

  dplyr::bind_rows(rows) %>%
    dplyr::arrange(match(.data$area, areas), .data$end_year, match(.data$code, codes))
}

planning_under75_years <- function(window = 3L) {
  years <- sort(intersect(death_totals_years(), life_death_years()))
  years[vapply(years, function(y) all(seq.int(y - window + 1L, y) %in% years), logical(1))]
}

planning_under75_note <- function(table) {
  if (nrow(table) == 0 || !any(nzchar(table$flag))) return(NULL)
  paste0(
    UNDER75_FLAG, " triénio que inclui 2014, ano em que o INE publicou a idade de apenas ",
    "80% dos óbitos por município. Nas áreas sem linha regional do INE, as quotas desse ",
    "triénio podem estar desviadas em cerca de 2 pontos percentuais."
  )
}
