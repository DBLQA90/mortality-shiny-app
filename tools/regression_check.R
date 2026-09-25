#!/usr/bin/env Rscript
# Regression check: recompute the numbers the app has been validated on and
# compare them with data/regression_baseline.csv.
#
#   Rscript tools/regression_check.R            # check; see the exit codes below
#   Rscript tools/regression_check.R --quick    # skip the two slow groups
#   Rscript tools/regression_check.R --update   # record the current values as the baseline
#   Rscript tools/regression_check.R --report=FILE   # also write the report there
#
# It runs after every scheduled refresh (tools/scheduled_refresh.sh), which is
# the point: INE revises published series quietly, and a revision that moves a
# municipality's death rate looks exactly like a bug in the app until someone
# compares. The baseline is committed, so `git diff` after --update shows what
# moved and by how much.
#
# Two kinds of measure:
#   invariante  something that must hold whatever the data says - a partition
#               that has to be exact, shares that have to add to 100, Portugal's
#               SMR against itself, an export that has to be written.
#   valor       a validated number, compared within a tolerance: the numbers the
#               audits checked against INE's own publications.
#
# Exit codes: 0 nothing moved, 1 a tracked value moved or is new/missing,
# 2 an invariant failed or a group could not be computed.

suppressMessages({library(dplyr); library(tibble); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
update_baseline <- "--update" %in% args
quick <- "--quick" %in% args
report_path <- sub("^--report=", "", grep("^--report=", args, value = TRUE))
baseline_path <- sub("^--baseline=", "", grep("^--baseline=", args, value = TRUE))
if (length(baseline_path) == 0) baseline_path <- "data/regression_baseline.csv"

script_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])))
setwd(normalizePath(file.path(script_dir, "..")))

app <- new.env(parent = globalenv())
assign("app_dir", getwd(), app)
for (f in c("R/config.R", "R/helpers.R", "R/regions.R", "R/metrics.R", "R/standardisation.R", "R/avoidable.R",
            "R/infant.R", "R/regional_rows.R", "R/planning_parish.R", "R/planning_indicators.R",
            "R/life_expectancy.R", "R/planning_under75.R", "R/planning_standardised.R", "R/planning_sns.R",
            "R/planning_weekly.R", "R/planning_charts.R", "R/planning_export.R", "R/planning_profile.R",
            "R/data_versions.R")) {
  sys.source(f, envir = app)
}
attach(app, warn.conflicts = FALSE)

# The years the baseline is measured on. They move only when the reference
# values are refetched (tools/fetch_validation_refs.R) and the baseline is
# re-recorded, so that a new year of data does not read as a regression.
REF_YEAR <- 2024L        # the year INE's municipal indicators cover
REF_TRIENNIUM <- 2024L   # the triennium ending here, for life expectancy and rates
REF_CENSUS <- 2021L

# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------
measures <- list()
measure <- function(key, value, tol = 0, kind = "valor", note = "") {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1) value <- NA_real_
  measures[[key]] <<- tibble(key = key, kind = kind, value = value, tol = tol, note = note)
  invisible(value)
}
invariant <- function(key, ok, note = "") measure(key, isTRUE(ok), tol = 0, kind = "invariante", note = note)

# A group that fails records the failure rather than stopping the run: one
# broken loader should not hide every other check. `prefixes` are the keys the
# group produces, so that --quick can report the ones it did not recompute as
# skipped rather than as missing.
skipped <- character(0)
group <- function(name, prefixes, expr, slow = FALSE) {
  if (slow && quick) {
    skipped <<- c(skipped, prefixes)
    message(sprintf("  %-26s ignorado (--quick)", name))
    return(invisible(NULL))
  }
  started <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    message("  ", name, ": ERRO - ", conditionMessage(e))
    measure(paste0(name, ".erro"), 1, kind = "invariante", note = conditionMessage(e))
    FALSE
  })
  message(sprintf("  %-26s %s (%.0f s)", name, if (ok) "ok" else "ERRO", as.numeric(difftime(Sys.time(), started, units = "secs"))))
}

lookup <- get_nuts_lookup("2024")
municipalities <- sort(unique(lookup$municipality))
dico <- tibble(dico = substr(lookup$municipality_code, 4, 7), municipality = lookup$municipality)

message("A recalcular (", if (quick) "rápido" else "completo", ")...")

# ---------------------------------------------------------------------------
# 1. Municipal values against INE's own municipal indicators
# ---------------------------------------------------------------------------
group("indicadores INE", "ine.", {
  refs <- readRDS("data/snapshots/validation/ine_municipal_refs.rds") %>%
    filter(nchar(code) == 7, indicator != "teen_5y") %>%
    mutate(municipality = dico$municipality[match(substr(code, 4, 7), dico$dico)]) %>%
    filter(!is.na(municipality)) %>%
    group_by(indicator, municipality) %>% slice(1) %>% ungroup()

  ours <- planning_indicator_table(
    municipalities, REF_YEAR,
    ids = c("ageing_index", "old_dependency", "youth_dependency", "birth_rate", "death_rate",
            "waste_per_capita", "rsi_rate", "pensioners_rate"), lookup = lookup) %>%
    transmute(area, indicator = recode(indicator, waste_per_capita = "waste", rsi_rate = "rsi_1564",
                                       pensioners_rate = "pens_1564"), value)
  pooled <- planning_components(municipalities, (REF_YEAR - 4L):REF_YEAR, lookup) %>%
    group_by(area) %>% summarise(infant_deaths = sum(infant_deaths), births = sum(births), .groups = "drop") %>%
    transmute(area, indicator = "infant_5y", value = infant_deaths / births * 1000)

  comparison <- bind_rows(ours, pooled) %>%
    inner_join(refs %>% select(indicator, area = municipality, ine = value), by = c("indicator", "area")) %>%
    mutate(diff = value - ine)
  summary <- comparison %>% group_by(indicator) %>%
    summarise(n = n(), max_abs = max(abs(diff), na.rm = TRUE),
              p95_abs = quantile(abs(diff), 0.95, na.rm = TRUE), .groups = "drop")

  invariant("ine.municipios_comparados", nrow(comparison) == 9 * 308,
            "9 indicadores x 308 municípios")
  for (i in seq_len(nrow(summary))) {
    id <- summary$indicator[[i]]
    # Rounding in INE's own table is +/-0.05; the two rates built from a
    # different denominator vintage are tracked with a wider tolerance.
    tol <- if (id %in% c("pens_1564", "rsi_1564")) 5 else 0.02
    measure(paste0("ine.", id, ".max_abs"), summary$max_abs[[i]], tol = tol,
            note = paste0("maior diferença face ao INE, ", REF_YEAR))
    measure(paste0("ine.", id, ".p95_abs"), summary$p95_abs[[i]], tol = tol / 2,
            note = "percentil 95 da diferença")
  }
})

# ---------------------------------------------------------------------------
# 2. Partitions: ULS and ARS cover the Continente exactly
# ---------------------------------------------------------------------------
group("partições ULS/ARS", "particao.", {
  health <- get_health_lookup()
  uls <- planning_uls_units()
  ars <- unique(health$unit[health$kind == "ARS"])
  continental <- sort(unique(lookup$municipality[lookup$nuts1 == "Continente"]))
  columns <- c("pop_total", "births", "deaths", "infant_deaths", "rsi", "pensioners",
               "waste_total", "employees_total", "perinatal_deaths", "births_weight_known")
  for (year in c(2014L, 2023L)) {
    components <- planning_components(c(uls, ars, continental), year, lookup)
    total <- function(areas) colSums(components[components$area %in% areas, columns], na.rm = TRUE)
    reference <- total(continental)
    mismatched <- names(which(abs(total(uls) - reference) > 1e-6 | abs(total(ars) - reference) > 1e-6))
    # 2023 waste is the known exception: INE withholds it for some municipalities.
    measure(paste0("particao.", year, ".componentes_divergentes"), length(mismatched), tol = 0,
            note = if (length(mismatched) == 0) "" else paste(mismatched, collapse = ","))
  }
})

# ---------------------------------------------------------------------------
# 3. Portugal's published row against the sum of its 308 municipalities
# ---------------------------------------------------------------------------
group("Portugal vs municípios", "pt_vs_municipios.", {
  columns <- c("pop_total", "births", "deaths", "infant_deaths", "waste_total",
               "employees_total", "pensioners", "rsi")
  for (year in c(2014L, 2021L, REF_YEAR)) {
    components <- planning_year_components(year)
    published <- components[components$area == "Portugal", columns, drop = FALSE]
    summed <- colSums(components[components$area %in% municipalities, columns, drop = FALSE], na.rm = TRUE)
    for (column in columns) {
      value <- as.numeric(published[[column]][[1]])
      if (is.na(value) || summed[[column]] == 0) next
      measure(paste0("pt_vs_municipios.", column, ".", year), (value / summed[[column]] - 1) * 100,
              tol = 0.02, note = "diferença percentual face à soma dos municípios")
    }
  }
})

# ---------------------------------------------------------------------------
# 4. Identities and ranges over every area and year (slow)
# ---------------------------------------------------------------------------
group("identidades (todas as áreas)", c("tabela.", "identidade.", "amplitude."), slow = TRUE, expr = {
  areas <- planning_area_levels(lookup)$area
  table <- planning_indicator_table(areas, 1995:2025, lookup = lookup)
  wide <- table %>% select(area, year, indicator, value) %>%
    pivot_wider(names_from = indicator, values_from = value)
  units <- PLANNING_INDICATORS %>% select(id, unit)
  typed <- table %>% left_join(units, by = c("indicator" = "id"))

  measure("tabela.linhas", nrow(table), tol = 0, note = "áreas x anos x indicadores")
  measure("tabela.valores", sum(!is.na(table$value)), tol = 0, note = "valores não vazios")
  invariant("tabela.finitos", all(is.na(table$value) | is.finite(table$value)))
  invariant("tabela.percentagens_ate_100", sum(typed$unit == "%" & typed$value > 100, na.rm = TRUE) == 0)
  # The only negatives are the census population change, which can be negative.
  negatives <- typed %>% filter(value < 0) %>% count(indicator)
  invariant("tabela.negativos_esperados", all(negatives$indicator == "census_population_change"),
            paste(negatives$indicator, collapse = ","))
  measure("tabela.negativos", sum(negatives$n), tol = 0, note = "variação da população censitária")

  identity <- wide %>% filter(!is.na(infant_rate), !is.na(neonatal_rate), !is.na(postneonatal_rate)) %>%
    mutate(gap = neonatal_rate + postneonatal_rate - infant_rate)
  invariant("identidade.neonatal_mais_pos_neonatal", max(abs(identity$gap)) < 1e-9)
  invariant("identidade.neonatal_precoce", sum(wide$early_neonatal_rate > wide$neonatal_rate + 1e-9, na.rm = TRUE) == 0)
  education <- wide %>% filter(!is.na(pct_education_none)) %>%
    mutate(total = pct_education_none + pct_education_basic + pct_education_secondary + pct_education_higher)
  invariant("identidade.escolaridade_100", max(abs(education$total - 100)) < 1e-6)
  sectors <- wide %>% filter(!is.na(pct_employees_primary)) %>%
    mutate(total = pct_employees_primary + pct_employees_secondary + pct_employees_tertiary)
  invariant("identidade.sectores_100", max(abs(sectors$total - 100)) < 1e-6)

  range_of <- function(id) range(table$value[table$indicator == id], na.rm = TRUE)
  for (id in c("life_expectancy", "life_expectancy_65", "fertility_index", "purchasing_power")) {
    bounds <- range_of(id)
    measure(paste0("amplitude.", id, ".min"), bounds[[1]], tol = 0.01)
    measure(paste0("amplitude.", id, ".max"), bounds[[2]], tol = 0.01)
  }
})

# ---------------------------------------------------------------------------
# 5. Standardised mortality and life expectancy
# ---------------------------------------------------------------------------
group("mortalidade padronizada", c("padronizadas.", "causas."), {
  areas <- c("Portugal", "Continente", "Norte", "Alentejo", "ULS Baixo Alentejo", "Matosinhos")
  table <- planning_indicator_table(areas, REF_TRIENNIUM, ids = c(standardised_ids, "life_expectancy"),
                                    lookup = lookup)
  value_of <- function(area, id) table$value[table$area == area & table$indicator == id]
  invariant("padronizadas.portugal_smr_100", abs(value_of("Portugal", "smr_all") - 100) < 1e-9)
  # INE publishes Norte's premature deaths as a row of its own: the app has to
  # reproduce it exactly (33.155 in the triennium ending 2024).
  measure("padronizadas.norte.premature_deaths", value_of("Norte", "premature_deaths"), tol = 0,
          note = "linha do INE para o Norte")
  for (area in areas) {
    for (id in c("dsr_all", "dsr_premature", "smr_all", "ypll_rate", "life_expectancy")) {
      measure(paste0("padronizadas.", gsub(" ", "_", tolower(area)), ".", id), value_of(area, id),
              tol = if (id == "life_expectancy") 0.01 else 0.5, note = paste("triénio", REF_TRIENNIUM))
    }
  }
  # The regions read their deaths by age from INE's own rows, so none of their
  # deaths is redistributed and none carries the ‡ flag.
  flagged <- table %>% filter(area %in% c("Norte", "Alentejo"), flag != "")
  invariant("padronizadas.regioes_sem_redistribuicao", nrow(flagged) == 0,
            paste(unique(flagged$area), collapse = ","))

  causes <- planning_cause_standardised(c("Portugal", "ULS Matosinhos"), REF_TRIENNIUM, lookup = lookup)
  benchmark <- causes %>% filter(area == "Portugal")
  invariant("causas.esperados_iguais_observados", max(abs(benchmark$observed - benchmark$expected)) < 1e-6)
  measure("causas.grupos", nrow(benchmark) - 1, tol = 0, note = "grupos de causas além do total")
  measure("causas.cobertura", sum(benchmark$observed[benchmark$code != "all"]) /
            benchmark$observed[benchmark$code == "all"], tol = 0.005,
          note = "parte dos óbitos coberta pelos grupos")
})

# ---------------------------------------------------------------------------
# 6. Parish split of the shared ULS
# ---------------------------------------------------------------------------
group("repartição por freguesia", "freguesias.", {
  shares <- planning_parish_shares()
  totals <- shares %>% group_by(municipality, basis) %>% summarise(total = sum(weight), .groups = "drop")
  invariant("freguesias.pesos_censo_somam_1", max(abs(totals$total - 1)) < 1e-9)
  for (kind in c("births", "deaths")) {
    actual <- planning_parish_actual(kind)
    sums <- actual %>% group_by(year, municipality) %>% summarise(total = sum(weight), .groups = "drop")
    invariant(paste0("freguesias.registos_somam_1.", kind), max(abs(sums$total - 1)) < 1e-9)
    measure(paste0("freguesias.registos_ultimo_ano.", kind), max(actual$year), tol = 0)
  }
  # In parish mode the ULS partition the Continente exactly. In whole-municipality
  # mode they overlap on purpose - a shared municipality counts for every ULS that
  # serves it - so the excess is tracked rather than required to be zero.
  all_uls <- planning_uls_units("units")
  for (mode in PLANNING_SPLIT_MODES) {
    table <- planning_indicator_table(c(all_uls, "Continente"), REF_YEAR,
                                      ids = c("pop_total", "deaths", "births"), lookup = lookup, mode = mode)
    for (id in c("pop_total", "deaths", "births")) {
      rows <- table[table$indicator == id, ]
      gap <- (sum(rows$value[rows$area != "Continente"]) / rows$value[rows$area == "Continente"] - 1) * 100
      if (identical(mode, "parish")) {
        invariant(paste0("freguesias.parish.", id, ".soma_exacta"), abs(gap) < 1e-9,
                  "soma das ULS igual ao Continente")
      } else {
        measure(paste0("freguesias.whole.", id, ".sobreposicao_pct"), gap, tol = 0.05,
                note = "municípios partilhados contados por cada ULS")
      }
    }
  }
})

# ---------------------------------------------------------------------------
# 7. SNS primary care and weekly excess mortality
# ---------------------------------------------------------------------------
group("SNS e óbitos semanais", c("sns.", "semanais."), {
  components <- planning_sns_components()
  measure("sns.unidades", n_distinct(components$unit), tol = 0)
  sns <- planning_sns_table(list(Continente = planning_uls_units("units")))
  measure("sns.indicadores", n_distinct(sns$indicator), tol = 0)
  measure("sns.ultimo_periodo", as.numeric(gsub("-", "", max(sns$period))), tol = 0)
  invariant("sns.valores_0_100", all(sns$value >= 0 & sns$value <= 100, na.rm = TRUE))

  weekly <- planning_weekly_data()
  measure("semanais.regioes", n_distinct(weekly$region), tol = 0)
  measure("semanais.ultimo_ano", max(weekly$year), tol = 0)
  excess <- planning_weekly_summary(planning_weekly_excess("Portugal", (REF_YEAR + 1L):(REF_YEAR + 2L)))
  for (i in seq_len(nrow(excess))) {
    measure(paste0("semanais.excesso_pct.", excess$year[[i]]), excess$excess_pct[[i]], tol = 0.5,
            note = paste("semanas", excess$weeks[[i]]))
  }
})

# ---------------------------------------------------------------------------
# 8. Exports, both modes (slow)
# ---------------------------------------------------------------------------
group("exportações", "exportacao.", slow = TRUE, expr = {
  areas <- bind_rows(tibble(area = "ULS São José", level = "Local"),
                     planning_comparators("ULS São José", lookup)[, c("area", "level")])
  for (mode in PLANNING_SPLIT_MODES) {
    workbook <- tempfile(fileext = ".xlsx"); profile <- tempfile(fileext = ".docx")
    on.exit(unlink(c(workbook, profile)), add = TRUE)
    written <- tryCatch({
      write_planning_workbook(workbook, areas, lookup = lookup, vintage = "2024",
                              data_date = format(Sys.Date()), focus = "ULS São José",
                              years = (REF_YEAR - 9L):(REF_YEAR + 1L), mode = mode)
      write_planning_profile(profile, areas, lookup = lookup, vintage = "2024",
                             data_date = format(Sys.Date()), mode = mode)
      TRUE
    }, error = function(e) { message("    ", mode, ": ", conditionMessage(e)); FALSE })
    invariant(paste0("exportacao.", mode), written && file.size(workbook) > 1e5 && file.size(profile) > 1e5)
  }
})

# ---------------------------------------------------------------------------
# Compare with the baseline
# ---------------------------------------------------------------------------
current <- bind_rows(measures) %>% arrange(key)
baseline <- if (file.exists(baseline_path)) {
  utils::read.csv(baseline_path, stringsAsFactors = FALSE, comment.char = "#", encoding = "UTF-8") %>% as_tibble()
} else {
  tibble(key = character(0), kind = character(0), value = numeric(0), tol = numeric(0), note = character(0))
}

comparison <- current %>%
  full_join(baseline %>% select(key, kind_old = kind, value_old = value, tol_old = tol), by = "key") %>%
  mutate(
    tol = ifelse(is.na(tol), tol_old, tol),
    kind = ifelse(is.na(kind), kind_old, kind),
    status = case_when(
      is.na(value_old) ~ "NOVO",
      is.na(value) & vapply(key, function(k) any(startsWith(k, skipped)), logical(1)) ~ "ignorado",
      is.na(value) ~ "EM FALTA",
      kind == "invariante" & value < 1 ~ "FALHA",
      abs(value - value_old) > tol + 1e-12 ~ "MUDOU",
      TRUE ~ "igual"
    )
  ) %>% arrange(factor(status, levels = c("FALHA", "MUDOU", "EM FALTA", "NOVO", "igual", "ignorado")), key)

log <- tryCatch(read_import_log("data"), error = function(e) NULL)
latest_import <- if (is.null(log) || nrow(log) == 0) "desconhecida" else max(log$imported_at, na.rm = TRUE)

lines <- c(
  sprintf("== Verificação de regressão  %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("   dados: última importação %s | método %s | base %s",
          latest_import, PLANNING_METHOD_VERSION, baseline_path),
  sprintf("   %d medidas: %s", nrow(comparison),
          paste(sprintf("%d %s", table(comparison$status), names(table(comparison$status))), collapse = ", "))
)
moved <- comparison %>% filter(!status %in% c("igual", "ignorado"))
if (nrow(moved) > 0) {
  lines <- c(lines, "", vapply(seq_len(nrow(moved)), function(i) {
    row <- moved[i, ]
    sprintf("%-9s %-46s %s -> %s  (tol %s)%s", row$status, row$key,
            if (is.na(row$value_old)) "-" else format(signif(row$value_old, 8)),
            if (is.na(row$value)) "-" else format(signif(row$value, 8)),
            format(row$tol), if (is.na(row$note) || row$note == "") "" else paste0("  ", row$note))
  }, character(1)))
}
report <- paste(lines, collapse = "\n")
cat(report, "\n", sep = "")
if (length(report_path) > 0) writeLines(lines, report_path[[1]])

if (update_baseline) {
  header <- c(
    "# Valores de referência de tools/regression_check.R.",
    "# Gerado com --update; cada linha é uma medida e a tolerância com que é comparada.",
    sprintf("# Gerado em %s | método %s | última importação %s",
            format(Sys.Date()), PLANNING_METHOD_VERSION, latest_import),
    sprintf("# Anos fixados: indicadores %d, triénio %d, censo %d", REF_YEAR, REF_TRIENNIUM, REF_CENSUS)
  )
  connection <- file(baseline_path, open = "w", encoding = "UTF-8")
  writeLines(header, connection)
  utils::write.csv(current, connection, row.names = FALSE)
  close(connection)
  message("Base de referência escrita em ", baseline_path, " (", nrow(current), " medidas).")
  quit(status = 0)
}

status <- 0L
if (any(comparison$status %in% c("MUDOU", "NOVO", "EM FALTA"))) status <- 1L
if (any(comparison$status == "FALHA")) status <- 2L
quit(status = status)
