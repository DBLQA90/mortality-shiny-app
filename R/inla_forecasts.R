# =========================================================
# Small-area forecasts of the INLA model, precomputed
# =========================================================
# tools/precompute_inla.R fits the 308 municipalities together (a Bayesian
# spatio-temporal Poisson model) and sums the posterior samples to the ULS,
# NUTS III, NUTS II and Portugal; the app only reads the results. A municipality
# borrows strength from its neighbours and from the national trend, so it is
# projected even where its own series has years with fewer than 10 deaths.
#
# Two switches (Previsão Avançada, "Quebras e Estrutura"):
#   B  level shifts at the breaks of INE's population series (2014, 2021;
#      with C only 2014, tools/precompute_inla.R), off by default: the 2014
#      step made the backtests worse;
#   C  an effect for each COVID year (2020-2022), on by default: it removed
#      most of the bias of the forecasts made after COVID.
# Each of the four combinations is precomputed.

# The causes tested in the backtests (8 of them, 2008-2022 origins).
INLA_CAUSES <- c(
  "Acidentes de transporte", "Suicídios e lesões autoprovocadas voluntariamente",
  "Tumor (neoplasma) maligno do estômago", "Tumor (neoplasma) maligno do cólon", "Diabetes mellitus",
  "Tumor (neoplasma) maligno da laringe, da traqueia, dos brônquios e dos pulmões", "Pneumonia",
  "Doenças isquémicas do coração"
)
INLA_COMBOS <- c("B1C1", "B0C1", "B1C0", "B0C0")
# The defaults of the switches, set by the backtests (see the 'Quebras e
# Estrutura' text, INLA_SWITCHES_TEXT).
INLA_DEFAULT_STEPS <- FALSE
INLA_DEFAULT_COVID <- TRUE
INLA_FIRST_YEAR <- 1999L
INLA_MODEL_ID <- "inla"

# The model's 95% predictive intervals held the observed rate less often than
# 95% in the backtests; they are widened on the log scale by these factors,
# measured there for each combination and level (every projected year 1-4
# ahead, 9 origins 2008-2022; forecast_smallarea_2026/38_widening.R).
INLA_WIDENING <- list(
  B0C0 = c(municipio = 1, ULS = 1.09, `NUTS III` = 1.11, `NUTS II` = 1.08, Portugal = 1),
  B0C1 = c(municipio = 1, ULS = 1.06, `NUTS III` = 1.08, `NUTS II` = 1.08, Portugal = 1.08),
  B1C0 = c(municipio = 1, ULS = 1.10, `NUTS III` = 1.13, `NUTS II` = 1.14, Portugal = 1.15),
  B1C1 = c(municipio = 1, ULS = 1.09, `NUTS III` = 1.13, `NUTS II` = 1.17, Portugal = 1.28)
)

# What the backtests found, for the summaries (decisions_C_all.txt there).
INLA_BACKTEST_TEXT <- paste0(
  "Nos testes com dados passados (previsões feitas a partir de 2008-2022, 4 anos à frente, 8 causas), ",
  "o viés nas ULS e NUTS ficou entre +0,5% e +2%, contra cerca de +10% da média de 6 modelos, com erro também menor; ",
  "os intervalos estão alargados para conter 95% dos valores observados."
)
# The backtests reach 4 years ahead: the widening is measured there, and the
# coverage still fell with the horizon (Portugal: 97% one year ahead, 93% four).
INLA_TESTED_HORIZON <- 4L

# The warning for a forecast that reaches past the tested years, or "".
inla_horizon_note <- function(last_year, origin = inla_origin_year()) {
  if (!is.finite(origin) || !is.finite(last_year) || last_year <= origin + INLA_TESTED_HORIZON) return("")
  paste0("Os intervalos foram testados até ", INLA_TESTED_HORIZON, " anos à frente (", origin + INLA_TESTED_HORIZON,
         "); para lá disso não há anos com que os testar, e são provavelmente estreitos de mais.")
}

# Why B is off and C on by default, for the 'Quebras e Estrutura' tab.
INLA_SWITCHES_TEXT <- paste0(
  "Por omissão, B desligado e C ligado, pelos testes com dados passados (9 anos de partida, 2008-2022). ",
  "A quebra de 2014 (B) aumentou o viés e o erro e baixou a cobertura dos intervalos (ULS: 89% contra 91% sem ela); ",
  "a de 2021 não pôde ainda ser testada e, com C ligado, fica de fora: não se distingue dos anos de COVID-19. ",
  "Separar os anos de COVID-19 (C) corrigiu quase todo o viés das previsões feitas depois da pandemia ",
  "(ULS, a partir de 2021: de -7,7% para -0,2%; cobertura de 86% para 94%), à custa de um erro um pouco maior ",
  "em cada previsão. Os intervalos de cada combinação estão alargados pelo factor medido para ela."
)

# The last observed year the forecasts start from.
inla_origin_year <- function() {
  m <- inla_manifest()
  if (is.null(m)) NA_integer_ else as.integer(m$origin)
}

inla_forecast_dir <- function() {
  Sys.getenv("MORTALITY_INLA_DIR", unset = file.path(app_dir_or_wd(), "data", "inla_forecasts"))
}

inla_key <- function(cause, sex) paste0(planning_cause_file_token(cause), "__", sex)

inla_combo <- function(steps = INLA_DEFAULT_STEPS, covid = INLA_DEFAULT_COVID) {
  paste0("B", as.integer(isTRUE(steps)), "C", as.integer(isTRUE(covid)))
}

inla_manifest <- function() {
  path <- file.path(inla_forecast_dir(), "manifest.rds")
  if (file.exists(path)) readRDS(path) else NULL
}

# Read once per session: a file is a cause x sex x combination, all areas.
inla_file_cache <- new.env(parent = emptyenv())
inla_read <- function(cause, sex, combo) {
  path <- file.path(inla_forecast_dir(), paste0(inla_key(cause, sex), "__", combo, ".rds"))
  if (!file.exists(path)) return(NULL)
  stamp <- paste(path, file.mtime(path))
  if (!exists(stamp, envir = inla_file_cache, inherits = FALSE)) assign(stamp, readRDS(path), envir = inla_file_cache)
  get(stamp, envir = inla_file_cache)
}

# Why the INLA forecast cannot serve a selection, or NULL when it can. The
# model projects the age-standardised rate of the whole population, on the
# corrected data, for one municipality, ULS of whole municipalities, NUTS of
# 2024 or Portugal, from the last year of causes of death.
inla_unavailable_reason <- function(areas, cause, sex, rate_type = "dsr", population = "Total", combo = inla_combo(),
                                    mode = data_mode(), vintage = "2024", last_year = NULL) {
  if (length(areas) != 1) return("a projecção espacial existe para um local de cada vez")
  if (!identical(rate_type, "dsr")) return("a projecção espacial é da taxa padronizada")
  if (!identical(population, "Total")) return("a projecção espacial é da população total")
  if (!identical(mode, "corrected")) return("a projecção espacial usa os dados corrigidos")
  if (length(cause) != 1 || !cause %in% INLA_CAUSES) return("a projecção espacial existe para 8 causas, testadas com dados passados")
  x <- inla_read(cause, sex, combo)
  if (is.null(x)) return("a projecção espacial ainda não foi calculada (tools/precompute_inla.R)")
  level <- x$forecast$level[match(areas, x$forecast$area)]
  if (is.na(level)) return("a projecção espacial existe para municípios, ULS de municípios inteiros, NUTS de 2024 e Portugal")
  if (level %in% c("NUTS III", "NUTS II") && !identical(as.character(vintage), "2024")) {
    return("a projecção espacial usa as NUTS de 2024")
  }
  if (!is.null(last_year) && !identical(as.integer(last_year), as.integer(x$origin))) {
    return(paste0("a projecção espacial parte de ", x$origin, ", o último ano da janela de ajuste tem de ser esse"))
  }
  NULL
}

# The forecast as the other models give it: year, mean (the posterior median),
# lower and upper of the predictive interval at `conf_level`. The stored
# quantiles are the 95% predictive interval; other levels scale its log-scale
# half-widths by the normal quantile, and both sides are widened by the factor
# measured in the backtests.
inla_forecast <- function(area, cause, sex, combo = inla_combo(), conf_level = 95, horizon = Inf, widening = INLA_WIDENING[[combo]]) {
  x <- inla_read(cause, sex, combo)
  if (is.null(x)) return(NULL)
  f <- x$forecast[x$forecast$area == area & x$forecast$year > x$origin & x$forecast$year <= x$origin + horizon, , drop = FALSE]
  if (nrow(f) == 0) return(NULL)
  k <- unname(widening[f$level[[1]]])
  if (length(k) != 1 || !is.finite(k)) k <- 1
  k <- k * stats::qnorm(1 - (1 - conf_level / 100) / 2) / stats::qnorm(0.975)
  centre <- f$q500
  side <- function(bound) ifelse(bound > 0 & centre > 0, centre * (bound / centre)^k, 0)
  tibble::tibble(year = as.integer(f$year), model = INLA_MODEL_ID, mean = centre, lower = side(f$p25),
                 upper = ifelse(centre > 0, centre * (f$p975 / centre)^k, f$p975))
}

# The estimated step and COVID effects of a fit, as rate ratios.
inla_effects <- function(cause, sex, combo = inla_combo()) {
  x <- inla_read(cause, sex, combo)
  if (is.null(x) || is.null(x$fixed) || nrow(x$fixed) == 0) return(NULL)
  labels <- c(step2014 = "Quebra da série de população, 2014", step2021 = "Quebra da série de população, 2021",
              covid2020 = "COVID-19, 2020", covid2021 = "COVID-19, 2021", covid2022 = "COVID-19, 2022")
  tibble::tibble(Termo = unname(labels[rownames(x$fixed)]), `Razão de taxas` = exp(x$fixed$mean),
                 `IC 95% inf.` = exp(x$fixed$`0.025quant`), `IC 95% sup.` = exp(x$fixed$`0.975quant`))
}

# The maps' projection (Mapas, "Projecção") for the areas INLA covers: the
# projected rate in `target` and the level the model estimates for the last
# observed year, which the change is measured against - for a small area the
# observed rate of one year is too noisy to compare with. NULL when the model
# cannot serve the map.
inla_map_rows <- function(areas, cause, sex, target, combo = inla_combo(), under75 = FALSE, mode = data_mode(), vintage = "2024") {
  if (isTRUE(under75) || !identical(mode, "corrected") || !cause %in% INLA_CAUSES) return(NULL)
  x <- inla_read(cause, sex, combo)
  if (is.null(x) || target <= x$origin || target > x$origin + x$horizon) return(NULL)
  covered <- unique(x$forecast[, c("area", "level")])
  if (!identical(as.character(vintage), "2024")) covered <- covered[!covered$level %in% c("NUTS III", "NUTS II"), , drop = FALSE]
  areas <- intersect(areas, covered$area)
  if (length(areas) == 0) return(NULL)
  rows <- lapply(areas, function(a) {
    f <- inla_forecast(a, cause, sex, combo, horizon = target - x$origin)
    base <- x$forecast$q500[x$forecast$area == a & x$forecast$year == x$origin]
    if (is.null(f) || length(base) != 1) return(NULL)
    f <- f[f$year == target, , drop = FALSE]
    tibble::tibble(area = a, last_year = as.integer(x$origin), last_value = base, fit_from = INLA_FIRST_YEAR,
                   value = f$mean, lower = f$lower, upper = f$upper, members = NA_integer_, method = INLA_MODEL_ID)
  })
  dplyr::bind_rows(rows)
}
