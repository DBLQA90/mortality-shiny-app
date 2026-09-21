# =========================================================
# Planning indicators: Excel export
# =========================================================
# Two workbooks from the same writer:
#
#   selection  the chosen location and its comparators, every indicator and
#              year, with a summary sheet and a long sheet carrying intervals,
#              numerators and denominators
#   full       every area the app can build (Portugal, NUTS I-III, ARS, ULS,
#              the 308 municipalities), every indicator and year - the
#              automated counterpart of the PLS support workbook
#
# Each indicator gets its own sheet, areas as rows and years (or triennia) as
# columns, like the workbook, with the bounds of the 95% interval in two blocks
# below the values. Values carrying a caveat (* thin denominator,
# † incomplete municipal source) are written in grey italics, and the read-me
# sheet states the import date of the data, so two files for the same area can
# be told apart.

PLANNING_SHEET_NAMES <- c(
  pop_total = "I1 População", census_population = "I2 População nos Censos",
  census_population_change = "I2 Variação da população", pct_education_none = "I24 Sem escolaridade",
  pct_education_basic = "I24 Ensino básico", pct_education_secondary = "I24 Ensino secundário",
  pct_education_higher = "I24 Ensino superior", illiteracy_rate = "I26 Analfabetismo", pct_0_14 = "I1 Jovens", pct_65_plus = "I1 Idosos",
  pct_75_plus = "I1 75 e mais anos", ageing_index = "I4 Envelhecimento",
  youth_dependency = "I5 Dependência jovens", old_dependency = "I6 Dependência idosos",
  births = "I7 Nados-vivos", birth_rate = "I8 Natalidade", fertility_index = "I9 Fecundidade",
  teen_births_pct = "I32 Mães com menos de 20", older_births_pct = "I33 Mães com 35 e mais",
  preterm_pct = "I35 Pré-termo", low_birth_weight_pct = "I36 Baixo peso", rsi_beneficiaries = "I13 RSI", rsi_rate = "I14 RSI por 1000",
  pensioners = "I15 Pensionistas", pensioners_rate = "I16 Pensionistas por 1000",
  pension_mean = "I17 Pensão média", earnings_mean = "I27 Ganho médio mensal",
  employees = "I12 Trabalhadores", pct_employees_primary = "I12 Sector primário",
  pct_employees_secondary = "I12 Sector secundário", pct_employees_tertiary = "I12 Sector terciário", purchasing_power = "I28 Poder de compra",
  waste_per_capita = "I64 Resíduos", waste_selective_per_capita = "I65 Resíduos selectivos",
  life_expectancy = "I10 Esperança de vida", life_expectancy_men = "I10 EV homens",
  life_expectancy_women = "I10 EV mulheres", life_expectancy_65 = "I10 EV aos 65 anos",
  life_expectancy_65_men = "I10 EV aos 65, homens", life_expectancy_65_women = "I10 EV aos 65, mulheres", deaths = "I37 Óbitos", death_rate = "I38 Mortalidade", infant_rate = "I39 Mortalidade infantil",
  neonatal_rate = "I40 Mortalidade neonatal", early_neonatal_rate = "I41 Neonatal precoce",
  postneonatal_rate = "I42 Pós-neonatal", late_fetal_rate = "I43 Fetal tardia",
  perinatal_rate = "I44 Perinatal", smr_all = "SMR todas as causas", dsr_all = "Taxa padronizada",
  dsr_premature = "Mortalidade prematura", premature_deaths = "Óbitos prematuros", dsr_preventable = "Evitável por prevenção",
  dsr_treatable = "Evitável por cuidados", avoidable_deaths = "Óbitos evitáveis", ypll_rate = "Anos de vida perdidos"
)

planning_number_format <- function(digits) {
  if (digits <= 0) "#,##0" else paste0("#,##0.", strrep("0", digits))
}

# Census years plus the latest, as in the workbook's pyramid sheet (I3).
planning_pyramid_years <- function(available) {
  sort(unique(c(intersect(c(2011L, 2021L), available), max(available))))
}

# Build the workbook for `areas` (a tibble of area and level) and write it to
# `path`. `focus` names the selected location for the summary sheet; NULL for the
# full file.
write_planning_workbook <- function(path,
                                    areas,
                                    lookup = get_nuts_lookup(),
                                    vintage = default_nuts_vintage(),
                                    data_date = NA_character_,
                                    focus = NULL,
                                    years = NULL,
                                    include_long = TRUE,
                                    education_min_age = 0L,
                                    benchmark = "Portugal",
                                    progress = function(value, detail = NULL) invisible(NULL)) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("The openxlsx package is required to write Excel files.", call. = FALSE)
  }
  areas <- dplyr::distinct(areas, area, .keep_all = TRUE)
  all_years <- sort(unique(unlist(lapply(PLANNING_INDICATORS$id, planning_indicator_years))))
  years <- if (is.null(years)) all_years else intersect(as.integer(years), all_years)

  progress(0.05, "indicadores")
  table <- planning_indicator_table(union(areas$area, benchmark), years, lookup = lookup, education_min_age = education_min_age, benchmark = benchmark) %>%
    planning_add_significance(benchmark) %>%
    dplyr::filter(.data$area %in% areas$area) %>%
    dplyr::left_join(areas, by = "area")

  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#e1e0d9", border = "bottom")
  title_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 13)
  note_style <- openxlsx::createStyle(fontColour = "#52514e", wrapText = FALSE)
  flagged_style <- openxlsx::createStyle(fontColour = "#898781", textDecoration = "italic")

  # --- Read-me ---------------------------------------------------------------
  openxlsx::addWorksheet(wb, "Leia-me")
  readme <- c(
    "Indicadores de apoio aos Planos Locais de Saúde",
    "",
    paste0("Portugal de referência (comparadores e significância): ",
           if (identical(benchmark, PLANNING_PORTUGAL_MUNICIPAL)) "soma dos 308 municípios, sem os acontecimentos de residência desconhecida." else "total publicado pelo INE, que inclui os acontecimentos de residência desconhecida."),
    if (education_min_age > 0) paste0("Escolaridade [I24]: população com ", education_min_age, " e mais anos (Censos de 2011 e 2021; o INE não publica a escolaridade por idade e município em 1991 e 2001).") else "Escolaridade [I24]: toda a população, como no ficheiro de apoio.",
    paste0("Significância face a ", benchmark, " (coluna «Face a Portugal» no Resumo e nos Dados): Superior ou Inferior quando o intervalo de confiança de 95% do local fica inteiramente acima ou abaixo do valor de referência no mesmo período; Semelhante quando o inclui. Só para indicadores com intervalo (taxas, proporções, esperança de vida). Não indica se a diferença é boa ou má."),
    "",
    "Notas assinaladas com * junto ao nome do indicador",
    paste("* Ganho médio e trabalhadores por sector [I27, I12]:", PLANNING_INDICATOR_NOTES[["earnings_mean"]]),
    paste("* Esperança de vida [I10]:", PLANNING_INDICATOR_NOTES[["life_expectancy"]]),
    "",
    paste0("Dados importados do INE até: ", ifelse(is.na(data_date), "desconhecido", data_date)),
    paste0("Ficheiro gerado em: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    paste0("Definição das regiões: NUTS ", vintage),
    paste0("Âmbito: ", if (is.null(focus)) paste0("todas as áreas (", nrow(areas), ")") else paste0(focus, " e comparadores")),
    paste0("Anos: ", min(years), "-", max(years)),
    paste0("Versão do método: ", PLANNING_METHOD_VERSION, " (nota metodológica)"),
    "",
    "Como são calculados",
    "Cada área é a soma dos seus municípios, e cada indicador é a razão dessas somas, nunca a média de valores municipais.",
    "Portugal e o Continente usam as linhas publicadas pelo INE, que incluem acontecimentos de residência desconhecida.",
    "As ULS e ARS usam a composição actual em municípios, aplicada a todos os anos. As cinco ULS que partilham municípios ao nível da freguesia aparecem em dois agrupamentos exactos.",
    "Indicadores de triénio (mortalidade infantil e neonatal, nascimentos por idade da mãe, pré-termo) somam três anos e são identificados pelo último.",
    "",
    "Marcas (valores a cinzento e itálico)",
    "* taxa sobre menos de 1.000 nados-vivos no triénio: exacta, mas instável.",
    "† triénio com anos (1995-2001) em que os óbitos com menos de 1 ano por município estão incompletos no INE: valor subestimado.",
    "‡ esperança de vida em que mais de 2% dos óbitos do triénio não tinham idade publicada por município e foram distribuídos pelas idades na proporção dos restantes.",
    "§ mortalidade proporcional com menos de 75 anos, num triénio que inclui 2014, numa área sem linha regional do INE: as quotas podem estar desviadas em cerca de 2 pontos percentuais.",
    "≈ trabalhadores por sector em que mais de 1% dos trabalhadores da área estão em sectores ocultados pelo INE (segredo estatístico) e foram repartidos na proporção do resto da NUTS III.",
    "",
    "Lacunas nos dados do INE",
    "Taxas brutas, RSI e resíduos por habitante usam a população média do ano (média das estimativas a 31 de Dezembro do ano anterior e do próprio), como o INE; pensionistas usam a de 31 de Dezembro.",
    "Uma célula em branco não é um zero: nos resíduos, pensões, trabalhadores e ganho médio, um município sem valor deixa sem valor as áreas que o contêm.",
    "Odivelas, Trofa e Vizela (criados em 1998): até 1998 os acontecimentos estão no município de origem (Loures, Santo Tirso, Guimarães); só há valores para áreas que contenham os dois. Os resíduos de Loures incluem sempre Odivelas (serviço conjunto SIMAR).",
    "Óbitos de todas as causas em branco no INE para um município (Vimioso 2024, por exemplo) são substituídos pela soma dos grupos etários publicados.",
    "",
    "Esperança de vida à nascença",
    "Tábua de mortalidade abreviada (Chiang II, grupos quinquenais até 85 e mais anos), por triénio, com o método e a variância de PHEindicatormethods. Omitida para populações até 5.000 e quando o intervalo de confiança excede 20 anos.",
    "Reproduz os valores do Eurostat para Portugal (2017-2019: 81,9 anos na aplicação; 82,0 no Eurostat em 2019), mas fica cerca de 0,8-0,9 anos acima dos valores publicados pelo INE (Metodologia 2007) para Portugal e NUTS III. A ordenação das regiões coincide com a do INE (correlação 0,97). Compare valores da aplicação entre si, não com os do INE.",
    "",
    "Mudanças de série",
    paste(unique(PLANNING_SERIES_BREAKS$note), collapse = " "),
    "População: série revista do INE (0012918) a partir de 2021; degrau de cerca de 1,7% em 2020/2021.",
    "",
    "Folhas",
    if (!is.null(focus)) "Resumo: todos os indicadores no último ano disponível, para o local e os comparadores.",
    "Uma folha por indicador: áreas em linhas, anos (ou triénios) em colunas, com os limites do intervalo de confiança de 95% em dois blocos por baixo, quando o indicador tem intervalo.",
    if (include_long) "Dados: formato longo, com intervalos de confiança de 95%, numerador e denominador.",
    "Pirâmide etária: população por grupo etário e sexo.",
    "Mortalidade proporcional: óbitos por grande grupo de causas, por triénio, para todas as idades [I45] e para as idades abaixo de 75 [I46].",
    "",
    "Fontes: INE (população 0003182/0008273/0012918; óbitos 0008206/0013166; nados-vivos 0000003/0008084/0012434; e os indicadores listados no manual da aplicação)."
  )
  openxlsx::writeData(wb, "Leia-me", data.frame(x = readme), colNames = FALSE)
  openxlsx::addStyle(wb, "Leia-me", title_style, rows = 1, cols = 1)
  openxlsx::addStyle(wb, "Leia-me", openxlsx::createStyle(textDecoration = "bold"),
                     rows = which(readme %in% c("Como são calculados", "Marcas (valores a cinzento e itálico)", "Mudanças de série", "Folhas", "Esperança de vida à nascença", "Notas assinaladas com * junto ao nome do indicador", "Lacunas nos dados do INE")), cols = 1, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Leia-me", cols = 1, widths = 140)

  # --- Summary (selection) ----------------------------------------------------
  if (!is.null(focus)) {
    openxlsx::addWorksheet(wb, "Resumo")
    summary <- planning_latest_summary(table, areas, education_min_age)
    openxlsx::writeData(wb, "Resumo", summary, headerStyle = header_style)
    openxlsx::setColWidths(wb, "Resumo", cols = seq_len(ncol(summary)), widths = c(16, 55, 8, 22, 8, 14, rep(18, ncol(summary) - 6)))
    openxlsx::freezePane(wb, "Resumo", firstRow = TRUE)
  }

  # --- One sheet per indicator -------------------------------------------------
  n <- nrow(PLANNING_INDICATORS)
  for (i in seq_len(n)) {
    spec <- as.list(PLANNING_INDICATORS[i, ])
    progress(0.1 + 0.6 * i / n, spec$label)
    sheet <- PLANNING_SHEET_NAMES[[spec$id]]
    rows <- table[table$indicator == spec$id, , drop = FALSE]
    available <- sort(unique(rows$year[!is.na(rows$value)]))
    if (length(available) == 0) next

    periods <- vapply(available, function(y) planning_period_label(spec$id, y), character(1))
    wide_values <- tidyr::pivot_wider(
      rows[rows$year %in% available, c("level", "area", "year", "value")],
      names_from = year, values_from = value
    )
    wide_flags <- tidyr::pivot_wider(
      rows[rows$year %in% available, c("area", "year", "flag")],
      names_from = year, values_from = flag, values_fill = ""
    )
    order <- match(wide_values$area, areas$area)
    wide_values <- wide_values[order(order), c("level", "area", as.character(available)), drop = FALSE]
    wide_flags <- wide_flags[match(wide_values$area, wide_flags$area), as.character(available), drop = FALSE]
    names(wide_values) <- c("Nível", "Local", periods)

    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, planning_indicator_label(spec$id, education_min_age), startRow = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1)
    notes <- paste0(
      "Unidade: ", spec$unit,
      if (spec$window > 1) ". Triénios identificados pelos três anos." else ". Anual.",
      if (!isTRUE(spec$comparable)) " Contagem absoluta: depende do tamanho da área." else "",
      if (any(nzchar(unlist(wide_flags)))) " Valores a cinzento e itálico têm uma marca (ver Leia-me)." else ""
    )
    breaks <- c(PLANNING_SERIES_BREAKS$note[PLANNING_SERIES_BREAKS$indicator == spec$id],
                if (spec$id %in% names(PLANNING_INDICATOR_NOTES)) paste("*", PLANNING_INDICATOR_NOTES[[spec$id]]))
    openxlsx::writeData(wb, sheet, paste(c(notes, breaks), collapse = " "), startRow = 2)
    openxlsx::addStyle(wb, sheet, note_style, rows = 2, cols = 1)
    openxlsx::writeData(wb, sheet, wide_values, startRow = 4, headerStyle = header_style)
    openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = planning_number_format(spec$digits)),
                       rows = 4 + seq_len(nrow(wide_values)), cols = 2 + seq_along(available), gridExpand = TRUE)
    flagged <- which(as.matrix(wide_flags) != "", arr.ind = TRUE)
    if (nrow(flagged) > 0) {
      openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = planning_number_format(spec$digits), fontColour = "#898781", textDecoration = "italic"),
                         rows = 4 + flagged[, 1], cols = 2 + flagged[, 2], gridExpand = FALSE, stack = FALSE)
    }
    # Below the values, the two bounds of the 95% interval, so every sheet
    # carries its uncertainty even in the file that has no long data sheet.
    bounds <- list(lower = "Limite inferior do intervalo de confiança de 95%",
                   upper = "Limite superior do intervalo de confiança de 95%")
    next_row <- 5 + nrow(wide_values)
    for (bound in names(bounds)) {
      if (all(is.na(rows[[bound]]))) next
      wide_bound <- tidyr::pivot_wider(
        rows[rows$year %in% available, c("level", "area", "year", bound)],
        names_from = year, values_from = dplyr::all_of(bound)
      )
      wide_bound <- wide_bound[match(wide_values$Local, wide_bound$area), c("level", "area", as.character(available)), drop = FALSE]
      names(wide_bound) <- c("Nível", "Local", periods)
      openxlsx::writeData(wb, sheet, bounds[[bound]], startRow = next_row + 1)
      openxlsx::addStyle(wb, sheet, note_style, rows = next_row + 1, cols = 1)
      openxlsx::writeData(wb, sheet, wide_bound, startRow = next_row + 2, headerStyle = header_style)
      openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = planning_number_format(spec$digits)),
                         rows = next_row + 2 + seq_len(nrow(wide_bound)), cols = 2 + seq_along(available), gridExpand = TRUE)
      next_row <- next_row + 2 + nrow(wide_bound)
    }

    openxlsx::setColWidths(wb, sheet, cols = 1:2, widths = c(11, 42))
    openxlsx::setColWidths(wb, sheet, cols = 2 + seq_along(available), widths = if (spec$window > 1) 11 else 9)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 3)
  }

  # --- Long data ------------------------------------------------------------------
  if (include_long) {
    progress(0.75, "dados")
    long <- table %>%
      dplyr::filter(!is.na(.data$value)) %>%
      dplyr::left_join(PLANNING_INDICATORS[, c("id", "label", "ref", "unit")], by = c("indicator" = "id")) %>%
      dplyr::transmute(
        `Nível` = .data$level, Local = .data$area, Indicador = .data$label, `Referência` = .data$ref,
        Ano = .data$year,
        `Período` = vapply(seq_along(.data$year), function(k) planning_period_label(.data$indicator[[k]], .data$year[[k]]), character(1)),
        Unidade = .data$unit, Valor = .data$value, `IC 95% inferior` = .data$lower, `IC 95% superior` = .data$upper,
        Numerador = .data$numerator, Denominador = .data$denominator, Marca = .data$flag,
        `Face a Portugal` = .data$significance
      )
    openxlsx::addWorksheet(wb, "Dados")
    openxlsx::writeData(wb, "Dados", long, headerStyle = header_style)
    openxlsx::freezePane(wb, "Dados", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Dados", cols = seq_len(ncol(long)), widths = c(11, 36, 50, 10, 6, 11, 22, 12, 14, 14, 12, 12, 7, 14))
  }

  # --- Pyramid ---------------------------------------------------------------------
  progress(0.85, "pirâmide etária")
  pyramid_years <- planning_pyramid_years(snapshot_years_for("population"))
  pyramid <- dplyr::bind_rows(lapply(pyramid_years, function(year) {
    dplyr::bind_rows(lapply(areas$area, planning_pyramid, year = year, lookup = lookup))
  }))
  if (nrow(pyramid) > 0) {
    pyramid <- pyramid %>%
      dplyr::left_join(areas, by = "area") %>%
      dplyr::transmute(
        `Nível` = .data$level, Local = .data$area, Ano = .data$year, `Grupo etário` = .data$age_band,
        Sexo = ifelse(.data$sex == "H", "Homens", "Mulheres"), `População` = .data$pop,
        `% da população do local` = round(.data$share, 3)
      )
    openxlsx::addWorksheet(wb, "Pirâmide etária")
    openxlsx::writeData(wb, "Pirâmide etária", pyramid, headerStyle = header_style)
    openxlsx::freezePane(wb, "Pirâmide etária", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Pirâmide etária", cols = 1:7, widths = c(11, 42, 6, 16, 10, 12, 20))
  }

  # --- Proportional mortality -----------------------------------------------------
  progress(0.92, "mortalidade proporcional")
  proportional_years <- planning_proportional_years()
  proportional_years <- proportional_years[proportional_years <= max(years) & proportional_years >= min(years)]
  if (length(proportional_years) > 0) {
    proportional <- planning_proportional_table(areas$area, proportional_years, lookup = lookup) %>%
      dplyr::left_join(areas, by = "area") %>%
      dplyr::transmute(
        `Nível` = .data$level, Local = .data$area, `Triénio` = .data$period, `Código` = .data$code,
        `Grupo de causas` = .data$group, `Óbitos` = .data$deaths, `%` = .data$share,
        `IC 95% inferior` = .data$lower, `IC 95% superior` = .data$upper
      )
    openxlsx::addWorksheet(wb, "Mortalidade proporcional")
    openxlsx::writeData(wb, "Mortalidade proporcional", proportional, headerStyle = header_style)
    openxlsx::freezePane(wb, "Mortalidade proporcional", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Mortalidade proporcional", cols = 1:9, widths = c(11, 42, 11, 8, 50, 10, 8, 14, 14))
  }

  # --- Proportional mortality under 75 ---------------------------------------------
  under75_years <- planning_under75_years()
  under75_years <- under75_years[under75_years <= max(years) & under75_years >= min(years)]
  if (length(under75_years) > 0) {
    progress(0.95, "mortalidade proporcional < 75 anos")
    under75 <- planning_under75_table(areas$area, under75_years, lookup = lookup, vintage = vintage) %>%
      dplyr::left_join(areas, by = "area") %>%
      dplyr::transmute(
        `Nível` = .data$level, Local = .data$area, `Triénio` = .data$period, `Código` = .data$code,
        `Grupo de causas` = .data$group, `Óbitos com menos de 75 anos` = .data$deaths, `%` = .data$share,
        `IC 95% inferior` = .data$lower, `IC 95% superior` = .data$upper, Marca = .data$flag
      )
    openxlsx::addWorksheet(wb, "Mortalidade proporcional <75")
    openxlsx::writeData(wb, "Mortalidade proporcional <75", under75, headerStyle = header_style)
    openxlsx::freezePane(wb, "Mortalidade proporcional <75", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Mortalidade proporcional <75", cols = 1:10, widths = c(11, 42, 11, 8, 50, 22, 8, 14, 14, 7))
  }

  # --- Mortality by cause group, standardised --------------------------------------
  cause_years <- planning_standardised_years()
  cause_years <- cause_years[cause_years <= max(years) & cause_years >= min(years)]
  # Every triennium for a selection; the latest one for the file of all areas.
  if (!include_long) cause_years <- utils::tail(cause_years, 1)
  if (length(cause_years) > 0) {
    progress(0.96, "mortalidade por causa")
    causes <- dplyr::bind_rows(lapply(cause_years, function(y) {
      planning_cause_standardised(union(areas$area, benchmark), y, lookup = lookup, benchmark = benchmark)
    })) %>%
      dplyr::filter(.data$area %in% areas$area) %>%
      dplyr::left_join(areas, by = "area") %>%
      dplyr::transmute(
        `Nível` = .data$level, Local = .data$area, `Triénio` = .data$period, `Código` = .data$code,
        `Grupo de causas` = .data$group, `Óbitos` = round(.data$observed, 1), `Óbitos esperados` = round(.data$expected, 1),
        SMR = round(.data$smr, 1), `SMR IC 95% inferior` = round(.data$smr_lower, 1), `SMR IC 95% superior` = round(.data$smr_upper, 1),
        `Face a Portugal` = .data$significance,
        `Taxa padronizada` = round(.data$dsr, 1), `Taxa padronizada IC inferior` = round(.data$dsr_lower, 1), `Taxa padronizada IC superior` = round(.data$dsr_upper, 1),
        `Taxa padronizada < 75` = round(.data$dsr75, 1), `< 75 IC inferior` = round(.data$dsr75_lower, 1), `< 75 IC superior` = round(.data$dsr75_upper, 1),
        Marca = .data$flag
      )
    openxlsx::addWorksheet(wb, "Mortalidade por causa (SMR)")
    openxlsx::writeData(wb, "Mortalidade por causa (SMR)", causes, headerStyle = header_style)
    openxlsx::freezePane(wb, "Mortalidade por causa (SMR)", firstRow = TRUE)
    openxlsx::setColWidths(wb, "Mortalidade por causa (SMR)", cols = seq_len(ncol(causes)), widths = c(11, 42, 11, 8, 50, rep(12, ncol(causes) - 5)))
  }

  progress(0.97, "a gravar")
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

# The latest year each indicator has for the focus area, with every area's value
# in that year side by side. Absolute counts are shown for the local area only.
planning_latest_summary <- function(table, areas, education_min_age = 0L) {
  focus <- areas$area[[1]]
  rows <- lapply(PLANNING_INDICATORS$id, function(id) {
    spec <- planning_indicator_spec(id)
    own <- table[table$indicator == id & table$area == focus & !is.na(table$value), , drop = FALSE]
    if (nrow(own) == 0) return(NULL)
    year <- max(own$year)
    values <- table[table$indicator == id & table$year == year, , drop = FALSE]
    cells <- stats::setNames(lapply(areas$area, function(a) {
      if (!identical(a, focus) && !isTRUE(spec$comparable)) return(NA_real_)
      v <- values$value[values$area == a]
      if (length(v) == 0) NA_real_ else round(v[[1]], spec$digits)
    }), paste0(areas$area, " (", areas$level, ")"))
    tibble::as_tibble(c(
      list(Tema = spec$theme, Indicador = planning_indicator_label(id, education_min_age),
           Unidade = spec$unit, `Período` = planning_period_label(id, year),
           Marca = values$flag[values$area == focus][[1]],
           `Face a Portugal` = if ("significance" %in% names(values)) values$significance[values$area == focus][[1]] else NA_character_),
      cells
    ))
  })
  dplyr::bind_rows(rows)
}

# Every area with its level, in the workbook's order, for the full export.
planning_full_export_areas <- function(lookup = get_nuts_lookup(), health = get_health_lookup()) {
  levels <- planning_area_levels(lookup, health)
  levels %>%
    # Portugal both as INE's total and as the sum of its municipalities.
    dplyr::bind_rows(tibble::tibble(area = PLANNING_PORTUGAL_MUNICIPAL, level = "Portugal")) %>%
    dplyr::mutate(order = match(.data$level, c("Portugal", "NUTS I", "NUTS II", "NUTS III", "ARS", "ULS", "Município"))) %>%
    dplyr::arrange(.data$order) %>%
    dplyr::select(area, level)
}
