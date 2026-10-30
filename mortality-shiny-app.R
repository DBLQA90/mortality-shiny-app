# =========================================================
# App bootstrap
# =========================================================

get_app_dir <- function() {
  ofiles <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]

  if (length(ofiles) > 0) {
    return(dirname(normalizePath(utils::tail(ofiles, 1), mustWork = FALSE)))
  }

  getwd()
}

app_dir <- get_app_dir()
bootstrap_env <- environment(get_app_dir)

sys.source(file.path(app_dir, "R/dependencies.R"), envir = bootstrap_env)
install_project_dependencies()

for (app_file in c(
  "R/config.R",
  "R/helpers.R",
  "R/cache.R",
  "R/snapshots.R",
  "R/ine_client.R",
  # regions before metadata: the area vocabulary is derived from the NUTS
  # lookup, so the lookup readers have to exist first.
  "R/regions.R",
  # infant before metadata too: the selectable year range includes years only
  # the infant datasets reach, so metadata has to be able to ask them.
  "R/infant.R",
  "R/regional_rows.R",
  "R/metadata.R",
  "R/metrics.R",
  "R/standardisation.R",
  "R/avoidable.R",
  "R/planning_parish.R",
  "R/planning_indicators.R",
  "R/life_expectancy.R",
  "R/planning_under75.R",
  "R/planning_standardised.R",
  "R/planning_sns.R",
  "R/planning_weekly.R",
  "R/planning_export.R",
  "R/planning_charts.R",
  "R/planning_profile.R",
  "R/planning_maps.R",
  "R/cause_maps.R",
  "R/data_versions.R",
  "R/forecast_helpers.R",
  "R/trend.R",
  "R/data_access.R",
  "R/ui_helpers.R"
)) {
  sys.source(file.path(app_dir, app_file), envir = bootstrap_env)
}

# =========================================================
# UI
# =========================================================

# The empty world topojson the plotly maps ask for (R/cause_maps.R).
shiny::addResourcePath("map-topojson", file.path(app_dir, "data", "boundaries", "topojson"))

ui <- navbarPage(
  title = "PNS Monitorização",
  id = "main_nav",
  header = nuts_vintage_control(),

  intro_tab_ui(),
  observed_mortality_tab_ui(),
  beginner_forecasting_tab_ui(),
  advanced_forecasting_tab_ui(),
  annual_metrics_tab_ui(),
  avoidable_tab_ui(),
  maps_tab_ui(),
  planning_tab_ui(),
  data_availability_tab_ui(),
  glossary_tab_ui()
)

# =========================================================
# Server
# =========================================================

server <- function(input, output, session) {
  cancel_seq <- reactiveValues(rates = 0L, beginner = 0L, forecast = 0L, annual = 0L, avoidable = 0L, maps = 0L)

  # -------------------------
  # NUTS vintage
  # -------------------------
  # One app-wide choice of which region definition groups the municipalities.
  # Every region expansion reads the lookup through these, so the whole app
  # switches together and no tab can be left grouping by the other vintage.
  value_or_default_input <- function(x, default) {
    if (is.null(x) || length(x) == 0 || !nzchar(x[[1]])) default else x[[1]]
  }

  active_nuts_vintage <- reactive({
    chosen <- input$nuts_vintage
    normalize_nuts_vintage(if (is.null(chosen) || !nzchar(chosen)) default_nuts_vintage() else chosen)
  })

  active_nuts_lookup <- reactive({
    get_nuts_lookup(active_nuts_vintage())
  })

  # How a region's deaths are built: INE's regional rows where they exist, or
  # always the sum of municipalities. App-wide, like the vintage, because it
  # changes what a regional figure means. See R/regional_rows.R.
  active_region_source <- reactive({
    normalize_region_source(value_or_default_input(input$region_source, default_region_source()))
  })

  # The denominator of every rate per inhabitant: the mean of two consecutive
  # end-of-year estimates (INE) or the end-of-year estimate alone. App-wide,
  # like the two above, so no tab divides by a different population.
  active_population_basis <- reactive({
    # Every computation that reads the denominator also reads the data mode, so
    # the planning tab recomputes when either changes.
    active_data_mode()
    normalize_population_basis(value_or_default_input(input$population_basis, default_population_basis()))
  })

  # Corrected data or INE's as published (R/helpers.R, DADOS_INE_LIMITACOES.md).
  # The mode is a process-wide option the data functions read, so it is set
  # here, as the first thing any computation that depends on it does.
  active_data_mode <- reactive({
    mode <- normalize_data_mode(value_or_default_input(input$data_mode, default_data_mode()))
    options(mortality.data_mode = mode)
    mode
  })

  observeEvent(input$data_mode, {
    active_data_mode()
    showNotification(
      paste0(
        "Dados: ", names(data_mode_choices)[match(active_data_mode(), data_mode_choices)], ". ",
        if (identical(active_data_mode(), "original")) {
          paste("Sem as correções da aplicação: os óbitos sem idade não são distribuídos, os totais confidenciais e o emprego",
                "ocultado ficam em falta, e o último ano sem causas de morte não é reconstruído. ")
        } else "",
        "O separador de planeamento atualiza-se; nos outros separadores, volte a carregar os dados para recalcular."
      ),
      type = "message", duration = 15
    )
  }, ignoreInit = TRUE)

  # The planning tab recomputes on its own; the tabs that load on a button keep
  # the rates they computed until the user loads again.
  observeEvent(input$population_basis, {
    showNotification(
      paste0(
        "População no denominador: ", population_basis_label(active_population_basis()), ". ",
        "O separador de planeamento atualiza-se; nos separadores de mortalidade observada, previsões, ",
        "métricas anuais e mortalidade evitável, volte a carregar os dados para recalcular."
      ),
      type = "message", duration = 12
    )
  }, ignoreInit = TRUE)

  # Applied immediately after every death load, so each tab gets the same
  # substitution without repeating it. Population is untouched.
  apply_region_source <- function(df, expanded_regions) {
    substitute_regional_deaths(
      df = df,
      expanded_regions = expanded_regions,
      vintage = active_nuts_vintage(),
      lookup = active_nuts_lookup(),
      source = active_region_source(),
      # Rows for the municipalities a composition subtracts, which are not
      # members of the region and so are never in the loaded frame. Read from
      # the snapshot archive, like the regional rows themselves.
      load_municipal_deaths = function(areas, years, causes) {
        tryCatch(
          get_death_data_for(area = areas, cause = causes, years = years, data_source = "snapshot"),
          error = function(e) NULL
        )
      }
    )
  }

  # Seam and municipal warnings shared by every tab that loads regional data.
  notify_region_source_warnings <- function(selected_areas, expanded_regions, causes, years) {
    # The rebuilt latest year is read for whole regions only, exactly; it is
    # neither a seam between sources nor a municipal sum with ages missing.
    years <- setdiff(as.integer(years), calendar_death_years)
    if (length(years) == 0) return(invisible())
    seam <- region_source_seam_warning(expanded_regions, active_nuts_vintage(), years, active_region_source())
    if (!is.null(seam)) notify_data(seam, type = "warning", duration = 20, id = "region-seam")

    municipal <- municipal_age_detail_warning(
      selected_areas, active_nuts_lookup(), causes, years,
      vintage = active_nuts_vintage(), source = active_region_source()
    )
    if (!is.null(municipal)) notify_data(municipal, type = "warning", duration = 20, id = "municipal-age-detail")
  }

  # Swap the region entries of every area selector when the vintage changes.
  #
  # Municipality choices are identical across vintages, so a municipal selection
  # survives untouched. A selected region may not exist in the new vintage -
  # "Grande Lisboa" has no NUTS-2013 counterpart - and those are dropped and
  # named, rather than silently left in place to be read as an unknown area.
  observeEvent(input$nuts_vintage, {
    vintage <- active_nuts_vintage()
    lookup <- active_nuts_lookup()

    if (nrow(lookup) == 0) {
      return()
    }

    area_inputs <- c("area", "area2", "beginner_area", "availability_area", "annual_area", "avoidable_area", "planning_area")
    dropped <- character(0)

    for (input_id in area_inputs) {
      reconciled <- reconcile_area_selection(
        selected = input[[input_id]],
        vintage = vintage,
        lookup = lookup,
        exclude = character(0)
      )

      dropped <- c(dropped, reconciled$dropped)
      updateSelectInput(
        session,
        input_id,
        choices = reconciled$choices,
        selected = reconciled$selected
      )
    }

    # The annual tab's comparators carry the automatic entry besides the areas.
    comparator_choices <- annual_comparator_choices(vintage, lookup)
    kept <- intersect(as.character(input$annual_comparators), comparator_choices)
    dropped <- c(dropped, setdiff(as.character(input$annual_comparators), comparator_choices))
    updateSelectizeInput(session, "annual_comparators", choices = comparator_choices, selected = kept)

    dropped <- unique(dropped)
    if (length(dropped) > 0) {
      showNotification(
        glue::glue(
          "Estas regiões não existem na definição NUTS {vintage} e foram ",
          "removidas da selecção: {paste(dropped, collapse = ', ')}."
        ),
        type = "warning",
        duration = 20
      )
    }
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_rates, {
    cancel_seq$rates <- cancel_seq$rates + 1L
    showNotification("Pedido de interrupção recebido (Taxas).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_forecast, {
    cancel_seq$forecast <- cancel_seq$forecast + 1L
    showNotification("Pedido de interrupção recebido (Projecções).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_beginner_forecast, {
    cancel_seq$beginner <- cancel_seq$beginner + 1L
    showNotification("Pedido de interrupção recebido (Previsão Guiada).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_annual_metrics, {
    cancel_seq$annual <- cancel_seq$annual + 1L
    showNotification("Pedido de interrupção recebido (Métricas anuais).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_avoidable, {
    cancel_seq$avoidable <- cancel_seq$avoidable + 1L
    showNotification("Pedido de interrupção recebido (Mortalidade evitável).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  # -------------------------
  # Warnings inside the tab
  # -------------------------
  # Warnings about the data an analysis read go to the top of the tab that is
  # computing it - the one on screen, as the outputs of hidden tabs do not run -
  # instead of floating over its charts; each new run of the tab clears them.
  # Elsewhere (no tab id in a test server, or a tab without a slot) they stay
  # notifications. Replies to the header controls and to a cancel stay
  # notifications too: they answer an action, not describe the data.
  TAB_NOTES <- c("Mortalidade Observada" = "observed", "Previsão Guiada" = "beginner", "Previsão Avançada" = "forecast",
                 "Métricas Anuais" = "annual", "Mortalidade Evitável" = "avoidable", "Mapas" = "maps")
  tab_notes <- do.call(reactiveValues, lapply(stats::setNames(nm = unname(TAB_NOTES)), function(tab) list()))

  notify_data <- function(ui, type = "message", duration = 5, id = NULL, ...) {
    tab <- unname(TAB_NOTES[isolate(input$main_nav) %||% ""])
    if (length(tab) != 1 || is.na(tab)) return(showNotification(ui, type = type, duration = duration, id = id, ...))
    text <- paste(as.character(ui), collapse = " ")
    notes <- isolate(tab_notes[[tab]])
    notes[[id %||% text]] <- list(text = text, type = type)
    tab_notes[[tab]] <- notes
    invisible(NULL)
  }

  clear_tab_notes <- function(tab) tab_notes[[tab]] <- list()
  for (pair in list(c("go_rates", "observed"), c("go_beginner_forecast", "beginner"), c("go_forecast", "forecast"),
                    c("go_annual_metrics", "annual"), c("go_avoidable", "avoidable"))) local({
    button <- pair[[1]]; tab <- pair[[2]]
    observeEvent(input[[button]], clear_tab_notes(tab), priority = 100, ignoreInit = TRUE)
  })

  render_tab_notes <- function(tab) renderUI({
    notes <- tab_notes[[tab]]
    if (length(notes) == 0) return(NULL)
    types <- vapply(notes, function(n) n$type, character(1))
    box <- function(items, class, heading) {
      if (length(items) == 0) return(NULL)
      div(
        class = paste("alert alert-dismissible tab-notes", class), role = "alert",
        tags$button(type = "button", class = "close", `data-dismiss` = "alert", `aria-label` = "Fechar", HTML("&times;")),
        tags$strong(heading),
        tags$ul(lapply(items, function(n) tags$li(n$text)))
      )
    }
    tagList(
      box(notes[types %in% c("warning", "error")], "alert-warning", "Avisos sobre estes dados"),
      box(notes[!types %in% c("warning", "error")], "alert-info", "Notas")
    )
  })
  output$observedNotes <- render_tab_notes("observed")
  output$beginnerNotes <- render_tab_notes("beginner")
  output$forecastNotes <- render_tab_notes("forecast")
  output$annualNotes <- render_tab_notes("annual")
  output$avoidableNotes <- render_tab_notes("avoidable")
  output$mapsNotes <- render_tab_notes("maps")

  # Downloads stay locked until their analysis has been run: before that they
  # wrote an error page instead of a file.
  DOWNLOAD_LOCKS <- list(
    go_rates = c("downloadRateSummaryCSV", "downloadRatePlot", "downloadRateTableCSV"),
    go_beginner_forecast = c("downloadBeginnerMultiCSV", "downloadBeginnerForecastPlot", "downloadBeginnerForecastCSV"),
    go_forecast = c("downloadForecastSpecCSV", "downloadForecastSummaryCSV", "downloadForecastPlot", "downloadForecastCSV",
                    "downloadDiagnosticResidualPlot", "downloadDiagnosticAcfPlot", "downloadDiagnosticPacfPlot",
                    "downloadDiagnosticLjungCSV", "downloadDiagnosticSummaryTXT", "downloadComparisonRankingCSV",
                    "downloadAccuracyCSV", "downloadComparisonPlot", "downloadBreakPlot", "downloadBreakTableCSV"),
    go_annual_metrics = c("downloadAnnualMetricsCSV", "downloadAnnualMetricsPlot", "downloadAnnualSourcesCSV"),
    go_avoidable = c("downloadAvoidableCSV", "downloadAvoidablePlot", "downloadAvoidableCausesCSV")
  )
  for (button in names(DOWNLOAD_LOCKS)) local({
    trigger <- button
    observe({
      session$sendCustomMessage("lock-downloads", list(ids = DOWNLOAD_LOCKS[[trigger]], locked = !isTRUE((input[[trigger]] %||% 0) > 0)))
    })
  })

  # The horizon sliders say the year they reach: the last year with data for
  # the selection, within the fit window, plus the horizon. The rebuilt latest
  # year serves only all causes, both sexes and whole regions.
  selection_last_year <- function(areas, cause, sex, fit_end) {
    active_data_mode()
    cause <- as.character(cause %||% "")[1]
    last <- if (is_detailed_tumour_cause(cause)) max(DETAILED_TUMOUR_YEARS) else max(planning_cause_years())
    calendar <- active_calendar_years()
    if (length(calendar) > 0 && length(areas) > 0 && isTRUE(calendar_selection_supported(areas, cause, sex %||% "HM"))) {
      last <- max(last, calendar)
    }
    min(as.integer(fit_end), last)
  }
  horizon_label <- function(horizon, last) paste0("Horizonte de projecção (anos): até ", last + as.integer(horizon))
  observe({
    req(input$beginner_horizon, input$beginner_years_fit)
    last <- selection_last_year(input$beginner_area, input$beginner_cause, input$beginner_sex, input$beginner_years_fit[[2]])
    updateSliderInput(session, "beginner_horizon", label = horizon_label(input$beginner_horizon, last))
  })
  observe({
    req(input$horizon, input$years_fit)
    last <- selection_last_year(input$area2, input$cause2, input$sex2, input$years_fit[[2]])
    updateSliderInput(session, "horizon", label = horizon_label(input$horizon, last))
  })

  abort_if_cancelled <- function(kind, token) {
    if (!identical(cancel_seq[[kind]], token)) {
      validate(need(FALSE, "Operação interrompida. A cache foi preservada."))
    }
  }

  year_load_state <- reactiveValues(
    rates_range = NULL,
    rates_order = "desc",
    beginner_range = NULL,
    beginner_order = "desc",
    forecast_range = NULL,
    forecast_order = "desc"
  )

  clean_year_range <- function(year_range) {
    year_range <- suppressWarnings(as.integer(year_range))
    year_range <- year_range[!is.na(year_range)]

    if (length(year_range) < 2) {
      return(NULL)
    }

    range(year_range)
  }

  detect_year_load_order <- function(previous_range, current_range) {
    current_range <- clean_year_range(current_range)
    previous_range <- clean_year_range(previous_range)

    if (is.null(current_range) || is.null(previous_range)) {
      return("desc")
    }

    current_mid <- mean(current_range)
    previous_mid <- mean(previous_range)

    if (current_mid < previous_mid) {
      return("desc")
    }
    if (current_mid > previous_mid) {
      return("asc")
    }
    if (current_range[[2]] < previous_range[[2]]) {
      return("desc")
    }
    if (current_range[[1]] > previous_range[[1]]) {
      return("asc")
    }

    "desc"
  }

  observeEvent(input$years_import, {
    current_range <- clean_year_range(input$years_import)
    year_load_state$rates_order <- detect_year_load_order(year_load_state$rates_range, current_range)
    year_load_state$rates_range <- current_range
  }, ignoreInit = FALSE, ignoreNULL = TRUE)

  observeEvent(input$years_fit, {
    current_range <- clean_year_range(input$years_fit)
    year_load_state$forecast_order <- detect_year_load_order(year_load_state$forecast_range, current_range)
    year_load_state$forecast_range <- current_range
  }, ignoreInit = FALSE, ignoreNULL = TRUE)

  observeEvent(input$beginner_years_fit, {
    current_range <- clean_year_range(input$beginner_years_fit)
    year_load_state$beginner_order <- detect_year_load_order(year_load_state$beginner_range, current_range)
    year_load_state$beginner_range <- current_range
  }, ignoreInit = FALSE, ignoreNULL = TRUE)

  get_rate_mapping <- function(rate_type, population = "Total") {
    if (identical(rate_type, "crude")) {
      list(
        value_col = "crude_rate",
        lower_col = "crude_lower",
        upper_col = "crude_upper",
        y_label = "Taxa Bruta por 100.000",
        rate_label = "Bruta"
      )
    } else {
      # The under-75 scope is standardised to the ESP-2013 0-74 sub-population
      # (conventional premature mortality), a different standard from the
      # all-age rate, so the label names the standard used to avoid implying
      # the two standardised rates share one basis.
      under75 <- identical(population, "Menos de 75 anos")
      list(
        value_col = "dsr",
        lower_col = "dsr_lower",
        upper_col = "dsr_upper",
        y_label = if (under75) {
          "Taxa Padronizada por 100.000 (padrão ESP 0-74)"
        } else {
          "Taxa Padronizada por 100.000"
        },
        rate_label = if (under75) "Padronizada (padrão ESP 0-74)" else "Padronizada"
      )
    }
  }

  # Every CSV and PNG this session writes says which data version produced it;
  # the helpers themselves are in R/helpers.R.
  current_export_stamp <- function() {
    date <- tryCatch(latest_import_date(read_import_log(app_data_root(get_snapshot_dir()))), error = function(e) NA_character_)
    vintage <- tryCatch(active_nuts_vintage(), error = function(e) NULL)
    export_stamp(date, vintage)
  }

  write_csv_utf8 <- function(x, file, stamp = current_export_stamp()) {
    helpers_write_csv_utf8(x, file, stamp = stamp)
  }

  save_ggplot_png <- function(file, plot_obj, width = 1200, height = 800, res = 150, stamp = current_export_stamp()) {
    helpers_save_ggplot_png(file, plot_obj, width = width, height = height, res = res, stamp = stamp)
  }

  save_base_plot_png <- function(file, plot_expr, width = 1200, height = 800, res = 150, stamp = current_export_stamp()) {
    helpers_save_base_plot_png(file, plot_expr, width = width, height = height, res = res, stamp = stamp)
  }

  # Shared historical-series pipeline:
  # 1) Build a stable query spec with the inputs that determine which raw data
  #    need to be downloaded and aggregated.
  # 2) Load one metric bundle (both population scopes, both rate outputs) for
  #    that query.
  # 3) Apply the final population/rate/year filter to get a reusable series.
  # A ULS that shares a municipality is read by whole municipalities outside
  # the planning tab's parish weights; say what that means for each one.
  notify_uls_split <- function(areas) {
    note <- uls_split_note(areas, "whole")
    if (length(note) > 0) notify_data(paste(note, collapse = "\n"), type = "message", duration = 25, id = "uls-split")
  }

  make_query_spec <- function(area, area_label, cause, sex, data_source = "ine") {
    area_key <- sort(unique(area))

    # Selected areas are summed into one geography, so overlapping choices
    # (Portugal with anything, a region with its own municipalities) would be
    # counted twice. Warn rather than block: a user may knowingly want the
    # combination, but should not get it by accident.
    # A region belonging to the other vintage cannot be resolved and would
    # otherwise reach the loader as an unknown place name.
    validate(need(
      is.null(vintage_mismatch_message(area_key, active_nuts_vintage(), active_nuts_lookup())),
      vintage_mismatch_message(area_key, active_nuts_vintage(), active_nuts_lookup())
    ))

    overlap <- overlapping_selection_warning(area_key, lookup = active_nuts_lookup())
    if (!is.null(overlap)) {
      notify_data(overlap, type = "warning", duration = 20)
    }
    notify_uls_split(area_key)

    # Expand any selected region into its municipalities. This applies to every
    # tab that goes through a query spec - observed mortality and both forecast
    # tabs - which previously used INE's regional rows directly and so carried
    # the NUTS-2024 discontinuity into trends and projections.
    resolved <- resolve_region_areas(
      areas = area_key,
      region_mode = default_region_mode(),
      lookup = active_nuts_lookup(),
      available_areas = get_available_areas(data_source)
    )
    for (message_text in resolved$warnings) {
      notify_data(message_text, type = "warning", duration = 15)
    }
    # Label from what the user chose, not from the expansion: a region must stay
    # named "Alentejo", not become a list of 47 municipalities.
    selection_label <- get_selection_label(area_key, area_label)

    list(
      area_key = resolved$areas,
      selected_areas = area_key,
      area_label = selection_label,
      expanded_regions = resolved$expanded,
      cause = cause,
      sex = sex,
      data_source = normalize_data_source(data_source)
    )
  }


  get_years_in_selected_range <- function(year_range, year_order = "asc") {
    selected_bounds <- suppressWarnings(as.integer(year_range))
    selected_bounds <- selected_bounds[!is.na(selected_bounds)]

    validate(
      need(length(selected_bounds) >= 2, "Seleccione um intervalo de anos para importar.")
    )

    lower_year <- min(selected_bounds)
    upper_year <- max(selected_bounds)

    validate(
      need(lower_year <= upper_year, "O ano inicial não pode ser posterior ao ano final.")
    )

    selected_years <- year_of_interest[
      year_of_interest >= lower_year &
        year_of_interest <= upper_year
    ]

    validate(
      need(length(selected_years) > 0, "O intervalo seleccionado não contém anos disponíveis nos indicadores.")
    )

    order_years(selected_years, year_order)
  }

  format_year_selection <- function(years) {
    years <- sort(unique(as.integer(years)))

    if (length(years) == 0) {
      return("Nenhum")
    }

    if (identical(years, seq.int(min(years), max(years)))) {
      return(glue::glue("{min(years)} - {max(years)}"))
    }

    paste(years, collapse = ", ")
  }

  snapshot_inventory_data <- reactive({
    get_snapshot_inventory()
  })

  output$snapshotInventorySummary <- renderTable({
    build_snapshot_inventory_summary(snapshot_inventory_data())
  })

  output$downloadSnapshotInventorySummaryCSV <- downloadHandler(
    filename = function() paste0("resumo_disponibilidade_rds_", Sys.Date(), ".csv"),
    content = function(file) {
      write_csv_utf8(build_snapshot_inventory_summary(snapshot_inventory_data()), file)
    }
  )

  snapshot_availability_data <- reactive({
    validate(need(length(input$availability_area) > 0, "Seleccione pelo menos um local."))

    years <- get_years_in_selected_range(input$availability_years)
    causes <- if (identical(input$availability_dataset, "deaths")) {
      input$availability_cause
    } else {
      NULL
    }

    build_snapshot_availability_table(
      dataset = input$availability_dataset,
      years = years,
      areas = input$availability_area,
      causes = causes,
      show_missing = isTRUE(input$availability_show_missing),
      inventory = snapshot_inventory_data()
    )
  })

  output$snapshotAvailabilityTable <- renderTable({
    snapshot_availability_data()
  })

  # -------------------------
  # Data history
  # -------------------------
  # Read once per session: the log changes only when data is refreshed, which
  # happens outside the app.
  import_log_data <- reactive({
    tryCatch(read_import_log(app_data_root(get_snapshot_dir())), error = function(e) read_import_log(tempdir()))
  })

  output$dataImportDatasets <- renderTable({
    summary <- dataset_import_summary(import_log_data())
    validate(need(nrow(summary) > 0, "Sem registo de importações (data/import_log.csv)."))
    summary %>%
      dplyr::transmute(
        `Conjunto de dados` = dataset_label(.data$dataset),
        Anos = .data$years,
        Ficheiros = .data$files,
        `Primeira importação` = .data$first_import,
        `Última importação` = .data$last_import,
        `Última revisão de valores` = dplyr::coalesce(.data$last_change, "—")
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  import_runs_data <- reactive({
    import_runs_summary(import_log_data())
  })

  observe({
    runs <- import_runs_data()
    if (nrow(runs) == 0) return()
    with_changes <- runs[runs$replaced > 0, , drop = FALSE]
    choices <- stats::setNames(runs$run_id, paste0(runs$date, " - ", runs$added, " novos, ", runs$replaced, " revistos"))
    updateSelectInput(
      session, "data_import_run",
      choices = choices,
      selected = if (nrow(with_changes) > 0) with_changes$run_id[[1]] else runs$run_id[[1]]
    )
  })

  output$dataImportRuns <- renderTable({
    runs <- import_runs_data()
    validate(need(nrow(runs) > 0, "Sem registo de importações."))
    runs %>%
      dplyr::transmute(
        Data = .data$date,
        `Ficheiros novos` = .data$added,
        `Ficheiros revistos` = .data$replaced,
        `Conjuntos de dados` = vapply(strsplit(.data$datasets, ", ", fixed = TRUE), function(x) paste(dataset_label(x), collapse = "; "), character(1)),
        Nota = .data$note,
        Identificador = .data$run_id
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  import_run_changes <- reactive({
    req(input$data_import_run)
    run_changes(import_log_data(), input$data_import_run)
  })

  output$dataImportChanges <- renderTable({
    changes <- import_run_changes()
    validate(need(nrow(changes) > 0, "Esta importação não reviu valores já existentes: só acrescentou ficheiros, ou nenhum mudou."))
    changes %>%
      dplyr::transmute(
        `Conjunto de dados` = dataset_label(.data$dataset),
        Ano = .data$year,
        Ficheiros = .data$files,
        `Linhas com valor alterado` = format(.data$rows_changed, big.mark = ".", decimal.mark = ",", scientific = FALSE),
        `Total Portugal anterior` = format(.data$total_old, big.mark = ".", decimal.mark = ",", scientific = FALSE),
        `Total Portugal novo` = format(.data$total_new, big.mark = ".", decimal.mark = ",", scientific = FALSE),
        `Variação` = ifelse(is.na(.data$change_pct), "", sprintf("%+.2f%%", .data$change_pct))
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$downloadImportLogCSV <- downloadHandler(
    filename = function() paste0("registo_importacoes_", Sys.Date(), ".csv"),
    content = function(file) write_csv_utf8(import_log_data(), file)
  )

  output$downloadSnapshotAvailabilityCSV <- downloadHandler(
    filename = function() paste0("cobertura_rds_", Sys.Date(), ".csv"),
    content = function(file) {
      write_csv_utf8(snapshot_availability_data(), file)
    }
  )

  make_series_spec <- function(query_spec, population, rate_type, year_range = range(year_of_interest)) {
    selected_years <- get_years_in_selected_range(year_range)

    # Observed rates and forecasts are always rate series, so any window
    # crossing the population revision inherits its step.
    revision_warning <- population_revision_warning(selected_years, "crude")
    if (!is.null(revision_warning)) {
      notify_data(revision_warning, type = "warning", duration = 20, id = "population-revision")
    }

    c(
      query_spec,
      list(
        population = population,
        rate_type = rate_type,
        years = selected_years,
        year_range = range(selected_years)
      )
    )
  }

  notify_snapshot_request_warnings <- function(
    years,
    areas,
    causes,
    data_source,
    include_population = TRUE,
    include_deaths = TRUE
  ) {
    if (!identical(normalize_data_source(data_source), "snapshot")) {
      return(invisible(character(0)))
    }

    messages <- build_snapshot_request_warnings(
      years = years,
      areas = areas,
      causes = causes,
      include_population = include_population,
      include_deaths = include_deaths
    )

    if (length(messages) > 0) {
      notify_data(
        paste(messages, collapse = "\n"),
        type = "warning",
        duration = 12
      )
    }

    invisible(messages)
  }

  load_metric_bundle <- function(query_spec, kind, token, year_range = range(year_of_interest), year_order = "asc") {
    validate(need(length(query_spec$area_key) > 0, "Selecione pelo menos um local de residência."))

    abort_if_cancelled(kind, token)
    incProgress(0.1)

    year_order <- normalize_year_order(year_order)
    years_to_load <- get_years_in_selected_range(year_range, year_order = year_order)

    # Every consumer of a metric bundle computes a rate, so a year without a
    # population estimate cannot contribute and would abort the whole load. The
    # rate-bearing sliders are already bounded to `population_years`; this is the
    # backstop for a bookmarked URL or an out-of-range value reaching here.
    unusable <- years_without_population(years_to_load)
    if (length(unusable) > 0) {
      years_to_load <- setdiff(years_to_load, unusable)
      validate(need(
        length(years_to_load) > 0,
        population_gap_message(unusable, "as taxas de mortalidade")
      ))
      notify_data(
        population_gap_message(unusable, "as taxas de mortalidade"),
        type = "warning",
        duration = 15
      )
    }

    # The latest year has deaths by age only for whole regions and all causes;
    # elsewhere it is left out of the series, and said so.
    calendar <- intersect(years_to_load, calendar_death_years)
    if (length(calendar) > 0) {
      if (calendar_selection_supported(query_spec$area_key, query_spec$cause, query_spec$sex)) {
        notify_data(
          paste0(paste(calendar, collapse = ", "), ": óbitos por idade reconstruídos a partir dos óbitos semanais do INE (provisório)."),
          type = "message", duration = 10, id = "calendar-year"
        )
      } else {
        years_to_load <- setdiff(years_to_load, calendar)
        validate(need(length(years_to_load) > 0, calendar_gap_message(calendar)))
        notify_data(calendar_gap_message(calendar), type = "warning", duration = 15, id = "calendar-year")
      }
    }

    notify_snapshot_request_warnings(
      years = setdiff(years_to_load, calendar_death_years),
      areas = query_spec$area_key,
      causes = query_spec$cause,
      data_source = query_spec$data_source
    )

    dat <- with_data_load_cancel_checker(
      cancel_checker = function() !identical(cancel_seq[[kind]], token),
      get_data_for(
        query_spec$area_key,
        query_spec$cause,
        years_to_load,
        year_order = year_order,
        data_source = query_spec$data_source,
        population_basis = active_population_basis()
      )
    )

    if (isTRUE(attr(dat, "cancelled"))) {
      showNotification(
        "Carregamento interrompido. As fatias já concluídas foram preservadas; os resultados usam apenas os dados disponíveis.",
        type = "warning",
        duration = 6
      )
    }
    incProgress(0.5)

    notify_region_source_warnings(
      selected_areas = query_spec$selected_areas,
      expanded_regions = query_spec$expanded_regions,
      causes = query_spec$cause,
      years = years_to_load
    )

    df_full <- apply_region_source(dat$full, query_spec$expanded_regions) %>%
      dplyr::filter(sex == query_spec$sex)
    # Derived from the substituted frame, so both population scopes use the
    # same deaths.
    df_trunc <- df_full %>%
      dplyr::filter(!age_band %in% exclude_bands)

    validate(
      need(nrow(df_full) > 0, "Não existem dados carregados para a selecção actual.")
    )

    metrics <- dplyr::bind_rows(
      compute_metrics(df_full) %>% dplyr::mutate(População = "Total"),
      compute_metrics(df_trunc) %>% dplyr::mutate(População = "Menos de 75 anos")
    ) %>%
      dplyr::arrange(year)

    abort_if_cancelled(kind, token)
    incProgress(0.4)

    list(
      query_spec = query_spec,
      metrics = metrics,
      years = years_to_load,
      source_summary = get_loaded_source_summary(df_full)
    )
  }

  build_historical_series <- function(metric_bundle, series_spec) {
    rate_map <- get_rate_mapping(series_spec$rate_type, series_spec$population)

    series <- metric_bundle$metrics %>%
      dplyr::filter(
        População == series_spec$population,
        year %in% series_spec$years
      ) %>%
      dplyr::arrange(year) %>%
      dplyr::transmute(
        year = as.integer(year),
        value = .data[[rate_map$value_col]],
        lower = .data[[rate_map$lower_col]],
        upper = .data[[rate_map$upper_col]]
      )

    list(
      spec = series_spec,
      metric_bundle = metric_bundle,
      series = series,
      area_label = series_spec$area_label,
      y_label = rate_map$y_label,
      rate_label = rate_map$rate_label
    )
  }

  build_forecast_plot <- function(dat) {
    build_advanced_forecast_plot(dat, view_mode = "compare")
  }

  build_forecast_display_table <- function(dat, view_mode = "compare", selected_model = NULL) {
    obs <- dat$obs %>%
      dplyr::mutate(year = as.integer(year))
    fc <- get_forecast_output_fc(
      dat = dat,
      view_mode = view_mode,
      selected_model = selected_model
    ) %>%
      dplyr::mutate(year = as.integer(year))

    obs_fmt <- obs %>%
      dplyr::transmute(
        Ano = year,
        Observado = glue::glue("{round(value, 2)}")
      )

    if (nrow(fc) == 0) {
      return(obs_fmt)
    }

    fc_fmt <- fc %>%
      dplyr::mutate(
        model_label = get_model_labels(model),
        texto = glue::glue(
          "{round(mean, 2)} ({round(lower, 2)}; {round(upper, 2)})"
        )
      )

    if (identical(view_mode, "single")) {
      return(
        dplyr::full_join(
          obs_fmt,
          fc_fmt %>%
            dplyr::transmute(
              Ano = year,
              `Previsão (IC)` = texto
            ),
          by = "Ano"
        ) %>%
          dplyr::arrange(Ano) %>%
          dplyr::mutate(dplyr::across(-Ano, ~ ifelse(is.na(.x), "—", as.character(.x))))
      )
    }

    fc_fmt <- fc_fmt %>%
      dplyr::select(year, model_label, texto)

    fc_wide <- fc_fmt %>%
      tidyr::pivot_wider(
        id_cols = year,
        names_from = model_label,
        values_from = texto
      )

    dplyr::full_join(obs_fmt, fc_wide, by = c("Ano" = "year")) %>%
      dplyr::arrange(Ano) %>%
      dplyr::mutate(dplyr::across(-Ano, ~ ifelse(is.na(.x), "—", as.character(.x))))
  }

  build_forecast_download_table <- function(dat, view_mode = "compare", selected_model = NULL) {
    obs <- dat$obs %>%
      dplyr::mutate(year = as.integer(year))
    fc <- get_forecast_output_fc(
      dat = dat,
      view_mode = view_mode,
      selected_model = selected_model
    ) %>%
      dplyr::mutate(year = as.integer(year))

    obs_fmt <- obs %>%
      dplyr::transmute(
        Ano = year,
        Observado = round(value, 2)
      )

    if (nrow(fc) == 0) {
      return(obs_fmt)
    }

    fc_fmt <- fc %>%
      dplyr::mutate(
        model_label = get_model_labels(model),
        texto = glue::glue(
          "{round(mean, 2)} ({round(lower, 2)}; {round(upper, 2)})"
        )
      )

    if (identical(view_mode, "single")) {
      return(
        dplyr::full_join(
          obs_fmt,
          fc_fmt %>%
            dplyr::transmute(
              Ano = year,
              Previsão = texto
            ),
          by = "Ano"
        ) %>%
          dplyr::arrange(Ano) %>%
          dplyr::mutate(dplyr::across(-Ano, ~ ifelse(is.na(.x), "—", as.character(.x))))
      )
    }

    fc_fmt <- fc_fmt %>%
      dplyr::select(year, model_label, texto)

    fc_wide <- fc_fmt %>%
      tidyr::pivot_wider(
        id_cols = year,
        names_from = model_label,
        values_from = texto
      )

    dplyr::full_join(obs_fmt, fc_wide, by = c("Ano" = "year")) %>%
      dplyr::arrange(Ano) %>%
      dplyr::mutate(dplyr::across(-Ano, ~ ifelse(is.na(.x), "—", as.character(.x))))
  }

  get_diagnostic_fit <- function(dat, model_id = NULL) {
    successful_models <- get_successful_model_ids(dat)

    validate(need(length(successful_models) > 0, "Nenhum modelo foi estimado com sucesso."))

    chosen_model <- resolve_selected_successful_model(dat, model_id)

    list(
      model_id = chosen_model,
      model_label = get_model_labels(chosen_model),
      fit = dat$fits[[chosen_model]]
    )
  }

  get_model_residual_values <- function(fit) {
    as.numeric(stats::na.omit(residuals(fit)))
  }

  build_diagnostic_residual_df <- function(fit) {
    res <- residuals(fit)

    tibble(
      time = as.numeric(time(res)),
      resid = as.numeric(res)
    ) %>%
      tidyr::drop_na()
  }

  build_diagnostic_residual_plot <- function(dat, model_id = NULL) {
    diag_fit <- get_diagnostic_fit(dat, model_id)
    resid_df <- build_diagnostic_residual_df(diag_fit$fit)

    validate(need(nrow(resid_df) > 0, "Os resíduos não estão disponíveis para o modelo seleccionado."))

    ggplot(resid_df, aes(x = time, y = resid)) +
      geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
      geom_line(color = "#0b2e4f", linewidth = 0.8) +
      geom_point(color = "#0b2e4f", size = 1.5) +
      labs(
        title = paste("Gráfico Temporal dos Resíduos -", diag_fit$model_label),
        x = "Tempo",
        y = "Resíduo"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }

  plot_diagnostic_correlation <- function(dat, model_id = NULL, partial = FALSE) {
    diag_fit <- get_diagnostic_fit(dat, model_id)
    res <- get_model_residual_values(diag_fit$fit)

    validate(need(length(res) >= 3, "São necessárias pelo menos 3 observações residuais para os diagnósticos de correlação."))

    lag_cap <- max(1, min(20, length(res) - 1))
    title_prefix <- if (partial) "PACF dos Resíduos" else "ACF dos Resíduos"

    if (partial) {
      stats::pacf(
        res,
        lag.max = lag_cap,
        main = paste(title_prefix, "-", diag_fit$model_label)
      )
    } else {
      stats::acf(
        res,
        lag.max = lag_cap,
        main = paste(title_prefix, "-", diag_fit$model_label)
      )
    }
  }

  get_model_fitdf <- function(fit, residual_count) {
    coef_count <- tryCatch(
      {
        coefs <- stats::coef(fit)
        sum(is.finite(coefs))
      },
      error = function(e) 0L
    )

    as.integer(min(coef_count, max(residual_count - 2L, 0L)))
  }

  build_ljung_box_table <- function(dat, model_id = NULL) {
    diag_fit <- get_diagnostic_fit(dat, model_id)
    res <- get_model_residual_values(diag_fit$fit)

    validate(need(length(res) >= 3, "São necessárias pelo menos 3 observações residuais para o teste de Ljung-Box."))

    fitdf <- get_model_fitdf(diag_fit$fit, length(res))
    max_lag <- max(length(res) - 1L, 1L)
    candidate_lags <- unique(pmin(c(5L, 10L, 15L), max_lag))
    candidate_lags <- candidate_lags[candidate_lags > fitdf]

    if (length(candidate_lags) == 0) {
      candidate_lag <- min(max_lag, fitdf + 1L)
      candidate_lags <- if (candidate_lag >= 1L) candidate_lag else 1L
    }

    dplyr::bind_rows(lapply(candidate_lags, function(current_lag) {
      lb <- stats::Box.test(
        res,
        lag = current_lag,
        type = "Ljung-Box",
        fitdf = fitdf
      )

      tibble(
        Modelo = diag_fit$model_label,
        Desfasamento = current_lag,
        FitDF = fitdf,
        Estatística = unname(lb$statistic),
        Valor_p = lb$p.value
      )
    }))
  }

  build_diagnostic_model_summary <- function(dat, model_id = NULL) {
    diag_fit <- get_diagnostic_fit(dat, model_id)

    summary_lines <- tryCatch(
      capture.output(summary(diag_fit$fit)),
      error = function(e) character(0)
    )

    if (length(summary_lines) == 0) {
      summary_lines <- tryCatch(
        capture.output(print(diag_fit$fit)),
        error = function(e) paste("Resumo do modelo indisponível:", e$message)
      )
    }

    paste(
      c(
        paste("Modelo:", diag_fit$model_label),
        paste("Classe:", paste(class(diag_fit$fit), collapse = ", ")),
        "",
        summary_lines
      ),
      collapse = "\n"
    )
  }

  rebuild_history_with_year_range <- function(history, year_range) {
    # An infant series is already built over every year; the earlier years
    # are a subset of it.
    if (!is.null(history$infant)) {
      year_range <- range(as.integer(year_range))
      keep <- function(x) x[x$year >= year_range[[1]] & x$year <= year_range[[2]], , drop = FALSE]
      history$series <- keep(history$series)
      history$infant$annual <- history$infant$annual[history$infant$annual$year <= year_range[[2]], , drop = FALSE]
      return(history)
    }
    build_historical_series(
      metric_bundle = history$metric_bundle,
      series_spec = make_series_spec(
        query_spec = list(
          area_key = history$spec$area_key,
          area_label = history$spec$area_label,
          cause = history$spec$cause,
          sex = history$spec$sex
        ),
        population = history$spec$population,
        rate_type = history$spec$rate_type,
        year_range = year_range
      )
    )
  }

  build_holdout_split <- function(history, holdout_k) {
    df <- history$series %>%
      dplyr::arrange(year)

    validate(need(holdout_k >= 1, "O holdout tem de ter pelo menos 1 ano."))
    validate(
      need(
        nrow(df) >= holdout_k + 3,
        "A janela de holdout seleccionada deixa demasiado poucos anos para treino. Reduza k ou alargue a janela de ajuste."
      )
    )

    train_df <- df %>%
      dplyr::slice_head(n = nrow(df) - holdout_k)
    holdout_df <- df %>%
      dplyr::slice_tail(n = holdout_k)

    list(
      training_history = rebuild_history_with_year_range(
        history = history,
        year_range = range(train_df$year)
      ),
      holdout_actual = holdout_df,
      full_observed = df
    )
  }

  build_holdout_metric_table <- function(training_result, holdout_actual) {
    successful_models <- get_successful_model_ids(training_result)

    if (length(successful_models) == 0) {
      return(
        tibble(
          Model = character(0),
          ME = numeric(0),
          RMSE = numeric(0),
          MAE = numeric(0),
          MAPE = numeric(0),
          MASE = numeric(0)
        )
      )
    }

    scale_denom <- if (nrow(training_result$obs) >= 2) {
      mean(abs(diff(training_result$obs$value)), na.rm = TRUE)
    } else {
      NA_real_
    }

    dplyr::bind_rows(lapply(successful_models, function(model_id) {
      pred_df <- training_result$fc %>%
        dplyr::filter(model == model_id) %>%
        dplyr::select(year, predicted = mean)
      eval_df <- holdout_actual %>%
        dplyr::select(year, actual = value) %>%
        dplyr::inner_join(pred_df, by = "year")

      validate(
        need(
          nrow(eval_df) > 0,
          paste0("Não foram encontrados anos de holdout sobrepostos para o modelo ", get_model_labels(model_id), ".")
        )
      )

      compute_forecast_error_metrics(
        actual = eval_df$actual,
        predicted = eval_df$predicted,
        scale_denom = scale_denom
      ) %>%
        dplyr::mutate(Model = model_id, .before = 1)
    }))
  }

  # `recommended_model`: the forecast actually shown, when it is not simply the
  # lowest error (the guided forecast's mean of six models).
  rank_model_metrics <- function(metric_tbl, recommended_model = choose_recommended_model(metric_tbl)) {
    if (nrow(metric_tbl) == 0) {
      return(tibble(Classificação = integer(0), Model = character(0), Recomendação = character(0)))
    }

    ranking_order <- intersect(c("RMSE", "MAE", "MASE", "MAPE"), names(metric_tbl))

    metric_tbl %>%
      dplyr::arrange(dplyr::across(dplyr::all_of(ranking_order)), Model) %>%
      dplyr::mutate(
        Classificação = dplyr::row_number(),
        Recomendação = dplyr::if_else(Model == recommended_model, "Seleccionado pela lógica actual", "")
      ) %>%
      dplyr::select(Classificação, Model, Recomendação)
  }

  build_comparison_plot <- function(comparison_dat) {
    if (identical(comparison_dat$mode, "holdout")) {
      validate(need(nrow(comparison_dat$forecast_df) > 0, "Não existem previsões de modelos disponíveis para a comparação holdout."))

      return(
        ggplot() +
          geom_line(
            data = comparison_dat$training_obs,
            aes(x = year, y = value),
            color = "grey75",
            linewidth = 0.8
          ) +
          geom_line(
            data = comparison_dat$holdout_actual,
            aes(x = year, y = value),
            color = "black",
            linewidth = 1
          ) +
          geom_point(
            data = comparison_dat$holdout_actual,
            aes(x = year, y = value),
            color = "black",
            size = 2
          ) +
          geom_line(
            data = comparison_dat$forecast_df,
            aes(x = year, y = mean, color = model),
            linetype = "dashed",
            linewidth = 1
          ) +
          geom_point(
            data = comparison_dat$forecast_df,
            aes(x = year, y = mean, color = model),
            size = 1.8
          ) +
          labs(
            title = paste("Comparação por Validação - Últimos", comparison_dat$holdout_k, "Anos"),
            subtitle = "A linha preta mostra os valores observados no período de validação; as linhas tracejadas mostram as previsões dos modelos.",
            x = "Ano",
            y = comparison_dat$y_label
          ) +
          scale_color_brewer(palette = "Set1", name = "Modelo", labels = get_model_labels) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }

    if (identical(comparison_dat$mode, "rolling")) {
      validate(need(nrow(comparison_dat$predictions) > 0, "Não existem previsões de validação móvel para apresentar."))

      return(
        ggplot() +
          geom_line(
            data = comparison_dat$full_obs,
            aes(x = year, y = value),
            color = "grey70",
            linewidth = 0.8
          ) +
          geom_point(
            data = comparison_dat$full_obs,
            aes(x = year, y = value),
            color = "grey40",
            size = 1.6
          ) +
          geom_line(
            data = comparison_dat$predictions,
            aes(x = year, y = predicted, color = model),
            linetype = "dashed",
            linewidth = 0.9
          ) +
          geom_point(
            data = comparison_dat$predictions,
            aes(x = year, y = predicted, color = model),
            size = 2
          ) +
          labs(
            title = "Validação móvel - previsões a um passo",
            subtitle = "Pontos cinzentos: série observada; linhas tracejadas: previsão de cada modelo a partir de cada origem.",
            x = "Ano",
            y = comparison_dat$y_label
          ) +
          scale_color_brewer(palette = "Set1", name = "Modelo", labels = get_model_labels) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }

    inverse_fn <- comparison_dat$inverse
    if (!is.function(inverse_fn)) {
      inverse_fn <- function(x) as.numeric(x)
    }
    fitted_df <- build_fitted_values_df(
      obs = comparison_dat$obs,
      fits = comparison_dat$fits,
      inverse = inverse_fn
    )

    validate(need(nrow(fitted_df) > 0, "Não existem valores ajustados disponíveis para a comparação de modelos."))

    ggplot() +
      geom_line(
        data = comparison_dat$obs,
        aes(x = year, y = value),
        color = "black",
        linewidth = 1
      ) +
      geom_line(
        data = fitted_df,
        aes(x = year, y = fitted, color = model),
        alpha = 0.9
      ) +
      labs(
        title = "Comparação de Modelos no Ajuste",
        subtitle = "A série observada surge a preto; as linhas coloridas mostram os valores ajustados de cada modelo.",
        x = "Ano",
        y = comparison_dat$y_label
      ) +
      scale_color_brewer(palette = "Set1", name = "Modelo", labels = get_model_labels) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }

  analyze_structural_breaks <- function(history) {
    df <- history$series %>%
      dplyr::arrange(year)

    validate(
      need(nrow(df) >= 5, "Seleccione pelo menos 5 anos observados para a análise de quebras estruturais.")
    )

    ts_y <- stats::ts(df$value, start = min(df$year), frequency = 1)
    time_index <- seq_len(nrow(df))

    # Segmented *trend* model, not a mean-shift model.
    #
    # `ts_y ~ 1` asks where the mean level changes. Portuguese mortality has
    # fallen steadily - the national standardised rate is down about 42% since
    # 1991 - and a mean-only model has no way to express a slope, so it explains
    # a smooth decline by chopping it into a staircase of level shifts. Those
    # breakpoints are artefacts of the trend, not events.
    #
    # `ts_y ~ time_index` gives each segment its own intercept and slope, so a
    # break is reported only where the level or the rate of change actually
    # shifts. It costs one extra parameter per segment, so it needs a longer
    # series; where the trend model cannot be estimated the app falls back to
    # the mean-only model and records which was used.
    break_model <- "trend"
    bp_obj <- tryCatch(
      strucchange::breakpoints(ts_y ~ time_index),
      error = function(e) NULL
    )

    if (is.null(bp_obj)) {
      break_model <- "mean"
      bp_obj <- tryCatch(
        strucchange::breakpoints(ts_y ~ 1),
        error = function(e) NULL
      )
    }

    break_index <- integer(0)

    if (!is.null(bp_obj) && length(bp_obj$breakpoints) > 0) {
      break_index <- as.integer(bp_obj$breakpoints)
      break_index <- break_index[is.finite(break_index) & break_index >= 1 & break_index < nrow(df)]
    }

    break_years <- if (length(break_index) > 0) {
      as.integer(df$year[break_index])
    } else {
      integer(0)
    }

    segment_starts <- c(1L, break_index + 1L)
    segment_ends <- c(break_index, nrow(df))

    segment_tbl <- tibble(
      Segment = seq_along(segment_starts),
      `Ano Inicial` = df$year[segment_starts],
      `Ano Final` = df$year[segment_ends],
      Anos = segment_ends - segment_starts + 1L,
      `Taxa Média` = round(
        purrr::map2_dbl(segment_starts, segment_ends, ~ mean(df$value[.x:.y], na.rm = TRUE)),
        2
      ),
      # Slope per year within the segment. With a segmented trend model this is
      # the quantity a break is about, so showing it lets the reader see whether
      # a detected break is a change of direction or only of level.
      `Variação Anual` = round(
        purrr::map2_dbl(segment_starts, segment_ends, function(from, to) {
          if (to - from < 1) return(NA_real_)
          stats::coef(stats::lm(df$value[from:to] ~ seq_len(to - from + 1)))[[2]]
        }),
        2
      ),
      `Ano da Quebra` = c(break_years, NA_integer_)
    )

    list(
      history = history,
      series = df,
      breakpoints = bp_obj,
      break_index = break_index,
      break_years = break_years,
      break_model = break_model,
      segments = segment_tbl
    )
  }

  safely_analyze_structural_breaks <- function(history) {
    tryCatch(
      analyze_structural_breaks(history),
      error = function(e) NULL
    )
  }

  build_break_plot <- function(break_info) {
    df <- break_info$series

    p <- ggplot(df, aes(x = year, y = value)) +
      geom_line(color = "#0b2e4f", linewidth = 1) +
      geom_point(color = "#0b2e4f", size = 2) +
      labs(
        title = paste("Possíveis Quebras Estruturais -", break_info$history$area_label),
        x = "Ano",
        y = break_info$history$y_label
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

    if (length(break_info$break_years) > 0) {
      marker_df <- tibble(
        year = break_info$break_years,
        value = df$value[match(break_info$break_years, df$year)]
      )

      p <- p +
        geom_vline(
          xintercept = break_info$break_years,
          color = "#b22222",
          linetype = "dashed",
          linewidth = 0.7
        ) +
        geom_point(
          data = marker_df,
          aes(x = year, y = value),
          inherit.aes = FALSE,
          color = "#b22222",
          size = 2.5
        )
    }

    p
  }

  # A detected break that sits on a population seam is at least partly the
  # denominator changing basis. Naming it here is the whole point of the
  # exercise: the breaks tab is where the reader meets the year, and without
  # this they have no way to tell an artefact from an event.
  #
  # A breakpoint is the last year of a segment, so a seam at `year` shows up as
  # a break reported at `year - 1` or, when the fit places it one step later, at
  # `year` itself.
  describe_breaks_on_population_seams <- function(break_years) {
    if (length(break_years) == 0) {
      return(NULL)
    }

    hits <- Filter(
      function(seam) any(as.integer(break_years) %in% c(seam$year - 1L, seam$year)),
      POPULATION_SEAMS
    )

    if (length(hits) == 0) {
      return(NULL)
    }

    described <- vapply(
      hits,
      function(seam) as.character(glue::glue("{seam$year - 1}/{seam$year}")),
      character(1)
    )

    as.character(glue::glue(
      "Atenção: {if (length(hits) == 1) 'uma das quebras coincide' else 'algumas das quebras coincidem'} ",
      "com {if (length(hits) == 1) 'uma mudança' else 'mudanças'} de base da população do INE, em ",
      "{paste(described, collapse = ' e em ')}. Nesse ponto a fonte da população muda de indicador, ",
      "o que desloca as taxas mesmo sem qualquer alteração da mortalidade — na taxa padronizada de ",
      "Portugal, cerca de metade da descida de 2013 para 2014 vem daí. Uma quebra detectada aqui não ",
      "deve ser lida como um acontecimento epidemiológico sem confirmação por outra via."
    ))
  }

  build_break_interpretation_text <- function(break_info) {
    break_n <- length(break_info$break_years)
    segment_n <- nrow(break_info$segments)
    seam_note <- describe_breaks_on_population_seams(break_info$break_years)
    seam_note <- if (is.null(seam_note)) "" else paste0(" ", seam_note)

    # The wording has to match the model that was actually fitted: a break in a
    # segmented trend means the level or the slope changed, whereas a break in
    # a mean-only model means only the average level moved.
    trend_model <- identical(break_info$break_model, "trend")
    what_changed <- if (trend_model) "o nível ou a inclinação da tendência" else "o nível médio"
    model_note <- if (trend_model) {
      "O modelo ajusta uma tendência própria a cada segmento, pelo que uma descida contínua e regular não é assinalada como quebra."
    } else {
      "Atenção: a série é demasiado curta para ajustar tendências por segmento, pelo que foi usado um modelo apenas de nível médio; numa série com tendência marcada este modelo tende a assinalar quebras que são apenas o declive."
    }

    if (break_n == 0) {
      return(glue::glue(
        "A rotina de breakpoint não detectou qualquer potencial mudança estrutural, pelo que {what_changed} parece relativamente estável no histórico seleccionado. {model_note}"
      ))
    }

    if (break_n == 1) {
      return(glue::glue(
        "Foi detectada uma potencial mudança estrutural em torno de {break_info$break_years[[1]]}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que {what_changed} poderá ter mudado nessa altura. {model_note}{seam_note}"
      ))
    }

    glue::glue(
      "Foram detectadas {break_n} potenciais mudanças estruturais em torno de {paste(break_info$break_years, collapse = ', ')}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que {what_changed} poderá não ser estável ao longo do período observado. {model_note}{seam_note}"
    )
  }

  get_simple_structural_warning_text <- function(history) {
    break_info <- safely_analyze_structural_breaks(history)

    if (is.null(break_info) || length(break_info$break_years) == 0) {
      return(NULL)
    }

    # Infant mortality is per live birth: INE's population revisions do not touch it.
    seam_note <- if (is.null(history$infant)) describe_breaks_on_population_seams(break_info$break_years)
    seam_note <- if (is.null(seam_note)) "" else paste0(" ", seam_note)

    if (length(break_info$break_years) == 1) {
      return(glue::glue(
        "Foi detectada uma possível mudança estrutural em torno de {break_info$break_years[[1]]}; as previsões baseadas em todo o histórico devem ser lidas com cautela acrescida.{seam_note}"
      ))
    }

    glue::glue("Foi detectada uma possível mudança estrutural na série histórica; as previsões baseadas em todo o histórico devem ser lidas com cautela acrescida.{seam_note}")
  }

  value_or_default <- function(x, default) {
    if (is.null(x)) default else x
  }

  get_successful_model_ids <- function(dat) {
    if (is.null(dat)) {
      return(character(0))
    }

    value_or_default(dat$fitted_models, character(0))
  }

  has_successful_forecast <- function(dat) {
    length(get_successful_model_ids(dat)) > 0 &&
      !is.null(dat$fc) &&
      nrow(dat$fc) > 0
  }

  validate_successful_forecast <- function(dat) {
    validate(
      need(
        has_successful_forecast(dat),
        "Erro detectado na previsão. Nenhum modelo pôde ser estimado com sucesso para esta série; consulte os avisos dos modelos."
      )
    )
  }

  get_default_forecast_output_model <- function(dat) {
    if (is.null(dat)) {
      return(NULL)
    }

    successful_models <- get_successful_model_ids(dat)

    if (!is.null(dat$recommended_model) && dat$recommended_model %in% successful_models) {
      return(dat$recommended_model)
    }

    if (length(successful_models) > 0) {
      return(successful_models[[1]])
    }

    NULL
  }

  resolve_selected_successful_model <- function(dat, selected_model = NULL) {
    successful_models <- get_successful_model_ids(dat)
    default_model <- get_default_forecast_output_model(dat)

    if (length(successful_models) == 0 || is.null(default_model)) {
      return(NULL)
    }

    selected_model <- value_or_default(selected_model, default_model)

    if (!is.null(selected_model) && selected_model %in% successful_models) {
      return(selected_model)
    }

    default_model
  }

  get_named_model_choices <- function(model_ids) {
    stats::setNames(model_ids, get_model_labels(model_ids))
  }

  build_forecast_download_filename <- function(
    prefix,
    extension,
    dat = NULL,
    fallback_area_label = NULL,
    fallback_cause_label = NULL,
    view_mode = "compare",
    selected_model = NULL
  ) {
    area_label <- if (!is.null(dat)) dat$history$spec$area_label else fallback_area_label
    cause_label <- if (!is.null(dat)) dat$history$spec$cause else fallback_cause_label
    model_label <- get_model_labels(value_or_default(selected_model, "single"))
    suffix <- if (identical(view_mode, "single")) {
      paste0("_", safe_filename_token(model_label))
    } else {
      "_compare"
    }

    paste0(
      prefix,
      "_",
      safe_filename_token(area_label),
      "_",
      safe_filename_token(cause_label),
      suffix,
      "_",
      Sys.Date(),
      extension
    )
  }

  get_forecast_output_fc <- function(dat, view_mode = "compare", selected_model = NULL) {
    fc <- dat$fc

    if (!identical(view_mode, "single")) {
      return(fc)
    }

    model_id <- value_or_default(selected_model, get_default_forecast_output_model(dat))

    if (is.null(model_id)) {
      return(fc[0, , drop = FALSE])
    }

    fc %>%
      dplyr::filter(model == model_id)
  }

  build_advanced_forecast_plot <- function(dat, view_mode = "compare", selected_model = NULL) {
    obs <- dat$obs
    fc <- get_forecast_output_fc(dat, view_mode = view_mode, selected_model = selected_model)

    validate(need(nrow(fc) > 0, "Nenhum modelo de projecção foi estimado com sucesso. Consulte os avisos abaixo."))

    if (identical(view_mode, "single")) {
      model_id <- unique(fc$model)[[1]]
      model_label <- get_model_labels(model_id)

      return(
        ggplot() +
          geom_line(data = obs, aes(x = year, y = value), linewidth = 1) +
          geom_point(data = obs, aes(x = year, y = value), size = 2) +
          geom_ribbon(
            data = fc,
            aes(x = year, ymin = lower, ymax = upper),
            fill = "#8fb8de",
            alpha = 0.3
          ) +
          geom_line(
            data = fc,
            aes(x = year, y = mean),
            color = "#0b2e4f",
            linetype = "dashed",
            linewidth = 1
          ) +
          labs(
            title = paste("Resultados da Projecção -", dat$area_label),
            subtitle = paste("Vista de modelo único:", model_label),
            x = "Ano",
            y = dat$history$y_label
          ) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }

    ggplot() +
      geom_line(data = obs, aes(x = year, y = value), linewidth = 1) +
      geom_point(data = obs, aes(x = year, y = value), size = 2) +
      geom_ribbon(
        data = fc,
        aes(x = year, ymin = lower, ymax = upper, fill = model),
        alpha = 0.15
      ) +
      geom_line(
        data = fc,
        aes(x = year, y = mean, color = model),
        linetype = "dashed",
        linewidth = 1
      ) +
      labs(
        title = paste("Resultados da Projecção -", dat$area_label),
        subtitle = "Vista comparativa dos modelos ajustados",
        x = "Ano",
        y = dat$history$y_label
      ) +
      scale_color_brewer(palette = "Set1", name = "Modelo", labels = get_model_labels) +
      scale_fill_brewer(palette = "Set1", name = "Modelo", labels = get_model_labels) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }

  build_forecast_summary_table <- function(dat, view_mode = "compare", selected_model = NULL) {
    successful_models <- get_successful_model_ids(dat)
    validate(
      need(length(successful_models) > 0, "Não existe resumo da projecção porque nenhum modelo foi estimado com sucesso.")
    )

    focal_model <- if (identical(view_mode, "single")) {
      value_or_default(selected_model, get_default_forecast_output_model(dat))
    } else {
      get_default_forecast_output_model(dat)
    }

    obs_last <- dat$obs %>%
      dplyr::slice_tail(n = 1)
    focal_fc <- dat$fc %>%
      dplyr::filter(model == focal_model) %>%
      dplyr::slice_tail(n = 1)

    validate(need(nrow(focal_fc) > 0, "Não existe resumo da projecção porque nenhum modelo foi estimado com sucesso."))

    abs_change <- focal_fc$mean - obs_last$value
    pct_change <- if (isTRUE(all.equal(obs_last$value, 0))) {
      NA_real_
    } else {
      (abs_change / obs_last$value) * 100
    }

    tibble(
      Item = c(
        "Vista apresentada",
        "Modelo de resumo",
        "Modelos seleccionados",
        "Modelos estimados com sucesso",
        "Último ano observado",
        "Última taxa observada",
        "Último ano projectado",
        "Valor final projectado",
        "Intervalo da projecção",
        "Variação absoluta face ao último observado",
        "Variação percentual face ao último observado"
      ),
      Valor = c(
        if (identical(view_mode, "single")) "Modelo único" else "Comparação entre modelos",
        get_model_labels(focal_model),
        paste(get_model_labels(dat$selected_model_ids), collapse = ", "),
        if (length(successful_models) > 0) paste(get_model_labels(successful_models), collapse = ", ") else "Nenhum",
        as.character(obs_last$year),
        sprintf("%.2f", obs_last$value),
        as.character(focal_fc$year),
        sprintf("%.2f", focal_fc$mean),
        glue::glue("{round(focal_fc$lower, 2)} a {round(focal_fc$upper, 2)}"),
        sprintf("%+.2f", abs_change),
        if (is.na(pct_change)) "N/D" else sprintf("%+.2f%%", pct_change)
      )
    )
  }

  build_forecast_warning_ui <- function(dat) {
    if (is.null(dat)) {
      return(NULL)
    }

    failures <- value_or_default(
      dat$failures,
      tibble::tibble(Model = character(0), Message = character(0))
    )
    successful_models <- get_successful_model_ids(dat)

    if (nrow(failures) == 0 && length(successful_models) > 0) {
      return(NULL)
    }

    wellPanel(
      h4(if (length(successful_models) == 0) "Erro detectado na previsão" else "Avisos dos modelos"),
      p(
        if (length(successful_models) == 0) {
          "Nenhum dos modelos pedidos pôde ser estimado para esta especificação. A aplicação não apresenta a previsão como válida para esta selecção."
        } else {
          "Alguns dos modelos pedidos não puderam ser estimados. Os resultados abaixo utilizam os modelos que foram ajustados com sucesso."
        }
      ),
      if (nrow(failures) > 0) {
        tags$ul(
          lapply(seq_len(nrow(failures)), function(i) {
            tags$li(
              paste0(get_model_labels(failures$Model[[i]]), ": ", failures$Message[[i]])
            )
          })
        )
      } else {
        p("Não foi devolvida uma mensagem técnica pelo estimador.")
      }
    )
  }

  build_advanced_model_panel <- function(model_id) {
    model_label <- get_model_labels(model_id)

    switch(
      model_id,
      arima = wellPanel(
        h4(model_label),
        radioButtons(
          "arima_mode",
          "Especificação ARIMA:",
          choices = c("Automática" = "auto", "Manual" = "manual"),
          selected = "auto",
          inline = TRUE
        ),
        conditionalPanel(
          "input.arima_mode == 'manual'",
          fluidRow(
            column(4, numericInput("arima_p", "p", value = 0, min = 0, max = 5)),
            column(4, numericInput("arima_d", "d", value = 1, min = 0, max = 2)),
            column(4, numericInput("arima_q", "q", value = 0, min = 0, max = 5))
          ),
          tags$hr(),
          h5("Termos sazonais"),
          fluidRow(
            column(3, numericInput("arima_P", "P", value = 0, min = 0, max = 3)),
            column(3, numericInput("arima_D", "D", value = 0, min = 0, max = 2)),
            column(3, numericInput("arima_Q", "Q", value = 0, min = 0, max = 3)),
            column(3, numericInput("arima_period", "Período", value = 1, min = 1, max = 20))
          ),
          checkboxInput("arima_include_constant", "Incluir constante / drift quando admissível", value = TRUE),
          helpText("Para dados anuais, mantenha o período sazonal em 1, salvo se pretender deliberadamente uma estrutura sazonal de vários anos.")
        )
      ),
      ets = wellPanel(
        h4(model_label),
        radioButtons(
          "ets_mode",
          "Especificação ETS:",
          choices = c("Automática" = "auto", "Manual" = "manual"),
          selected = "auto",
          inline = TRUE
        ),
        conditionalPanel(
          "input.ets_mode == 'manual'",
          selectInput("ets_error", "Erro", choices = c("Automático" = "Z", "Aditivo" = "A", "Multiplicativo" = "M"), selected = "Z"),
          selectInput("ets_trend", "Tendência", choices = c("Automática" = "Z", "Nenhuma" = "N", "Aditiva" = "A", "Multiplicativa" = "M"), selected = "Z"),
          selectInput("ets_season", "Sazonalidade", choices = c("Automática" = "Z", "Nenhuma" = "N", "Aditiva" = "A", "Multiplicativa" = "M"), selected = "N"),
          radioButtons(
            "ets_damped",
            "Tendência amortecida:",
            choices = c("Automática" = "auto", "Não" = "no", "Sim" = "yes"),
            selected = "auto",
            inline = TRUE
          ),
          numericInput("ets_period", "Período sazonal", value = 1, min = 1, max = 20),
          helpText("Uma especificação ETS sazonal requer um período sazonal superior a 1.")
        )
      ),
      tbats = wellPanel(
        h4(model_label),
        checkboxInput("tbats_use_box_cox", "Permitir transformação Box-Cox interna", value = FALSE),
        checkboxInput("tbats_use_trend", "Permitir tendência", value = TRUE),
        checkboxInput("tbats_use_damped_trend", "Permitir tendência amortecida", value = TRUE),
        checkboxInput("tbats_use_arma_errors", "Permitir erros ARMA", value = TRUE),
        textInput("tbats_seasonal_periods", "Períodos sazonais (separados por vírgulas, opcional)", value = ""),
        helpText("Deixe os períodos sazonais em branco para a série anual actual, salvo se pretender explicitamente ciclos sazonais de vários anos.")
      ),
      rwf = wellPanel(
        h4(model_label),
        p("Este modelo utiliza um passeio aleatório com drift. Não são expostas aqui definições adicionais específicas desta família.")
      ),
      naive = wellPanel(
        h4(model_label),
        p("Este modelo projecta em frente o nível observado mais recente. Não são expostas aqui definições adicionais específicas desta família.")
      ),
      theta = wellPanel(
        h4(model_label),
        p("Este modelo utiliza a abordagem padrão Theta. Não são expostas aqui definições adicionais específicas desta família.")
      ),
      holt = wellPanel(
        h4(model_label),
        p("Este modelo utiliza o método de tendência linear de Holt com parâmetros de suavização automáticos.")
      ),
      holt_damped = wellPanel(
        h4(model_label),
        p("Este modelo utiliza o método de tendência amortecida de Holt com parâmetros de suavização automáticos.")
      )
    )
  }

  parse_tbats_periods <- function(x) {
    if (is.null(x) || !nzchar(trimws(x))) {
      return(NULL)
    }

    vals <- strsplit(x, ",", fixed = TRUE)[[1]] %>%
      trimws()

    vals <- vals[nzchar(vals)]

    nums <- suppressWarnings(as.numeric(vals))

    if (length(nums) == 0 || any(!is.finite(nums)) || any(nums <= 1)) {
      validate(
        need(FALSE, "Os períodos sazonais de TBATS têm de ser uma lista de números superiores a 1, separados por vírgulas.")
      )
    }

    nums
  }

  get_model_frequency <- function(model_id, model_spec) {
    if (identical(model_id, "arima") &&
        identical(model_spec$mode, "manual") &&
        isTRUE(model_spec$seasonal_period > 1) &&
        any(model_spec$seasonal_order > 0)) {
      return(model_spec$seasonal_period)
    }

    if (identical(model_id, "ets") &&
        identical(model_spec$mode, "manual") &&
        !identical(model_spec$season, "N") &&
        isTRUE(model_spec$seasonal_period > 1)) {
      return(model_spec$seasonal_period)
    }

    1
  }

  get_transformation_setup <- function(values, transform_method) {
    if (identical(transform_method, "none")) {
      return(list(
        method = "none",
        label = "Sem transformação",
        offset = 0,
        zero_count = 0L,
        forward = function(x) as.numeric(x),
        inverse = function(x) as.numeric(x),
        inverse_mean = function(x, sigma) as.numeric(x)
      ))
    }

    if (any(values < 0, na.rm = TRUE)) {
      validate(
        need(FALSE, "A série contém valores negativos e não pode usar a transformação log.")
      )
    }

    # The log offset is a data-dependent pseudo-count (half the smallest
    # positive rate). It matters most when the series contains zeros, so expose
    # its value in the label and flag the zero case for the user.
    min_positive <- suppressWarnings(min(values[values > 0], na.rm = TRUE))
    log_offset <- if (is.finite(min_positive)) min_positive / 2 else 1e-6
    zero_count <- sum(values == 0, na.rm = TRUE)

    label <- as.character(glue::glue("Transformação log com offset (offset = {signif(log_offset, 4)})"))
    if (zero_count > 0) {
      label <- as.character(glue::glue(
        "{label}; {zero_count} valor(es) zero na série — o offset influencia a projecção e os intervalos"
      ))
    }

    list(
      method = "log_offset",
      label = label,
      offset = log_offset,
      zero_count = as.integer(zero_count),
      forward = function(x) log(as.numeric(x) + log_offset),
      # Back-transforming a forecast of log(rate) with exp() returns the
      # *median*, not the mean: for a lognormal, E[Y] = exp(mu + sigma^2 / 2).
      # `inverse` is therefore correct for interval limits, which are quantiles
      # and map through any monotone transform unchanged.
      inverse = function(x) pmax(exp(as.numeric(x)) - log_offset, 0),
      # `inverse_mean` adds the variance term to recover the expected value.
      # The gap grows with the forecast variance, so it is widest exactly where
      # the app allows the longest horizons - the point forecast would otherwise
      # sit increasingly below the mean the further out it goes.
      inverse_mean = function(x, sigma) {
        sigma <- if (length(sigma) == 0) 0 else sigma
        pmax(exp(as.numeric(x) + sigma^2 / 2) - log_offset, 0)
      }
    )
  }

  build_advanced_model_specs <- reactive({
    selected_models <- value_or_default(input$models, character(0))

    lapply(selected_models, function(model_id) {
      switch(
        model_id,
        arima = list(
          mode = value_or_default(input$arima_mode, "auto"),
          order = c(value_or_default(input$arima_p, 0), value_or_default(input$arima_d, 1), value_or_default(input$arima_q, 0)),
          seasonal_order = c(value_or_default(input$arima_P, 0), value_or_default(input$arima_D, 0), value_or_default(input$arima_Q, 0)),
          seasonal_period = value_or_default(input$arima_period, 1),
          include_constant = value_or_default(input$arima_include_constant, TRUE)
        ),
        ets = list(
          mode = value_or_default(input$ets_mode, "auto"),
          error = value_or_default(input$ets_error, "Z"),
          trend = value_or_default(input$ets_trend, "Z"),
          season = value_or_default(input$ets_season, "N"),
          damped = value_or_default(input$ets_damped, "auto"),
          seasonal_period = value_or_default(input$ets_period, 1)
        ),
        tbats = list(
          use_box_cox = value_or_default(input$tbats_use_box_cox, FALSE),
          use_trend = value_or_default(input$tbats_use_trend, TRUE),
          use_damped_trend = value_or_default(input$tbats_use_damped_trend, TRUE),
          use_arma_errors = value_or_default(input$tbats_use_arma_errors, TRUE),
          seasonal_periods = parse_tbats_periods(input$tbats_seasonal_periods)
        ),
        rwf = list(),
        naive = list(),
        theta = list(),
        holt = list(),
        holt_damped = list()
      )
    }) %>%
      stats::setNames(selected_models)
  })

  summarize_model_spec <- function(model_id, model_spec) {
    switch(
      model_id,
      arima = if (identical(model_spec$mode, "auto")) {
        "ARIMA automática"
      } else {
        seasonal_txt <- if (isTRUE(model_spec$seasonal_period > 1) &&
                            any(model_spec$seasonal_order > 0)) {
          paste0(
            "; seasonal (",
            paste(model_spec$seasonal_order, collapse = ","),
            ")[",
            model_spec$seasonal_period,
            "]"
          )
        } else {
          "; sem termos sazonais"
        }
        paste0(
          "ARIMA manual(",
          paste(model_spec$order, collapse = ","),
          ")",
          seasonal_txt
        )
      },
      ets = if (identical(model_spec$mode, "auto")) {
        "ETS automática"
      } else {
        paste0(
          "ETS manual(",
          paste0(model_spec$error, model_spec$trend, model_spec$season),
          "), amortecida=",
          model_spec$damped,
          if (!identical(model_spec$season, "N")) paste0(", período=", model_spec$seasonal_period) else ""
        )
      },
      tbats = paste0(
        "use.box.cox=", model_spec$use_box_cox,
        ", tendência=", model_spec$use_trend,
        ", amortecida=", model_spec$use_damped_trend,
        ", erros.arma=", model_spec$use_arma_errors,
        if (length(model_spec$seasonal_periods) > 0) {
          paste0(", períodos.sazonais=", paste(model_spec$seasonal_periods, collapse = ","))
        } else {
          ", períodos.sazonais=nenhum"
        }
      ),
      rwf = "Passeio aleatório com drift padrão",
      naive = "Modelo naive padrão",
      theta = "Método Theta padrão",
      holt = "Tendência linear de Holt padrão",
      holt_damped = "Tendência amortecida de Holt padrão"
    )
  }

  output$advancedModelParameterPanels <- renderUI({
    selected_models <- value_or_default(input$models, character(0))

    if (length(selected_models) == 0) {
      return(NULL)
    }

    tagList(lapply(selected_models, build_advanced_model_panel))
  })

  choose_recommended_model <- function(accuracy_tbl) {
    if (nrow(accuracy_tbl) == 0) {
      return(NULL)
    }

    metric_priority <- c("RMSE", "MAE", "MASE", "MAPE")

    for (metric in metric_priority) {
      candidates <- accuracy_tbl %>%
        dplyr::filter(is.finite(.data[[metric]])) %>%
        dplyr::arrange(.data[[metric]], Model)

      if (nrow(candidates) > 0) {
        return(candidates$Model[[1]])
      }
    }

    accuracy_tbl$Model[[1]]
  }

  # Shared model runner used by both the guided and advanced forecasting paths.
  # It preserves the current model fitting logic and adds a recommended-model
  # pick based on the same in-sample comparison metrics already used elsewhere.
  run_forecast_models <- function(
    history,
    model_ids,
    horizon,
    kind = NULL,
    token = NULL,
    conf_level = 95,
    transform_method = "log_offset",
    model_specs = list(),
    bias_adjust = TRUE
  ) {
    abort_if_requested <- function() {
      if (!is.null(kind) && !is.null(token)) {
        abort_if_cancelled(kind, token)
      }
    }

    df_vals <- history$series

    if (nrow(df_vals) < 3) {
      validate(
        need(FALSE, "Selecione pelo menos 3 anos para ajustar o modelo de projecção.")
      )
    }

    if (length(model_ids) == 0) {
      validate(
        need(FALSE, "Selecione pelo menos um modelo de projecção.")
      )
    }

    if (any(!is.finite(df_vals$value))) {
      validate(
        need(FALSE, "A série de taxas contém valores não finitos e não pode ser modelada.")
      )
    }

    transform_setup <- get_transformation_setup(df_vals$value, transform_method)
    transformed_values <- transform_setup$forward(df_vals$value)

    fit_results <- lapply(model_ids, function(m) {
      abort_if_requested()
      model_spec <- model_specs[[m]]
      if (is.null(model_spec)) {
        model_spec <- list()
      }

      ts_data <- stats::ts(
        transformed_values,
        start = min(df_vals$year),
        frequency = get_model_frequency(m, model_spec)
      )

      tryCatch({
        fit <- switch(
          m,
          arima = {
            if (identical(model_spec$mode, "manual")) {
              if (any(model_spec$seasonal_order > 0) && model_spec$seasonal_period <= 1) {
              stop("Os termos sazonais do ARIMA manual requerem um período sazonal superior a 1.")
              }

              arima_args <- list(
                y = ts_data,
                order = model_spec$order,
                include.constant = isTRUE(model_spec$include_constant)
              )

              if (model_spec$seasonal_period > 1 && any(model_spec$seasonal_order > 0)) {
                arima_args$seasonal <- list(
                  order = model_spec$seasonal_order,
                  period = model_spec$seasonal_period
                )
              }

              do.call(forecast::Arima, arima_args)
            } else {
              forecast::auto.arima(ts_data)
            }
          },
          ets = {
            if (identical(model_spec$mode, "manual")) {
              if (!identical(model_spec$season, "N") && model_spec$seasonal_period <= 1) {
                stop("A sazonalidade ETS manual requer um período sazonal superior a 1.")
              }

              ets_args <- list(
                y = ts_data,
                model = paste0(model_spec$error, model_spec$trend, model_spec$season)
              )

              if (!identical(model_spec$damped, "auto")) {
                ets_args$damped <- identical(model_spec$damped, "yes")
              }

              do.call(forecast::ets, ets_args)
            } else {
              forecast::ets(ts_data)
            }
          },
          rwf = forecast::rwf(ts_data, drift = TRUE, h = horizon, level = conf_level),
          naive = forecast::naive(ts_data, h = horizon, level = conf_level),
          theta = forecast::thetaf(ts_data, h = horizon, level = conf_level),
          tbats = {
            if (!identical(transform_method, "none") && isTRUE(model_spec$use_box_cox)) {
              stop("A Box-Cox interna do TBATS só pode ser activada quando a transformação global está definida para 'Sem transformação'.")
            }

            forecast::tbats(
              y = ts_data,
              use.box.cox = isTRUE(model_spec$use_box_cox),
              use.trend = isTRUE(model_spec$use_trend),
              use.damped.trend = isTRUE(model_spec$use_damped_trend),
              use.arma.errors = isTRUE(model_spec$use_arma_errors),
              seasonal.periods = model_spec$seasonal_periods,
              use.parallel = FALSE
            )
          },
          holt = forecast::holt(ts_data, h = horizon, level = conf_level),
          holt_damped = forecast::holt(ts_data, h = horizon, damped = TRUE, level = conf_level),
          stop("Modelo desconhecido: ", m)
        )

        fc <- if (inherits(fit, "forecast")) {
          fit
        } else {
          forecast::forecast(fit, h = horizon, level = conf_level)
        }

        list(status = "ok", fit = fit, forecast = fc)
      }, error = function(e) {
        list(status = "failed", message = conditionMessage(e))
      })
    })
    names(fit_results) <- model_ids

    successful_ids <- names(fit_results)[vapply(fit_results, function(x) identical(x$status, "ok"), logical(1))]
    failed_ids <- names(fit_results)[vapply(fit_results, function(x) identical(x$status, "failed"), logical(1))]

    fits <- lapply(fit_results[successful_ids], `[[`, "fit")
    fc_list <- lapply(fit_results[successful_ids], `[[`, "forecast")

    fc_df <- dplyr::bind_rows(lapply(names(fc_list), function(m) {
      fc <- fc_list[[m]]
      back_transform <- transform_setup$inverse

      mean_t <- as.numeric(fc$mean)
      lower_t <- if (is.null(dim(fc$lower))) as.numeric(fc$lower) else as.numeric(fc$lower[, ncol(fc$lower)])
      upper_t <- if (is.null(dim(fc$upper))) as.numeric(fc$upper) else as.numeric(fc$upper[, ncol(fc$upper)])

      # Forecast standard deviation on the modelling scale, recovered from the
      # interval the model reported. Needed to turn the back-transformed median
      # into a mean; see `inverse_mean` in get_transformation_setup().
      level <- if (!is.null(fc$level)) utils::tail(as.numeric(fc$level), 1) else conf_level
      z <- stats::qnorm(1 - (1 - level / 100) / 2)
      sigma <- if (is.finite(z) && z > 0) pmax((upper_t - lower_t) / (2 * z), 0) else rep(0, length(mean_t))
      sigma[!is.finite(sigma)] <- 0

      tibble(
        year  = seq(max(df_vals$year) + 1, by = 1, length.out = horizon),
        model = m,
        # Interval limits are quantiles and back-transform directly; only the
        # point forecast needs the bias correction.
        mean  = if (isTRUE(bias_adjust)) transform_setup$inverse_mean(mean_t, sigma) else back_transform(mean_t),
        lower = back_transform(lower_t),
        upper = back_transform(upper_t)
      )
    }))

    failure_tbl <- tibble(
      Model = failed_ids,
      Message = vapply(fit_results[failed_ids], `[[`, character(1), "message")
    )

    accuracy_tbl <- build_accuracy_table(
      fits,
      obs = df_vals,
      inverse = transform_setup$inverse
    )

    list(
      history = history,
      obs = df_vals,
      fc = fc_df,
      fits = fits,
      fitted_models = successful_ids,
      selected_model_ids = model_ids,
      failures = failure_tbl,
      accuracy = accuracy_tbl,
      recommended_model = choose_recommended_model(accuracy_tbl),
      horizon = horizon,
      conf_level = conf_level,
      transform_method = transform_method,
      bias_adjust = isTRUE(bias_adjust),
      transform_label = transform_setup$label,
      transform_inverse = transform_setup$inverse,
      model_specs = model_specs,
      area_label = history$area_label
    )
  }

  # ---- Out-of-sample model selection --------------------------------------
  # The recommended model is chosen from out-of-sample forecast accuracy rather
  # than in-sample fit, so short annual series do not simply reward the most
  # flexible model. The most recent `test_fraction` of the selected series forms
  # the evaluation region:
  #   * "single"  fits once on the earlier years and scores one multi-step
  #               forecast over the whole region;
  #   * "rolling" re-forecasts one step ahead from every origin in the region
  #               (expanding training window) and pools the errors.
  # When the series is too short to leave >= MIN_VALIDATION_TRAIN training years
  # and >= 1 test year, selection falls back to the in-sample accuracy table.
  run_rolling_validation <- function(base_result, test_size) {
    df <- base_result$obs %>% dplyr::arrange(year)
    n <- nrow(df)
    origins <- seq.int(n - test_size + 1L, n)
    scale_denom <- if (n >= 2) mean(abs(diff(df$value)), na.rm = TRUE) else NA_real_

    steps <- dplyr::bind_rows(lapply(origins, function(pos) {
      train_years <- df$year[seq_len(pos - 1L)]
      origin_year <- df$year[[pos]]
      training_history <- rebuild_history_with_year_range(
        history = base_result$history,
        year_range = range(train_years)
      )
      fold <- run_forecast_models(
        history = training_history,
        model_ids = base_result$selected_model_ids,
        horizon = 1L,
        conf_level = base_result$conf_level,
        transform_method = base_result$transform_method,
        model_specs = base_result$model_specs,
        bias_adjust = isTRUE(base_result$bias_adjust)
      )
      fold$fc %>%
        dplyr::filter(year == origin_year) %>%
        dplyr::transmute(model, year, predicted = mean, actual = df$value[[pos]])
    }))

    if (nrow(steps) == 0) {
      return(list(metrics = base_result$accuracy[0, , drop = FALSE], predictions = steps))
    }

    metrics <- steps %>%
      dplyr::group_by(model) %>%
      dplyr::group_modify(~ compute_forecast_error_metrics(.x$actual, .x$predicted, scale_denom)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Model = model) %>%
      dplyr::select(Model, ME, RMSE, MAE, MAPE, MASE)

    list(metrics = metrics, predictions = steps)
  }

  evaluate_model_selection <- function(base_result, method = "rolling", test_fraction = 0.25) {
    n <- nrow(base_result$obs)
    fitted_models <- get_successful_model_ids(base_result)
    requested <- if (method %in% c("insample", "single", "rolling")) method else "rolling"

    finalize <- function(method_used, metric_tbl, test_size = 0L, note = NULL, plot_data = NULL) {
      usable <- metric_tbl %>% dplyr::filter(Model %in% fitted_models)
      recommended <- choose_recommended_model(usable)
      if (is.null(recommended)) {
        recommended <- get_default_forecast_output_model(base_result)
      }
      train_range <- NULL
      test_range <- NULL
      if (test_size > 0L && !identical(method_used, "insample")) {
        yrs <- sort(base_result$obs$year)
        train_range <- range(utils::head(yrs, n - test_size))
        test_range <- range(utils::tail(yrs, test_size))
      }
      list(
        requested_method = requested,
        method_used = method_used,
        test_size = test_size,
        train_range = train_range,
        test_range = test_range,
        metrics = metric_tbl,
        ranking = rank_model_metrics(metric_tbl),
        recommended_model = recommended,
        fallback = !identical(method_used, requested),
        note = note,
        plot_data = plot_data
      )
    }

    if (identical(requested, "insample")) {
      return(finalize("insample", base_result$accuracy))
    }

    test_size <- compute_validation_test_size(n, test_fraction)
    if (test_size < 1L) {
      return(finalize(
        "insample",
        base_result$accuracy,
        note = glue::glue(
          "Série demasiado curta ({n} ano{if (n == 1) '' else 's'}) para reservar treino e teste; a recomendação recorre ao ajuste dentro da amostra."
        )
      ))
    }

    if (identical(requested, "single")) {
      split <- build_holdout_split(base_result$history, holdout_k = test_size)
      training_result <- run_forecast_models(
        history = split$training_history,
        model_ids = base_result$selected_model_ids,
        horizon = nrow(split$holdout_actual),
        conf_level = base_result$conf_level,
        transform_method = base_result$transform_method,
        model_specs = base_result$model_specs,
        bias_adjust = isTRUE(base_result$bias_adjust)
      )
      metric_tbl <- build_holdout_metric_table(training_result, split$holdout_actual)
      plot_data <- list(
        training_obs = split$training_history$series,
        holdout_actual = split$holdout_actual,
        forecast_df = training_result$fc,
        holdout_k = nrow(split$holdout_actual)
      )
      return(finalize("single", metric_tbl, test_size = test_size, plot_data = plot_data))
    }

    roll <- run_rolling_validation(base_result, test_size)
    plot_data <- list(
      full_obs = base_result$obs,
      predictions = roll$predictions,
      test_size = test_size
    )
    finalize("rolling", roll$metrics, test_size = test_size, plot_data = plot_data)
  }

  describe_validation_selection <- function(validation) {
    if (is.null(validation)) {
      return(NULL)
    }

    method_label <- switch(
      validation$method_used,
      insample = "ajuste dentro da amostra",
      single = "divisão única treino/teste",
      rolling = "validação móvel (origem deslizante)",
      validation$method_used
    )

    parts <- if (identical(validation$recommended_model, FORECAST_MEAN_ID)) {
      glue::glue("Projecção: média de 6 modelos (ARIMA, ETS, Theta, TBATS, Holt e Holt amortecido); a {method_label} mostra o erro de cada modelo.")
    } else {
      glue::glue("Escolha do modelo: {method_label}.")
    }
    if (!is.null(validation$train_range) && !is.null(validation$test_range)) {
      parts <- c(parts, glue::glue(
        "Treino {validation$train_range[[1]]}-{validation$train_range[[2]]}; teste {validation$test_range[[1]]}-{validation$test_range[[2]]} ({validation$test_size} ano{if (validation$test_size == 1) '' else 's'})."
      ))
    }
    if (!is.null(validation$note)) {
      parts <- c(parts, validation$note)
    }

    paste(parts, collapse = " ")
  }

  get_beginner_training_range <- function(years, year_range) {
    selected_years <- get_years_in_selected_range(year_range)
    training_years <- sort(intersect(as.integer(years), selected_years))

    validate(
      need(
        length(training_years) >= 3,
        "Seleccione pelo menos 3 anos de ajuste que já estejam carregados em 'Mortalidade Observada'."
      )
    )

    range(training_years)
  }

  # The guided forecast's infant mortality series (R/forecast_helpers.R): the
  # yearly rates when the models are to be fitted to them, the window rates
  # when they are to be fitted to those. `infant` keeps what the window
  # aggregation and the re-fits of the validation need.
  build_infant_history <- function(query_spec, year_range, window, method) {
    window <- as.integer(window)
    years <- sort(intersect(get_years_in_selected_range(year_range), infant_mortality_years()))
    validate(need(length(years) >= 3 + window - 1L, "Seleccione pelo menos 3 períodos com nados-vivos e óbitos infantis."))
    undercount <- infant_undercount_message(years, query_spec$area_key)
    if (!is.null(undercount)) notify_data(undercount, type = "warning", duration = 20)
    annual <- infant_annual_series(years, query_spec$area_key)
    validate(need(nrow(annual) >= 3, "Sem nados-vivos e óbitos infantis suficientes para esta selecção."))
    fit_annual <- window == 1L || identical(method, "annual")
    series <- infant_window_rates(annual, if (fit_annual) 1L else window)
    period <- c("1" = "anual", "3" = "triénio", "5" = "quinquénio")[[as.character(window)]]
    list(
      spec = c(query_spec, list(population = "Total", rate_type = "infant", years = series$year, year_range = range(series$year))),
      metric_bundle = NULL,
      series = series[, c("year", "value", "lower", "upper")],
      area_label = query_spec$area_label,
      y_label = paste0("Mortalidade infantil (", period, "), por 1.000 nados-vivos"),
      rate_label = "Mortalidade infantil",
      unit_label = "por 1.000 nados-vivos",
      measure_label = "a mortalidade infantil",
      infant = list(annual = annual, window = window, method = if (fit_annual) "annual" else "rolling")
    )
  }

  # Window rates in place of the yearly ones, once the models have run on them.
  aggregate_infant_result <- function(result, history) {
    info <- history$infant
    if (is.null(info) || info$window == 1L || !identical(info$method, "annual")) return(result)
    windowed <- infant_window_rates(info$annual, info$window)[, c("year", "value", "lower", "upper")]
    result$fc <- infant_window_forecast(info$annual, result$fc, info$window)
    result$obs <- windowed
    result$full_history$series <- windowed
    result$history$series <- windowed
    result
  }

  get_beginner_training_label <- function(training_history) {
    paste0("os anos ", format_year_selection(training_history$series$year))
  }

  build_beginner_forecast_plot <- function(dat) {
    validate_successful_forecast(dat)

    full_obs <- dat$full_history$series
    train_obs <- dat$obs
    recommended_fc <- dat$fc %>%
      dplyr::filter(model == dat$recommended_model)

    p <- ggplot() +
      geom_line(
        data = full_obs,
        aes(x = year, y = value),
        color = "grey75",
        linewidth = 0.8
      ) +
      geom_line(
        data = train_obs,
        aes(x = year, y = value),
        color = "#1f4e79",
        linewidth = 1
      ) +
      geom_point(
        data = train_obs,
        aes(x = year, y = value),
        color = "#1f4e79",
        size = 2
      ) +
      geom_ribbon(
        data = recommended_fc,
        aes(x = year, ymin = lower, ymax = upper),
        fill = "#8fb8de",
        alpha = 0.3
      ) +
      geom_line(
        data = recommended_fc,
        aes(x = year, y = mean),
        color = "#0b2e4f",
        linetype = "dashed",
        linewidth = 1.1
      ) +
      labs(
        title = paste("Previsão Guiada -", dat$area_label),
        subtitle = if (identical(dat$mode, "compare")) {
          "As linhas cinzentas mostram trajectórias alternativas dos modelos."
        } else {
          NULL
        },
        x = "Ano",
        y = dat$history$y_label,
        caption = paste(
          "Azul: taxa observada usada no ajuste. Cinzento: restante histórico.",
          "Linha tracejada e zona sombreada: previsão e a sua incerteza",
          "(quanto mais larga, menos certa)."
        )
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0)
      )

    if (identical(dat$mode, "compare")) {
      p <- p +
        geom_line(
          data = dat$fc %>% dplyr::filter(model != dat$recommended_model),
          aes(x = year, y = mean, group = model),
          color = "grey60",
          linetype = "dotted",
          alpha = 0.8
        )
    }

    p
  }

  build_beginner_summary_ui <- function(dat) {
    validate_successful_forecast(dat)

    last_observed <- dat$full_history$series %>%
      dplyr::slice_tail(n = 1)
    final_forecast <- dat$fc %>%
      dplyr::filter(model == dat$recommended_model) %>%
      dplyr::slice_tail(n = 1)
    model_endpoints <- dat$fc %>%
      dplyr::group_by(model) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()

    change_abs <- final_forecast$mean - last_observed$value
    change_pct <- if (isTRUE(all.equal(last_observed$value, 0))) {
      NA_real_
    } else {
      (change_abs / last_observed$value) * 100
    }
    interval_ratio <- (final_forecast$upper - final_forecast$lower) / pmax(final_forecast$mean, 1e-6)
    agreement_ratio <- (max(model_endpoints$mean) - min(model_endpoints$mean)) / pmax(abs(final_forecast$mean), 1e-6)

    trend_text <- if (is.na(change_pct) || abs(change_pct) < 5) {
      "se manterá relativamente estável"
    } else if (change_pct > 0) {
      "aumentará"
    } else {
      "diminuirá"
    }

    uncertainty_text <- if (interval_ratio < 0.25) {
      "O intervalo da projecção mantém-se relativamente estreito no final do horizonte."
    } else if (interval_ratio < 0.6) {
      "O intervalo da projecção alarga-se ao longo do tempo, pelo que os anos mais distantes devem ser lidos com alguma cautela."
    } else {
      "O intervalo da projecção torna-se amplo no final do horizonte, pelo que a trajectória de longo prazo é incerta."
    }

    compare_text <- if (agreement_ratio < 0.15) {
      "Os diferentes modelos plausíveis contam uma história muito semelhante."
    } else if (agreement_ratio < 0.35) {
      "Os diferentes modelos plausíveis apontam na mesma direcção geral, mas não exactamente para o mesmo nível."
    } else {
      "Os diferentes modelos plausíveis divergem bastante quanto ao nível final."
    }

    unit <- dat$history$unit_label %||% "por 100.000"
    measure <- dat$history$measure_label %||% "a taxa de mortalidade"

    wellPanel(
      h4("Resumo em linguagem simples"),
      p(glue::glue(
        "Usando {dat$training_label}, a previsão recomendada sugere que {measure} {trend_text} ao longo dos próximos {dat$horizon} anos."
      )),
      p(glue::glue(
        "Em {final_forecast$year}, a projecção central é de {round(final_forecast$mean, 2)} {unit}, em comparação com {round(last_observed$value, 2)} em {last_observed$year}."
      )),
      p(glue::glue(
        "Isto corresponde a uma variação de {sprintf('%+.2f', change_abs)} {unit}{if (!is.na(change_pct)) paste0(' (', sprintf('%+.2f%%', change_pct), ')') else ''}."
      )),
      p(uncertainty_text),
      p(compare_text)
    )
  }

  build_beginner_reliability_ui <- function(dat) {
    validate_successful_forecast(dat)

    recommended_fc <- dat$fc %>%
      dplyr::filter(model == dat$recommended_model) %>%
      dplyr::slice_tail(n = 1)
    model_endpoints <- dat$fc %>%
      dplyr::group_by(model) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()

    training_years <- nrow(dat$history$series)
    interval_ratio <- (recommended_fc$upper - recommended_fc$lower) / pmax(recommended_fc$mean, 1e-6)
    agreement_ratio <- (max(model_endpoints$mean) - min(model_endpoints$mean)) / pmax(abs(recommended_fc$mean), 1e-6)

    training_score <- dplyr::case_when(
      training_years >= 15 ~ 2L,
      training_years >= 8 ~ 1L,
      TRUE ~ 0L
    )
    interval_score <- dplyr::case_when(
      interval_ratio < 0.25 ~ 2L,
      interval_ratio < 0.6 ~ 1L,
      TRUE ~ 0L
    )
    agreement_score <- dplyr::case_when(
      agreement_ratio < 0.15 ~ 2L,
      agreement_ratio < 0.35 ~ 1L,
      TRUE ~ 0L
    )

    overall_score <- training_score + interval_score + agreement_score
    overall_label <- dplyr::case_when(
      overall_score >= 5 ~ "Elevada",
      overall_score >= 3 ~ "Moderada",
      TRUE ~ "Mais baixa"
    )

    training_message <- dplyr::case_when(
      training_years >= 15 ~ "A previsão baseia-se num período relativamente longo de dados observados.",
      training_years >= 8 ~ "A previsão utiliza uma quantidade moderada de histórico recente.",
      TRUE ~ "A previsão baseia-se numa janela de ajuste curta, pelo que é mais sensível às flutuações recentes."
    )
    interval_message <- dplyr::case_when(
      interval_ratio < 0.25 ~ "A banda de incerteza é relativamente estreita no final da projecção.",
      interval_ratio < 0.6 ~ "A banda de incerteza é moderada e torna-se mais visível com o passar do tempo.",
      TRUE ~ "A banda de incerteza é ampla no ano final."
    )
    agreement_message <- dplyr::case_when(
      agreement_ratio < 0.15 ~ "Os diferentes modelos concordam de perto quanto à direcção e ao nível da variação.",
      agreement_ratio < 0.35 ~ "Os diferentes modelos concordam em termos gerais quanto à direcção, mas diferem no nível final.",
      TRUE ~ "Os diferentes modelos dão pontos finais bastante distintos, o que reduz a confiança."
    )
    structural_warning <- get_simple_structural_warning_text(dat$full_history)
    validation_text <- describe_validation_selection(dat$validation)
    zero_note <- if (
      identical(dat$transform_method, "log_offset") &&
        any(dat$obs$value == 0, na.rm = TRUE)
    ) {
      "A série de treino contém anos com taxa zero; a transformação logarítmica usa um pequeno valor de deslocamento (offset) que influencia a projecção e os respectivos intervalos."
    } else {
      NULL
    }

    wellPanel(
      h4("Fiabilidade"),
      p(tags$strong(paste("Leitura global:", overall_label))),
      p(training_message),
      p(interval_message),
      p(agreement_message),
      if (!is.null(validation_text)) {
        p(tags$em(validation_text))
      },
      if (!is.null(zero_note)) {
        p(tags$strong(zero_note))
      },
      if (!is.null(structural_warning)) {
        p(tags$strong(structural_warning))
      }
    )
  }

  # Areas the active data source can actually supply, used to tell an
  # incomplete municipal rebuild apart from a missing region.
  get_available_areas <- function(data_source = "ine") {
    if (!identical(normalize_data_source(data_source), "snapshot")) {
      return(NULL)
    }

    tryCatch(
      {
        inventory <- get_snapshot_inventory()
        if (is.null(inventory) || nrow(inventory) == 0) NULL else unique(unlist(inventory$areas, use.names = FALSE))
      },
      error = function(e) NULL
    )
  }

  # Build one area spec, expanding a region into its municipalities when the
  # municipal region mode is active. `warnings` carries any incomplete-coverage
  # message up to the caller so it can be shown rather than swallowed.
  make_annual_area_spec <- function(label, areas, region_mode, available_areas,
                                    lookup = active_nuts_lookup()) {
    resolved <- resolve_region_areas(
      areas = areas,
      region_mode = region_mode,
      lookup = lookup,
      available_areas = available_areas
    )

    list(
      label = label,
      areas = resolved$areas,
      expanded = resolved$expanded,
      warnings = resolved$warnings
    )
  }

  # The annual tab's columns: the comparators, in the order chosen (the
  # automatic entry resolved to the local area's NUTS II region), then the local
  # area. resolve_annual_comparators() is in R/regions.R.
  annual_comparator_areas <- function(selected_areas) {
    resolve_annual_comparators(
      value_or_default(input$annual_comparators, ANNUAL_DEFAULT_COMPARATORS),
      selected_areas,
      lookup = active_nuts_lookup(),
      members = function(area) planning_area_members(area, active_nuts_lookup())
    )
  }

  get_annual_area_specs <- function(selected_areas,
                                    custom_label = NULL,
                                    region_mode = "original",
                                    available_areas = NULL) {
    selected_areas <- sort(unique(as.character(selected_areas)))

    validate(
      need(length(selected_areas) > 0, "Seleccione o local a comparar.")
    )
    comparators <- annual_comparator_areas(selected_areas)
    mismatch <- vintage_mismatch_message(comparators, active_nuts_vintage(), active_nuts_lookup())
    validate(need(is.null(mismatch), mismatch))

    mismatch <- vintage_mismatch_message(
      selected_areas,
      active_nuts_vintage(),
      active_nuts_lookup()
    )
    validate(need(is.null(mismatch), mismatch))

    # The local column sums whatever is selected into one geography, so an
    # overlapping choice double-counts exactly as it does on the other tabs.
    # Only reachable since NUTS I regions became selectable: before that every
    # offered region sat at one level and so was disjoint from the others.
    notify_uls_split(selected_areas)
    overlap <- overlapping_selection_warning(selected_areas, lookup = active_nuts_lookup())
    if (!is.null(overlap)) {
      notify_data(overlap, type = "warning", duration = 20)
    }

    comparator_specs <- lapply(comparators, function(area) {
      # Portugal is a national total under every NUTS vintage, so it is never
      # rebuilt from municipalities.
      if (identical(area, "Portugal")) {
        return(list(label = "Portugal", areas = "Portugal", expanded = character(0), warnings = character(0)))
      }
      make_annual_area_spec(area, area, region_mode, available_areas)
    })
    c(comparator_specs, list(make_annual_area_spec(
      get_selection_label(selected_areas, custom_label),
      selected_areas,
      region_mode,
      available_areas
    )))
  }

  # Years contributing to a pooled value centred on `year`, clipped to the
  # years the app knows about so an edge window is truncated rather than
  # requesting data that cannot exist.
  get_pooled_years <- function(year, window, metric_id = NULL) {
    window <- normalize_pooling_window(window)
    half <- (window - 1L) %/% 2L
    candidate <- seq.int(as.integer(year) - half, as.integer(year) + half)
    years <- sort(intersect(candidate, as.integer(year_of_interest)))

    # A window is cut to the years that have what the metric needs, rather than
    # failing the whole request: deaths for every metric but the infant ones,
    # and population as well for a rate (deaths without person-years would
    # inflate the numerator). Either side can run ahead - population now
    # reaches 2025 while deaths by cause stop at 2024, so a window centred on
    # 2024 reached 2025 and every metric was refused. The period label reports
    # the span actually used.
    if (!is.null(metric_id) && !metric_id %in% infant_metric_ids) {
      # Counts reach the rebuilt latest year as rates do; whether the selection
      # can use it is decided afterwards (calendar_selection_supported()).
      have <- if (metric_id %in% annual_metrics_needing_population) rate_years else union(death_archive_years, calendar_death_years)
      usable <- intersect(years, as.integer(have))
      if (length(usable) > 0) years <- usable
    }

    # Same reasoning for the infant metrics, whose windows are bounded by their
    # own datasets rather than by population: the rate needs published live
    # births, the count only the under-1 death archive.
    if (!is.null(metric_id) && metric_id %in% infant_metric_ids) {
      usable <- intersect(years, infant_metric_years(metric_id))
      if (length(usable) > 0) years <- usable
    }

    years
  }

  get_annual_metric_label <- function(metric_id) {
    metric_label <- names(annual_metric_choices)[match(metric_id, unname(annual_metric_choices))]
    ifelse(is.na(metric_label), NA_character_, metric_label)
  }

  first_numeric_or_na <- function(x) {
    if (length(x) == 0) {
      return(NA_real_)
    }
    as.numeric(x[[1]])
  }

  # Metrics needing a population denominator: rates, direct standardisation and
  # both indirect measures.
  annual_metrics_needing_population <- c("crude", "dsr", "smr", "isr")

  load_annual_cause_data <- function(area_spec, cause, sex, years, metric_id, data_source = "ine") {
    # The infant metrics read their own parallel datasets and take nothing from
    # the main archive. Loading it anyway is not merely wasted work: 2025 has
    # infant deaths and live births but no by-cause deaths, so the loader would
    # throw about a missing snapshot file for a year the metric can answer.
    if (metric_id %in% infant_metric_ids) {
      return(tibble(
        year = integer(0), area = character(0), sex = character(0),
        cause = character(0), age_band = character(0), deaths = numeric(0)
      ))
    }

    data <- if (metric_id %in% annual_metrics_needing_population) {
      get_data_for(
        area = area_spec$areas,
        cause = cause,
        years = years,
        year_order = "desc",
        data_source = data_source,
        population_basis = active_population_basis()
      )$full
    } else {
      get_death_data_for(
        area = area_spec$areas,
        cause = cause,
        years = years,
        year_order = "desc",
        data_source = data_source
      )
    }

    apply_region_source(data, area_spec$expanded) %>%
      dplyr::filter(sex == .env$sex)
  }

  # Collapse a multi-year, multi-area selection to one row per age band, which
  # is the shape every metric below consumes. Deaths and population are summed,
  # so a pooled denominator is person-years and a multi-area selection is one
  # combined geography - the same convention the app already uses for areas.
  #
  # Summing population is correct only because the frame holds exactly one
  # cause: each (year, area, band) then appears once. With two causes loaded the
  # same population would be counted twice and every rate would be halved -
  # which is precisely the mistake the avoidable-mortality tab made in the other
  # direction, where population repeated once per cause and had to be
  # deduplicated. Fail loudly here rather than return a plausible wrong number.
  collapse_annual_cause_data <- function(cause_data) {
    if ("cause" %in% names(cause_data)) {
      loaded_causes <- unique(as.character(cause_data$cause))
      if (length(loaded_causes) > 1) {
        stop(
          "collapse_annual_cause_data() recebeu ", length(loaded_causes),
          " causas; a população seria contada uma vez por causa. Carregue uma causa de cada vez.",
          call. = FALSE
        )
      }
    }

    cause_data %>%
      dplyr::group_by(age_band) %>%
      dplyr::summarise(
        deaths = sum(deaths, na.rm = TRUE),
        pop = if ("pop" %in% names(cause_data)) sum(pop, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
  }

  calculate_annual_metric_values <- function(area_spec,
                                             causes,
                                             metric_id,
                                             sex,
                                             year,
                                             token,
                                             data_source = "ine",
                                             pooling_window = 1L,
                                             reference_areas = "Portugal") {
    abort_if_cancelled("annual", token)

    validate(
      need(!is.na(get_annual_metric_label(metric_id)), "Seleccione uma métrica válida.")
    )

    pooling_window <- normalize_pooling_window(pooling_window)
    years <- get_pooled_years(year, pooling_window, metric_id)

    # The latest year has deaths by age only for whole regions and all causes
    # (and the reference of an indirect rate must be one too). Outside that, a
    # window reaching it is cut back to the years before; the year itself is
    # refused.
    calendar <- intersect(years, calendar_death_years)
    calendar_ok <- FALSE
    if (length(calendar) > 0 && !(metric_id %in% infant_metric_ids)) {
      calendar_ok <- metric_id %in% c("deaths", "crude", "dsr", "smr", "isr", "ypll") &&
        calendar_selection_supported(area_spec$areas, causes, sex) &&
        (!metric_id %in% c("smr", "isr") || calendar_selection_supported(reference_areas, causes, sex))
      if (!calendar_ok) {
        validate(need(!as.integer(year) %in% calendar, calendar_gap_message(calendar)))
        years <- setdiff(years, calendar)
      }
    }
    period_label <- make_period_label(min(years), max(years))

    # A rate needs a denominator for every year in the window. INE publishes
    # deaths ahead of population, so the most recent year(s) can be selectable
    # for counts while being impossible for rates; say so plainly instead of
    # failing later with an empty join.
    if (metric_id %in% metrics_requiring_population) {
      validate(
        need(
          length(years_without_population(years)) == 0,
          population_gap_message(years, get_annual_metric_label(metric_id))
        )
      )
    }

    # Everything except the infant pair reads the by-cause death archive, which
    # now ends a year before population does. Say so, rather than failing deep
    # in the loader with a message about a missing snapshot file.
    if (!(metric_id %in% infant_metric_ids)) {
      validate(
        need(
          length(years_without_deaths(years, calendar_ok)) == 0,
          death_gap_message(years, get_annual_metric_label(metric_id))
        )
      )
    }

    all_cause_deaths <- NA_real_
    all_cause_source <- NA_character_
    if (identical(metric_id, "proportional")) {
      # Built exactly like the numerator, or the share would divide a regional
      # row by a municipal sum.
      all_cause_data <- apply_region_source(
        get_death_data_for(
          area = area_spec$areas,
          cause = "Todas as causas de morte",
          years = years,
          year_order = "desc",
          data_source = data_source
        ),
        area_spec$expanded
      ) %>%
        dplyr::filter(sex == .env$sex)

      all_cause_deaths <- sum(all_cause_data$deaths, na.rm = TRUE)
      all_cause_source <- get_loaded_source_summary(all_cause_data)$death_source

      validate(
        need(
          is.finite(all_cause_deaths) && all_cause_deaths > 0,
          glue::glue("Sem denominador válido de todas as causas para {area_spec$label}.")
        )
      )
    }

    metric_label <- get_annual_metric_label(metric_id)

    dplyr::bind_rows(lapply(causes, function(cause) {
      abort_if_cancelled("annual", token)

      cause_data <- load_annual_cause_data(
        area_spec = area_spec,
        cause = cause,
        sex = sex,
        years = years,
        metric_id = metric_id,
        data_source = data_source
      )

      validate(
        need(
          nrow(cause_data) > 0 || metric_id %in% infant_metric_ids,
          glue::glue("Sem dados para {area_spec$label} / {cause}.")
        )
      )

      deaths <- sum(cause_data$deaths, na.rm = TRUE)
      source_summary <- get_loaded_source_summary(cause_data)

      # Age-band totals across the pooled window and the selected geography,
      # reused by every metric that works from an age distribution.
      pooled_bands <- collapse_annual_cause_data(cause_data)

      metric_values <- switch(
        metric_id,
        deaths = {
          ci <- compute_count_interval(deaths)
          list(
            value = deaths,
            lower = unname(ci[[1]]),
            upper = unname(ci[[2]]),
            source_detail = source_summary$death_source
          )
        },
        crude = {
          # Computed from the pooled age bands, so a 3- or 5-year window gives
          # total deaths over total person-years rather than a mean of yearly
          # rates.
          metric_row <- compute_metrics(
            pooled_bands %>% dplyr::mutate(year = .env$year, sex = .env$sex, cause = .env$cause)
          ) %>% dplyr::slice(1)
          list(
            value = first_numeric_or_na(metric_row$crude_rate),
            lower = first_numeric_or_na(metric_row$crude_lower),
            upper = first_numeric_or_na(metric_row$crude_upper),
            source_detail = glue::glue("Pop.: {source_summary$population_source}; Óbitos: {source_summary$death_source}")
          )
        },
        dsr = {
          metric_row <- compute_metrics(
            pooled_bands %>% dplyr::mutate(year = .env$year, sex = .env$sex, cause = .env$cause)
          ) %>% dplyr::slice(1)
          list(
            value = first_numeric_or_na(metric_row$dsr),
            lower = first_numeric_or_na(metric_row$dsr_lower),
            upper = first_numeric_or_na(metric_row$dsr_upper),
            source_detail = glue::glue("Pop.: {source_summary$population_source}; Óbitos: {source_summary$death_source}")
          )
        },
        smr = ,
        isr = {
          # Indirect standardisation: the reference area's age-specific rates
          # are applied to this area's age structure. The reference is loaded
          # over the same pooled window, cause and sex, so both sides of the
          # comparison rest on identical data.
          # The reference must be built the same way as the areas it is compared
          # against. Regions are summed from their municipalities, so a reference
          # of "Norte" has to be expanded too - otherwise a municipal-sum area
          # would be standardised against INE's own regional row, mixing the two
          # conventions inside a single ratio.
          reference_resolved <- resolve_region_areas(
            areas = reference_areas,
            region_mode = default_region_mode(),
            lookup = active_nuts_lookup(),
            available_areas = get_available_areas(data_source)
          )

          reference_bands <- collapse_annual_cause_data(
            load_annual_cause_data(
              # Carries the expansion so a regional reference is built the same
              # way as the regional areas compared against it.
              area_spec = list(
                label = "referência",
                areas = reference_resolved$areas,
                expanded = reference_resolved$expanded
              ),
              cause = cause,
              sex = sex,
              years = years,
              metric_id = "smr",
              data_source = data_source
            )
          )

          smr_row <- compute_smr(df_area = pooled_bands, df_ref = reference_bands)

          if (identical(metric_id, "smr")) {
            list(
              value = first_numeric_or_na(smr_row$smr),
              lower = first_numeric_or_na(smr_row$smr_lower),
              upper = first_numeric_or_na(smr_row$smr_upper),
              source_detail = glue::glue(
                "Pop.: {source_summary$population_source}; Óbitos: {source_summary$death_source}; ",
                "Ref.: {paste(reference_areas, collapse = '+')} ",
                "(O={round(smr_row$observed)}, E={round(smr_row$expected, 1)})"
              )
            )
          } else {
            list(
              value = first_numeric_or_na(smr_row$isr),
              lower = first_numeric_or_na(smr_row$isr_lower),
              upper = first_numeric_or_na(smr_row$isr_upper),
              source_detail = glue::glue(
                "Pop.: {source_summary$population_source}; Óbitos: {source_summary$death_source}; ",
                "Ref.: {paste(reference_areas, collapse = '+')}"
              )
            )
          }
        },
        proportional = {
          ci <- compute_proportion_interval(deaths, all_cause_deaths)
          list(
            value = deaths / all_cause_deaths * 100,
            lower = unname(ci[[1]]),
            upper = unname(ci[[2]]),
            source_detail = glue::glue("Óbitos: {source_summary$death_source}; Denom.: {all_cause_source}")
          )
        },
        infant_deaths = {
          # The count on its own. Available wherever the under-1 archive
          # reaches, with no birth denominator involved, so it is readable at
          # municipal scale where the rate is not.
          count <- infant_deaths_total(years, area_spec$areas, cause = cause, sex = sex)
          ci <- compute_count_interval(count)
          list(
            value = count,
            lower = unname(ci[[1]]),
            upper = unname(ci[[2]]),
            source_detail = get_loaded_source_summary(
              get_infant_death_data(years, area_spec$areas, cause = cause, sex = sex)
            )$death_source
          )
        },
        infant = {
          # Reads the parallel under-1 and births datasets rather than the main
          # pipeline, which cannot separate infant deaths from deaths at ages
          # one to four.
          imr <- compute_infant_mortality(
            deaths = infant_deaths_total(years, area_spec$areas, cause = cause, sex = sex),
            births = get_births_total(years, area_spec$areas)
          )

          list(
            value = imr$value,
            lower = imr$lower,
            upper = imr$upper,
            # Marks a rate whose denominator is too small to be read as a rate.
            # The value is still shown; see infant_instability_footnote().
            flag = if (infant_rate_is_unstable(imr$births)) infant_instability_mark else "",
            source_detail = glue::glue(
              "Óbitos <1 ano: {imr$deaths}; Nados-vivos: {format(imr$births, big.mark = '')}"
            )
          )
        },
        ypll = {
          # Correct the weight of the `0 - 4 anos` band before summing years
          # lost: nearly all of its deaths are infants, who lose almost the
          # whole cutoff rather than the 67.5 years the band midpoint implies.
          # Returns the frame unchanged when the under-1 counts do not cover the
          # window, so AVPP is never half-corrected.
          under_one <- infant_deaths_total(years, area_spec$areas, cause = cause, sex = sex, complete = FALSE)
          ypll_bands <- split_infant_age_band(pooled_bands, under_one)

          ci <- compute_ypll_interval(ypll_bands, cutoff = 70)
          list(
            value = unname(ci[["estimate"]]),
            lower = unname(ci[["lower"]]),
            upper = unname(ci[["upper"]]),
            source_detail = if (is.finite(under_one)) {
              glue::glue("{source_summary$death_source}; <1 ano separado ({round(under_one)} óbitos)")
            } else {
              glue::glue("{source_summary$death_source}; <1 ano não separado")
            }
          )
        },
        list(value = NA_real_, lower = NA_real_, upper = NA_real_, source_detail = "N/D")
      )

      # Pooling is presented as a moving average, and rates already behave that
      # way: they divide by pooled person-years, so they stay per-year. Counts
      # do not - summing three years of deaths triples the number, which reads
      # as a tripling of mortality rather than as smoothing, and cannot be
      # compared with an unpooled value. Averaging them over the contributing
      # years puts every metric on the same per-year footing. Intervals are
      # scaled with the estimate, so they remain the interval of the average.
      #
      # Infant deaths are the exception and stay a window total - see
      # annualised_metric_ids for why.
      annualise <- if (metric_id %in% annualised_metric_ids && length(years) > 1L) {
        length(years)
      } else {
        1L
      }

      tibble(
        location = area_spec$label,
        cause = cause,
        metric_id = metric_id,
        metric = metric_label,
        period = period_label,
        n_years = length(years),
        value = metric_values$value / annualise,
        lower = metric_values$lower / annualise,
        upper = metric_values$upper / annualise,
        # Marks a value that should not be read at face value. Only the infant
        # rate sets one today; every other metric carries an empty string so the
        # column is always present and the table can paste it unconditionally.
        flag = if (is.null(metric_values$flag)) "" else as.character(metric_values$flag),
        source_detail = as.character(metric_values$source_detail)
      )
    }))
  }

  round_annual_metric_value <- function(metric_id, value) {
    dplyr::case_when(
      metric_id %in% count_metric_ids ~ round(value, 0),
      metric_id %in% c("proportional", "smr", "isr") ~ round(value, 1),
      TRUE ~ round(value, 2)
    )
  }

  format_annual_metric_interval <- function(metric_id, value, lower, upper) {
    if (!is.finite(value)) {
      return("N/D")
    }

    value <- round_annual_metric_value(metric_id, value)
    lower <- round_annual_metric_value(metric_id, lower)
    upper <- round_annual_metric_value(metric_id, upper)

    if (!is.finite(lower) || !is.finite(upper)) {
      return(as.character(value))
    }

    glue::glue("{value} ({lower}; {upper})")
  }

  get_annual_sort_location <- function(metric_long) {
    unique(metric_long$sort_location)[[1]]
  }

  get_annual_cause_order <- function(metric_long) {
    sort_location <- get_annual_sort_location(metric_long)

    metric_long %>%
      dplyr::filter(location == .env$sort_location) %>%
      dplyr::arrange(dplyr::desc(value)) %>%
      dplyr::pull(cause)
  }

  # Rows carry a `flag` column marking values that should not be read at face
  # value. Defaulted here so a frame built without one still renders.
  with_annual_flags <- function(metric_long) {
    if (!("flag" %in% names(metric_long))) {
      metric_long$flag <- ""
    }
    metric_long$flag <- ifelse(is.na(metric_long$flag), "", metric_long$flag)
    metric_long
  }

  # Which footnotes the current table needs. Only rendered when a flagged row is
  # actually on screen, so the note never appears next to a table it explains
  # nothing about.
  annual_metrics_footnotes <- function(metric_long, selected_metric) {
    shown <- with_annual_flags(metric_long) %>%
      dplyr::filter(metric_id == selected_metric)

    if (!any(nzchar(shown$flag))) {
      return(character(0))
    }

    infant_instability_footnote()
  }

  build_annual_metrics_table <- function(metric_long, selected_metric) {
    cause_order <- get_annual_cause_order(metric_long)
    sort_location <- get_annual_sort_location(metric_long)
    location_order <- unique(metric_long$location)

    with_annual_flags(metric_long) %>%
      dplyr::filter(metric_id == selected_metric) %>%
      dplyr::mutate(
        cause = factor(cause, levels = cause_order),
        # The mark is appended to the formatted value so it travels into the
        # CSV too: the caveat should follow the number rather than stay behind
        # in the browser.
        display_value = paste0(
          purrr::pmap_chr(
            list(metric_id, value, lower, upper),
            format_annual_metric_interval
          ),
          flag
        )
      ) %>%
      dplyr::arrange(cause) %>%
      dplyr::select(`Causa de Morte` = cause, location, display_value) %>%
      tidyr::pivot_wider(names_from = location, values_from = display_value) %>%
      dplyr::arrange(`Causa de Morte`) %>%
      dplyr::mutate(`Causa de Morte` = as.character(`Causa de Morte`)) %>%
      dplyr::select(`Causa de Morte`, dplyr::all_of(location_order))
  }

  build_annual_sources_table <- function(metric_long, selected_metric) {
    cause_order <- get_annual_cause_order(metric_long)
    location_order <- unique(metric_long$location)

    metric_long %>%
      dplyr::filter(metric_id == selected_metric) %>%
      dplyr::mutate(cause = factor(cause, levels = cause_order)) %>%
      dplyr::arrange(cause) %>%
      dplyr::select(`Causa de Morte` = cause, location, source_detail) %>%
      tidyr::pivot_wider(names_from = location, values_from = source_detail) %>%
      dplyr::arrange(`Causa de Morte`) %>%
      dplyr::mutate(`Causa de Morte` = as.character(`Causa de Morte`)) %>%
      dplyr::select(`Causa de Morte`, dplyr::all_of(location_order))
  }

  build_annual_metrics_plot <- function(metric_long, selected_metric) {
    cause_order <- get_annual_cause_order(metric_long)
    plot_df <- with_annual_flags(metric_long) %>%
      dplyr::filter(metric_id == selected_metric) %>%
      dplyr::mutate(
        cause = factor(cause, levels = rev(cause_order)),
        value = round_annual_metric_value(metric_id, value),
        lower = round_annual_metric_value(metric_id, lower),
        upper = round_annual_metric_value(metric_id, upper),
        # Bars carry the mark on the location label, since a chart has no cell
        # to append it to. Flagging the whole series is right: the denominator
        # is a property of the place and period, not of one cause.
        location = paste0(location, ifelse(nzchar(flag), paste0(" ", flag), ""))
      )

    validate(need(nrow(plot_df) > 0, "Seleccione pelo menos uma causa de morte para apresentar."))

    period_text <- if ("period" %in% names(plot_df)) {
      periods <- unique(as.character(plot_df$period))
      n_years <- if ("n_years" %in% names(plot_df)) max(plot_df$n_years, na.rm = TRUE) else 1L
      if (length(periods) == 1 && n_years > 1) {
        basis <- if (selected_metric %in% annualised_metric_ids) {
          "média anual"
        } else if (identical(selected_metric, infant_count_metric_id)) {
          "total do período"
        } else {
          "pessoas-ano"
        }
        glue::glue("{periods} ({n_years} anos agregados, {basis})")
      } else {
        paste(periods, collapse = ", ")
      }
    } else {
      NULL
    }

    plot <- ggplot(plot_df, aes(x = cause, y = value, fill = location)) +
      geom_col(position = position_dodge(width = 0.72), width = 0.65) +
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.72),
        width = 0.2,
        linewidth = 0.4,
        na.rm = TRUE
      )

    # An SMR is read against its reference, so draw the reference line rather
    # than leaving the reader to judge 100 by eye.
    if (identical(selected_metric, "smr")) {
      plot <- plot +
        geom_hline(yintercept = 100, linetype = "dashed", colour = "#0b2e4f", linewidth = 0.5)
    }

    # Downloaded as a PNG, the chart leaves the app and loses the footnote under
    # the table, so the caveat is drawn onto it.
    caption_text <- if (any(nzchar(plot_df$flag))) {
      paste(strwrap(infant_instability_footnote(), width = 95), collapse = "\n")
    } else {
      NULL
    }

    plot +
      coord_flip() +
      labs(
        title = paste(unique(plot_df$metric), "por causa de morte"),
        subtitle = period_text,
        caption = caption_text,
        x = NULL,
        y = unique(plot_df$metric)
      ) +
      scale_fill_brewer(palette = "Set2", name = "Local") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0, size = 8, colour = "#4a4a4a"),
        legend.position = "bottom"
      )
  }
  
  # -------------------------
  # Observed Mortality
  # -------------------------
  observed_metric_bundle <- eventReactive(input$go_rates, {
    token <- isolate(cancel_seq$rates)
    year_order <- isolate(year_load_state$rates_order)
    query_spec <- make_query_spec(
      input$area,
      input$area_label,
      input$cause,
      input$sex,
      input$data_source
    )

    shiny::withProgress(message = get_data_source_progress_message(query_spec$data_source), value = 0, {
      load_metric_bundle(
        query_spec,
        "rates",
        token,
        year_range = input$years_import,
        year_order = year_order
      )
    })
  })

  # Keep the final population/rate filter reactive so the observed tab can reuse
  # the already-prepared metric bundle without reloading INE data.
  observed_series_spec <- reactive({
    req(input$go_rates > 0)
    make_series_spec(
      query_spec = observed_metric_bundle()$query_spec,
      population = input$population,
      rate_type = input$rate_type,
      year_range = input$years_import
    )
  })

  # `observed_history()` is the final shared series object used by the observed
  # plot/table. It is rebuilt from the cached metric bundle, not from raw data.
  observed_history <- reactive({
    req(input$go_rates > 0)
    build_historical_series(observed_metric_bundle(), observed_series_spec())
  })

  observed_summary <- reactive({
    dat <- observed_history()
    df <- dat$series
    req(nrow(df) > 0)
    first_row <- dplyr::slice(df, 1)
    last_row <- dplyr::slice(df, nrow(df))
    source_summary <- dat$metric_bundle$source_summary

    absolute_change <- last_row$value - first_row$value
    percent_change <- if (isTRUE(all.equal(first_row$value, 0))) {
      NA_real_
    } else {
      (absolute_change / first_row$value) * 100
    }

    tibble(
      Métrica = c(
        "Fonte de dados",
        "Fonte de população",
        "Fonte de óbitos",
        "Último ano",
        "Última taxa (IC 95%)",
        "Variação absoluta no período observado",
        "Variação percentual no período observado"
      ),
      Valor = c(
        get_data_source_label(dat$spec$data_source),
        source_summary$population_source,
        source_summary$death_source,
        as.character(last_row$year),
        glue::glue("{round(last_row$value, 2)} ({round(last_row$lower, 2)}; {round(last_row$upper, 2)})"),
        sprintf("%+.2f", absolute_change),
        if (is.na(percent_change)) "N/D" else sprintf("%+.2f%%", percent_change)
      )
    )
  })

  build_observed_rate_plot <- function(dat) {
    df <- dat$series

    ggplot(df, aes(x = year, group = 1)) +
      geom_ribbon(
        aes(ymin = lower, ymax = upper),
        fill = "grey80", alpha = 0.4
      ) +
      geom_line(aes(y = value), linewidth = 1) +
      geom_point(aes(y = value), size = 2) +
      labs(
        title = paste(dat$area_label, "-", dat$spec$cause, dat$spec$sex, "(", dat$spec$population, ")"),
        x = "Ano",
        y = dat$y_label,
        caption = "A zona sombreada mostra o intervalo de confiança de 95% (a incerteza em torno da taxa)."
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0)
      )
  }

  build_observed_rate_table <- function(dat) {
    dat$series %>%
      dplyr::transmute(
        Ano = year,
        `Taxa (Intervalo de Confiança 95%)` = glue::glue(
          "{round(value, 2)} ({round(lower, 2)}; {round(upper, 2)})"
        )
      )
  }

  output$rateSummaryTable <- renderTable({
    observed_summary()
  }, bordered = TRUE, spacing = "s")

  output$downloadRateSummaryCSV <- downloadHandler(
    filename = function() {
      paste0("observed_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(observed_summary(), file)
    }
  )
  
  # Wrap a ggplot as an interactive plotly widget: hovering a point shows the
  # year and value, plus zoom/pan. The ggplot object is unchanged, so the PNG
  # download handlers keep using it directly. Note that plotly ignores
  # labs(caption), so the plain-language captions are repeated as static text
  # under these plots in the UI.
  to_interactive_plot <- function(p, tooltip = c("x", "y")) {
    gp <- plotly::ggplotly(p, tooltip = tooltip)
    gp <- plotly::layout(gp, hovermode = "closest", margin = list(t = 60))
    plotly::config(gp, displaylogo = FALSE, modeBarButtonsToRemove = list("lasso2d", "select2d"))
  }

  # Annual percentage change over the chosen span of the observed series.
  observed_trend <- reactive({
    dat <- observed_history()
    series <- dat$series[is.finite(dat$series$value), , drop = FALSE]
    span <- input$observed_trend_span %||% "10"
    if (!identical(span, "all")) series <- series[series$year > max(series$year) - as.integer(span), , drop = FALSE]
    trend_apc(series$year, series$value, series$lower, series$upper,
              max_joinpoints = as.integer(input$observed_trend_joinpoints %||% "0"))
  })

  build_trend_table <- function(trend) {
    rows <- trend$segments
    if (nrow(rows) > 1) rows <- dplyr::bind_rows(rows, trend$aapc)
    tibble::tibble(
      `Período` = paste0(rows$from, "-", rows$to, ifelse(seq_len(nrow(rows)) > nrow(trend$segments), " (média, AAPC)", "")),
      `Variação anual (IC 95%)` = vapply(seq_len(nrow(rows)), function(i) trend_label(rows[i, ]), character(1)),
      Leitura = vapply(seq_len(nrow(rows)), function(i) trend_reading(rows[i, ]), character(1))
    )
  }

  output$rateTrendTable <- renderTable({
    trend <- observed_trend()
    validate(need(!is.null(trend), "São precisos pelo menos 5 anos para estimar a tendência."))
    build_trend_table(trend)
  }, bordered = TRUE, spacing = "s")

  output$rateTrendNote <- renderUI({
    trend <- tryCatch(observed_trend(), error = function(e) NULL)
    req(trend)
    note <- trend_break_note(trend$fitted$year, trend$joinpoints)
    if (is.null(note)) return(NULL)
    tags$p(class = "text-warning", note)
  })

  output$ratePlot <- plotly::renderPlotly({
    p <- build_observed_rate_plot(observed_history())
    trend <- tryCatch(observed_trend(), error = function(e) NULL)
    if (!is.null(trend)) {
      p <- p + geom_line(data = trend$fitted, aes(x = year, y = fitted), inherit.aes = FALSE,
                         linetype = "dashed", colour = "#eb6834", linewidth = 0.8)
    }
    to_interactive_plot(p)
  })

  output$downloadRatePlot <- downloadHandler(
    filename = function() {
      paste0("observed_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(file, build_observed_rate_plot(observed_history()))
    }
  )
  
  output$rateTable <- renderTable({
    build_observed_rate_table(observed_history())
  }, sanitize.text.function = identity)

  output$downloadRateTableCSV <- downloadHandler(
    filename = function() {
      paste0("observed_series_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(build_observed_rate_table(observed_history()), file)
    }
  )
  
  # -------------------------
  # Beginner Forecasting
  # -------------------------
  beginner_training_history <- eventReactive(input$go_beginner_forecast, {
    query_spec <- make_query_spec(
      input$beginner_area,
      input$beginner_area_label,
      input$beginner_cause,
      input$beginner_sex,
      input$beginner_data_source
    )
    build_guided_history(query_spec, isolate(cancel_seq$beginner), isolate(year_load_state$beginner_order))
  }, ignoreNULL = TRUE)

  # The guided forecast's training series for one selection.
  build_guided_history <- function(query_spec, token, year_order) {
    if (identical(input$beginner_indicator, "infant")) {
      return(build_infant_history(query_spec, input$beginner_years_fit,
                                  value_or_default(input$beginner_infant_window, "5"),
                                  value_or_default(input$beginner_infant_method, "annual")))
    }

    shiny::withProgress(message = get_data_source_progress_message(query_spec$data_source), value = 0, {
      metric_bundle <- load_metric_bundle(
        query_spec,
        "beginner",
        token,
        year_range = input$beginner_years_fit,
        year_order = year_order
      )
      training_range <- get_beginner_training_range(
        years = metric_bundle$years,
        year_range = input$beginner_years_fit
      )

      build_historical_series(
        metric_bundle = metric_bundle,
        series_spec = make_series_spec(
          query_spec = metric_bundle$query_spec,
          population = input$beginner_population,
          rate_type = input$beginner_rate_type,
          year_range = training_range
        )
      )
    })
  }

  # Adds the two averages of models (R/forecast_helpers.R) as further
  # candidates, and picks the forecast shown:
  #   - the combination of the three best, its members chosen at each
  #     validation step from the earlier steps only (needs the rolling
  #     validation);
  #   - the mean of six fixed models, the default projection (`method` "mean6").
  # Both are scored on the same one-step errors as the single models, so the
  # validation table compares them all. With `method` "selected" the forecast
  # shown is the candidate with the lowest error, as before October 2026.
  add_forecast_averages <- function(result, validation, method = "mean6") {
    steps <- validation$plot_data$predictions
    fitted <- get_successful_model_ids(result)
    rolling <- identical(validation$method_used, "rolling") && !is.null(steps) && nrow(steps) > 0
    setup <- get_transformation_setup(result$obs$value, result$transform_method)
    scale_denom <- mean(abs(diff(result$obs$value)), na.rm = TRUE)
    score <- function(pred, id) {
      compute_forecast_error_metrics(pred$actual, pred$predicted, scale_denom) %>%
        dplyr::mutate(Model = id) %>%
        dplyr::select(Model, ME, RMSE, MAE, MAPE, MASE)
    }
    new_fc <- list(); new_steps <- list(); new_metrics <- list()
    if (rolling && length(fitted) >= FORECAST_COMBINATION_SIZE) {
      initial <- result$accuracy %>% dplyr::filter(is.finite(RMSE)) %>% dplyr::arrange(RMSE)
      combo_steps <- combination_rolling_predictions(steps, initial$Model)
      members <- validation$metrics %>% dplyr::filter(Model %in% fitted, is.finite(RMSE)) %>% dplyr::arrange(RMSE)
      members <- utils::head(members$Model, FORECAST_COMBINATION_SIZE)
      new_fc$combo <- combine_model_forecasts(result$fc, members, setup$forward, setup$inverse, result$conf_level)
      new_steps$combo <- combo_steps
      new_metrics$combo <- score(combo_steps, FORECAST_COMBINATION_ID)
      result$combination_members <- members
    }
    mean_members <- intersect(FORECAST_MEAN_MEMBERS, fitted)
    if (length(mean_members) >= FORECAST_MEAN_MIN_MEMBERS) {
      new_fc$mean <- combine_model_forecasts(result$fc, mean_members, setup$forward, setup$inverse, result$conf_level,
                                             id = FORECAST_MEAN_ID)
      result$mean_members <- mean_members
      if (rolling) {
        mean_steps <- mean_rolling_predictions(steps, mean_members)
        new_steps$mean <- mean_steps
        new_metrics$mean <- score(mean_steps, FORECAST_MEAN_ID)
      }
    }
    result$fc <- dplyr::bind_rows(result$fc, dplyr::bind_rows(new_fc))
    result$fitted_models <- c(fitted, unique(dplyr::bind_rows(new_fc)$model))
    if (length(new_metrics) > 0) {
      validation$metrics <- dplyr::bind_rows(validation$metrics, dplyr::bind_rows(new_metrics))
      validation$plot_data$predictions <- dplyr::bind_rows(steps, dplyr::bind_rows(new_steps))
      validation$recommended_model <- choose_recommended_model(validation$metrics)
    }
    # The candidate with the lowest validation error, whatever is shown.
    result$validated_best <- validation$recommended_model
    use_mean <- identical(method, "mean6") && FORECAST_MEAN_ID %in% result$fitted_models
    result$projection_method <- if (use_mean) "mean6" else "selected"
    if (use_mean) validation$recommended_model <- FORECAST_MEAN_ID
    validation$ranking <- rank_model_metrics(validation$metrics, validation$recommended_model)
    list(result = result, validation = validation)
  }

  # Fit, validate and combine for one training series.
  run_guided_forecast <- function(training_history, token) {
    guided_result <- run_forecast_models(
      history = training_history,
      model_ids = unname(forecast_model_choices),
      horizon = input$beginner_horizon,
      kind = "beginner",
      token = token
    )
    validation <- evaluate_model_selection(
      base_result = guided_result,
      method = value_or_default(input$beginner_validation, "rolling"),
      test_fraction = value_or_default(input$beginner_test_pct, 25) / 100
    )
    combined <- add_forecast_averages(guided_result, validation, value_or_default(input$beginner_method, "mean6"))
    guided_result <- combined$result
    validation <- combined$validation
    if (!is.null(validation$recommended_model)) {
      guided_result$recommended_model <- validation$recommended_model
    }
    aggregate_infant_result(c(
      guided_result,
      list(
        full_history = training_history,
        horizon = input$beginner_horizon,
        mode = input$beginner_mode,
        training_label = get_beginner_training_label(training_history),
        validation = validation
      )
    ), training_history)
  }

  beginner_forecast <- eventReactive(input$go_beginner_forecast, {
    training_history <- beginner_training_history()
    shiny::withProgress(message = "A ajustar e validar modelos...", value = 0.3, {
      run_guided_forecast(training_history, isolate(cancel_seq$beginner))
    })
  }, ignoreNULL = TRUE)

  # Several areas, each forecast on its own, under one model policy:
  #   "mean6"  the mean of six models for every area (the default) - the same
  #            fixed method everywhere, so differences come from the data;
  #   "combo"  the combination of each area's three best models;
  #   "common" the single model with the lowest mean relative error across
  #            the areas;
  #   "best"   each area's own recommended model.
  GUIDED_MAX_AREAS <- 12L
  beginner_multi_forecast <- eventReactive(input$go_beginner_forecast, {
    areas <- unique(as.character(input$beginner_area))
    validate(need(length(areas) >= 2, "Seleccione pelo menos dois locais para os comparar."))
    validate(need(length(areas) <= GUIDED_MAX_AREAS, paste("Seleccione no máximo", GUIDED_MAX_AREAS, "locais.")))
    token <- isolate(cancel_seq$beginner)
    order <- isolate(year_load_state$beginner_order)
    results <- list()
    shiny::withProgress(message = "A ajustar e validar modelos, área a área...", value = 0, {
      for (area in areas) {
        query_spec <- make_query_spec(area, "", input$beginner_cause, input$beginner_sex, input$beginner_data_source)
        history <- build_guided_history(query_spec, token, order)
        results[[area]] <- tryCatch(run_guided_forecast(history, token), error = function(e) {
          if (inherits(e, "shiny.silent.error")) stop(e)
          e
        })
        incProgress(1 / length(areas))
      }
    })
    # An area whose models all failed validation has no recommended model; it
    # is reported with the failures rather than stopping the comparison.
    failed <- vapply(results, function(r) inherits(r, "error") || is.null(r$recommended_model), logical(1))
    if (any(failed)) {
      notify_data(paste("Sem previsão para:", paste(names(results)[failed], collapse = ", ")), type = "warning", duration = 15)
    }
    results <- results[!failed]
    validate(need(length(results) > 0, "Nenhuma área pôde ser projectada."))

    policy <- value_or_default(input$beginner_model_policy, "mean6")
    chosen <- vapply(results, function(r) r$validated_best %||% r$recommended_model, character(1))
    if (identical(policy, "mean6")) {
      has_mean <- vapply(results, function(r) FORECAST_MEAN_ID %in% r$fc$model, logical(1))
      chosen[has_mean] <- FORECAST_MEAN_ID
    } else if (identical(policy, "combo")) {
      has_combo <- vapply(results, function(r) FORECAST_COMBINATION_ID %in% r$fc$model, logical(1))
      chosen[has_combo] <- FORECAST_COMBINATION_ID
    } else if (identical(policy, "common")) {
      common <- choose_common_model(
        lapply(results, function(r) r$validation$metrics),
        vapply(results, function(r) mean(r$obs$value, na.rm = TRUE), numeric(1))
      )
      if (!is.null(common)) chosen[] <- common
    }
    list(results = results, chosen = chosen, policy = policy)
  }, ignoreNULL = TRUE)

  beginner_multi_table <- reactive({
    multi <- beginner_multi_forecast()
    dplyr::bind_rows(lapply(names(multi$results), function(area) {
      r <- multi$results[[area]]
      model <- multi$chosen[[area]]
      last <- utils::tail(r$full_history$series, 1)
      end <- utils::tail(r$fc[r$fc$model == model, , drop = FALSE], 1)
      tibble::tibble(
        Local = area,
        Modelo = get_model_labels(model),
        `Último observado` = paste0(round(last$value, 2), " (", last$year, ")"),
        `Projecção (IC)` = paste0(round(end$mean, 2), " (", round(end$lower, 2), "; ", round(end$upper, 2), ") em ", end$year),
        `Variação` = sprintf("%+.1f%%", (end$mean / last$value - 1) * 100),
        `Menor erro no local` = get_model_labels(r$validated_best %||% r$recommended_model)
      )
    }))
  })

  output$beginnerMultiTable <- renderTable({
    req(input$go_beginner_forecast > 0)
    beginner_multi_table()
  }, bordered = TRUE, spacing = "s")

  output$downloadBeginnerMultiCSV <- downloadHandler(
    filename = function() paste0("previsao_por_local_", Sys.Date(), ".csv"),
    content = function(file) write_csv_utf8(beginner_multi_table(), file)
  )

  output$beginnerMultiNote <- renderUI({
    req(input$go_beginner_forecast > 0)
    multi <- beginner_multi_forecast()
    policy <- c(mean6 = "a média de 6 modelos (ARIMA, ETS, Theta, TBATS, Holt e Holt amortecido), igual em todos",
                combo = "a combinação dos três modelos que melhor previram os anos de teste em cada local",
                common = paste0("o mesmo modelo em todos (", get_model_labels(unique(multi$chosen))[[1]], ", o de menor erro relativo médio)"),
                best = "o modelo com menor erro de validação em cada local")[[multi$policy]]
    helpText(paste0(
      "Cada local é projectado em separado com ", policy, ". Linha cheia: observado; tracejada e faixa: projecção e ",
      "intervalo. A última coluna diz qual foi o modelo com menor erro de validação nesse local."
    ))
  })

  output$beginnerMultiPlot <- plotly::renderPlotly({
    req(input$go_beginner_forecast > 0)
    multi <- beginner_multi_forecast()
    observed <- dplyr::bind_rows(lapply(names(multi$results), function(a) dplyr::mutate(multi$results[[a]]$full_history$series, local = a)))
    projected <- dplyr::bind_rows(lapply(names(multi$results), function(a) {
      r <- multi$results[[a]]
      dplyr::mutate(r$fc[r$fc$model == multi$chosen[[a]], , drop = FALSE], local = a)
    }))
    p <- ggplot() +
      geom_ribbon(data = projected, aes(x = year, ymin = lower, ymax = upper), fill = "#8fb8de", alpha = 0.35) +
      geom_line(data = observed, aes(x = year, y = value), colour = "#1f4e79", linewidth = 0.8) +
      geom_line(data = projected, aes(x = year, y = mean), colour = "#0b2e4f", linetype = "dashed", linewidth = 0.9) +
      facet_wrap(~factor(local, levels = names(multi$results))) +
      labs(x = "Ano", y = multi$results[[1]]$history$y_label) +
      theme_minimal()
    to_interactive_plot(p)
  })

  output$beginnerForecastPlot <- plotly::renderPlotly({
    req(input$go_beginner_forecast > 0)
    to_interactive_plot(build_beginner_forecast_plot(beginner_forecast()))
  })

  output$downloadBeginnerForecastPlot <- downloadHandler(
    filename = function() {
      paste0("guided_forecast_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(file, build_beginner_forecast_plot(beginner_forecast()))
    }
  )

  output$beginnerForecastWarnings <- renderUI({
    req(input$go_beginner_forecast > 0)
    build_forecast_warning_ui(beginner_forecast())
  })

  output$beginnerForecastSummary <- renderUI({
    req(input$go_beginner_forecast > 0)
    build_beginner_summary_ui(beginner_forecast())
  })

  output$beginnerForecastReliability <- renderUI({
    req(input$go_beginner_forecast > 0)
    build_beginner_reliability_ui(beginner_forecast())
  })

  beginner_forecast_table <- reactive({
    dat <- beginner_forecast()
    validate_successful_forecast(dat)

    if (identical(dat$mode, "recommended")) {
      build_forecast_display_table(
        dat,
        view_mode = "single",
        selected_model = dat$recommended_model
      )
    } else {
      build_forecast_display_table(dat, view_mode = "compare")
    }
  })

  # On screen the last observed years and the forecast; the CSV keeps the series.
  recent_forecast_rows <- function(table, dat, years = FORECAST_TABLE_RECENT_YEARS) {
    last <- suppressWarnings(max(as.integer(dat$obs$year), na.rm = TRUE))
    if (!is.finite(last)) return(table)
    table[table$Ano > last - years, , drop = FALSE]
  }

  output$beginnerForecastTable <- renderTable({
    req(input$go_beginner_forecast > 0)
    recent_forecast_rows(beginner_forecast_table(), beginner_forecast())
  }, sanitize.text.function = identity)

  output$downloadBeginnerForecastCSV <- downloadHandler(
    filename = function() {
      paste0("guided_forecast_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(beginner_forecast_table(), file)
    }
  )

  # -------------------------
  # Annual Metrics
  # -------------------------
  annual_metrics_long <- eventReactive(input$go_annual_metrics, {
    token <- isolate(cancel_seq$annual)
    year <- as.integer(input$annual_year)
    selected_causes <- input$annual_cause
    selected_metric <- input$annual_metric

    validate(
      need(length(selected_causes) > 0, "Seleccione pelo menos uma causa de morte."),
      need(length(selected_metric) == 1, "Seleccione uma métrica.")
    )

    pooling_window <- normalize_pooling_window(input$annual_pooling)
    region_mode <- default_region_mode()
    reference_areas <- value_or_default(input$annual_smr_reference, "Portugal")
    pooled_years <- get_pooled_years(year, pooling_window, selected_metric)

    area_specs <- get_annual_area_specs(
      selected_areas = input$annual_area,
      custom_label = input$annual_area_label,
      region_mode = region_mode,
      available_areas = get_available_areas(input$annual_data_source)
    )
    annual_causes_needed <- unique(c(
      selected_causes,
      if (selected_metric %in% "proportional") "Todas as causas de morte" else character(0)
    ))
    annual_areas_needed <- sort(unique(unlist(
      c(lapply(area_specs, `[[`, "areas"),
        if (selected_metric %in% indirect_metric_ids) list(reference_areas) else NULL),
      use.names = FALSE
    )))

    # Surface any incomplete municipal rebuild, and the NUTS-2024 break when
    # the user chose to keep INE's own regional rows.
    for (message_text in unique(unlist(lapply(area_specs, `[[`, "warnings"), use.names = FALSE))) {
      notify_data(message_text, type = "warning", duration = 15)
    }
    compared <- setdiff(annual_comparator_areas(input$annual_area), "Portugal")
    vintage_warning <- region_vintage_warning(
      areas = unique(c(compared, as.character(input$annual_area))),
      years = pooled_years,
      region_mode = region_mode
    )
    if (!is.null(vintage_warning)) {
      notify_data(vintage_warning, type = "warning", duration = 20)
    }

    if (!selected_metric %in% infant_metric_ids) {
      notify_region_source_warnings(
        selected_areas = unique(c(compared, as.character(input$annual_area))),
        expanded_regions = unique(unlist(lapply(area_specs, `[[`, "expanded"), use.names = FALSE)),
        causes = selected_causes,
        years = pooled_years
      )
    }

    # A pooled window straddling 2020/2021 mixes the two population series.
    revision_warning <- population_revision_warning(pooled_years, selected_metric)
    if (!is.null(revision_warning)) {
      notify_data(revision_warning, type = "warning", duration = 20, id = "population-revision")
    }

    # The infant metrics are bounded by their own datasets, not by population:
    # the rate needs published live births (1995-), the count only the under-1
    # death archive. Years outside those are refused for that metric alone.
    if (selected_metric %in% infant_metric_ids) {
      infant_gap <- infant_gap_message(pooled_years, selected_metric)
      validate(need(
        length(intersect(pooled_years, infant_metric_years(selected_metric))) > 0,
        if (is.null(infant_gap)) "Sem dados de mortalidade infantil para o período seleccionado." else infant_gap
      ))
      if (!is.null(infant_gap)) {
        notify_data(infant_gap, type = "warning", duration = 15)
      }

      # A year published only as an all-cause, both-sexes total cannot answer a
      # cause- or sex-specific request; say which, rather than reporting zero.
      detail_cause <- if (length(selected_causes) == 1) selected_causes else "Todas as causas de morte"
      detail_gap <- infant_detail_message(
        pooled_years,
        cause = detail_cause,
        sex = input$annual_sex
      )
      if (!is.null(detail_gap)) {
        # Refuse outright when no year in the window can answer; warn when only
        # part of it cannot, since the rest still produces a usable figure.
        validate(need(
          length(infant_detail_years(pooled_years, detail_cause, input$annual_sex)) > 0,
          detail_gap
        ))
        notify_data(detail_gap, type = "warning", duration = 20)
      }

      undercount <- infant_undercount_message(
        pooled_years,
        areas = c(compared, as.character(input$annual_area)),
        cause = detail_cause,
        sex = input$annual_sex,
        municipalities = active_nuts_lookup()$municipality
      )
      if (!is.null(undercount)) {
        notify_data(undercount, type = "warning", duration = 20)
      }
    }

    notify_snapshot_request_warnings(
      years = pooled_years,
      areas = annual_areas_needed,
      causes = annual_causes_needed,
      data_source = input$annual_data_source,
      include_population = selected_metric %in% annual_metrics_needing_population,
      include_deaths = TRUE
    )

    shiny::withProgress(message = get_data_source_progress_message(input$annual_data_source, context = "annual"), value = 0, {
      with_data_load_cancel_checker(
        cancel_checker = function() !identical(cancel_seq$annual, token),
        {
          out <- dplyr::bind_rows(lapply(seq_along(area_specs), function(i) {
            incProgress(0.2 / length(area_specs))
            calculate_annual_metric_values(
              area_spec = area_specs[[i]],
              causes = selected_causes,
              metric_id = selected_metric,
              sex = input$annual_sex,
              year = year,
              token = token,
              data_source = input$annual_data_source,
              pooling_window = pooling_window,
              reference_areas = reference_areas
            )
          })) %>%
            # Causes are sorted by the local area, the last column.
            dplyr::mutate(sort_location = area_specs[[length(area_specs)]]$label)

          incProgress(0.8)
          out
        }
      )
    })
  }, ignoreNULL = TRUE)

  annual_metrics_table <- reactive({
    req(input$go_annual_metrics > 0)
    build_annual_metrics_table(annual_metrics_long(), input$annual_metric)
  })

  annual_sources_table <- reactive({
    req(input$go_annual_metrics > 0)
    build_annual_sources_table(annual_metrics_long(), input$annual_metric)
  })

  annual_metrics_plot <- reactive({
    req(input$go_annual_metrics > 0)
    build_annual_metrics_plot(annual_metrics_long(), input$annual_metric)
  })

  output$annualMetricsTable <- renderTable({
    annual_metrics_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 2)

  output$annualMetricsFootnote <- renderUI({
    req(input$go_annual_metrics > 0)
    notes <- annual_metrics_footnotes(annual_metrics_long(), input$annual_metric)

    if (length(notes) == 0) {
      return(NULL)
    }

    tagList(lapply(notes, function(note) helpText(note)))
  })

  output$annualSourcesTable <- renderTable({
    annual_sources_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$annualMetricsPlot <- renderPlot({
    annual_metrics_plot()
  })

  output$downloadAnnualMetricsCSV <- downloadHandler(
    filename = function() {
      paste0("annual_", input$annual_metric, "_", input$annual_year, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(annual_metrics_table(), file)
    }
  )

  output$downloadAnnualSourcesCSV <- downloadHandler(
    filename = function() {
      paste0("annual_sources_", input$annual_metric, "_", input$annual_year, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(annual_sources_table(), file)
    }
  )

  output$downloadAnnualMetricsPlot <- downloadHandler(
    filename = function() {
      paste0("annual_", input$annual_metric, "_", input$annual_year, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(file, annual_metrics_plot())
    }
  )

  # -------------------------
  # Avoidable mortality
  # -------------------------
  # One load serves the whole tab: the two lists, the deliberately unresolved
  # set and the all-cause denominator, all under 75. Loading them together is
  # what lets the reconciliation add up to the whole - computing the groups
  # separately would leave no way to show what is missing from them.
  avoidable_data <- eventReactive(input$go_avoidable, {
    token <- isolate(cancel_seq$avoidable)
    year <- as.integer(input$avoidable_year)
    pooling_window <- normalize_pooling_window(input$avoidable_pooling)
    years <- get_pooled_years(year, pooling_window, "dsr")
    # Avoidable mortality needs the causes, which the latest year does not have.
    validate(need(!year %in% calendar_death_years, calendar_gap_message(year)))
    years <- setdiff(years, calendar_death_years)

    validate(need(length(input$avoidable_area) > 0, "Seleccione pelo menos um local."))
    validate(need(
      length(years_without_population(years)) == 0,
      population_gap_message(years, "a mortalidade evitável")
    ))

    # Not get_annual_area_specs(): that one adds the comparator columns. Here
    # every area is an ordinary selection, and Portugal is the default.
    selected_areas <- sort(unique(as.character(input$avoidable_area)))

    mismatch <- vintage_mismatch_message(
      selected_areas,
      active_nuts_vintage(),
      active_nuts_lookup()
    )
    validate(need(is.null(mismatch), mismatch))

    notify_uls_split(selected_areas)
    overlap <- overlapping_selection_warning(selected_areas, lookup = active_nuts_lookup())
    if (!is.null(overlap)) {
      notify_data(overlap, type = "warning", duration = 20)
    }

    area_spec <- make_annual_area_spec(
      label = get_selection_label(selected_areas, input$avoidable_area_label),
      areas = selected_areas,
      region_mode = default_region_mode(),
      available_areas = get_available_areas(input$avoidable_data_source)
    )

    for (message_text in area_spec$warnings) {
      notify_data(message_text, type = "warning", duration = 15)
    }
    revision_warning <- population_revision_warning(years, "dsr")
    if (!is.null(revision_warning)) {
      showNotification(revision_warning, type = "warning", duration = 20, id = "population-revision")
    }

    causes <- avoidable_required_causes()
    notify_region_source_warnings(selected_areas, area_spec$expanded, causes, years)

    shiny::withProgress(message = "A carregar mortalidade evitável...", value = 0, {
      loaded <- dplyr::bind_rows(lapply(seq_along(causes), function(i) {
        abort_if_cancelled("avoidable", token)
        incProgress(1 / length(causes))
        apply_region_source(
          get_data_for(
            area = area_spec$areas,
            cause = causes[[i]],
            years = years,
            year_order = "desc",
            data_source = input$avoidable_data_source,
            population_basis = active_population_basis()
          )$full,
          area_spec$expanded
        ) %>%
          dplyr::filter(sex == input$avoidable_sex)
      }))

      validate(need(nrow(loaded) > 0, "Sem dados para a selecção."))

      list(
        data = loaded,
        years = years,
        period = make_period_label(min(years), max(years)),
        label = area_spec$label,
        n_years = length(years)
      )
    })
  }, ignoreNULL = TRUE)

  # Deaths, share of under-75 mortality, and the under-75 standardised rate for
  # each part. Counts are annualised over a pooled window, as elsewhere; the
  # standardised rate divides by pooled person-years and is already per year.
  avoidable_table <- reactive({
    req(input$go_avoidable > 0)
    dat <- avoidable_data()
    breakdown <- avoidable_breakdown(dat$data)

    rate_for <- function(group) {
      avoidable_group_dsr(
        df = dat$data,
        causes = switch(
          group,
          preventable = AVOIDABLE_PREVENTABLE,
          treatable = AVOIDABLE_TREATABLE,
          avoidable = avoidable_group_causes("avoidable"),
          unresolved = AVOIDABLE_UNRESOLVED$cause,
          total = "Todas as causas de morte",
          character(0)
        ),
        year = dat$years[[1]]
      )
    }

    annualise <- if (dat$n_years > 1L) dat$n_years else 1L

    breakdown %>%
      dplyr::mutate(
        dsr = vapply(group, rate_for, numeric(1)),
        Óbitos = round(.data$deaths / annualise, 0),
        `% dos óbitos < 75 anos` = round(.data$share, 1),
        `Taxa padronizada < 75 (por 100.000)` = round(.data$dsr, 1)
      ) %>%
      dplyr::select(
        Grupo = label,
        Óbitos,
        `% dos óbitos < 75 anos`,
        `Taxa padronizada < 75 (por 100.000)`
      )
  })

  avoidable_causes_table <- reactive({
    dplyr::bind_rows(
      tibble(Grupo = "Prevenível", `Causa de morte` = AVOIDABLE_PREVENTABLE, Nota = ""),
      tibble(Grupo = "Tratável", `Causa de morte` = AVOIDABLE_TREATABLE, Nota = ""),
      tibble(
        Grupo = "Excluída *",
        `Causa de morte` = AVOIDABLE_UNRESOLVED$cause,
        Nota = AVOIDABLE_UNRESOLVED$reason
      )
    )
  })

  avoidable_plot <- reactive({
    req(input$go_avoidable > 0)
    dat <- avoidable_data()
    plot_df <- avoidable_table() %>%
      dplyr::filter(Grupo %in% c("Prevenível", "Tratável", "Não classificada (por resolver) *")) %>%
      dplyr::mutate(Grupo = factor(Grupo, levels = rev(unique(Grupo))))

    ggplot(plot_df, aes(x = Grupo, y = Óbitos, fill = Grupo)) +
      geom_col(width = 0.6) +
      coord_flip() +
      labs(
        title = glue::glue("Mortalidade evitável < 75 anos - {dat$label}"),
        subtitle = glue::glue("{dat$period}{if (dat$n_years > 1) ' (media anual)' else ''}"),
        x = NULL,
        y = "Óbitos"
      ) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })

  output$avoidableTable <- renderTable({
    avoidable_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 1)

  output$avoidableNote <- renderUI({
    req(input$go_avoidable > 0)
    helpText(avoidable_unresolved_note())
  })

  output$avoidableCausesTable <- renderTable({
    avoidable_causes_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$avoidablePlot <- renderPlot({
    avoidable_plot()
  })

  output$downloadAvoidableCSV <- downloadHandler(
    filename = function() paste0("evitavel_", input$avoidable_year, "_", Sys.Date(), ".csv"),
    content = function(file) write_csv_utf8(avoidable_table(), file)
  )

  output$downloadAvoidableCausesCSV <- downloadHandler(
    filename = function() paste0("evitavel_causas_", Sys.Date(), ".csv"),
    content = function(file) write_csv_utf8(avoidable_causes_table(), file)
  )

  output$downloadAvoidablePlot <- downloadHandler(
    filename = function() paste0("evitavel_", input$avoidable_year, "_", Sys.Date(), ".png"),
    content = function(file) save_ggplot_png(file, avoidable_plot())
  )

  # -------------------------
  # Planning indicators
  # -------------------------
  # Built around one location. Its comparators are the areas that contain it,
  # offered only for indicators that do not depend on the size of an area.
  planning_data_date <- reactive({
    tryCatch(latest_import_date(read_import_log(app_data_root(get_snapshot_dir()))), error = function(e) NA_character_)
  })

  # The education age restricts I24 everywhere in the tab; the label says so.
  planning_education_age <- reactive({
    age <- suppressWarnings(as.integer(input$planning_education_age))
    if (length(age) == 0 || is.na(age)) 0L else age
  })

  planning_spec <- reactive({
    req(input$planning_indicator)
    spec <- planning_indicator_spec(input$planning_indicator)
    spec$label <- planning_indicator_label(spec$id, planning_education_age(), ref = FALSE, star = FALSE)
    spec
  })

  # Whole municipality or parish weights for the six ULS that share one.
  planning_split_mode <- reactive({
    mode <- input$planning_split_mode %||% PLANNING_DEFAULT_SPLIT_MODE
    if (identical(mode, "parish") && !planning_parish_available()) PLANNING_DEFAULT_SPLIT_MODE else mode
  })

  # Portugal as INE's published total or as the sum of its municipalities, for
  # every comparator, the significance marks, the ranking and the funnel.
  planning_benchmark <- reactive(planning_portugal_area(input$planning_portugal %||% "published"))

  planning_swap_portugal <- function(areas) {
    areas$area[areas$area == "Portugal"] <- planning_benchmark()
    areas
  }

  planning_location <- reactive({
    req(input$planning_area)
    area <- as.character(input$planning_area)[[1]]
    mismatch <- vintage_mismatch_message(area, active_nuts_vintage(), active_nuts_lookup())
    validate(need(is.null(mismatch), mismatch))
    area
  })

  planning_available_comparators <- reactive({
    planning_comparators(planning_location(), active_nuts_lookup())
  })

  # Comparators start all selected for a new location and keep the user's choice
  # while only the indicator or years change.
  planning_comparator_state <- reactiveValues(area = NULL, selected = NULL)

  observeEvent(input$planning_comparators, {
    planning_comparator_state$selected <- input$planning_comparators
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output$planningComparatorsUI <- renderUI({
    spec <- planning_spec()
    comparators <- planning_available_comparators()
    area <- planning_location()

    if (!identical(isolate(planning_comparator_state$area), area)) {
      planning_comparator_state$area <- area
      planning_comparator_state$selected <- comparators$area
    }
    if (!isTRUE(spec$comparable)) {
      return(helpText(
        tags$b("Sem comparadores. "),
        spec$label, " é uma contagem absoluta: depende do tamanho da área, ",
        "por isso não se compara com áreas maiores."
      ))
    }
    if (nrow(comparators) == 0) {
      return(helpText("Portugal não tem níveis acima para comparar."))
    }
    checkboxGroupInput(
      "planning_comparators", "Comparar com:",
      choices = stats::setNames(comparators$area, paste0(comparators$level, ": ", comparators$area)),
      selected = isolate(intersect(planning_comparator_state$selected, comparators$area))
    )
  })

  planning_areas <- reactive({
    local <- planning_location()
    comparators <- planning_available_comparators()
    chosen <- if (isTRUE(planning_spec()$comparable)) {
      intersect(input$planning_comparators, comparators$area)
    } else {
      character(0)
    }
    planning_swap_portugal(dplyr::bind_rows(
      tibble::tibble(area = local, level = "Local"),
      comparators[comparators$area %in% chosen, c("area", "level")]
    ))
  })

  planning_year_range <- reactive({
    req(input$planning_years)
    seq.int(as.integer(input$planning_years[[1]]), as.integer(input$planning_years[[2]]))
  })

  planning_series <- reactive({
    spec <- planning_spec()
    areas <- planning_areas()
    years <- intersect(planning_year_range(), planning_indicator_years(spec$id))
    validate(need(length(years) > 0, paste0(
      "Sem dados de «", spec$label, "» nos anos escolhidos. Disponível: ",
      paste(range(planning_indicator_years(spec$id)), collapse = "-"), "."
    )))
    # The benchmark is computed even when not shown, for the significance marks.
    benchmark <- planning_benchmark()
    table <- planning_indicator_table(union(areas$area, benchmark), years, ids = spec$id, lookup = active_nuts_lookup(), benchmark = benchmark,
                                      education_min_age = planning_education_age(), mode = planning_split_mode(),
                                      basis = active_population_basis())
    # The rebuilt latest year exists for whole regions only; elsewhere it is a
    # gap at the end of the line, not a point.
    table <- table[!(table$year %in% planning_calendar_years() & is.na(table$value)), , drop = FALSE]
    planning_add_significance(table, benchmark)
  })

  output$planningIndicatorHeader <- renderUI({
    spec <- planning_spec()
    tagList(
      h4(planning_indicator_label(spec$id, planning_education_age())),
      helpText(paste0(
        "Unidade: ", spec$unit, ". ",
        if (spec$window > 1) "Cada valor soma três anos e é identificado pelo último. " else "",
        if (!isTRUE(spec$comparable)) "Contagem absoluta, mostrada só para o local." else "Linhas: o local e os comparadores; faixa: intervalo de confiança de 95% do local, quando existe."
      ))
    )
  })

  output$planningIndicatorPlot <- plotly::renderPlotly({
    spec <- planning_spec()
    series <- planning_series()
    validate(need(any(!is.na(series$value)), "Sem valores para este local nos anos escolhidos."))
    series <- series[series$area %in% planning_areas()$area, , drop = FALSE]
    if (isTRUE(spec$comparable)) {
      planning_trend_chart(series, spec, planning_areas())
    } else {
      planning_count_chart(series, spec, planning_location())
    }
  })

  output$planningIndicatorNotes <- renderUI({
    spec <- planning_spec()
    series <- planning_series()
    notes <- character(0)
    if (any(grepl("\\*", series$flag))) notes <- c(notes, "* taxa sobre menos de 1.000 nados-vivos no triénio: exacta, mas instável.")
    if (any(grepl("\u2020", series$flag))) notes <- c(notes, "\u2020 triénio com anos em que os óbitos com menos de 1 ano por município estão incompletos no INE (1995-2001): valor subestimado.")
    if (any(grepl("\u2021", series$flag))) notes <- c(notes, "\u2021 mais de 2% dos óbitos do triénio não tinham idade publicada por município e foram distribuídos pelas idades na proporção dos restantes.")
    if (any(grepl(PLANNING_ESTIMATED_FLAG, series$flag, fixed = TRUE))) notes <- c(notes, if (spec$id %in% c("employees", "earnings_mean")) {
      "\u2248 inclui valores que o INE ocultou num município (segredo estatístico), recuperados do total publicado do Continente ou de Portugal menos os municípios publicados."
    } else {
      "\u2248 mais de 1% dos trabalhadores da área estão em sectores ocultados pelo INE (segredo estatístico) e foram repartidos na proporção do resto da NUTS III."
    })
    if (any(grepl(PLANNING_PROVISIONAL_FLAG, series$flag, fixed = TRUE))) notes <- c(notes, "p provisório: período com um ano que o INE ainda não publicou por causa de morte; óbitos por freguesia e, por idade, reconstruídos dos óbitos semanais (ver Notas).")
    if (spec$id %in% names(PLANNING_INDICATOR_NOTES)) notes <- c(notes, paste("*", PLANNING_INDICATOR_NOTES[[spec$id]]))
    if (spec$id %in% PLANNING_EDUCATION_IDS && planning_education_age() > 0) notes <- c(notes, PLANNING_EDUCATION_NOTE)
    if (any(planning_areas()$area %in% planning_parish_units())) {
      notes <- c(notes, uls_split_note(planning_areas()$area, planning_split_mode()))
      if (identical(planning_split_mode(), "parish")) {
        notes <- c(notes, paste0(
          "Ponderação por freguesias: as partes somam o município; os nascimentos e óbitos vêm dos registos por freguesia ",
          "do INE desde 2014, e o resto da população das freguesias nos Censos de 2021, por grupo etário."
        ))
      }
    }
    if (any(!is.na(series$significance))) {
      notes <- c(notes, paste0(
        "\u25b2 / \u25bc / = : intervalo de confiança de 95% inteiramente acima, abaixo ou a incluir o valor de ",
        planning_benchmark(), " no mesmo período (sem juízo sobre qual o sentido desejável)."
      ))
    }
    breaks <- PLANNING_SERIES_BREAKS[PLANNING_SERIES_BREAKS$indicator == spec$id & PLANNING_SERIES_BREAKS$year %in% planning_year_range(), , drop = FALSE]
    if (nrow(breaks) > 0) notes <- c(notes, paste0("Linha pontilhada em ", breaks$year, ": ", breaks$note))
    if (length(notes) == 0) return(NULL)
    helpText(HTML(paste(htmltools::htmlEscape(notes), collapse = "<br>")))
  })

  # Annual percentage change of each area over the chosen years
  # (R/trend.R). A triennial value shares two years with its neighbour, so
  # only non-overlapping periods ending in the last year are used; otherwise
  # the interval would be far too narrow.
  planning_trend <- reactive({
    spec <- planning_spec()
    series <- planning_series()
    series <- series[is.finite(series$value) & series$value > 0, , drop = FALSE]
    validate(need(nrow(series) > 0, "Sem valores para estimar a tendência."))
    window <- as.integer(spec$window)
    areas <- planning_areas()
    # One period for all: the last year every area has, so the rows compare
    # (the rebuilt latest year exists for whole regions only).
    shown <- series[series$area %in% areas$area, , drop = FALSE]
    last <- min(tapply(shown$year, shown$area, max))
    series <- series[series$year <= last, , drop = FALSE]
    rows <- lapply(seq_len(nrow(areas)), function(i) {
      own <- series[series$area == areas$area[[i]], , drop = FALSE]
      if (nrow(own) == 0) return(NULL)
      if (window > 1) own <- own[(last - own$year) %% window == 0, , drop = FALSE]
      trend <- trend_apc(own$year, own$value, own$lower, own$upper, min_points = 4L)
      if (is.null(trend)) return(NULL)
      tibble::tibble(
        Área = planning_series_name(areas$area[[i]], areas$level[[i]]),
        `Período` = paste0(if (window > 1) paste0(trend$aapc$from - window + 1, "-") else "", trend$aapc$to,
                           " (", trend$n, if (window > 1) " triénios)" else " anos)"),
        `Variação anual (IC 95%)` = trend_label(trend$aapc),
        Leitura = trend_reading(trend$aapc)
      )
    })
    dplyr::bind_rows(rows)
  })

  output$planningTrendTable <- renderTable({
    table <- planning_trend()
    validate(need(nrow(table) > 0, "São precisos pelo menos 4 períodos para estimar a tendência."))
    table
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$planningTrendNote <- renderUI({
    spec <- planning_spec()
    helpText(paste0(
      "Variação percentual anual (APC) de uma regressão log-linear ao longo dos anos escolhidos, cada ponto pesado ",
      "pela sua precisão (método Joinpoint, sem mudanças de tendência). ",
      if (spec$window > 1) "Valores trienais: só triénios que não se sobrepõem, terminando no último ano. " else "",
      "Significativa quando o intervalo não inclui 0; não diz se a subida ou a descida é desejável. ",
      "Com mudanças de tendência, veja o separador Mortalidade Observada. ",
      trend_break_note(planning_year_range()) %||% ""
    ))
  })

  output$planningIndicatorTable <- renderTable({
    planning_series_display(planning_series(), planning_spec(), planning_areas())
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  # Every indicator, up to the last year chosen, for the summary tab.
  planning_all_indicators <- reactive({
    areas <- planning_areas_all_comparators()
    last <- max(planning_year_range())
    # Back to the last census, so the census indicators appear in the summary.
    years <- seq.int(max(1991L, min(last - 4L, 10L * (last %/% 10L) + 1L)), last)
    benchmark <- planning_benchmark()
    table <- planning_indicator_table(union(areas$area, benchmark), years, lookup = active_nuts_lookup(), benchmark = benchmark,
                                      education_min_age = planning_education_age(), mode = planning_split_mode(),
                                      basis = active_population_basis())
    planning_add_significance(table, benchmark)
  })

  # The summary shows every checked comparator whatever the indicator in view.
  planning_areas_all_comparators <- reactive({
    comparators <- planning_available_comparators()
    chosen <- intersect(planning_comparator_state$selected, comparators$area)
    planning_swap_portugal(dplyr::bind_rows(
      tibble::tibble(area = planning_location(), level = "Local"),
      comparators[comparators$area %in% chosen, c("area", "level")]
    ))
  })

  output$planningSummaryTable <- renderTable({
    planning_summary_display(planning_all_indicators(), planning_areas_all_comparators(), planning_education_age())
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$planningSummaryNotes <- renderUI({
    table <- planning_all_indicators()
    shown <- unique(table$indicator[table$area == planning_location() & !is.na(table$value)])
    notes <- unique(PLANNING_INDICATOR_NOTES[intersect(names(PLANNING_INDICATOR_NOTES), shown)])
    helpText(HTML(paste(htmltools::htmlEscape(c(
      paste0("\u25b2 / \u25bc / = : intervalo de confiança de 95% acima, abaixo ou a incluir o valor de ", planning_benchmark(),
             ". Só para os indicadores com intervalo (taxas, proporções, esperança de vida)."),
      paste("*", notes)
    )), collapse = "<br>")))
  })

  planning_ranking_year <- reactive({
    spec <- planning_spec()
    # ULS, municipalities and NUTS III have no deaths by age in the rebuilt
    # latest year: rank them on the last year INE published by cause.
    available <- planning_small_area_years(spec$id)
    available <- available[available <= max(planning_year_range())]
    validate(need(length(available) > 0, "Sem anos disponíveis para este indicador."))
    max(available)
  })

  output$planningRankingNote <- renderUI({
    spec <- planning_spec()
    helpText(if (isTRUE(spec$comparable)) {
      paste0(
        "Todas as ULS do Continente em ", planning_period_label(spec$id, planning_ranking_year()),
        " (o último período disponível até ao ano escolhido). Contorno escuro: a ULS do local; a tracejado, ", planning_benchmark(), ". ",
        "Laranja acima e azul abaixo com significância, cinzento sem diferença. ",
        if (identical(planning_split_mode(), "whole")) {
          "As seis ULS que partilham Lisboa, Loures ou Porto levam o município inteiro: os seus valores sobrepõem-se e não devem ser somados."
        } else {
          "As seis ULS que partilham Lisboa, Loures ou Porto estão repartidas pelas freguesias."
        }
      )
    } else {
      "As contagens absolutas não se comparam entre ULS de tamanhos diferentes. Escolha uma taxa, proporção ou índice."
    })
  })

  output$planningRankingPlot <- plotly::renderPlotly({
    spec <- planning_spec()
    validate(need(isTRUE(spec$comparable), "Sem comparação para contagens absolutas."))
    year <- planning_ranking_year()
    units <- planning_uls_units("units")
    benchmark <- planning_benchmark()
    ranking <- planning_indicator_table(c(benchmark, units), year, ids = spec$id, lookup = active_nuts_lookup(), benchmark = benchmark,
                                        education_min_age = planning_education_age(), mode = planning_split_mode(),
                                      basis = active_population_basis())
    local <- planning_location()
    highlight <- c(local, planning_available_comparators()$area[planning_available_comparators()$level == "ULS"])
    planning_ranking_chart(ranking, spec, year, highlight = highlight, benchmark_area = benchmark)
  })

  # Maps of Portugal (R/planning_maps.R): up to four indicators side by side,
  # every area of one level, the latest period of each up to the last year
  # chosen.
  planning_map_view <- reactive({
    ids <- utils::head(input$planning_map_indicators, PLANNING_MAP_MAX_INDICATORS)
    validate(need(length(ids) > 0, "Escolha pelo menos um indicador."))
    level <- input$planning_map_level %||% "municipio"
    lookup <- active_nuts_lookup()
    boundaries <- planning_map_boundaries(level, active_nuts_vintage())
    validate(need(!is.null(boundaries), "Faltam os limites dos mapas (data/boundaries): corra tools/build_boundaries.R."))
    benchmark <- planning_benchmark()
    values <- shiny::withProgress(message = "A calcular os indicadores de todas as áreas...", value = 0.3, {
      planning_map_values(ids, level, lookup, max(planning_year_range()), benchmark = benchmark,
                          smoothed = isTRUE(input$planning_map_smoothed), mode = planning_split_mode(),
                          basis = active_population_basis())
    })
    validate(need(nrow(values) > 0, "Sem valores para estes indicadores nos anos escolhidos."))
    list(values = values, boundaries = boundaries, level = level, benchmark = benchmark)
  })

  # The outline of the chosen location, apart from the values, so that moving
  # it does not recompute every area of every indicator.
  planning_map_highlight <- reactive({
    level <- input$planning_map_level %||% "municipio"
    lookup <- active_nuts_lookup()
    units <- planning_map_units(level, lookup)
    location <- planning_location()
    if (location %in% units) {
      location
    } else if (identical(level, "municipio") && !location %in% c("Portugal", PLANNING_PORTUGAL_MUNICIPAL, "Continente")) {
      intersect(planning_area_members(location, lookup), units)
    } else {
      character(0)
    }
  })

  planning_map_plot <- function() {
    view <- planning_map_view()
    planning_map_chart(view$values, view$boundaries, input$planning_map_scale %||% "index",
                       highlight = planning_map_highlight(), benchmark = view$benchmark)
  }

  # On screen, one interactive map per indicator (R/cause_maps.R), two per row;
  # the PNG keeps the static drawing.
  output$planningMapWidgets <- renderUI({
    n <- length(unique(planning_map_view()$values$panel))
    width <- if (n == 1) 12 else 6
    do.call(fluidRow, lapply(seq_len(n), function(i) {
      column(width, plotly::plotlyOutput(paste0("planningMap", i), height = if (n == 1) "720px" else "580px"))
    }))
  })
  for (i in seq_len(PLANNING_MAP_MAX_INDICATORS)) local({
    k <- i
    output[[paste0("planningMap", k)]] <- plotly::renderPlotly({
      view <- planning_map_view()
      panels <- unique(view$values$panel)
      req(k <= length(panels))
      v <- view$values[view$values$panel == panels[[k]], , drop = FALSE]
      v$hover <- planning_map_hover(v, view$benchmark)
      direction <- v$direction[[1]]
      colours <- switch(direction,
                        worse = PLANNING_MAP_DIVERGING,
                        better = c(low = PLANNING_MAP_DIVERGING[["high"]], mid = PLANNING_MAP_DIVERGING[["mid"]], high = PLANNING_MAP_DIVERGING[["low"]]),
                        PLANNING_MAP_NEUTRAL)
      meaning <- switch(direction, worse = "vermelho: pior", better = "vermelho: pior", "sem melhor nem pior")
      scale <- input$planning_map_scale %||% "index"
      # An indicator without an interval has no classes: it keeps the colours
      # against the benchmark, and its legend says so.
      without_interval <- identical(scale, "significance") && all(is.na(v$significance))
      if (without_interval) scale <- "index"
      significance <- identical(scale, "significance")
      interactive_choropleth(v, view$boundaries, scale = if (significance) "classes" else scale, colours = colours,
                             classes = map_class_colours(unname(MAP_SIGNIFICANCE_LABELS), direction),
                             title = gsub("\n", "<br>", panels[[k]]),
                             legend_title = if (significance) paste0("Face a ", view$benchmark) else
                               paste0(view$benchmark, " = 100", if (without_interval) " (sem intervalo)" else "", "<br>", meaning),
                             source = "planning_map", highlight = planning_map_highlight(),
                             legend = if (length(panels) > 1) "bottom" else "right")
    })
  })

  output$downloadPlanningMapPNG <- downloadHandler(
    filename = function() paste0("mapas_", input$planning_map_level %||% "municipio", "_", Sys.Date(), ".png"),
    content = function(file) {
      panels <- length(unique(planning_map_view()$values$panel))
      grDevices::png(file, width = 10, height = if (panels <= 2) 6.5 else 12, units = "in", res = 150, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE)
      planning_map_draw(planning_map_plot())
    }
  )

  output$downloadPlanningMapCSV <- downloadHandler(
    filename = function() paste0("mapas_", input$planning_map_level %||% "municipio", "_", Sys.Date(), ".csv"),
    content = function(file) {
      v <- planning_map_view()$values
      write_csv_utf8(tibble::tibble(
        indicador = gsub("\n", " ", v$panel), area = v$area, valor = v$value, valor_observado = v$observed,
        marca = v$flag, referencia = v$benchmark_value, indice_referencia_100 = round(v$index, 1),
        face_a_referencia = v$significance, quintil = v$quintile
      ), file)
    }
  )

  output$planningMapNote <- renderUI({
    view <- planning_map_view()
    scale <- input$planning_map_scale %||% "index"
    helpText(paste0(
      if (identical(scale, "index")) {
        paste0("Cada indicador face a ", view$benchmark, " (= 100), de 67 a 150 (os valores fora ficam na cor do extremo), ",
               "colorido pelo que a diferença significa: vermelho onde a área está pior do que ", view$benchmark, ", azul onde está melhor ",
               "- uma mortalidade alta e uma esperança de vida baixa ficam ambas a vermelho. Os indicadores sem melhor nem pior ",
               "(estrutura etária, natalidade, sectores de emprego, taxa bruta de mortalidade, que depende sobretudo da idade da população) ",
               "ficam a violeta abaixo e amarelo acima. Assim, indicadores em unidades diferentes lêem-se juntos. ")
      } else if (identical(scale, "significance")) {
        paste0("Três cores por indicador: \u25b2 acima e \u25bc abaixo de ", view$benchmark, " quando o intervalo de 95% o exclui, com as cores do que ",
               "a diferença significa (vermelho pior, azul melhor; violeta e amarelo sem melhor nem pior). Só os indicadores com intervalo têm classes. ")
      } else {
        "Cada indicador dividido em cinco grupos com o mesmo número de áreas: mostra a ordenação dentro do indicador, não a distância a Portugal nem se é melhor ou pior. "
      },
      "Cada mapa usa o último período disponível até ao último ano escolhido, indicado no título. ",
      if (isTRUE(input$planning_map_smoothed)) "A suavização só se aplica às taxas com modelo de contagem (as do funil); os títulos dizem quais. " else "",
      if (identical(view$level, "uls")) "As ULS que partilham Lisboa, Loures, Porto ou Gondomar têm os valores das suas freguesias (a Lomba na ULS Vila Nova de Gaia/Espinho). As ilhas não têm ULS. " else "",
      "Cinzento: sem valor. Contorno escuro: o local escolhido. Passe o rato sobre uma área para ver o valor, use a roda do rato para aproximar ",
      "e clique numa área para ver os seus valores em todos os indicadores."
    ))
  })

  output$planningMapClicked <- renderUI({
    # Read only once there is a click: before the maps are drawn, plotly warns
    # that the event is not registered.
    req(input[["plotly_click-planning_map"]])
    click <- plotly::event_data("plotly_click", source = "planning_map")
    req(click)
    view <- planning_map_view()
    area <- as.character(click$customdata %||% click$location)[1]
    req(area, area %in% view$values$area)
    rows <- view$values[view$values$area == area, , drop = FALSE]
    table <- tibble::tibble(
      Indicador = gsub("\n", " · ", rows$panel),
      Valor = paste0(planning_format_value(rows$value, rows$digits), ifelse(is.na(rows$flag), "", rows$flag)),
      !!view$benchmark := planning_format_value(rows$benchmark_value, rows$digits),
      `Índice` = planning_format_value(rows$index, 0),
      Quintil = ifelse(is.na(rows$quintile), "—", as.character(rows$quintile))
    )
    tagList(h5(area), renderTable(table, striped = TRUE, bordered = TRUE, spacing = "s", align = "l"))
  })

  # -------------------------
  # Maps of any cause (R/cause_maps.R)
  # -------------------------
  # The settings of the map on screen. They are captured when the map is asked
  # for - by its button or from another tab - rather than read from the inputs,
  # which a tab switch updates only after this round.
  map_settings <- reactiveVal(NULL)
  # The areas chosen in the tab the map was opened from, outlined on the map.
  map_focus <- reactiveVal(character(0))

  map_end_year_choices <- function(cause) {
    years <- rev(cause_map_end_years(cause))
    stats::setNames(as.character(years), paste0(years - 2L, "-", years))
  }
  # Earlier triennia with no year in common with the chosen one; by default the
  # one ten years before, or the oldest.
  map_base_year_choices <- function(cause, end_year) {
    choices <- map_end_year_choices(cause)
    choices[as.integer(choices) <= as.integer(end_year) - 3L]
  }
  observeEvent(input$map_cause, {
    choices <- map_end_year_choices(input$map_cause)
    keep <- if (isTRUE(input$map_end_year %in% choices)) input$map_end_year else utils::head(choices, 1)
    updateSelectInput(session, "map_end_year", choices = choices, selected = keep)
  })
  observeEvent(list(input$map_cause, input$map_end_year), {
    req(input$map_end_year)
    choices <- map_base_year_choices(input$map_cause, input$map_end_year)
    target <- as.character(as.integer(input$map_end_year) - 9L)
    keep <- if (isTRUE(input$map_base_year %in% choices)) input$map_base_year else if (target %in% choices) target else utils::tail(choices, 1)
    updateSelectInput(session, "map_base_year", choices = choices, selected = keep)
  })

  request_map <- function(cause = input$map_cause, sex = input$map_sex, end_year = input$map_end_year) {
    clear_tab_notes("maps")
    view <- input$map_view %||% "period"
    map_settings(list(
      cause = as.character(cause)[1], sex = sex %||% "HM", metric = input$map_metric %||% "smr",
      end_year = suppressWarnings(as.integer(end_year)), level = input$map_level %||% "municipio",
      scale = input$map_scale %||% "index", view = view,
      trend_n = suppressWarnings(as.integer(input$map_trend_n %||% "5")),
      base_year = suppressWarnings(as.integer(input$map_base_year)),
      target = suppressWarnings(as.integer(input$map_target_year %||% "2030")),
      smoothed = isTRUE(input$map_smoothed) && identical(input$map_metric, "smr") && identical(view, "period")
    ))
  }
  observeEvent(input$go_map, request_map())

  # "Ver esta causa no mapa" in the mortality tabs: the tab's cause, sex and
  # areas (outlined), the latest triennium, the map's other settings as they are.
  open_map_for <- function(cause, sex, areas = character(0)) {
    cause <- as.character(cause %||% "")[1]
    if (!nzchar(cause)) cause <- "Todas as causas de morte"
    choices <- map_end_year_choices(cause)
    updateSelectInput(session, "map_cause", selected = cause)
    updateRadioButtons(session, "map_sex", selected = sex %||% "HM")
    updateSelectInput(session, "map_end_year", choices = choices, selected = utils::head(choices, 1))
    updateNavbarPage(session, "main_nav", selected = "Mapas")
    map_focus(as.character(areas))
    request_map(cause, sex, utils::head(choices, 1))
  }
  observeEvent(input$map_from_observed, open_map_for(input$cause, input$sex, input$area))
  observeEvent(input$map_from_beginner, open_map_for(input$beginner_cause, input$beginner_sex, input$beginner_area))
  observeEvent(input$map_from_annual, open_map_for(utils::head(input$annual_cause, 1), input$annual_sex, input$annual_area))

  # The outline: the focus areas at this level, or, on the municipal map, the
  # municipalities of a focus region.
  map_highlight <- function(level, lookup) {
    focus <- setdiff(map_focus(), c("Portugal", PLANNING_PORTUGAL_MUNICIPAL, "Continente"))
    if (length(focus) == 0) return(character(0))
    units <- planning_map_units(level, lookup)
    direct <- intersect(focus, units)
    if (identical(level, "municipio")) {
      direct <- union(direct, intersect(unlist(lapply(setdiff(focus, units), planning_area_members, lookup = lookup)), units))
    }
    direct
  }

  map_period_label <- function(end_year) paste0(end_year - 2L, "-", end_year)

  cause_map <- reactive({
    settings <- map_settings()
    validate(need(!is.null(settings), "Escolha a causa, o sexo e a métrica e carregue em «Desenhar mapa»."))
    lookup <- active_nuts_lookup()
    boundaries <- planning_map_boundaries(settings$level, active_nuts_vintage())
    validate(need(!is.null(boundaries), "Faltam os limites dos mapas (data/boundaries): corra tools/build_boundaries.R."))
    years <- cause_map_end_years(settings$cause)
    validate(need(length(years) > 0, "Esta causa não tem três anos seguidos com óbitos por causa."))
    end_year <- if (isTRUE(settings$end_year %in% years)) settings$end_year else max(years)
    basis <- active_population_basis()
    out <- list(boundaries = boundaries, settings = settings, end_year = end_year,
                highlight = map_highlight(settings$level, lookup))

    if (identical(settings$view, "projection")) {
      validate(need(settings$level %in% PROJECTION_LEVELS,
                    "As projecções fazem-se por NUTS III, ULS ou NUTS II: por município as séries são curtas e instáveis demais."))
      units <- planning_map_units(settings$level, lookup)
      values <- shiny::withProgress(message = paste0("A projectar ", length(units), " áreas até ", settings$target, "..."), value = 0, {
        projection_map_values(settings$cause, settings$sex, identical(settings$metric, "dsr75"), settings$target, settings$level, lookup)
      })
      validate(need(nrow(values) > 0 && any(is.finite(values$value)), "Nenhuma área pôde ser projectada para esta causa."))
      values$hover <- cause_map_projection_hover(values, settings$target)
      # The change to the target year in %: a few points either way, which the
      # 67-150 index scale would leave nearly white.
      values$z <- values$index - 100
      out$target <- settings$target
      missing <- sum(!is.finite(values$value))
      if (missing > 0) notify_data(paste0(
        missing, " de ", nrow(values), " áreas sem projecção (a cinzento): menos de ", PROJECTION_MIN_YEARS,
        " anos seguidos com taxa no fim da série, ou nenhum dos modelos se ajustou."
      ), type = "message", id = "map-missing")
      if (identical(settings$metric, "smr")) notify_data(
        "A projecção usa a taxa padronizada (a da Previsão Guiada); a SMR não se projecta.", type = "message", id = "map-projection-metric")
    } else if (identical(settings$view, "trend")) {
      n <- settings$trend_n %||% 5L
      validate(need(length(cause_map_triennia(settings$cause, end_year, n)) >= 4,
                    "A tendência precisa de pelo menos 4 triénios com óbitos por causa: escolha um triénio mais recente ou outra causa."))
      values <- shiny::withProgress(message = "A calcular a causa em cada triénio...", value = 0, {
        cause_map_trend(settings$cause, settings$sex, end_year, n, settings$metric, settings$level, lookup, basis = basis,
                        progress = function(i, total) shiny::setProgress(i / total))
      })
      validate(need(nrow(values) > 0 && any(is.finite(values$z)), "Sem tendência para esta causa nestes triénios."))
      values$hover <- cause_map_trend_hover(values, settings$metric)
      out$periods <- attr(values, "periods")
      missing <- sum(!is.finite(values$z))
      if (missing > 0) notify_data(paste0(missing, " de ", nrow(values), " áreas sem tendência (a cinzento): menos de 4 triénios com valor."),
                                   type = "message", id = "map-missing")
    } else if (identical(settings$view, "change")) {
      base <- settings$base_year
      validate(need(isTRUE(base %in% years) && base <= end_year - 3L, "Escolha um triénio de comparação anterior, sem anos em comum."))
      values <- shiny::withProgress(message = "A calcular a causa nos dois triénios...", value = 0.3, {
        cause_map_change(settings$cause, settings$sex, end_year, base, settings$metric, settings$level, lookup, basis = basis)
      })
      validate(need(nrow(values) > 0 && any(is.finite(values$index)), "Sem valores para comparar estes triénios."))
      values$hover <- cause_map_change_hover(values, settings$metric, map_period_label(base), map_period_label(end_year))
      out$base_year <- base
      missing <- sum(!is.finite(values$index))
      if (missing > 0) notify_data(paste0(missing, " de ", nrow(values), " áreas sem valor num dos triénios (a cinzento)."),
                                   type = "message", id = "map-missing")
    } else {
      values <- shiny::withProgress(message = "A calcular a causa em todas as áreas...", value = 0.3, {
        cause_map_values(settings$cause, settings$sex, end_year, settings$metric, settings$level, lookup,
                         basis = basis, smoothed = settings$smoothed)
      })
      validate(need(nrow(values) > 0 && any(is.finite(values$value)), "Sem valores para esta causa neste triénio."))
      missing <- sum(!is.finite(values$value))
      if (missing > 0) {
        notify_data(paste0(
          missing, " de ", nrow(values), " áreas sem valor (a cinzento): sem óbitos por idade publicados nesses anos, ",
          "ou mais de 25% dos óbitos sem idade (o INE oculta as idades em células pequenas)."
        ), type = "message", id = "map-missing")
      }
      values$hover <- cause_map_hover(values, settings$metric)
    }
    out$values <- values
    out
  })

  # Projections by area: the guided forecast's default projection, the mean of
  # its six models, fitted to each area's annual standardised rate. The annual
  # series of every area come at once from the maps' engine with one-year
  # windows (planning_cause_standardised(span = 1)): the guided tab's own
  # series, built area by area, took 15-20 s per ULS, and the two agree to
  # 0.00% in the ULS and NUTS II tested (to 3% in a year of a NUTS III built
  # from INE's rows and municipalities). As in the guided tab, a year with
  # fewer than 10 deaths has no rate. The validation is left out: it does not
  # change the mean of six, only the error shown beside it. A series with
  # years without a rate keeps its most recent run of consecutive years with
  # one, at least PROJECTION_MIN_YEARS long, as decided for the 82_MRT
  # projections; shorter series are not projected.
  PROJECTION_MIN_YEARS <- 8L
  PROJECTION_MIN_DEATHS <- 10
  PROJECTION_LEVELS <- c("nuts3", "uls", "nuts2")
  projection_cache <- new.env(parent = emptyenv())

  projection_state_key <- function(...) paste(..., active_nuts_vintage(), active_data_mode(), active_population_basis(),
                                              active_region_source(), sep = "|")

  # Annual standardised rates of Portugal and every unit of `level`.
  projection_series <- function(cause, sex, under75, level, lookup) {
    key <- projection_state_key("series", cause, sex, under75, level)
    if (exists(key, envir = projection_cache, inherits = FALSE)) return(get(key, envir = projection_cache))
    units <- planning_map_units(level, lookup)
    mode <- if (identical(level, "uls")) "parish" else PLANNING_DEFAULT_SPLIT_MODE
    years <- seq.int(min(rate_years), max(rate_years))
    rows <- lapply(seq_along(years), function(i) {
      shiny::setProgress(0.7 * i / length(years), detail = paste("série anual,", years[[i]]))
      tab <- planning_cause_standardised(c(units, "Portugal"), years[[i]], lookup = lookup, sex = sex, mode = mode,
                                         basis = active_population_basis(), causes = cause, span = 1L)
      if (nrow(tab) == 0) return(NULL)
      value <- if (under75) tab$dsr75 else tab$dsr
      deaths <- if (under75) tab$observed75 else tab$observed
      tibble::tibble(area = tab$area, year = years[[i]],
                     value = ifelse(is.finite(deaths) & deaths >= PROJECTION_MIN_DEATHS, value, NA_real_),
                     lower = if (under75) tab$dsr75_lower else tab$dsr_lower,
                     upper = if (under75) tab$dsr75_upper else tab$dsr_upper)
    })
    out <- dplyr::bind_rows(rows)
    assign(key, out, envir = projection_cache)
    out
  }

  projection_one <- function(series, target, token) {
    # The rebuilt latest year serves only whole regions, all causes and both
    # sexes: elsewhere it is empty, and the series ends with INE's last year of
    # causes of death, as in the guided tab.
    series <- series[is.finite(series$value) | series$year <= max(planning_cause_years()), , drop = FALSE]
    series <- series[order(series$year), , drop = FALSE]
    ok <- is.finite(series$value) & series$value > 0
    if (length(ok) == 0 || !ok[[length(ok)]]) return(NULL)
    last <- max(series$year)
    start <- last
    while ((start - 1L) %in% series$year[ok]) start <- start - 1L
    if (last - start + 1L < PROJECTION_MIN_YEARS || target <= last) return(NULL)
    history <- list(series = series[series$year >= start, c("year", "value", "lower", "upper"), drop = FALSE])
    result <- run_forecast_models(history, FORECAST_MEAN_MEMBERS, horizon = target - last, kind = "maps", token = token)
    members <- intersect(FORECAST_MEAN_MEMBERS, get_successful_model_ids(result))
    if (length(members) < FORECAST_MEAN_MIN_MEMBERS) return(NULL)
    setup <- get_transformation_setup(result$obs$value, result$transform_method)
    fc <- combine_model_forecasts(result$fc, members, setup$forward, setup$inverse, result$conf_level, id = FORECAST_MEAN_ID)
    row <- fc[fc$year == target, , drop = FALSE]
    if (nrow(row) != 1) return(NULL)
    tibble::tibble(last_year = last, last_value = series$value[series$year == last], fit_from = start,
                   value = row$mean, lower = row$lower, upper = row$upper, members = length(members))
  }

  projection_map_values <- function(cause, sex, under75, target, level, lookup) {
    key <- projection_state_key("projection", cause, sex, under75, target, level)
    if (exists(key, envir = projection_cache, inherits = FALSE)) return(get(key, envir = projection_cache))
    series <- projection_series(cause, sex, under75, level, lookup)
    if (nrow(series) == 0) return(tibble::tibble())
    token <- isolate(cancel_seq$maps)
    areas <- c("Portugal", planning_map_units(level, lookup))
    rows <- lapply(seq_along(areas), function(i) {
      shiny::setProgress(0.7 + 0.3 * i / length(areas), detail = areas[[i]])
      one <- tryCatch(projection_one(series[series$area == areas[[i]], , drop = FALSE], target, token), error = function(e) NULL)
      if (is.null(one)) NULL else dplyr::mutate(one, area = areas[[i]], .before = 1)
    })
    projected <- dplyr::bind_rows(rows)
    if (nrow(projected) == 0) return(tibble::tibble())
    reference <- projected$value[projected$area == "Portugal"][1]
    out <- dplyr::left_join(tibble::tibble(area = areas[-1]), projected, by = "area")
    out$benchmark_value <- if (length(reference) == 1) reference else NA_real_
    out$index <- ifelse(is.finite(out$value) & is.finite(out$last_value) & out$last_value > 0, out$value / out$last_value * 100, NA_real_)
    out$class <- ifelse(!is.finite(out$index), NA_character_,
                        ifelse(out$lower > out$last_value, MAP_CHANGE_LABELS[[3]],
                               ifelse(out$upper < out$last_value, MAP_CHANGE_LABELS[[1]], MAP_CHANGE_LABELS[[2]])))
    out$quintile <- map_quintile(out$index)
    assign(key, out, envir = projection_cache)
    out
  }

  # Title, colours and legend of the map on screen.
  cause_map_drawing <- function(m) {
    s <- m$settings
    sex_label <- c(HM = "ambos os sexos", H = "homens", M = "mulheres")[[s$sex]]
    metric <- cause_map_metric_label(s$metric)
    if (identical(s$view, "projection")) {
      metric <- if (identical(s$metric, "dsr75")) cause_map_metric_label("dsr75") else cause_map_metric_label("dsr")
      list(title = paste0(s$cause, " · ", sex_label, " · projecção para ", m$target,
                          "<br><sup>", metric, ", média de 6 modelos, face ao último ano observado</sup>"),
           scale = switch(s$scale, significance = "classes", quintile = "quintile", "linear"),
           classes = map_class_colours(MAP_CHANGE_LABELS),
           legend = if (identical(s$scale, "significance")) "Projecção" else "% face ao último ano<br>vermelho: sobe",
           limits = c(-MAP_PROJECTION_LIMIT, MAP_PROJECTION_LIMIT), ticks = c(-20, -10, 0, 10, 20),
           tick_labels = c("-20", "-10", "0", "+10", "+20"))
    } else if (identical(s$view, "trend")) {
      first <- min(m$periods) - 2L
      list(title = paste0(s$cause, " · ", sex_label, " · ", first, "-", m$end_year,
                          "<br><sup>Variação anual: ", metric, ", ", length(m$periods), " triénios</sup>"),
           scale = switch(s$scale, significance = "classes", quintile = "quintile", "linear"),
           classes = map_class_colours(MAP_TREND_LABELS), legend = if (identical(s$scale, "significance")) "Tendência" else "% por ano<br>vermelho: a subir",
           limits = c(-MAP_TREND_LIMIT, MAP_TREND_LIMIT), ticks = c(-5, -2.5, 0, 2.5, 5), tick_labels = c("-5", "-2,5", "0", "+2,5", "+5"))
    } else if (identical(s$view, "change")) {
      list(title = paste0(s$cause, " · ", sex_label, "<br><sup>", metric, ": ", map_period_label(m$base_year), " → ",
                          map_period_label(m$end_year), "</sup>"),
           scale = switch(s$scale, significance = "classes", quintile = "quintile", "index"),
           classes = map_class_colours(MAP_CHANGE_LABELS),
           legend = if (identical(s$scale, "significance")) "Variação" else paste0(map_period_label(m$base_year), " = 100<br>vermelho: subiu"))
    } else {
      list(title = cause_map_title(s$cause, s$sex, m$end_year, s$metric, s$smoothed),
           scale = switch(s$scale, significance = "classes", quintile = "quintile", "index"),
           classes = map_class_colours(unname(MAP_SIGNIFICANCE_LABELS)),
           legend = if (identical(s$scale, "significance")) "Face a Portugal" else "Portugal = 100<br>vermelho: acima")
    }
  }

  output$causeMap <- plotly::renderPlotly({
    m <- cause_map()
    d <- cause_map_drawing(m)
    interactive_choropleth(m$values, m$boundaries, scale = d$scale, title = d$title, legend_title = d$legend,
                           source = "cause_map", highlight = m$highlight, classes = d$classes,
                           limits = d$limits, ticks = d$ticks, tick_labels = d$tick_labels)
  })

  output$causeMapNote <- renderUI({
    m <- cause_map()
    s <- m$settings
    smr <- identical(s$metric, "smr")
    helpText(paste0(
      switch(
        s$view,
        projection = paste0(
          "A projecção de cada área para ", m$target, " com o método da Previsão Guiada: a média de 6 modelos (ARIMA, ETS, Theta, TBATS, Holt e ",
          "Holt amortecido) ajustados à taxa padronizada anual, face ao último ano observado. ",
          if (identical(s$scale, "significance")) "\u25b2 e \u25bc quando todo o intervalo de previsão de 95% fica acima ou abaixo do último valor. " else
            paste0("Escala de -", MAP_PROJECTION_LIMIT, "% a +", MAP_PROJECTION_LIMIT, "%: vermelho sobe, azul desce. "),
          "Uma projecção não é uma meta: diz para onde a série vai se o passado se repetir. O rato mostra o intervalo de previsão. "
        ),
        trend = paste0(
          "A variação percentual anual de cada área, numa recta log-linear sobre os ", length(m$periods), " triénios, cada um pesado pela ",
          "sua precisão (como a tendência dos Indicadores de Planeamento). ",
          if (identical(s$scale, "significance")) "▲ a subir e ▼ a descer quando o intervalo de 95% exclui 0. " else
            paste0("Escala de -", MAP_TREND_LIMIT, "% a +", MAP_TREND_LIMIT, "% por ano: vermelho a subir, azul a descer. "),
          if (smr) "Com a SMR, a tendência é face a Portugal: 0% quer dizer que a área acompanha Portugal. " else "",
          "Os triénios de 2020-2022 incluem a COVID-19 e, a partir de 2021, a população revista do INE. "
        ),
        change = paste0(
          "O valor no triénio mais recente face ao do triénio de comparação (= 100), com um intervalo dos dois intervalos. ",
          if (identical(s$scale, "significance")) "▲ subiu e ▼ desceu quando o intervalo de 95% exclui 100. " else
            "Escala de 67 a 150: vermelho subiu, azul desceu. ",
          if (smr) "Com a SMR, a variação é da posição face a Portugal, não da mortalidade em si. " else ""
        ),
        paste0(
          if (identical(s$scale, "index")) {
            paste0("Cada área face a Portugal (= 100), de 67 a 150 (os valores fora ficam na cor do extremo): vermelho acima de Portugal, ",
                   "azul abaixo. ", if (smr) "A SMR já é esse índice. " else "O índice é a taxa da área a dividir pela de Portugal. ")
          } else if (identical(s$scale, "significance")) {
            "Três cores: ▲ acima e ▼ abaixo de Portugal quando o intervalo de 95% o exclui; as outras áreas são semelhantes. "
          } else {
            "Cinco grupos com o mesmo número de áreas: mostra a ordenação, não a distância a Portugal. "
          },
          "Óbitos e pessoas-ano somados nos três anos do triénio. ",
          if (isTRUE(s$smoothed)) "A SMR suavizada aproxima as áreas pequenas da média; o rato mostra também a observada. " else
            "Em municípios pequenos uma diferença grande pode ser só acaso: o rato diz se o intervalo de 95% exclui Portugal (▲ ▼). "
        )
      ),
      if (identical(s$scale, "quintile")) "Quintis: cinco grupos com o mesmo número de áreas. " else "",
      if (identical(s$level, "uls")) "As ULS que partilham Lisboa, Loures, Porto ou Gondomar têm os valores das suas freguesias. As ilhas não têm ULS. " else "",
      if (length(m$highlight) > 0) "Contorno escuro: o local escolhido no separador de onde veio. " else "",
      "Passe o rato para ver os valores, use a roda do rato para aproximar e clique numa área para abrir a sua série."
    ))
  })

  # As for the planning maps, read only once there is a click.
  cause_map_area <- reactive({
    req(input[["plotly_click-cause_map"]])
    click <- plotly::event_data("plotly_click", source = "cause_map")
    req(click)
    as.character(click$customdata %||% click$location)[1]
  })

  output$causeMapClicked <- renderUI({
    m <- cause_map()
    area <- cause_map_area()
    row <- m$values[m$values$area == area, , drop = FALSE]
    req(nrow(row) == 1)
    digits <- if (identical(m$settings$metric, "smr")) 0 else 1
    pct <- function(x, d = 1) formatC(x, format = "f", digits = d, decimal.mark = ",", flag = "+")
    table <- switch(
      m$settings$view,
      projection = tibble::tibble(
        `Último ano` = if (is.finite(row$last_year)) paste0(planning_format_value(row$last_value, digits), " (", row$last_year, ")") else "—",
        !!paste0("Projecção ", m$target) := planning_format_value(row$value, 1),
        `Intervalo de previsão 95%` = paste0(planning_format_value(row$lower, 1), " a ", planning_format_value(row$upper, 1)),
        `Variação` = if (is.finite(row$index)) paste0(pct(row$index - 100, 0), "%") else "—",
        `Portugal ` = planning_format_value(row$benchmark_value, 1),
        `Ajuste` = if (is.finite(row$fit_from)) paste0(row$fit_from, "-", row$last_year) else "—"
      ),
      trend = tibble::tibble(
        `Métrica` = cause_map_metric_label(m$settings$metric),
        `Variação anual` = if (is.finite(row$z)) paste0(pct(row$z), "%") else "—",
        `IC 95%` = if (is.finite(row$z)) paste0(pct(row$apc_lower), "% a ", pct(row$apc_upper), "%") else "—",
        `Leitura` = if (is.na(row$class)) "—" else row$class,
        `Primeiro triénio` = planning_format_value(row$first, digits),
        `Último triénio` = planning_format_value(row$last, digits),
        `Triénios` = as.character(row$n)
      ),
      change = tibble::tibble(
        `Métrica` = cause_map_metric_label(m$settings$metric),
        !!map_period_label(m$base_year) := planning_format_value(row$value_base, digits),
        !!map_period_label(m$end_year) := planning_format_value(row$value, digits),
        `Variação` = if (is.finite(row$index)) paste0(pct(row$index - 100, 0), "%") else "—",
        `IC 95%` = if (is.finite(row$change_lower)) paste0(pct(row$change_lower - 100, 0), "% a ", pct(row$change_upper - 100, 0), "%") else "—",
        `Leitura` = if (is.na(row$class)) "—" else row$class
      ),
      tibble::tibble(
        `Métrica` = cause_map_metric_label(m$settings$metric),
        Valor = paste0(planning_format_value(row$value, digits), row$flag),
        `IC 95%` = paste0(planning_format_value(row$lower, digits), " a ", planning_format_value(row$upper, digits)),
        `Óbitos` = planning_format_value(row$observed, 0),
        Esperados = planning_format_value(row$expected, 1),
        Portugal = planning_format_value(row$benchmark_value, digits),
        `Face a Portugal` = if (is.na(row$significance)) "—" else row$class
      )
    )
    openable <- area %in% area_choices_for(active_nuts_vintage(), active_nuts_lookup())
    tagList(
      h5(area),
      renderTable(table, striped = TRUE, bordered = TRUE, spacing = "s", align = "l"),
      if (openable) {
        actionButton("map_open_series", "Abrir a série na Mortalidade Observada", icon = icon("chart-line"))
      } else {
        helpText("Esta área não se escolhe na Mortalidade Observada (as ULS partilhadas só existem nos Indicadores de Planeamento).")
      }
    )
  })

  # The clicked area's series in the observed tab: same cause and sex, the
  # standardised rate (under 75 for that metric), loaded at once.
  observeEvent(input$map_open_series, {
    m <- cause_map()
    area <- cause_map_area()
    updateSelectInput(session, "area", selected = area)
    updateSelectInput(session, "cause", selected = m$settings$cause)
    updateSelectInput(session, "sex", selected = m$settings$sex)
    updateRadioButtons(session, "rate_type", selected = "dsr")
    updateRadioButtons(session, "population", selected = if (identical(m$settings$metric, "dsr75")) "Menos de 75 anos" else "Total")
    updateNavbarPage(session, "main_nav", selected = "Mortalidade Observada")
    session$sendCustomMessage("click-button", "go_rates")
  })

  observe({
    session$sendCustomMessage("lock-downloads", list(ids = c("downloadCauseMapCSV", "downloadCauseMapPNG"), locked = is.null(map_settings())))
  })

  cause_map_file <- function(extension) {
    s <- map_settings()
    paste0("mapa_", planning_cause_file_token(s$cause), "_", s$sex, "_", s$metric, "_", s$view, "_", s$level, "_", Sys.Date(), extension)
  }
  output$downloadCauseMapCSV <- downloadHandler(
    filename = function() cause_map_file(".csv"),
    content = function(file) {
      m <- cause_map()
      v <- m$values
      s <- m$settings
      head <- tibble::tibble(causa = s$cause, sexo = s$sex, metrica = cause_map_metric_label(s$metric), area = v$area)
      table <- switch(
        s$view,
        projection = dplyr::bind_cols(head, tibble::tibble(
          ano_projectado = m$target, ultimo_ano = v$last_year, valor_ultimo_ano = v$last_value, projeccao = v$value,
          ip95_inferior = v$lower, ip95_superior = v$upper, indice_ultimo_ano_100 = round(v$index, 1), leitura = v$class,
          ajuste_desde = v$fit_from, modelos = v$members, portugal_projeccao = v$benchmark_value, quintil = v$quintile)),
        trend = dplyr::bind_cols(head, tibble::tibble(
          periodo = paste0(min(m$periods) - 2L, "-", m$end_year), triennios = v$n, variacao_anual_pct = v$z,
          ic95_inferior = v$apc_lower, ic95_superior = v$apc_upper, leitura = v$class,
          primeiro_triennio = v$first, ultimo_triennio = v$last, quintil = v$quintile)),
        change = dplyr::bind_cols(head, tibble::tibble(
          triennio_base = map_period_label(m$base_year), triennio = map_period_label(m$end_year),
          valor_base = v$value_base, valor = v$value, indice_base_100 = round(v$index, 1),
          ic95_inferior = v$change_lower, ic95_superior = v$change_upper, leitura = v$class, quintil = v$quintile)),
        dplyr::bind_cols(head, tibble::tibble(
          triennio = map_period_label(m$end_year), valor = v$value, valor_observado = v$observed_value,
          ic95_inferior = v$lower, ic95_superior = v$upper, obitos = v$observed, esperados = v$expected, marca = v$flag,
          portugal = v$benchmark_value, indice_portugal_100 = round(v$index, 1), face_a_portugal = v$significance, quintil = v$quintile))
      )
      write_csv_utf8(table, file)
    }
  )
  output$downloadCauseMapPNG <- downloadHandler(
    filename = function() cause_map_file(".png"),
    content = function(file) {
      m <- cause_map()
      d <- cause_map_drawing(m)
      title <- trimws(gsub("<br><sup>|</sup>", "\n", d$title))
      values <- dplyr::mutate(m$values, panel = title, direction = "worse")
      grDevices::png(file, width = 8, height = 9, units = "in", res = 150, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE)
      planning_map_draw(planning_map_chart(values, m$boundaries, d$scale, highlight = m$highlight,
                                           classes = d$classes, limits = d$limits, ticks = d$ticks, tick_labels = d$tick_labels,
                                           legend_title = gsub("<br>", "\n", d$legend)))
    }
  )

  planning_funnel <- reactive({
    spec <- planning_spec()
    model <- PLANNING_FUNNEL_MODELS[PLANNING_FUNNEL_MODELS$id == spec$id, , drop = FALSE]
    validate(need(nrow(model) == 1, paste0(
      "O funil aplica-se às taxas de acontecimentos e às proporções de nascimentos, que têm um modelo de contagem ",
      "(Poisson ou binomial). Escolha, por exemplo, a mortalidade infantil, a natalidade ou os nascimentos pré-termo."
    )))
    year <- planning_ranking_year()
    kind <- input$planning_funnel_units %||% "ULS"
    lookup <- active_nuts_lookup()
    units <- if (identical(kind, "ULS")) {
      planning_uls_units("units")
    } else {
      levels <- planning_area_levels(lookup)
      levels$area[levels$level == kind]
    }
    benchmark <- planning_benchmark()
    table <- planning_indicator_table(c(benchmark, units), year, ids = spec$id, lookup = lookup, benchmark = benchmark, mode = planning_split_mode(),
                                      basis = active_population_basis())
    data <- planning_funnel_data(table, spec, benchmark)
    validate(need(!is.null(data) && nrow(data) > 0, "Sem valores para o funil neste período."))
    list(data = planning_eb_smooth(data, model), year = year, kind = kind, model = model)
  })

  planning_funnel_smoothed <- reactive(identical(input$planning_funnel_values %||% "observed", "smoothed"))

  output$planningFunnelNote <- renderUI({
    funnel <- planning_funnel()
    if (planning_funnel_smoothed()) {
      w <- funnel$data$shrinkage_weight
      return(helpText(paste0(
        "Suavização por Bayes empírico (Marshall, 1991): cada valor é aproximado da média de todas as unidades (",
        planning_format_value(attr(funnel$data, "smoothing_mean"), planning_spec()$digits), ") tanto mais quanto menor ",
        "o seu denominador, porque aí o acaso pesa mais. Círculo vazio: observado; cheio: suavizado; o traço é a deslocação. ",
        "Nesta selecção, o valor próprio pesa de ", planning_format_value(100 * min(w, na.rm = TRUE), 0), "% a ",
        planning_format_value(100 * max(w, na.rm = TRUE), 0), "%. Serve para ordenar ou mapear áreas pequenas com valores ",
        "mais estáveis; não diz se uma unidade difere de ", planning_benchmark(), " - para isso, use o funil (valores observados). ",
        "Os valores suavizados são estimativas e não substituem os observados em relatórios."
      )))
    }
    counts <- table(funnel$data$position)
    outside <- sum(counts[names(counts) != "Dentro dos limites"])
    helpText(paste0(
      "Cada ponto é uma unidade (", nrow(funnel$data), ") em ", planning_period_label(planning_spec()$id, funnel$year),
      ": o valor contra a dimensão do denominador. As linhas marcam o que uma unidade desse tamanho mostraria só por acaso ",
      "à volta de ", planning_benchmark(), " (95%: tracejado; 99,8%: pontilhado). Fora dos limites: ", outside,
      ". Com 20 unidades, uma fora do limite de 95% é esperada por acaso; fora do de 99,8%, quase nunca. ",
      "Contorno escuro: o local e a sua ULS.",
      if (identical(planning_split_mode(), "whole") && identical(input$planning_funnel_units %||% "ULS", "ULS")) " As seis ULS que partilham Lisboa, Loures ou Porto levam o município inteiro, pelo que se sobrepõem." else "",
      if (planning_spec()$id %in% c("birth_rate", "death_rate")) " As taxas brutas dependem da estrutura etária: uma população mais velha tem mais óbitos e menos nascimentos, e a dispersão no funil reflecte sobretudo isso." else ""
    ))
  })

  output$planningFunnelPlot <- plotly::renderPlotly({
    funnel <- planning_funnel()
    local <- planning_location()
    highlight <- c(local, planning_available_comparators()$area[planning_available_comparators()$level %in% c("ULS", funnel$kind)])
    if (planning_funnel_smoothed()) {
      return(planning_eb_chart(funnel$data, planning_spec(), funnel$year, planning_benchmark(), highlight = highlight,
                               denominator_label = funnel$model$denominator_label))
    }
    planning_funnel_chart(funnel$data, planning_spec(), funnel$year, planning_benchmark(), highlight = highlight,
                          denominator_label = funnel$model$denominator_label)
  })

  planning_pyramid_year <- reactive({
    available <- snapshot_years_for("population")
    available <- available[available <= max(planning_year_range())]
    validate(need(length(available) > 0, "Sem população para os anos escolhidos."))
    max(available)
  })

  output$planningPyramidPlot <- plotly::renderPlotly({
    year <- planning_pyramid_year()
    lookup <- active_nuts_lookup()
    local <- planning_location()
    own <- planning_pyramid(local, year, lookup)
    validate(need(nrow(own) > 0, "Sem população para este local."))
    comparators <- planning_areas_all_comparators()
    comparator <- if (nrow(comparators) > 1) comparators$area[[2]] else if (!identical(local, "Portugal")) planning_benchmark() else NULL
    other <- if (is.null(comparator)) NULL else planning_pyramid(comparator, year, lookup)
    planning_pyramid_chart(own, other, paste0(local, " ", year), if (is.null(comparator)) NULL else paste0(comparator, " ", year))
  })

  planning_proportional_view <- reactive({
    under75 <- identical(input$planning_proportional_ages, "under75")
    available <- if (under75) planning_under75_years() else planning_proportional_years()
    available <- available[available <= max(planning_year_range())]
    validate(need(length(available) > 0, "A mortalidade proporcional precisa de três anos de óbitos até ao ano escolhido."))
    end_year <- max(available)
    areas <- planning_areas_all_comparators()$area
    table <- if (under75) {
      planning_under75_table(areas, end_year, lookup = active_nuts_lookup(), vintage = active_nuts_vintage(), mode = planning_split_mode())
    } else {
      planning_proportional_table(areas, end_year, lookup = active_nuts_lookup(), mode = planning_split_mode())
    }
    list(end_year = end_year, under75 = under75, table = table)
  })

  output$planningProportionalNote <- renderUI({
    view <- planning_proportional_view()
    note <- paste0(
      if (view$under75) "Óbitos antes dos 75 anos" else "Todas as idades",
      ", ambos os sexos, triénio ", view$end_year - 2L, "-", view$end_year,
      if (view$under75) " [I46]. " else " [I45]. ",
      "Barras: o local; pontos: os comparadores. «Restantes causas» reúne o que os 13 grandes ",
      "grupos não cobrem, para que o total feche em 100%."
    )
    flagged <- if (view$under75) planning_under75_note(view$table) else NULL
    tagList(helpText(note), if (!is.null(flagged)) helpText(flagged))
  })

  output$planningProportionalPlot <- plotly::renderPlotly({
    view <- planning_proportional_view()
    planning_proportional_chart(view$table, planning_areas_all_comparators())
  })

  output$planningProportionalTable <- renderTable({
    view <- planning_proportional_view()
    areas <- planning_areas_all_comparators()
    columns <- planning_series_name(areas$area, areas$level)
    view$table %>%
      dplyr::mutate(
        column = columns[match(.data$area, areas$area)],
        cell = paste0(planning_format_value(.data$share, 1), "%",
                      if ("flag" %in% names(view$table)) .data$flag else "",
                      " (", planning_format_value(.data$deaths, 0), ")"),
        order = match(.data$code, unique(.data$code))
      ) %>%
      dplyr::select(order, `Grupo de causas` = group, column, cell) %>%
      tidyr::pivot_wider(names_from = column, values_from = cell) %>%
      dplyr::arrange(order) %>%
      dplyr::select(-order)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  planning_cause_view <- reactive({
    # By cause: not the rebuilt latest year, which has all causes only.
    available <- planning_cause_years()
    available <- available[available <= max(planning_year_range())]
    validate(need(length(available) > 0, "A mortalidade por causa precisa de três anos de óbitos e população até ao ano escolhido."))
    end_year <- max(available)
    areas <- planning_areas_all_comparators()
    benchmark <- planning_benchmark()
    sex <- input$planning_cause_sex %||% "HM"
    table <- planning_cause_standardised(areas$area, end_year, lookup = active_nuts_lookup(), benchmark = benchmark, sex = sex, mode = planning_split_mode(),
                                         basis = active_population_basis())
    validate(need(nrow(table) > 0 && any(is.finite(table$smr[table$area == areas$area[[1]]])), "Sem valores para este local neste triénio."))
    list(table = table, end_year = end_year, areas = areas, benchmark = benchmark, sex = sex)
  })

  output$planningCauseNote <- renderUI({
    view <- planning_cause_view()
    sex_label <- c(HM = "ambos os sexos", H = "homens", M = "mulheres")[[view$sex]]
    helpText(HTML(paste(htmltools::htmlEscape(c(
      paste0(
        "Triénio ", view$end_year - 2L, "-", view$end_year, ", ", sex_label, ". A razão padronizada de mortalidade (SMR) compara os óbitos ",
        "observados com os que o local teria com as taxas por idade de ", view$benchmark, " no mesmo triénio (", view$benchmark, " = 100). ",
        "Ponto cheio: o local, laranja acima e azul abaixo de ", view$benchmark, " com significância, cinzento sem diferença; círculos: os comparadores. Escala logarítmica, para que 50 e 200 fiquem à mesma distância de 100."
      ),
      "Taxas padronizadas: População Padrão Europeia de 2013, por 100.000 habitantes; a de menos de 75 anos é a mortalidade prematura. \u2021 mais de 2% dos óbitos sem idade publicada por município foram redistribuídos; as regiões usam as linhas por idade do próprio INE."
    )), collapse = "<br>")))
  })

  output$planningCausePlot <- plotly::renderPlotly({
    view <- planning_cause_view()
    planning_cause_smr_chart(view$table, view$areas, view$benchmark)
  })

  output$planningCauseTable <- renderTable({
    view <- planning_cause_view()
    planning_cause_display(view$table, view$areas)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  # Primary care, from the SNS Transparency portal: the location's ULS (or the
  # ULS that make it up), its ARS and the Continente.
  planning_sns_view <- reactive({
    validate(need(planning_sns_available(), "Sem dados do Portal da Transparência do SNS. Corra tools/fetch_sns.R."))
    areas <- planning_areas_all_comparators()
    resolved <- planning_sns_sets(areas, active_nuts_lookup())
    validate(need(!is.null(planning_sns_area_units(areas$area[[1]], active_nuts_lookup())), paste0(
      "Os dados do SNS são por ULS e ", areas$area[[1]], " não coincide com um conjunto de ULS (nem está numa só). ",
      "Escolha um município, uma ULS, uma ARS ou Portugal."
    )))
    spec <- as.list(SNS_INDICATORS[SNS_INDICATORS$id == (input$planning_sns_indicator %||% "sns_no_gp"), ])
    table <- planning_sns_table(resolved$sets, spec$id)
    every_uls <- planning_uls_units("units")
    ranking <- planning_sns_table(c(stats::setNames(as.list(every_uls), every_uls), list(Continente = every_uls)), spec$id)
    list(spec = spec, table = table, ranking = ranking, levels = resolved$levels, notes = resolved$notes)
  })

  output$planningSnsNote <- renderUI({
    view <- planning_sns_view()
    last <- max(view$table$period)
    helpText(HTML(paste(htmltools::htmlEscape(c(
      paste0(view$spec$label, ": ", SNS_CYCLE_LABELS[[view$spec$cycle]], ".",
             if (view$spec$cycle != "month") paste0(" A linha sobe ao longo de cada ciclo e recomeça; os pontos marcam os fins de ciclo, os únicos valores comparáveis entre si. O ciclo em curso vai até ", planning_sns_period_label(last), ": compare-o com o mesmo mês do ano anterior (tabela).") else ""),
      "Base: utentes inscritos nos cuidados de saúde primários, não residentes. ULS desde Janeiro de 2024; as cinco ULS de Lisboa e do Porto aparecem nos dois agrupamentos exactos. O último mês é provisório.",
      view$notes
    )), collapse = "<br>")))
  })

  output$planningSnsPlot <- plotly::renderPlotly({
    view <- planning_sns_view()
    planning_sns_chart(view$table, view$spec, view$levels)
  })

  output$planningSnsTable <- renderTable({
    view <- planning_sns_view()
    t <- view$table
    complete <- t[t$complete, , drop = FALSE]
    last_complete <- max(complete$period)
    previous <- paste0(as.integer(substr(last_complete, 1, 4)) - 1L, substr(last_complete, 5, 7))
    latest <- max(t$period)
    latest_previous <- paste0(as.integer(substr(latest, 1, 4)) - 1L, substr(latest, 5, 7))
    reference <- complete$value[complete$area == "Continente" & complete$period == last_complete]
    cell <- function(area, period) {
      r <- t[t$area == area & t$period == period, , drop = FALSE]
      if (nrow(r) == 0) return("\u2014")
      paste0(planning_format_value(r$value, 1), "% (", planning_format_value(r$lower, 1), " - ", planning_format_value(r$upper, 1), ")")
    }
    rows <- lapply(names(view$levels), function(area) {
      r <- complete[complete$area == area & complete$period == last_complete, , drop = FALSE]
      sig <- if (area != "Continente" && nrow(r) == 1) planning_significance_mark(planning_significance(r$lower, r$upper, reference)) else ""
      out <- list(Local = area)
      out[[planning_sns_period_label(last_complete)]] <- paste0(cell(area, last_complete), sig)
      out[[planning_sns_period_label(previous)]] <- cell(area, previous)
      if (latest != last_complete) {
        out[[paste0(planning_sns_period_label(latest), " (acumulado)")]] <- cell(area, latest)
        out[[paste0(planning_sns_period_label(latest_previous), " (acumulado)")]] <- cell(area, latest_previous)
      }
      tibble::as_tibble(out)
    })
    dplyr::bind_rows(rows)
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$planningSnsRanking <- plotly::renderPlotly({
    view <- planning_sns_view()
    planning_sns_ranking_chart(view$ranking, view$spec, highlight = names(view$levels)[view$levels %in% c("Local", "ULS")])
  })

  # The indicator on a map: the ACES of 2015-2023 (data/sns_aces_crosswalk.csv)
  # or the ULS from 2024, at the last end of its cycle in the year chosen.
  observe({
    req(planning_sns_available())
    years <- rev(planning_sns_map_years())
    keep <- if (isTRUE(input$planning_sns_map_year %in% years)) input$planning_sns_map_year else years[[1]]
    updateSelectInput(session, "planning_sns_map_year", choices = years, selected = keep)
  })

  planning_sns_map <- reactive({
    validate(need(planning_sns_available(), "Sem dados do Portal da Transparência do SNS. Corra tools/fetch_sns.R."))
    id <- input$planning_sns_indicator %||% sns_ids[[1]]
    year <- suppressWarnings(as.integer(input$planning_sns_map_year))
    req(is.finite(year))
    map <- planning_sns_map_values(id, year)
    validate(need(!is.null(map), "Sem valores deste indicador nesse ano."))
    boundaries <- planning_map_boundaries(map$kind, active_nuts_vintage())
    validate(need(!is.null(boundaries), "Faltam os limites dos ACES (data/boundaries/aces.rds): corra tools/build_boundaries.R levels=aces."))
    map$values$hover <- planning_sns_map_hover(map$values, map$period)
    map$boundaries <- boundaries
    map$spec <- SNS_INDICATORS[SNS_INDICATORS$id == id, ]
    map
  })

  output$planningSnsMap <- plotly::renderPlotly({
    map <- planning_sns_map()
    scale <- input$planning_sns_map_scale %||% "index"
    colours <- if (identical(map$direction, "worse")) PLANNING_MAP_DIVERGING else
      c(low = PLANNING_MAP_DIVERGING[["high"]], mid = PLANNING_MAP_DIVERGING[["mid"]], high = PLANNING_MAP_DIVERGING[["low"]])
    interactive_choropleth(
      map$values, map$boundaries, scale = if (identical(scale, "significance")) "classes" else scale, colours = colours,
      classes = map_class_colours(planning_sns_significance_labels(), map$direction),
      title = paste0(map$spec$label, "<br><sup>", if (identical(map$kind, "aces")) "ACES" else "ULS", ", ", map$period, "</sup>"),
      legend_title = if (identical(scale, "significance")) "Face ao Continente" else "Continente = 100<br>vermelho: pior",
      source = "planning_sns_map"
    )
  })

  output$planningSnsMapNote <- renderUI({
    map <- planning_sns_map()
    helpText(paste0(
      if (identical(map$kind, "aces")) {
        paste0("Até 2023 o portal publica os dados por ACES (e pelas sete ULS que já existiam). ", map$units, " de 55 unidades reportaram em ",
               map$period, ". Os ACES não coincidem com as ULS de 2024: uma ULS pode juntar vários ACES. ")
      } else {
        paste0("Desde 2024, as 39 ULS do Continente (", map$units, " reportaram em ", map$period, "). ")
      },
      "O valor é o do fim do ciclo do indicador nesse ano (", SNS_CYCLE_LABELS[[map$spec$cycle]], "), com a base de inscritos, não de residentes. ",
      "Vermelho onde a unidade está pior do que o Continente, azul onde está melhor. As ilhas não têm dados no portal."
    ))
  })

  # Weekly deaths, the location's smallest NUTS region.
  planning_weekly_view <- reactive({
    validate(need(planning_weekly_available(), "Sem óbitos semanais. Corra tools/fetch_weekly_deaths.R."))
    region <- planning_weekly_region(planning_location(), active_nuts_lookup())
    age <- input$planning_weekly_age %||% "all"
    data <- planning_weekly_data()
    years <- sort(unique(data$year[data$region == region$region]))
    years <- years[years >= max(years) - 2L]
    excess <- planning_weekly_excess(region$region, years, age, active_nuts_lookup())
    validate(need(nrow(excess) > 0, "Sem óbitos semanais para esta região."))
    list(region = region, age = age, excess = excess, summary = planning_weekly_summary(excess))
  })

  output$planningWeeklyNote <- renderUI({
    view <- planning_weekly_view()
    last <- view$excess[view$excess$year == max(view$excess$year), ]
    helpText(HTML(paste(htmltools::htmlEscape(c(
      paste0("Óbitos semanais do INE (0012100) em ", view$region$region, ", até à semana ", max(last$week), " de ", max(last$year),
             ". Os dados das últimas semanas são provisórios e crescem com os registos em atraso."),
      paste0("Esperados: taxas de mortalidade semanais por idade (<65, 65-74, 75-84, 85+) dos anos de base, aplicadas à população do ano, ",
             "com intervalo de previsão de 95%. Anos de base: até cinco anteriores, desde 2023 - depois do excesso da COVID-19 (2020-2022) e na série de ",
             "população revista pelo INE a partir de 2021. 2024 não tem anos de base, e com menos de três anos de base ",
             "mostra-se o valor esperado sem intervalo (dois anos não chegam para estimar a variação entre anos)."),
      view$region$note
    )), collapse = "<br>")))
  })

  output$planningWeeklyPlot <- plotly::renderPlotly({
    view <- planning_weekly_view()
    planning_weekly_chart(view$excess, view$region$region)
  })

  output$planningWeeklyTable <- renderTable({
    s <- planning_weekly_view()$summary
    validate(need(nrow(s) > 0, "Ainda sem anos de base para calcular óbitos esperados."))
    s %>% dplyr::transmute(
      Ano = as.character(.data$year),
      Semanas = paste0("1-", .data$last_week),
      `Óbitos` = planning_format_value(.data$observed, 0),
      Esperados = planning_format_value(.data$expected, 0),
      Excesso = ifelse(is.na(.data$multiplier), planning_format_value(.data$excess, 0),
                       paste0(planning_format_value(.data$excess, 0), " (", planning_format_value(.data$excess_lower, 0), " a ", planning_format_value(.data$excess_upper, 0), ")")),
      `Excesso (%)` = planning_format_value(.data$excess_pct, 1),
      `Anos de base` = .data$baseline
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

  output$planningWeeklyCumulative <- plotly::renderPlotly({
    planning_weekly_cumulative_chart(planning_weekly_view()$excess)
  })

  output$planningDataDate <- renderUI({
    date <- planning_data_date()
    helpText(paste0(
      "Os ficheiros indicam a data de importação dos dados",
      if (is.na(date)) "." else paste0(" (", date, ")."),
      " Definição das regiões: NUTS ", active_nuts_vintage(), "."
    ))
  })

  output$downloadPlanningSelectionXLSX <- downloadHandler(
    filename = function() {
      paste0("indicadores_planeamento_", gsub("[^A-Za-z0-9]+", "_", iconv(planning_location(), to = "ASCII//TRANSLIT")), "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      areas <- planning_areas_all_comparators()
      shiny::withProgress(message = "A preparar o Excel...", value = 0, {
        write_planning_workbook(
          file, areas,
          lookup = active_nuts_lookup(), vintage = active_nuts_vintage(), data_date = planning_data_date(),
          focus = planning_location(), years = planning_year_range(),
          education_min_age = planning_education_age(), benchmark = planning_benchmark(), mode = planning_split_mode(),
          basis = active_population_basis(),
          progress = function(value, detail = NULL) shiny::setProgress(value, detail = detail)
        )
      })
    }
  )

  output$downloadPlanningProfile <- downloadHandler(
    filename = function() {
      paste0("perfil_", gsub("[^A-Za-z0-9]+", "_", iconv(planning_location(), to = "ASCII//TRANSLIT")), "_", Sys.Date(), ".docx")
    },
    content = function(file) {
      shiny::withProgress(message = "A preparar o perfil...", value = 0, {
        write_planning_profile(
          file, planning_areas_all_comparators(),
          lookup = active_nuts_lookup(), vintage = active_nuts_vintage(), data_date = planning_data_date(),
          benchmark = planning_benchmark(), education_min_age = planning_education_age(), mode = planning_split_mode(),
          basis = active_population_basis(),
          last_year = max(planning_year_range()),
          progress = function(value, detail = NULL) shiny::setProgress(value, detail = detail)
        )
      })
    }
  )

  # The full file depends only on the data and the NUTS vintage, so it is built
  # once per import and kept in the cache.
  output$downloadPlanningFullXLSX <- downloadHandler(
    filename = function() paste0("indicadores_planeamento_todas_as_areas_NUTS", active_nuts_vintage(),
                                 if (identical(active_data_mode(), "original")) "_INE_originais" else "", "_", Sys.Date(), ".xlsx"),
    content = function(file) {
      date <- planning_data_date()
      cached <- file.path(
        get_app_dir(), ".mortality-shiny-cache", "exports",
        paste0("indicadores_planeamento_NUTS", active_nuts_vintage(), "_dados_", ifelse(is.na(date), "sem-data", date),
               "_metodo_", PLANNING_METHOD_VERSION, "_populacao_", active_population_basis(),
               if (identical(active_data_mode(), "original")) "_INE_originais" else "", ".xlsx")
      )
      if (!file.exists(cached)) {
        dir.create(dirname(cached), recursive = TRUE, showWarnings = FALSE)
        shiny::withProgress(message = "A gerar o Excel de todas as áreas...", value = 0, {
          write_planning_workbook(
            paste0(cached, ".tmp.xlsx"), planning_full_export_areas(active_nuts_lookup()),
            lookup = active_nuts_lookup(), vintage = active_nuts_vintage(), data_date = date,
            include_long = FALSE, basis = active_population_basis(),
            progress = function(value, detail = NULL) shiny::setProgress(value, detail = detail)
          )
          file.rename(paste0(cached, ".tmp.xlsx"), cached)
        })
      }
      file.copy(cached, file, overwrite = TRUE)
    }
  )

  # -------------------------
  # Advanced Forecasting: Shared Frozen Specification
  # -------------------------
  # Freeze the filtered historical series at click time so the advanced tabs all
  # reflect the same technical specification until the user reruns the model.
  forecast_history <- eventReactive(input$go_forecast, {
    token <- isolate(cancel_seq$forecast)
    year_order <- isolate(year_load_state$forecast_order)
    query_spec <- make_query_spec(
      input$area2,
      input$area_label2,
      input$cause2,
      input$sex2,
      input$data_source2
    )

    shiny::withProgress(message = get_data_source_progress_message(query_spec$data_source), value = 0, {
      metric_bundle <- load_metric_bundle(
        query_spec,
        "forecast",
        token,
        year_range = input$years_fit,
        year_order = year_order
      )
      build_historical_series(
        metric_bundle = metric_bundle,
        series_spec = make_series_spec(
          query_spec = metric_bundle$query_spec,
          population = input$population2,
          rate_type = input$rate_type2,
          year_range = input$years_fit
        )
      )
    })
  }, ignoreNULL = TRUE)

  forecast_sel <- eventReactive(input$go_forecast, {
    run_forecast_models(
      history = forecast_history(),
      model_ids = input$models,
      horizon = input$horizon,
      kind = "forecast",
      token = isolate(cancel_seq$forecast),
      conf_level = input$conf_level2,
      transform_method = input$transform2,
      bias_adjust = value_or_default(input$bias_adjust2, TRUE),
      model_specs = build_advanced_model_specs()
    )
  }, ignoreNULL = TRUE)

  # These thin wrappers keep the downstream renderers simple and make it clear
  # that all advanced sub-tabs are reading from one frozen result object.
  # The frozen fitting result (models fitted on the full selected series).
  advanced_forecast_base <- reactive({
    req(input$go_forecast > 0)
    forecast_sel()
  })

  # Out-of-sample selection runs live off the frozen base so users can switch
  # validation method/test size without refitting the forward forecast.
  advanced_validation <- reactive({
    evaluate_model_selection(
      base_result = advanced_forecast_base(),
      method = value_or_default(input$comparison_validation_mode, "rolling"),
      test_fraction = value_or_default(input$comparison_test_pct, 25) / 100
    )
  })

  # Downstream tabs read the recommended model from here; it now reflects the
  # out-of-sample choice rather than the in-sample fit.
  advanced_forecast_result <- reactive({
    base <- advanced_forecast_base()
    # Keep the other advanced tabs working even if out-of-sample validation
    # fails for this selection; fall back to the in-sample recommendation.
    validation <- tryCatch(advanced_validation(), error = function(e) NULL)
    if (!is.null(validation) && !is.null(validation$recommended_model)) {
      base$recommended_model <- validation$recommended_model
    }
    base
  })

  advanced_forecast_history <- reactive({
    req(input$go_forecast > 0)
    forecast_history()
  })

  # Keep each tab's model picker separate from model specification itself. The
  # selected model here only controls which already-fitted model is foregrounded
  # in the relevant output tab.
  advanced_forecast_focus_model <- reactive({
    resolve_selected_successful_model(
      dat = advanced_forecast_result(),
      selected_model = input$forecast_output_model
    )
  })

  advanced_diagnostic_model <- reactive({
    resolve_selected_successful_model(
      dat = advanced_forecast_result(),
      selected_model = input$diagnostic_model
    )
  })

  # -------------------------
  # Advanced Forecasting: Model Specification
  # -------------------------
  forecast_spec_table <- reactive({
    dat <- advanced_forecast_result()
    spec <- dat$history$spec
    model_rows <- tibble(
      Item = paste("Definições do modelo -", get_model_labels(names(dat$model_specs))),
      Valor = vapply(
        names(dat$model_specs),
        function(model_id) summarize_model_spec(model_id, dat$model_specs[[model_id]]),
        character(1)
      )
    )

    dplyr::bind_rows(
      tibble(
        Item = c(
          "Local de residência",
          "Causa de morte",
          "Sexo",
          "Fonte de dados",
          "População",
          "Taxa",
          "Anos para ajuste",
          "Nível de confiança",
          "Transformação",
          "Modelos seleccionados",
          "Modelos estimados com sucesso",
          "Modelos com falha",
          "Horizonte temporal"
        ),
        Valor = c(
          spec$area_label,
          spec$cause,
          spec$sex,
          get_data_source_label(spec$data_source),
          spec$population,
          dat$history$rate_label,
          format_year_selection(spec$years),
          glue::glue("{dat$conf_level}%"),
          dat$transform_label,
          paste(get_model_labels(dat$selected_model_ids), collapse = ", "),
          if (length(dat$fitted_models) > 0) paste(get_model_labels(dat$fitted_models), collapse = ", ") else "Nenhum",
          if (nrow(dat$failures) > 0) paste(get_model_labels(dat$failures$Model), collapse = ", ") else "Nenhum",
          glue::glue("{dat$horizon} anos")
        )
      ),
      model_rows
    )
  })

  output$forecastSpecTable <- renderTable({
    forecast_spec_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$downloadForecastSpecCSV <- downloadHandler(
    filename = function() {
      paste0("forecast_spec_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(forecast_spec_table(), file)
    }
  )

  # -------------------------
  # Advanced Forecasting: Forecast Output
  # -------------------------
  output$forecastOutputModelSelector <- renderUI({
    dat <- advanced_forecast_result()
    successful_models <- get_successful_model_ids(dat)

    if (length(successful_models) == 0) {
      return(p("Nenhum modelo foi estimado com sucesso para esta especificação."))
    }

    if (!identical(input$forecast_output_view, "single")) {
      return(
        p(glue::glue(
          "A vista comparativa sobrepõe {length(successful_models)} modelo{if (length(successful_models) == 1) '' else 's'} estimado{if (length(successful_models) == 1) '' else 's'} com sucesso."
        ))
      )
    }

    selectInput(
      "forecast_output_model",
      "Modelo:",
      choices = get_named_model_choices(successful_models),
      selected = advanced_forecast_focus_model()
    )
  })

  output$forecastWarnings <- renderUI({
    build_forecast_warning_ui(advanced_forecast_result())
  })

  output$forecastSummaryTable <- renderTable({
    build_forecast_summary_table(
      dat = advanced_forecast_result(),
      view_mode = input$forecast_output_view,
      selected_model = advanced_forecast_focus_model()
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$downloadForecastSummaryCSV <- downloadHandler(
    filename = function() {
      paste0("forecast_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(
        build_forecast_summary_table(
          dat = advanced_forecast_result(),
          view_mode = input$forecast_output_view,
          selected_model = advanced_forecast_focus_model()
        ),
        file
      )
    }
  )

  output$forecastPlot <- plotly::renderPlotly({
    to_interactive_plot(build_advanced_forecast_plot(
      dat = advanced_forecast_result(),
      view_mode = input$forecast_output_view,
      selected_model = advanced_forecast_focus_model()
    ))
  })
  
  output$downloadForecastPlot <- downloadHandler(
    filename = function() {
      dat <- if (input$go_forecast > 0) advanced_forecast_result() else NULL
      build_forecast_download_filename(
        prefix = "forecast_plot",
        extension = ".png",
        dat = dat,
        fallback_area_label = get_selection_label(input$area2, input$area_label2),
        fallback_cause_label = input$cause2,
        view_mode = input$forecast_output_view,
        selected_model = if (input$go_forecast > 0) advanced_forecast_focus_model() else NULL
      )
    },
    content = function(file) {
      dat <- advanced_forecast_result()
      validate(need(length(dat$fits) > 0, "Nenhum modelo foi estimado com sucesso."))

      grDevices::png(file, width = 1200, height = 800, res = 150)
      print(
        build_advanced_forecast_plot(
          dat = dat,
          view_mode = input$forecast_output_view,
          selected_model = advanced_forecast_focus_model()
        )
      )
      grDevices::dev.off()
    }
  )
  
  output$forecastTable <- renderTable({
    dat <- advanced_forecast_result()
    validate(need(length(dat$fits) > 0, "Nenhum modelo foi estimado com sucesso."))
    recent_forecast_rows(build_forecast_display_table(
      dat = dat,
      view_mode = input$forecast_output_view,
      selected_model = advanced_forecast_focus_model()
    ), dat)
  }, sanitize.text.function = identity)
  
  output$downloadForecastCSV <- downloadHandler(
    filename = function() {
      dat <- if (input$go_forecast > 0) advanced_forecast_result() else NULL
      build_forecast_download_filename(
        prefix = "forecast",
        extension = ".csv",
        dat = dat,
        fallback_area_label = get_selection_label(input$area2, input$area_label2),
        fallback_cause_label = input$cause2,
        view_mode = input$forecast_output_view,
        selected_model = if (input$go_forecast > 0) advanced_forecast_focus_model() else NULL
      )
    },
    content = function(file) {
      dat <- advanced_forecast_result()
      validate(need(length(dat$fits) > 0, "Nenhum modelo foi estimado com sucesso."))
      utils::write.csv(
        build_forecast_download_table(
          dat = dat,
          view_mode = input$forecast_output_view,
          selected_model = advanced_forecast_focus_model()
        ),
        file,
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  )

  # -------------------------
  # Advanced Forecasting: Diagnostics
  # -------------------------
  output$diagnosticModelSelector <- renderUI({
    dat <- advanced_forecast_result()
    successful_models <- get_successful_model_ids(dat)

    if (length(successful_models) == 0) {
      return(p("Nenhum modelo foi estimado com sucesso para diagnóstico."))
    }

    selectInput(
      "diagnostic_model",
      "Modelo a inspeccionar:",
      choices = get_named_model_choices(successful_models),
      selected = advanced_diagnostic_model()
    )
  })

  output$diagnosticResidualPlot <- renderPlot({
    build_diagnostic_residual_plot(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model()
    )
  })

  output$downloadDiagnosticResidualPlot <- downloadHandler(
    filename = function() {
      paste0("diagnostic_residuals_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(
        file,
        build_diagnostic_residual_plot(
          dat = advanced_forecast_result(),
          model_id = advanced_diagnostic_model()
        )
      )
    }
  )

  output$diagnosticAcfPlot <- renderPlot({
    plot_diagnostic_correlation(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model(),
      partial = FALSE
    )
  })

  output$downloadDiagnosticAcfPlot <- downloadHandler(
    filename = function() {
      paste0("diagnostic_acf_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_base_plot_png(
        file,
        plot_diagnostic_correlation(
          dat = advanced_forecast_result(),
          model_id = advanced_diagnostic_model(),
          partial = FALSE
        )
      )
    }
  )

  output$diagnosticPacfPlot <- renderPlot({
    plot_diagnostic_correlation(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model(),
      partial = TRUE
    )
  })

  output$downloadDiagnosticPacfPlot <- downloadHandler(
    filename = function() {
      paste0("diagnostic_pacf_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_base_plot_png(
        file,
        plot_diagnostic_correlation(
          dat = advanced_forecast_result(),
          model_id = advanced_diagnostic_model(),
          partial = TRUE
        )
      )
    }
  )

  diagnostic_ljung_table <- reactive({
    build_ljung_box_table(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model()
    )
  })

  output$diagnosticLjungBoxTable <- renderTable({
    diagnostic_ljung_table()
  }, digits = 4, striped = TRUE, bordered = TRUE, spacing = "s")

  output$downloadDiagnosticLjungCSV <- downloadHandler(
    filename = function() {
      paste0("diagnostic_ljung_box_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(diagnostic_ljung_table(), file)
    }
  )

  diagnostic_summary_text <- reactive({
    build_diagnostic_model_summary(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model()
    )
  })

  output$diagnosticModelSummary <- renderText({
    diagnostic_summary_text()
  })

  output$downloadDiagnosticSummaryTXT <- downloadHandler(
    filename = function() {
      paste0("diagnostic_model_summary_", Sys.Date(), ".txt")
    },
    content = function(file) {
      writeLines(diagnostic_summary_text(), file, useBytes = TRUE)
    }
  )

  # -------------------------
  # Advanced Forecasting: Backtesting & Comparison
  # -------------------------
  output$comparisonHoldoutControl <- renderUI({
    if (identical(input$comparison_validation_mode, "insample")) {
      return(
        p("As métricas do ajuste actual usam toda a série seleccionada, sem reservar anos para teste.")
      )
    }

    series_n <- nrow(advanced_forecast_history()$series)
    if (series_n < MIN_VALIDATION_TRAIN + 1L) {
      return(
        p("A validação fora da amostra requer pelo menos 4 anos observados no histórico de ajuste; a recomendação recorre ao ajuste dentro da amostra.")
      )
    }

    sliderInput(
      "comparison_test_pct",
      "Tamanho do teste (% dos anos):",
      min = 10,
      max = 40,
      value = 25,
      step = 5,
      post = "%"
    )
  })

  output$comparisonValidationInfo <- renderText({
    describe_validation_selection(advanced_validation())
  })

  advanced_comparison <- reactive({
    base <- advanced_forecast_base()
    validation <- advanced_validation()
    method <- validation$method_used

    common <- list(
      metrics = validation$metrics,
      ranking = validation$ranking,
      failures = base$failures,
      fitted_models = base$fitted_models,
      selected_model_ids = base$selected_model_ids,
      validation = validation,
      y_label = base$history$y_label
    )

    if (identical(method, "single")) {
      pd <- validation$plot_data
      return(c(common, list(
        mode = "holdout",
        training_obs = pd$training_obs,
        holdout_actual = pd$holdout_actual,
        forecast_df = pd$forecast_df,
        holdout_k = pd$holdout_k
      )))
    }

    if (identical(method, "rolling")) {
      pd <- validation$plot_data
      return(c(common, list(
        mode = "rolling",
        full_obs = pd$full_obs,
        predictions = pd$predictions
      )))
    }

    c(common, list(
      mode = "insample",
      obs = base$obs,
      fits = base$fits,
      inverse = base$transform_inverse
    ))
  })

  output$comparisonWarnings <- renderUI({
    build_forecast_warning_ui(advanced_comparison())
  })

  comparison_ranking_table <- reactive({
    comparison_dat <- advanced_comparison()
    validate(need(nrow(comparison_dat$ranking) > 0, "Nenhum modelo foi estimado com sucesso para a configuração de comparação seleccionada."))

    comparison_dat$ranking %>%
      dplyr::mutate(Model = get_model_labels(Model)) %>%
      dplyr::rename(Modelo = Model)
  })

  output$comparisonRankingTable <- renderTable({
    comparison_ranking_table()
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$downloadComparisonRankingCSV <- downloadHandler(
    filename = function() {
      paste0("comparison_ranking_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(comparison_ranking_table(), file)
    }
  )

  accuracy_display_table <- reactive({
    comparison_dat <- advanced_comparison()
    validate(need(nrow(comparison_dat$metrics) > 0, "Nenhum modelo foi estimado com sucesso para a configuração de comparação seleccionada."))

    comparison_dat$metrics %>%
      dplyr::left_join(comparison_dat$ranking %>% dplyr::select(Model, Classificação), by = "Model") %>%
      dplyr::mutate(Model = get_model_labels(Model)) %>%
      dplyr::rename(Modelo = Model) %>%
      dplyr::select(Classificação, Modelo, dplyr::everything()) %>%
      dplyr::arrange(Classificação)
  })

  output$accuracyTable <- renderTable({
    accuracy_display_table()
  }, digits = 3)

  output$downloadAccuracyCSV <- downloadHandler(
    filename = function() {
      paste0("comparison_metrics_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(accuracy_display_table(), file)
    }
  )

  output$comparisonPlot <- renderPlot({
    build_comparison_plot(advanced_comparison())
  })

  output$downloadComparisonPlot <- downloadHandler(
    filename = function() {
      paste0("comparison_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(file, build_comparison_plot(advanced_comparison()))
    }
  )
  
  # -------------------------
  # Advanced Forecasting: Breaks & Structure
  # -------------------------
  # The breaks tab now reuses the frozen advanced historical series instead of
  # defining a separate set of controls. This keeps structural-break analysis
  # aligned with the exact history chosen in "Model Specification".
  advanced_break_analysis <- reactive({
    analyze_structural_breaks(advanced_forecast_history())
  })

  output$breakInterpretation <- renderUI({
    break_info <- advanced_break_analysis()

    wellPanel(
      h4("Interpretação"),
      p(build_break_interpretation_text(break_info))
    )
  })

  output$breakPlot <- renderPlot({
    build_break_plot(advanced_break_analysis())
  })

  output$downloadBreakPlot <- downloadHandler(
    filename = function() {
      paste0("break_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      save_ggplot_png(file, build_break_plot(advanced_break_analysis()))
    }
  )

  output$breakTable <- renderTable({
    advanced_break_analysis()$segments
  }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 2)

  output$downloadBreakTableCSV <- downloadHandler(
    filename = function() {
      paste0("break_segments_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv_utf8(advanced_break_analysis()$segments, file)
    }
  )
}

# =========================================================
# Run app
# =========================================================

shinyApp(ui, server)
