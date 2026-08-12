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
  "R/metadata.R",
  "R/metrics.R",
  "R/standardisation.R",
  "R/regions.R",
  "R/forecast_helpers.R",
  "R/data_access.R",
  "R/ui_helpers.R"
)) {
  sys.source(file.path(app_dir, app_file), envir = bootstrap_env)
}

# =========================================================
# UI
# =========================================================

ui <- navbarPage(
  title = "PNS Monitorização",

  intro_tab_ui(),
  observed_mortality_tab_ui(),
  beginner_forecasting_tab_ui(),
  advanced_forecasting_tab_ui(),
  annual_metrics_tab_ui(),
  data_availability_tab_ui(),
  glossary_tab_ui()
)

# =========================================================
# Server
# =========================================================

server <- function(input, output, session) {
  cancel_seq <- reactiveValues(rates = 0L, beginner = 0L, forecast = 0L, annual = 0L)

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

  write_csv_utf8 <- function(x, file) {
    utils::write.csv(
      x,
      file,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  save_ggplot_png <- function(file, plot_obj, width = 1200, height = 800, res = 150) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot_obj)
  }

  save_base_plot_png <- function(file, plot_expr, width = 1200, height = 800, res = 150) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(plot_expr)
  }

  # Shared historical-series pipeline:
  # 1) Build a stable query spec with the inputs that determine which raw data
  #    need to be downloaded and aggregated.
  # 2) Load one metric bundle (both population scopes, both rate outputs) for
  #    that query.
  # 3) Apply the final population/rate/year filter to get a reusable series.
  make_query_spec <- function(area, area_label, cause, sex, data_source = "ine") {
    area_key <- sort(unique(area))

    # Selected areas are summed into one geography, so overlapping choices
    # (Portugal with anything, a region with its own municipalities) would be
    # counted twice. Warn rather than block: a user may knowingly want the
    # combination, but should not get it by accident.
    overlap <- overlapping_selection_warning(area_key)
    if (!is.null(overlap)) {
      showNotification(overlap, type = "warning", duration = 20)
    }

    # Expand any selected region into its municipalities. This applies to every
    # tab that goes through a query spec - observed mortality and both forecast
    # tabs - which previously used INE's regional rows directly and so carried
    # the NUTS-2024 discontinuity into trends and projections.
    resolved <- resolve_region_areas(
      areas = area_key,
      region_mode = default_region_mode(),
      available_areas = get_available_areas(data_source)
    )
    for (message_text in resolved$warnings) {
      showNotification(message_text, type = "warning", duration = 15)
    }
    # Label from what the user chose, not from the expansion: a region must stay
    # named "Alentejo", not become a list of 47 municipalities.
    selection_label <- get_selection_label(area_key, area_label)

    list(
      area_key = resolved$areas,
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

  output$downloadSnapshotAvailabilityCSV <- downloadHandler(
    filename = function() paste0("cobertura_rds_", Sys.Date(), ".csv"),
    content = function(file) {
      write_csv_utf8(snapshot_availability_data(), file)
    }
  )

  make_series_spec <- function(query_spec, population, rate_type, year_range = range(year_of_interest)) {
    selected_years <- get_years_in_selected_range(year_range)

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
      showNotification(
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
      showNotification(
        population_gap_message(unusable, "as taxas de mortalidade"),
        type = "warning",
        duration = 15
      )
    }

    notify_snapshot_request_warnings(
      years = years_to_load,
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
        data_source = query_spec$data_source
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

    df_full <- dat$full %>%
      dplyr::filter(sex == query_spec$sex)
    df_trunc <- dat$trunc %>%
      dplyr::filter(sex == query_spec$sex)

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
          dplyr::arrange(Ano)
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
      dplyr::arrange(Ano)
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
          dplyr::arrange(Ano)
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
      dplyr::arrange(Ano)
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

  rank_model_metrics <- function(metric_tbl) {
    if (nrow(metric_tbl) == 0) {
      return(tibble(Classificação = integer(0), Model = character(0), Recomendação = character(0)))
    }

    recommended_model <- choose_recommended_model(metric_tbl)
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

  build_break_interpretation_text <- function(break_info) {
    break_n <- length(break_info$break_years)
    segment_n <- nrow(break_info$segments)

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
        "Foi detectada uma potencial mudança estrutural em torno de {break_info$break_years[[1]]}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que {what_changed} poderá ter mudado nessa altura. {model_note}"
      ))
    }

    glue::glue(
      "Foram detectadas {break_n} potenciais mudanças estruturais em torno de {paste(break_info$break_years, collapse = ', ')}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que {what_changed} poderá não ser estável ao longo do período observado. {model_note}"
    )
  }

  get_simple_structural_warning_text <- function(history) {
    break_info <- safely_analyze_structural_breaks(history)

    if (is.null(break_info) || length(break_info$break_years) == 0) {
      return(NULL)
    }

    if (length(break_info$break_years) == 1) {
      return(glue::glue(
        "Foi detectada uma possível mudança estrutural em torno de {break_info$break_years[[1]]}; as previsões baseadas em todo o histórico devem ser lidas com cautela acrescida."
      ))
    }

    "Foi detectada uma possível mudança estrutural na série histórica; as previsões baseadas em todo o histórico devem ser lidas com cautela acrescida."
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

    parts <- glue::glue("Escolha do modelo: {method_label}.")
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

    wellPanel(
      h4("Resumo em linguagem simples"),
      p(glue::glue(
        "Usando {dat$training_label}, a previsão recomendada sugere que a taxa de mortalidade {trend_text} ao longo dos próximos {dat$horizon} anos."
      )),
      p(glue::glue(
        "Em {final_forecast$year}, a projecção central é de {round(final_forecast$mean, 2)} por 100.000, em comparação com {round(last_observed$value, 2)} em {last_observed$year}."
      )),
      p(glue::glue(
        "Isto corresponde a uma variação de {sprintf('%+.2f', change_abs)} por 100.000{if (!is.na(change_pct)) paste0(' (', sprintf('%+.2f%%', change_pct), ')') else ''}."
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
  make_annual_area_spec <- function(label, areas, region_mode, available_areas) {
    resolved <- resolve_region_areas(
      areas = areas,
      region_mode = region_mode,
      available_areas = available_areas
    )

    list(
      label = label,
      areas = resolved$areas,
      expanded = resolved$expanded,
      warnings = resolved$warnings
    )
  }

  get_annual_area_specs <- function(selected_areas,
                                    custom_label = NULL,
                                    region_mode = "original",
                                    available_areas = NULL) {
    selected_areas <- setdiff(sort(unique(as.character(selected_areas))), c("Portugal", "Norte"))

    validate(
      need(length(selected_areas) > 0, "Seleccione pelo menos um local adicional para a terceira coluna.")
    )

    list(
      # Portugal is a national total under every NUTS vintage, so it is never
      # rebuilt from municipalities.
      list(label = "Portugal", areas = "Portugal", expanded = character(0), warnings = character(0)),
      make_annual_area_spec("Norte", "Norte", region_mode, available_areas),
      make_annual_area_spec(
        get_selection_label(selected_areas, custom_label),
        selected_areas,
        region_mode,
        available_areas
      )
    )
  }

  # Years contributing to a pooled value centred on `year`, clipped to the
  # years the app knows about so an edge window is truncated rather than
  # requesting data that cannot exist.
  get_pooled_years <- function(year, window, metric_id = NULL) {
    window <- normalize_pooling_window(window)
    half <- (window - 1L) %/% 2L
    candidate <- seq.int(as.integer(year) - half, as.integer(year) + half)
    years <- sort(intersect(candidate, as.integer(year_of_interest)))

    # A rate pools deaths against person-years, so a year contributing deaths
    # without a population estimate would inflate the numerator against a
    # denominator that never included it. Those years are dropped from the
    # window rather than failing the whole request: a window centred on 2022
    # reaches 2024, which has deaths but no published population, and refusing
    # outright would disable pooling exactly where it is most useful. The
    # period label reports the span actually used.
    if (!is.null(metric_id) && metric_id %in% annual_metrics_needing_population) {
      usable <- intersect(years, as.integer(population_years))
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
    data <- if (metric_id %in% annual_metrics_needing_population) {
      get_data_for(
        area = area_spec$areas,
        cause = cause,
        years = years,
        year_order = "desc",
        data_source = data_source
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

    data %>%
      dplyr::filter(sex == .env$sex)
  }

  # Collapse a multi-year, multi-area selection to one row per age band, which
  # is the shape every metric below consumes. Deaths and population are summed,
  # so a pooled denominator is person-years and a multi-area selection is one
  # combined geography - the same convention the app already uses for areas.
  collapse_annual_cause_data <- function(cause_data) {
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

    all_cause_deaths <- NA_real_
    all_cause_source <- NA_character_
    if (identical(metric_id, "proportional")) {
      all_cause_data <- get_death_data_for(
        area = area_spec$areas,
        cause = "Todas as causas de morte",
        years = years,
        year_order = "desc",
        data_source = data_source
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
        need(nrow(cause_data) > 0, glue::glue("Sem dados para {area_spec$label} / {cause}."))
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
            available_areas = get_available_areas(data_source)
          )

          reference_bands <- collapse_annual_cause_data(
            load_annual_cause_data(
              area_spec = list(label = "referência", areas = reference_resolved$areas),
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
        ypll = {
          ci <- compute_ypll_interval(cause_data, cutoff = 70)
          list(
            value = unname(ci[["estimate"]]),
            lower = unname(ci[["lower"]]),
            upper = unname(ci[["upper"]]),
            source_detail = source_summary$death_source
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
      count_metric <- metric_id %in% c("deaths", "ypll")
      annualise <- if (count_metric && length(years) > 1L) length(years) else 1L

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
        source_detail = as.character(metric_values$source_detail)
      )
    }))
  }

  round_annual_metric_value <- function(metric_id, value) {
    dplyr::case_when(
      metric_id %in% c("deaths", "ypll") ~ round(value, 0),
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

  build_annual_metrics_table <- function(metric_long, selected_metric) {
    cause_order <- get_annual_cause_order(metric_long)
    sort_location <- get_annual_sort_location(metric_long)
    location_order <- unique(metric_long$location)

    metric_long %>%
      dplyr::filter(metric_id == selected_metric) %>%
      dplyr::mutate(
        cause = factor(cause, levels = cause_order),
        display_value = purrr::pmap_chr(
          list(metric_id, value, lower, upper),
          format_annual_metric_interval
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
    plot_df <- metric_long %>%
      dplyr::filter(metric_id == selected_metric) %>%
      dplyr::mutate(
        cause = factor(cause, levels = rev(cause_order)),
        value = round_annual_metric_value(metric_id, value),
        lower = round_annual_metric_value(metric_id, lower),
        upper = round_annual_metric_value(metric_id, upper)
      )

    validate(need(nrow(plot_df) > 0, "Seleccione pelo menos uma causa de morte para apresentar."))

    period_text <- if ("period" %in% names(plot_df)) {
      periods <- unique(as.character(plot_df$period))
      n_years <- if ("n_years" %in% names(plot_df)) max(plot_df$n_years, na.rm = TRUE) else 1L
      if (length(periods) == 1 && n_years > 1) {
        basis <- if (identical(selected_metric, "deaths") || identical(selected_metric, "ypll")) {
          "média anual"
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

    plot +
      coord_flip() +
      labs(
        title = paste(unique(plot_df$metric), "por causa de morte"),
        subtitle = period_text,
        x = NULL,
        y = unique(plot_df$metric)
      ) +
      scale_fill_brewer(palette = "Set2", name = "Local") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
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

  output$ratePlot <- plotly::renderPlotly({
    to_interactive_plot(build_observed_rate_plot(observed_history()))
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
    token <- isolate(cancel_seq$beginner)
    year_order <- isolate(year_load_state$beginner_order)
    query_spec <- make_query_spec(
      input$beginner_area,
      input$beginner_area_label,
      input$beginner_cause,
      input$beginner_sex,
      input$beginner_data_source
    )

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
  }, ignoreNULL = TRUE)

  beginner_forecast <- eventReactive(input$go_beginner_forecast, {
    training_history <- beginner_training_history()

    shiny::withProgress(message = "A ajustar e validar modelos...", value = 0, {
      guided_result <- run_forecast_models(
        history = training_history,
        model_ids = unname(forecast_model_choices),
        horizon = input$beginner_horizon,
        kind = "beginner",
        token = isolate(cancel_seq$beginner)
      )
      incProgress(0.5)

      validation <- evaluate_model_selection(
        base_result = guided_result,
        method = value_or_default(input$beginner_validation, "rolling"),
        test_fraction = value_or_default(input$beginner_test_pct, 25) / 100
      )
      if (!is.null(validation$recommended_model)) {
        guided_result$recommended_model <- validation$recommended_model
      }
      incProgress(0.5)

      c(
        guided_result,
        list(
          full_history = training_history,
          horizon = input$beginner_horizon,
          mode = input$beginner_mode,
          training_label = get_beginner_training_label(training_history),
          validation = validation
        )
      )
    })
  }, ignoreNULL = TRUE)

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

  output$beginnerForecastTable <- renderTable({
    req(input$go_beginner_forecast > 0)
    beginner_forecast_table()
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
      showNotification(message_text, type = "warning", duration = 15)
    }
    vintage_warning <- region_vintage_warning(
      areas = unique(c("Norte", as.character(input$annual_area))),
      years = pooled_years,
      region_mode = region_mode
    )
    if (!is.null(vintage_warning)) {
      showNotification(vintage_warning, type = "warning", duration = 20)
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
            dplyr::mutate(sort_location = area_specs[[3]]$label)

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
    build_forecast_display_table(
      dat = dat,
      view_mode = input$forecast_output_view,
      selected_model = advanced_forecast_focus_model()
    )
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
