# =========================================================
# UI Helpers
# =========================================================

year_range_slider <- function(input_id, label, value = range(year_of_interest)) {
  sliderInput(
    input_id,
    label,
    min   = min(year_of_interest),
    max   = max(year_of_interest),
    value = value,
    step  = 1,
    sep   = ""
  )
}

data_source_input <- function(input_id) {
  selectInput(
    input_id,
    "Fonte de dados:",
    choices = data_source_choices,
    selected = normalize_data_source(Sys.getenv("MORTALITY_DEFAULT_DATA_SOURCE", "ine"))
  )
}

# Shared control panels ----------------------------------------------------

forecast_controls_panel <- function() {
  tagList(
    selectInput("area2", "Local de residência:", choices = local_area, multiple = TRUE, selected = get_default_area_selection()),
    textInput("area_label2", "Nome da selecção (opcional):", placeholder = "Ex.: AML"),
    selectInput("cause2", "Causa de Morte:", choices = diseases),
    selectInput("sex2", "Sexo:", choices = sex_levels, selected = "HM"),
    radioButtons(
      "population2",
      "População:",
      choices = c("Total", "Menos de 75 anos")
    ),
    radioButtons(
      "rate_type2",
      "Taxa:",
      choices = c("Bruta" = "crude", "Padronizada" = "dsr")
    ),
    year_range_slider(
      "years_fit",
      "Anos a importar / ajustar:"
    ),
    data_source_input("data_source2"),
    checkboxGroupInput(
      "models",
      "Famílias de modelos:",
      choices = forecast_model_choices,
      selected = c("arima", "ets")
    ),
    sliderInput(
      "horizon",
      "Horizonte de projecção (anos):",
      min   = 1,
      max   = 8,
      value = 7
    ),
    sliderInput(
      "conf_level2",
      "Nível de confiança (%):",
      min = 80,
      max = 99,
      value = 95,
      step = 1
    ),
    selectInput(
      "transform2",
      "Transformação:",
      choices = c(
        "Transformação log com offset" = "log_offset",
        "Sem transformação" = "none"
      ),
      selected = "log_offset"
    ),
    uiOutput("advancedModelParameterPanels"),
    actionButton("go_forecast", "Carregar projecções"),
    br(), br(),
    actionButton("cancel_forecast", "Interromper carregamento")
  )
}

beginner_forecast_controls_panel <- function() {
  tagList(
    p("Esta previsão guiada utiliza a série actualmente carregada em 'Mortalidade Observada'."),
    sliderInput(
      "beginner_horizon",
      "Horizonte de projecção (anos):",
      min = 1,
      max = 8,
      value = 5
    ),
    year_range_slider(
      "beginner_years_fit",
      "Janela de ajuste:"
    ),
    radioButtons(
      "beginner_mode",
      "Modo:",
      choices = c(
        "Previsão recomendada" = "recommended",
        "Comparar modelos" = "compare"
      ),
      selected = "recommended"
    ),
    actionButton("go_beginner_forecast", "Gerar previsão")
  )
}

forecast_selection_note_ui <- function() {
  wellPanel(
    p("Estas vistas avançadas utilizam a especificação técnica definida em 'Especificação do Modelo'."),
    p("Execute 'Carregar projecções' nesse separador para preencher os resultados abaixo.")
  )
}

data_availability_tab_ui <- function() {
  tabPanel(
    "Disponibilidade de Dados",
    sidebarLayout(
      sidebarPanel(
        radioButtons(
          "availability_dataset",
          "Conjunto:",
          choices = c("Óbitos" = "deaths", "População" = "population"),
          selected = "deaths",
          inline = TRUE
        ),
        year_range_slider(
          "availability_years",
          "Anos:"
        ),
        conditionalPanel(
          "input.availability_dataset == 'deaths'",
          selectInput(
            "availability_cause",
            "Causa de Morte:",
            choices = diseases,
            selected = if ("Todas as causas de morte" %in% diseases) "Todas as causas de morte" else utils::head(diseases, 1),
            multiple = TRUE
          )
        ),
        selectInput(
          "availability_area",
          "Local:",
          choices = local_area,
          multiple = TRUE,
          selected = c("Portugal", "Norte")
        ),
        checkboxInput(
          "availability_show_missing",
          "Mostrar indisponíveis",
          value = TRUE
        )
      ),
      mainPanel(
        h4("Resumo RDS"),
        tableOutput("snapshotInventorySummary"),
        downloadButton("downloadSnapshotInventorySummaryCSV", "Descarregar resumo (CSV)"),
        br(), br(),
        h4("Cobertura Seleccionada"),
        tableOutput("snapshotAvailabilityTable"),
        downloadButton("downloadSnapshotAvailabilityCSV", "Descarregar cobertura (CSV)")
      )
    )
  )
}

observed_mortality_tab_ui <- function() {
  tabPanel(
    "Mortalidade Observada",
    sidebarLayout(
      sidebarPanel(
        selectInput("area", "Local de residência:", choices = local_area, multiple = TRUE, selected = get_default_area_selection()),
        textInput("area_label", "Nome da selecção (opcional):", placeholder = "Ex.: AML"),
        selectInput("cause", "Causa de Morte:", choices = diseases),
        selectInput("sex", "Sexo:", choices = sex_levels, selected = "HM"),
        radioButtons(
          "population",
          "População:",
          choices = c("Total", "Menos de 75 anos")
        ),
        radioButtons(
          "rate_type",
          "Taxa:",
          choices = c("Bruta" = "crude", "Padronizada" = "dsr")
        ),
        year_range_slider(
          "years_import",
          "Anos a importar:"
        ),
        data_source_input("data_source"),
        actionButton("go_rates", "Carregar dados"),
        br(), br(),
        actionButton("cancel_rates", "Interromper carregamento")
      ),
      mainPanel(
        h4("Resumo"),
        tableOutput("rateSummaryTable"),
        downloadButton("downloadRateSummaryCSV", "Descarregar resumo (CSV)"),
        br(),
        plotOutput("ratePlot", height = "400px"),
        br(),
        downloadButton("downloadRatePlot", "Descarregar gráfico (PNG)"),
        br(),
        h4("Série anual observada"),
        tableOutput("rateTable"),
        downloadButton("downloadRateTableCSV", "Descarregar tabela (CSV)")
      )
    )
  )
}

beginner_forecasting_tab_ui <- function() {
  tabPanel(
    "Previsão Guiada",
    sidebarLayout(
      sidebarPanel(
        beginner_forecast_controls_panel()
      ),
      mainPanel(
        plotOutput("beginnerForecastPlot", height = "400px"),
        br(),
        downloadButton("downloadBeginnerForecastPlot", "Descarregar gráfico (PNG)"),
        br(),
        uiOutput("beginnerForecastSummary"),
        br(),
        uiOutput("beginnerForecastReliability"),
        br(),
        h4("Tabela da previsão"),
        tableOutput("beginnerForecastTable"),
        downloadButton("downloadBeginnerForecastCSV", "Descarregar tabela (CSV)")
      )
    )
  )
}

annual_metrics_tab_ui <- function() {
  tabPanel(
    "Métricas Anuais",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "annual_year",
          "Ano:",
          choices = sort(year_of_interest),
          selected = max(year_of_interest)
        ),
        selectInput(
          "annual_cause",
          "Causa de Morte:",
          choices = diseases,
          selected = if ("Todas as causas de morte" %in% diseases) "Todas as causas de morte" else utils::head(diseases, 1),
          multiple = TRUE
        ),
        selectInput("annual_sex", "Sexo:", choices = sex_levels, selected = "HM"),
        selectInput(
          "annual_area",
          "Local adicional:",
          choices = setdiff(local_area, c("Portugal", "Norte")),
          multiple = TRUE,
          selected = character(0)
        ),
        textInput("annual_area_label", "Nome da selecção (opcional):", placeholder = "Ex.: ACES / ULS"),
        selectInput(
          "annual_metric",
          "Métrica:",
          choices = annual_metric_choices,
          selected = "deaths"
        ),
        data_source_input("annual_data_source"),
        actionButton("go_annual_metrics", "Carregar métricas"),
        br(), br(),
        actionButton("cancel_annual_metrics", "Interromper carregamento")
      ),
      mainPanel(
        tableOutput("annualMetricsTable"),
        br(),
        downloadButton("downloadAnnualMetricsCSV", "Descarregar tabela (CSV)"),
        br(), br(),
        plotOutput("annualMetricsPlot", height = "420px"),
        br(),
        downloadButton("downloadAnnualMetricsPlot", "Descarregar gráfico (PNG)")
      )
    )
  )
}

advanced_model_spec_tab_ui <- function() {
  tabPanel(
    "Especificação do Modelo",
    fluidRow(
      column(
        width = 4,
        wellPanel(
          forecast_controls_panel()
        )
      ),
      column(
        width = 8,
        forecast_selection_note_ui(),
        tableOutput("forecastSpecTable"),
        downloadButton("downloadForecastSpecCSV", "Descarregar especificação (CSV)")
      )
    )
  )
}

advanced_forecast_output_tab_ui <- function() {
  tabPanel(
    "Resultados da Previsão",
    forecast_selection_note_ui(),
    fluidRow(
      column(
        width = 4,
        radioButtons(
          "forecast_output_view",
          "Vista:",
          choices = c(
            "Modelo único" = "single",
            "Comparar modelos ajustados" = "compare"
          ),
          selected = "single",
          inline = TRUE
        )
      ),
      column(
        width = 8,
        uiOutput("forecastOutputModelSelector")
      )
    ),
    uiOutput("forecastWarnings"),
    br(),
    tableOutput("forecastSummaryTable"),
    downloadButton("downloadForecastSummaryCSV", "Descarregar resumo (CSV)"),
    br(),
    plotOutput("forecastPlot", height = "400px"),
    br(),
    downloadButton("downloadForecastPlot", "Descarregar gráfico (PNG)"),
    br(), br(),
    tableOutput("forecastTable"),
    br(),
    downloadButton("downloadForecastCSV", "Descarregar tabela (CSV)")
  )
}

advanced_diagnostics_tab_ui <- function() {
  tabPanel(
    "Diagnóstico",
    forecast_selection_note_ui(),
    uiOutput("diagnosticModelSelector"),
    br(),
    plotOutput("diagnosticResidualPlot", height = "260px"),
    downloadButton("downloadDiagnosticResidualPlot", "Descarregar resíduos (PNG)"),
    br(),
    fluidRow(
      column(
        width = 6,
        plotOutput("diagnosticAcfPlot", height = "260px"),
        downloadButton("downloadDiagnosticAcfPlot", "Descarregar ACF (PNG)")
      ),
      column(
        width = 6,
        plotOutput("diagnosticPacfPlot", height = "260px"),
        downloadButton("downloadDiagnosticPacfPlot", "Descarregar PACF (PNG)")
      )
    ),
    br(),
    h4("Teste de Ljung-Box"),
    tableOutput("diagnosticLjungBoxTable"),
    downloadButton("downloadDiagnosticLjungCSV", "Descarregar Ljung-Box (CSV)"),
    br(),
    h4("Resumo do Modelo"),
    verbatimTextOutput("diagnosticModelSummary"),
    downloadButton("downloadDiagnosticSummaryTXT", "Descarregar resumo (TXT)")
  )
}

advanced_backtesting_tab_ui <- function() {
  tabPanel(
    "Retroteste e Comparação",
    forecast_selection_note_ui(),
    fluidRow(
      column(
        width = 6,
        radioButtons(
          "comparison_validation_mode",
          "Abordagem de validação:",
          choices = c(
            "Métricas do ajuste actual" = "insample",
            "Validação nos últimos k anos" = "holdout"
          ),
          selected = "insample",
          inline = TRUE
        )
      ),
      column(
        width = 6,
        uiOutput("comparisonHoldoutControl")
      )
    ),
    uiOutput("comparisonWarnings"),
    br(),
    h4("Classificação"),
    tableOutput("comparisonRankingTable"),
    downloadButton("downloadComparisonRankingCSV", "Descarregar classificação (CSV)"),
    br(),
    h4("Valores das Métricas"),
    tableOutput("accuracyTable"),
    downloadButton("downloadAccuracyCSV", "Descarregar métricas (CSV)"),
    br(),
    plotOutput("comparisonPlot", height = "360px"),
    downloadButton("downloadComparisonPlot", "Descarregar gráfico (PNG)")
  )
}

advanced_breaks_tab_ui <- function() {
  tabPanel(
    "Quebras e Estrutura",
    forecast_selection_note_ui(),
    uiOutput("breakInterpretation"),
    br(),
    plotOutput("breakPlot", height = "400px"),
    downloadButton("downloadBreakPlot", "Descarregar gráfico (PNG)"),
    br(),
    tableOutput("breakTable"),
    downloadButton("downloadBreakTableCSV", "Descarregar tabela (CSV)")
  )
}

advanced_forecasting_tab_ui <- function() {
  tabPanel(
    "Previsão Avançada",
    tabsetPanel(
      advanced_model_spec_tab_ui(),
      advanced_forecast_output_tab_ui(),
      advanced_diagnostics_tab_ui(),
      advanced_backtesting_tab_ui(),
      advanced_breaks_tab_ui()
    )
  )
}
