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
    selected = get_default_data_source()
  )
}

# Plain-language help notes shown under controls, for users unfamiliar with the
# underlying statistics. Kept as small helpers so the wording stays consistent
# across the observed, guided and advanced panels.
rate_type_help <- function() {
  helpText(
    "Bruta: mortes por 100.000 habitantes. Padronizada: ajustada à idade, ",
    "para comparar de forma justa locais com estruturas etárias diferentes."
  )
}
population_help <- function() {
  helpText(
    "'Menos de 75 anos' foca a mortalidade prematura (potencialmente evitável) ",
    "e não representa a mortalidade total."
  )
}
horizon_help <- function() {
  helpText("Número de anos a projectar para o futuro. Quanto maior o horizonte, maior a incerteza.")
}
models_help <- function() {
  helpText(
    "Cada família é um método estatístico diferente. Em caso de dúvida, use a ",
    "'Previsão Guiada', que escolhe um método por si."
  )
}
confidence_help <- function() {
  helpText("Largura do intervalo de incerteza apresentado. 95% é o valor habitual.")
}
transform_help <- function() {
  helpText(
    "A transformação log estabiliza séries positivas e evita previsões negativas; ",
    "'Sem transformação' modela a taxa directamente."
  )
}
beginner_validation_help <- function() {
  helpText(
    "A aplicação escolhe o melhor método testando a previsão em anos recentes ",
    "reservados para avaliação. 'Validação móvel' é a opção mais fiável. ",
    "Se não souber, mantenha as predefinições."
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
    population_help(),
    radioButtons(
      "rate_type2",
      "Taxa:",
      choices = c("Bruta" = "crude", "Padronizada" = "dsr")
    ),
    rate_type_help(),
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
    models_help(),
    sliderInput(
      "horizon",
      "Horizonte de projecção (anos):",
      min   = 1,
      max   = 30,
      value = 7
    ),
    horizon_help(),
    sliderInput(
      "conf_level2",
      "Nível de confiança (%):",
      min = 80,
      max = 99,
      value = 95,
      step = 1
    ),
    confidence_help(),
    selectInput(
      "transform2",
      "Transformação:",
      choices = c(
        "Transformação log com offset" = "log_offset",
        "Sem transformação" = "none"
      ),
      selected = "log_offset"
    ),
    transform_help(),
    uiOutput("advancedModelParameterPanels"),
    actionButton("go_forecast", "Carregar projecções"),
    br(), br(),
    actionButton("cancel_forecast", "Interromper carregamento")
  )
}

beginner_forecast_controls_panel <- function() {
  tagList(
    selectInput(
      "beginner_area",
      "Local de residência:",
      choices = local_area,
      multiple = TRUE,
      selected = get_default_area_selection()
    ),
    textInput("beginner_area_label", "Nome da selecção (opcional):", placeholder = "Ex.: AML"),
    selectInput("beginner_cause", "Causa de Morte:", choices = diseases),
    selectInput("beginner_sex", "Sexo:", choices = sex_levels, selected = "HM"),
    radioButtons(
      "beginner_population",
      "População:",
      choices = c("Total", "Menos de 75 anos")
    ),
    population_help(),
    radioButtons(
      "beginner_rate_type",
      "Taxa:",
      choices = c("Bruta" = "crude", "Padronizada" = "dsr")
    ),
    rate_type_help(),
    data_source_input("beginner_data_source"),
    sliderInput(
      "beginner_horizon",
      "Horizonte de projecção (anos):",
      min = 1,
      max = 30,
      value = 5
    ),
    horizon_help(),
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
    # Keep the default path simple: the model-selection controls are hidden
    # until the user opts in. Their defaults (rolling, 25%) still apply.
    checkboxInput("beginner_show_advanced", "Mostrar opções avançadas", value = FALSE),
    conditionalPanel(
      "input.beginner_show_advanced == true",
      radioButtons(
        "beginner_validation",
        "Como escolher o modelo recomendado:",
        choices = c(
          "Validação móvel (recomendada)" = "rolling",
          "Divisão única treino/teste" = "single",
          "Ajuste dentro da amostra" = "insample"
        ),
        selected = "rolling"
      ),
      conditionalPanel(
        "input.beginner_validation != 'insample'",
        sliderInput(
          "beginner_test_pct",
          "Tamanho do teste (% dos anos):",
          min = 10,
          max = 40,
          value = 25,
          step = 5,
          post = "%"
        )
      ),
      beginner_validation_help()
    ),
    actionButton("go_beginner_forecast", "Gerar previsão"),
    br(), br(),
    actionButton("cancel_beginner_forecast", "Interromper carregamento")
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

intro_tab_ui <- function() {
  tabPanel(
    "Introdução",
    fluidRow(
      column(
        width = 10, offset = 1,
        h2("Bem-vindo"),
        p(
          "Esta aplicação ajuda a explorar dados de mortalidade em Portugal, a partir dos ",
          "indicadores do Instituto Nacional de Estatística (INE). Pode ver como a mortalidade ",
          "evoluiu ao longo dos anos, comparar territórios e obter previsões simples para os ",
          "próximos anos — sem precisar de saber estatística ou programação."
        ),

        h3("O que é uma previsão (e o que não é)"),
        p(
          "Uma previsão, ou projecção, é uma estimativa de como uma taxa de mortalidade poderá ",
          "evoluir no futuro, com base no padrão dos anos anteriores."
        ),
        tags$ul(
          tags$li("Não é uma certeza: é um cenário possível, não o que vai necessariamente acontecer."),
          tags$li("Não é uma meta nem um número oficial."),
          tags$li(
            "A incerteza aumenta com o tempo: os primeiros anos são mais fiáveis do que, ",
            "por exemplo, daqui a 20 ou 30 anos."
          )
        ),
        p(tags$strong(
          "Interprete os resultados como apoio à exploração, e não como conclusões definitivas."
        )),

        h3("Começar em 3 passos"),
        tags$ol(
          tags$li(
            tags$strong("Escolha o que quer ver. "),
            "No painel à esquerda, seleccione um local de residência (por exemplo, um concelho ou ",
            tags$em("Portugal"), ") e uma causa de morte."
          ),
          tags$li(
            tags$strong("Veja a evolução histórica. "),
            "No separador ", tags$em("Mortalidade Observada"), ", clique em ",
            tags$em("Carregar dados"), " para ver o gráfico dos últimos anos."
          ),
          tags$li(
            tags$strong("Peça uma previsão simples. "),
            "No separador ", tags$em("Previsão Guiada"), ", clique em ", tags$em("Gerar previsão"),
            ". A aplicação escolhe automaticamente um método adequado e explica o resultado em ",
            "linguagem simples, incluindo o grau de fiabilidade."
          )
        ),

        h3("Que separador usar"),
        tags$ul(
          tags$li(tags$strong("Mortalidade Observada"), " — ver a evolução histórica de uma taxa e exportar tabelas e gráficos."),
          tags$li(tags$strong("Previsão Guiada"), " — obter uma previsão simples, com explicação. Recomendado para a maioria dos utilizadores."),
          tags$li(tags$strong("Previsão Avançada"), " — controlar modelos, diagnósticos e testes. Destinado a utilizadores técnicos; não é necessário para uma previsão simples."),
          tags$li(tags$strong("Métricas Anuais"), " — comparar Portugal, Norte e um local à escolha num único ano."),
          tags$li(tags$strong("Disponibilidade de Dados"), " — verificar que anos, locais e causas existem nos ficheiros antes de carregar.")
        ),

        h3("Alguns termos que vai encontrar"),
        tags$ul(
          tags$li(tags$strong("Taxa bruta:"), " número de mortes por 100.000 habitantes."),
          tags$li(tags$strong("Taxa padronizada:"), " taxa ajustada à idade, que permite comparar de forma justa locais com populações mais jovens ou mais envelhecidas."),
          tags$li(tags$strong("Intervalo de confiança:"), " a margem de incerteza à volta de um valor (a zona sombreada nos gráficos)."),
          tags$li(tags$strong("Horizonte:"), " quantos anos para o futuro a previsão vai.")
        ),
        p("Para mais termos explicados em linguagem simples, consulte o separador ", tags$strong("Glossário"), "."),

        h3("Uma dica sobre a fonte de dados"),
        p(
          "No início de cada análise pode escolher a ", tags$em("Fonte de dados"), ". Comece com ",
          tags$em("Ficheiros RDS"), ", que lê dados já preparados e é muito mais rápido. Use ",
          tags$em("INE em directo"), " apenas se precisar de dados que os ficheiros não contêm — pode ser lento."
        ),

        wellPanel(
          p(tags$em(
            "Ferramenta não oficial e exploratória. Não substitui validação epidemiológica, ",
            "análise clínica ou produtos estatísticos oficiais. Para um guia detalhado, ",
            "separador a separador, consulte o Manual do Utilizador (USER_MANUAL.md)."
          ))
        )
      )
    )
  )
}

glossary_tab_ui <- function() {
  # term/definition pairs grouped into sections; rendered as a definition list.
  glossary_section <- function(title, entries) {
    tagList(
      h3(title),
      tags$dl(
        do.call(tagList, lapply(entries, function(e) {
          tagList(tags$dt(e[[1]]), tags$dd(e[[2]]))
        }))
      )
    )
  }

  mortality_terms <- list(
    list("Óbitos", "Número absoluto de mortes na selecção (ano, local, causa, sexo e idade)."),
    list("Taxa bruta", "Número de mortes por 100.000 habitantes. É simples, mas depende muito da idade da população."),
    list("Taxa padronizada", "Taxa ajustada à idade, que permite comparar de forma justa locais com populações mais jovens ou mais envelhecidas. Usa a População Padrão Europeia de 2013."),
    list("População padrão (ESP 2013)", "Uma estrutura etária de referência comum, aplicada na padronização para que as comparações não sejam distorcidas pela idade."),
    list("Mortalidade proporcional", "Percentagem das mortes de uma causa face ao total de mortes, no mesmo ano, sexo e local."),
    list("AVPP (anos de vida potencialmente perdidos)", "Medida do impacto da morte prematura: soma os anos que faltavam até aos 70 em cada morte antes dessa idade. Dá mais peso às mortes em idades jovens."),
    list("Mortalidade prematura", "Mortes antes de uma certa idade (aqui, antes dos 75 anos), muitas vezes consideradas potencialmente evitáveis."),
    list("Intervalo de confiança", "A margem de incerteza à volta de um valor estimado. Um intervalo de 95% indica uma gama de valores plausíveis; é a zona sombreada nos gráficos.")
  )

  forecast_terms <- list(
    list("Previsão (projecção)", "Estimativa de como uma taxa poderá evoluir no futuro, a partir do padrão dos anos anteriores. Não é uma certeza nem uma meta."),
    list("Horizonte", "Quantos anos para o futuro a previsão vai. Quanto maior, maior a incerteza."),
    list("Janela de ajuste (treino)", "Os anos usados para o modelo aprender o padrão da série."),
    list("Teste / validação", "Anos recentes reservados para avaliar quão bem o modelo prevê, antes de confiar na projecção futura."),
    list("Validação móvel", "Forma de validação que repete a previsão a partir de várias origens e combina os erros. É a mais fiável em séries curtas."),
    list("Divisão única (treino/teste)", "Forma de validação que reserva os últimos anos uma só vez para testar o modelo."),
    list("Ajuste dentro da amostra", "Avaliação usando o ajuste à série completa, sem reservar anos. É menos exigente e serve apenas como referência."),
    list("Retroteste (backtesting)", "Testar a previsão contra anos que realmente já aconteceram."),
    list("Modelo", "Um método matemático que descreve o padrão da série para o projectar (por exemplo ARIMA, ETS, Holt, Naive). Na Previsão Guiada, a aplicação escolhe um por si."),
    list("Transformação log", "Um passo opcional que estabiliza séries positivas e evita previsões negativas; a previsão é feita na escala transformada e depois reconvertida."),
    list("Métricas de erro (RMSE, MAE, MAPE, MASE)", "Números que medem quão longe as previsões ficam dos valores reais; servem para comparar modelos. Valores mais baixos são melhores."),
    list("Quebra estrutural", "Uma mudança no padrão da série (por exemplo no nível ou na tendência), que pode dever-se a alterações reais, de codificação ou de registo."),
    list("Resíduos e diagnósticos", "Os resíduos são as diferenças entre o observado e o ajustado; os diagnósticos (ACF, PACF, Ljung-Box) ajudam a verificar se o modelo captou bem o padrão.")
  )

  data_terms <- list(
    list("INE", "Instituto Nacional de Estatística, a fonte oficial dos dados de mortalidade e população."),
    list("Indicador", "Um conjunto de dados específico do INE (por exemplo, óbitos por causa), identificado por um código."),
    list("Ficheiros RDS", "Dados já preparados e guardados no repositório, que a aplicação lê rapidamente sem consultar o INE em directo."),
    list("Fonte de dados", "A escolha entre ler os Ficheiros RDS (rápido) ou consultar o INE em directo (mais lento, para dados não incluídos).")
  )

  tabPanel(
    "Glossário",
    fluidRow(
      column(
        width = 10, offset = 1,
        h2("Glossário"),
        p("Explicações simples dos termos usados na aplicação. Não é preciso conhecê-los todos para começar."),
        glossary_section("Conceitos de mortalidade", mortality_terms),
        glossary_section("Conceitos de previsão", forecast_terms),
        glossary_section("Dados", data_terms)
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
        population_help(),
        radioButtons(
          "rate_type",
          "Taxa:",
          choices = c("Bruta" = "crude", "Padronizada" = "dsr")
        ),
        rate_type_help(),
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
        plotly::plotlyOutput("ratePlot", height = "400px"),
        helpText(
          "A zona sombreada mostra o intervalo de confiança de 95% (a incerteza ",
          "em torno da taxa). Passe o rato sobre um ponto para ver os valores."
        ),
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
        plotly::plotlyOutput("beginnerForecastPlot", height = "400px"),
        helpText(
          "Azul: taxa observada usada no ajuste. Cinzento: restante histórico. ",
          "Linha tracejada e zona sombreada: previsão e a sua incerteza. ",
          "Passe o rato sobre um ponto para ver os valores."
        ),
        br(),
        uiOutput("beginnerForecastWarnings"),
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
        downloadButton("downloadAnnualMetricsPlot", "Descarregar gráfico (PNG)"),
        br(), br(),
        h4("Fontes usadas"),
        tableOutput("annualSourcesTable"),
        downloadButton("downloadAnnualSourcesCSV", "Descarregar fontes (CSV)")
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
    plotly::plotlyOutput("forecastPlot", height = "400px"),
    helpText("Passe o rato sobre um ponto para ver os valores."),
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
            "Divisão única (últimos %)" = "single",
            "Validação móvel (últimos %)" = "rolling"
          ),
          selected = "rolling",
          inline = TRUE
        )
      ),
      column(
        width = 6,
        uiOutput("comparisonHoldoutControl")
      )
    ),
    textOutput("comparisonValidationInfo"),
    br(),
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
    wellPanel(
      p(tags$strong("Para uma previsão simples, use o separador 'Previsão Guiada'.")),
      p(
        "Este separador destina-se a utilizadores técnicos: permite escolher modelos, ",
        "ajustar parâmetros e ver diagnósticos. Não é necessário para obter uma previsão."
      )
    ),
    tabsetPanel(
      advanced_model_spec_tab_ui(),
      advanced_forecast_output_tab_ui(),
      advanced_diagnostics_tab_ui(),
      advanced_backtesting_tab_ui(),
      advanced_breaks_tab_ui()
    )
  )
}
