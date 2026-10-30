# =========================================================
# UI Helpers
# =========================================================

# `years` bounds the slider. Tabs that can only show rates pass `rate_years`,
# the years that have both a numerator and a denominator: with the full range
# the observed tab defaulted to a span whose last year had no denominator and
# failed on load. It has to be the intersection, not either side - population
# now runs a year further than deaths by cause, having previously run a year
# behind. The data-availability tab keeps the full range, as its job is to show
# what exists.
year_range_slider <- function(input_id, label, years = year_of_interest, value = range(years)) {
  sliderInput(
    input_id,
    label,
    min   = min(years),
    max   = max(years),
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
trend_help <- function() {
  helpText(
    "Variação percentual anual (APC): o ritmo a que a taxa sobe ou desce por ano, com o intervalo de confiança de 95%, ",
    "de uma regressão log-linear em que cada ano pesa conforme a sua precisão (método Joinpoint do National Cancer ",
    "Institute). Com mudanças de tendência, a recta pode quebrar em até um ou dois anos, escolhidos pelo critério BIC, ",
    "e a última linha dá a variação média anual (AAPC) de todo o período. Uma tendência é significativa quando o ",
    "intervalo não inclui 0. A linha tracejada no gráfico é a tendência ajustada."
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
annual_metric_help <- function() {
  helpText(
    "Padronizada (directa): compara com a população-padrão europeia; instável ",
    "em concelhos pequenos. SMR (indirecta): compara os óbitos observados com ",
    "os esperados se o local tivesse as taxas da referência — 100 = igual à ",
    "referência. É a métrica indicada para concelhos com poucos óbitos. ",
    "Mortalidade infantil: disponível como taxa (por 1.000 nados-vivos) e como ",
    "contagem de óbitos. Num concelho com poucos nascimentos a taxa oscila ",
    "centenas de pontos com um único óbito; nesses casos o valor é assinalado ",
    "com * e a contagem é mais informativa."
  )
}
smr_reference_help <- function() {
  helpText(
    "As taxas por idade desta referência são aplicadas à estrutura etária de ",
    "cada local para calcular os óbitos esperados."
  )
}
pooling_help <- function() {
  helpText(
    "Soma óbitos e população de vários anos (denominador em pessoas-ano) para ",
    "estabilizar concelhos pequenos e causas raras. Reduz a largura dos ",
    "intervalos, mas suaviza variações anuais reais."
  )
}
# What a regional figure is made of: the vintage groups the municipalities, the
# region source decides whether deaths come from INE's regional rows or from the
# municipal sum. Both are chosen in the header; this note says what they mean.
region_aggregation_note <- function() {
  helpText(
    "As regiões seguem a definição escolhida no topo da página (NUTS 2013 ou ",
    "2024), aplicada a todos os anos. Os óbitos de cada região vêm, por ",
    "predefinição, das linhas regionais do INE sempre que existem, porque a ",
    "soma dos municípios perde parte da repartição por idade dos óbitos por ",
    "causa. A população é sempre a soma dos municípios."
  )
}

# When the data in use was last imported, on every page: two exports of the same
# analysis with different figures can then be told apart by this date.
data_version_badge <- function() {
  date <- tryCatch(
    latest_import_date(read_import_log(app_data_root(get_snapshot_dir()))),
    error = function(e) NA_character_
  )
  if (is.na(date)) return(NULL)
  tags$span(
    style = "font-size:0.85em; opacity:0.75; margin-left:auto; white-space:nowrap;",
    title = "Data da importação mais recente dos dados em uso. O histórico está no separador Disponibilidade de Dados.",
    paste0("Dados importados até ", date)
  )
}

# A single app-wide control, in the navbar header rather than in each tab: the
# vintage is a definition, not a per-analysis parameter, and six region names
# mean different things under the two vintages, so it must be on screen
# wherever a regional figure is read. The other three sit beside it because
# they too change what every tab's figures mean.
#
# One compact line: short labels, native selects as wide as their text (the
# selectize boxes cut it off), and each explanation in the control's tooltip.
header_select <- function(input_id, label, choices, selected, tooltip) {
  tags$span(
    class = "app-header-item",
    title = tooltip,
    tags$label(label, `for` = input_id),
    selectInput(input_id, label = NULL, choices = choices, selected = selected, selectize = FALSE, width = NULL)
  )
}

nuts_vintage_control <- function() {
  tags$div(
    class = "nuts-vintage-header app-header",
    tags$style(HTML(paste(
      ".app-header { display:flex; align-items:center; gap:0.35rem 1.1rem; flex-wrap:wrap;",
      "  padding:0.3rem 0.9rem; margin:0 0 0.5rem 0; border-bottom:1px solid rgba(128,128,128,0.25); font-size:0.9em; }",
      ".app-header-item { display:inline-flex; align-items:center; gap:0.4rem; }",
      ".app-header-item label { margin:0; font-weight:600; white-space:nowrap; }",
      ".app-header-item .form-group { margin:0; display:inline-block; width:auto; }",
      ".app-header-item select.form-control { width:auto; height:28px; padding:2px 6px; font-size:0.95em; }",
      ".navbar-nav > li > a { padding-left:9px; padding-right:9px; }",
      ".download-locked { pointer-events:none; opacity:0.45; }",
      ".tab-notes { margin-bottom:0.8rem; padding-top:0.6rem; padding-bottom:0.6rem; }",
      ".tab-notes li { white-space:pre-line; }",
      ".tab-notes ul { margin:0.3rem 0 0 0; padding-left:1.2rem; }",
      sep = "\n"
    ))),
    header_select(
      "nuts_vintage", "Regiões:", nuts_vintage_choices, default_nuts_vintage(),
      "Definição das regiões. Agrupa os mesmos 308 municípios de outra forma; não altera os dados lidos."
    ),
    header_select(
      "region_source", "Óbitos das regiões:", region_source_choices, default_region_source(),
      paste(
        "Como se obtêm os óbitos de uma região: das linhas regionais que o INE publica (com a repartição por idade completa)",
        "ou da soma dos seus municípios."
      )
    ),
    # The denominator of every rate per inhabitant, in every tab. Life
    # expectancy and the fertility index keep mid-year person-years.
    header_select(
      "population_basis", "Denominador:", population_basis_choices, default_population_basis(),
      paste(
        "População no denominador. As estimativas de população do INE referem-se a 31 de dezembro.",
        "O INE divide os acontecimentos de um ano pela média das estimativas no fim do ano anterior e no fim desse ano.",
        "A população a 31 de dezembro sozinha dá taxas mais baixas quando a população cresce e mais altas quando diminui (nas regiões, em 2019-2025, até cerca de 4%).",
        "Não se aplica à esperança de vida nem ao índice de fecundidade."
      )
    ),
    # Whether the app's corrections of INE's data apply (DADOS_INE_LIMITACOES.md).
    header_select(
      "data_mode", "Dados:", data_mode_choices, default_data_mode(),
      paste(
        "Corrigidos: a aplicação corrige erros de publicação do INE, distribui pelas idades os óbitos sem idade,",
        "reconstrói totais confidenciais, recupera emprego e ganhos ocultados e reconstrói o último ano sem causas de morte.",
        "Originais do INE: os dados tal como o INE os publica, sem nada disto.",
        "A lista completa e as razões estão em DADOS_INE_LIMITACOES.md."
      )
    ),
    data_version_badge(),
    # Downloads stay locked until their analysis has data (lock-downloads), and
    # the map can open a series in another tab (click-button).
    tags$script(HTML(paste(
      "Shiny.addCustomMessageHandler('lock-downloads', function(m) {",
      "  (m.ids || []).forEach(function(id) {",
      "    var el = document.getElementById(id); if (!el) return;",
      "    el.classList.toggle('download-locked', !!m.locked);",
      "    el.setAttribute('aria-disabled', m.locked ? 'true' : 'false');",
      "    el.title = m.locked ? 'Carregue primeiro os dados desta análise.' : '';",
      "  });",
      "});",
      "Shiny.addCustomMessageHandler('click-button', function(id) {",
      "  setTimeout(function() { var el = document.getElementById(id); if (el) el.click(); }, 400);",
      "});",
      sep = "\n"
    )))
  )
}
bias_adjust_help <- function() {
  helpText(
    "Com transformação log, a retransformação directa devolve a mediana, não a ",
    "média. Esta correcção devolve o valor esperado; a diferença cresce com o ",
    "horizonte. Os intervalos não mudam, porque são quantis."
  )
}
transform_help <- function() {
  helpText(
    "A transformação log estabiliza séries positivas e evita previsões negativas; ",
    "'Sem transformação' modela a taxa directamente."
  )
}
model_policy_help <- function() {
  helpText(
    "Para comparar locais ou fixar metas, use o mesmo método em todos: com «o melhor de cada local», um local pode ",
    "ficar com uma tendência e o vizinho com uma linha plana só porque o acaso favoreceu modelos diferentes. A ",
    "média de 6 modelos é igual em todos os locais e foi a que menos errou a 4 anos em testes com dados passados. ",
    "A combinação usa os três modelos que melhor previram os anos de teste de cada local; «o mesmo modelo» usa o ",
    "que teve, em média, o menor erro relativo em todos os locais."
  )
}
infant_method_help <- function() {
  helpText(
    "Série anual: os modelos são ajustados às taxas de cada ano, com a variação real de ano para ano, e cada ",
    "período projectado é depois a média, pesada pelos nascimentos, dos anos que contém (os anos futuros com os ",
    "nascimentos do último ano observado). O intervalo trata os erros dos anos como totalmente correlacionados, ",
    "pelo que é conservador. ",
    "Médias móveis: os modelos são ajustados directamente à série de taxas de 3 ou 5 anos. Cada valor partilha ",
    "quase todos os anos com o anterior, pelo que a série é suave por construção: os modelos tendem a ver ",
    "tendências onde há acaso e dão intervalos demasiado estreitos. Fica disponível para comparar com valores ",
    "calculados dessa forma. Em ambos, a escolha do modelo usa os erros da série ajustada."
  )
}
beginner_method_help <- function() {
  helpText(
    "Por omissão, a projecção é a média de 6 modelos (ARIMA, ETS, Theta, TBATS, Holt e Holt amortecido), com um ",
    "intervalo que junta a incerteza de cada um e o desacordo entre eles. Em testes com dados passados, a 4 anos, ",
    "errou menos do que escolher o modelo com menor erro e foi mais estável: a escolha muda muitas vezes de modelo ",
    "quando chega mais um ano de dados, e os intervalos de um só modelo são em geral estreitos demais."
  )
}
beginner_validation_help <- function() {
  helpText(
    "A validação testa cada modelo, e as duas médias, a prever anos recentes reservados para avaliação. ",
    "Mostra o erro de cada um e, na opção «modelo com menor erro», decide a projecção apresentada. ",
    "'Validação móvel' é a opção mais fiável. Se não souber, mantenha as predefinições."
  )
}

# Shared control panels ----------------------------------------------------

forecast_controls_panel <- function() {
  tagList(
    selectInput("area2", "Local de residência:", choices = local_area, multiple = TRUE, selected = get_default_area_selection()),
    textInput("area_label2", "Nome da selecção (opcional):", placeholder = "Ex.: AML"),
    selectInput("cause2", "Causa de Morte:", choices = cause_choices),
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
    region_aggregation_note(),
    year_range_slider(
      years = rate_years,
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
    conditionalPanel(
      "input.transform2 != 'none'",
      checkboxInput(
        "bias_adjust2",
        "Corrigir enviesamento da retransformação",
        value = TRUE
      ),
      bias_adjust_help()
    ),
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
    radioButtons(
      "beginner_areas_mode",
      "Com vários locais:",
      choices = c("Somá-los numa só selecção" = "pooled", "Projectar cada um em separado, para comparar" = "separate"),
      selected = "pooled"
    ),
    conditionalPanel(
      "input.beginner_areas_mode == 'separate'",
      radioButtons(
        "beginner_model_policy",
        "Modelo em cada local:",
        choices = c(
          "Média de 6 modelos (recomendado)" = "mean6",
          "Combinação dos 3 melhores modelos de cada local" = "combo",
          "O mesmo modelo para todos" = "common",
          "O melhor modelo de cada local" = "best"
        ),
        selected = "mean6"
      ),
      model_policy_help()
    ),
    radioButtons(
      "beginner_indicator",
      "Indicador:",
      choices = c("Taxa de mortalidade" = "mortality", "Mortalidade infantil (por 1.000 nados-vivos)" = "infant"),
      selected = "mortality"
    ),
    conditionalPanel(
      "input.beginner_indicator == 'mortality'",
      selectInput("beginner_cause", "Causa de Morte:", choices = cause_choices),
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
      region_aggregation_note()
    ),
    conditionalPanel(
      "input.beginner_indicator == 'infant'",
      selectInput(
        "beginner_infant_window",
        "Período de cada valor:",
        choices = c("Anual" = "1", "Triénio" = "3", "Quinquénio" = "5"),
        selected = "5"
      ),
      conditionalPanel(
        "input.beginner_infant_window != '1'",
        radioButtons(
          "beginner_infant_method",
          "Como projectar:",
          choices = c(
            "Projectar a série anual e agregar no fim (recomendado)" = "annual",
            "Projectar a série de médias móveis" = "rolling"
          ),
          selected = "annual"
        ),
        infant_method_help()
      )
    ),
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
      years = rate_years,
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
        "beginner_method",
        "Projecção apresentada:",
        choices = c(
          "Média de 6 modelos (recomendada)" = "mean6",
          "O modelo com menor erro na validação" = "selected"
        ),
        selected = "mean6"
      ),
      beginner_method_help(),
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
    actionButton("cancel_beginner_forecast", "Interromper carregamento"),
    map_link_button("map_from_beginner")
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
            choices = cause_choices,
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
        downloadButton("downloadSnapshotAvailabilityCSV", "Descarregar cobertura (CSV)"),
        br(), br(),
        h4("Histórico dos dados"),
        helpText(
          "O INE revê dados já publicados: anos provisórios passam a definitivos, e ",
          "séries de população são re-estimadas. Cada importação fica registada com a ",
          "sua data, e a versão anterior de qualquer ficheiro que mude é guardada. Se ",
          "uma análise der hoje valores diferentes dos de uma análise anterior, veja ",
          "aqui se os dados foram revistos entretanto."
        ),
        h5("Conjuntos de dados"),
        tableOutput("dataImportDatasets"),
        h5("Importações"),
        tableOutput("dataImportRuns"),
        selectInput("data_import_run", "Ver os valores revistos numa importação:", choices = NULL, width = "100%"),
        tableOutput("dataImportChanges"),
        downloadButton("downloadImportLogCSV", "Descarregar registo completo (CSV)"),
        br(), br(),
        helpText(
          "Para repetir uma análise com os dados tal como estavam numa data: ",
          tags$code("Rscript tools/data_as_of.R date=AAAA-MM-DD"),
          " e depois abrir a aplicação com ",
          tags$code("MORTALITY_SNAPSHOT_DIR=<pasta indicada>/snapshots"),
          ". Ver o manual, secção sobre o histórico dos dados."
        )
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
          tags$li(tags$strong("Métricas Anuais"), " — comparar um local com Portugal e a sua região, ou com outras áreas à escolha, num único ano."),
          tags$li(tags$strong("Mortalidade Evitável"), " — os óbitos antes dos 75 anos que se poderiam ter evitado, por prevenção ou por tratamento."),
          tags$li(tags$strong("Mapas"), " — qualquer causa de morte em todos os municípios, NUTS III, ULS ou regiões, num triénio. Passe o rato para ver os valores; clique numa área para abrir a sua série. Os separadores de mortalidade têm um botão ", tags$em("Ver esta causa no mapa"), "."),
          tags$li(tags$strong("Indicadores de Planeamento"), " — os indicadores dos Planos Locais de Saúde (população, natalidade, mortalidade, esperança de vida, condições sociais) de um local face aos seus comparadores, com mapas, funil e descarga em Excel."),
          tags$li(tags$strong("Disponibilidade de Dados"), " — verificar que anos, locais e causas existem nos ficheiros antes de carregar."),
          tags$li(tags$strong("Glossário"), " — os termos usados, em linguagem simples.")
        ),

        h3("Os controlos do topo da página"),
        p("Aplicam-se a todos os separadores. As predefinições servem para quase todos os usos; passe o rato sobre cada um para ver o que faz."),
        tags$ul(
          tags$li(tags$strong("Regiões:"), " a definição das regiões, NUTS 2024 (a actual) ou NUTS 2013. Agrupa os mesmos municípios de outra forma."),
          tags$li(tags$strong("Óbitos das regiões:"), " das linhas regionais que o INE publica (recomendado) ou da soma dos municípios."),
          tags$li(tags$strong("Denominador:"), " a população por que se dividem os óbitos, como o INE (média de dois fins de ano) ou a de 31 de dezembro."),
          tags$li(tags$strong("Dados:"), " os dados corrigidos pela aplicação (recomendado) ou tal como o INE os publica. As correcções estão explicadas em DADOS_INE_LIMITACOES.md.")
        ),
        p("Os avisos sobre os dados de uma análise aparecem no topo do seu separador. As descargas ficam activas quando a análise tem dados."),

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
            "separador a separador, consulte o Manual do Utilizador (USER_MANUAL.md); para as escolhas ",
            "metodológicas e as suas razões, em linguagem simples, ESCOLHAS_METODOLOGICAS.md; para os problemas ",
            "encontrados nos dados do INE e o que a aplicação faz com cada um, DADOS_INE_LIMITACOES.md."
          ))
        )
      )
    )
  )
}

# Opens the Maps tab with the cause and sex chosen in this tab.
map_link_button <- function(input_id) {
  tagList(br(), br(), actionButton(input_id, "Ver esta causa no mapa", icon = icon("map")))
}

FORECAST_TABLE_RECENT_YEARS <- 5L
forecast_table_note <- function() {
  helpText(paste0("Os últimos ", FORECAST_TABLE_RECENT_YEARS, " anos observados e a previsão. O CSV tem a série completa."))
}

maps_tab_ui <- function() {
  tabPanel(
    "Mapas",
    sidebarLayout(
      sidebarPanel(
        selectInput("map_cause", "Causa de Morte:", choices = cause_choices, selected = "Todas as causas de morte"),
        radioButtons("map_sex", "Sexo:", choices = c("Ambos" = "HM", "Homens" = "H", "Mulheres" = "M"), selected = "HM", inline = TRUE),
        selectInput("map_metric", "Métrica:", choices = CAUSE_MAP_METRICS, selected = "smr"),
        helpText(
          "SMR: óbitos observados a dividir pelos esperados com as taxas por idade de Portugal (100 = Portugal). ",
          "É a métrica indicada para comparar municípios. A taxa padronizada usa a População Padrão Europeia de 2013."
        ),
        radioButtons(
          "map_view", "Mostrar:",
          choices = c("Um triénio" = "period", "Tendência ao longo dos triénios" = "trend", "Variação entre dois triénios" = "change",
                      "Projecção (NUTS III, ULS, NUTS II)" = "projection"),
          selected = "period"
        ),
        conditionalPanel("input.map_view != 'projection'", selectInput("map_end_year", "Triénio:", choices = NULL)),
        conditionalPanel(
          "input.map_view == 'trend'",
          selectInput("map_trend_n", "Número de triénios:", choices = c("4 (12 anos)" = "4", "5 (15 anos)" = "5", "6 (18 anos)" = "6"), selected = "5"),
          helpText("Triénios seguidos, sem anos em comum, a terminar no triénio escolhido.")
        ),
        conditionalPanel(
          "input.map_view == 'change'",
          selectInput("map_base_year", "Comparar com o triénio:", choices = NULL)
        ),
        conditionalPanel(
          "input.map_view == 'projection'",
          selectInput("map_target_year", "Projectar até:", choices = as.character(2026:2035), selected = "2030"),
          helpText("Taxa padronizada anual (todas as idades, ou menos de 75 anos com essa métrica), com o método da Previsão Guiada. ",
                   "Da primeira vez demora cerca de um minuto; depois fica guardada.")
        ),
        radioButtons("map_level", "Unidades:", choices = PLANNING_MAP_LEVELS, selected = "municipio"),
        radioButtons(
          "map_scale", "Cores:",
          choices = c("Escala contínua" = "index", "Significância (\u25b2 = \u25bc)" = "significance", "Quintis" = "quintile"),
          selected = "index"
        ),
        helpText(
          "Escala contínua: face a Portugal (= 100) num triénio, em % por ano na tendência, face ao triénio anterior (= 100) na variação. ",
          "Significância: só três cores, conforme o intervalo de 95% exclui ou não Portugal, a ausência de tendência ou de mudança."
        ),
        conditionalPanel(
          "input.map_metric == 'smr' && input.map_view == 'period'",
          checkboxInput("map_smoothed", "Suavizar a SMR das áreas pequenas (Bayes empírico)", value = FALSE)
        ),
        actionButton("go_map", "Desenhar mapa")
      ),
      mainPanel(
        uiOutput("mapsNotes"),
        plotly::plotlyOutput("causeMap", height = "720px"),
        uiOutput("causeMapNote"),
        uiOutput("causeMapClicked"),
        br(),
        downloadButton("downloadCauseMapCSV", "Descarregar valores (CSV)"),
        downloadButton("downloadCauseMapPNG", "Descarregar mapa (PNG)")
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
    list("Padronização directa", "Aplica as taxas por idade do local a uma população padrão externa. Dá a Taxa padronizada. Precisa de uma taxa estimável em cada idade, por isso é instável em locais pequenos."),
    list("Padronização indirecta", "Aplica as taxas por idade da referência à estrutura etária do local, para calcular quantos óbitos seriam de esperar. Dá o SMR. É estável em locais pequenos porque só precisa do total de óbitos observados."),
    list("SMR", "Óbitos observados a dividir pelos esperados, vezes 100. A referência vale 100; 120 são 20% mais óbitos do que o esperado. Compare cada SMR com 100, nunca com outro SMR."),
    list("Óbitos esperados", "Os óbitos que o local teria tido com as taxas por idade da referência e a sua própria estrutura etária."),
    list("Agregação plurianual", "Calcular a métrica sobre 3 ou 5 anos em vez de 1, para estabilizar locais pequenos e causas raras."),
    list("Pessoas-ano", "O denominador de uma taxa agregada. Cinco anos de um concelho com 10.000 habitantes são 50.000 pessoas-ano, o que mantém a taxa por ano e comparável com um ano isolado."),
    list("Mortalidade infantil", "Óbitos antes do primeiro ano de vida por 1.000 nados-vivos."),
    list("Nados-vivos", "Nascimentos com vida. São o denominador da mortalidade infantil, em vez da população, porque nenhum indicador de população tem uma banda 'menos de 1 ano' e porque correspondem melhor ao grupo em risco."),
    list("Asterisco (*)", "Marca uma taxa de mortalidade infantil calculada sobre menos de 1.000 nados-vivos, em que um único óbito desloca o valor em mais de uma unidade por 1.000. O valor é exacto; a marca avisa que não é comparável."),
    list("Índice de envelhecimento", "Pessoas com 65 e mais anos por cada 100 com 0-14 anos. Acima de 100 há mais idosos do que jovens."),
    list("Índices de dependência", "Jovens (0-14) ou idosos (65+) por cada 100 pessoas em idade activa (15-64 anos). Medem o peso das idades dependentes sobre as activas."),
    list("Taxa bruta de natalidade / mortalidade", "Nados-vivos ou óbitos por 1.000 habitantes. No separador Indicadores de Planeamento usa-se ‰, como nos Planos Locais de Saúde; nos separadores de mortalidade as taxas são por 100.000."),
    list("Intervalo de confiança", "A margem de incerteza à volta de um valor estimado. Um intervalo de 95% indica uma gama de valores plausíveis; é a zona sombreada nos gráficos. Um intervalo largo não é um defeito: é o que há a dizer quando os acontecimentos são poucos.")
  )

  geography_terms <- list(
    list("ULS (Unidade Local de Saúde)", "A unidade de organização do SNS a que corresponde a população de um conjunto de municípios. Na aplicação é somada a partir dos seus municípios. Seis ULS partilham Lisboa, Loures ou Porto ao nível da freguesia e aparecem também em dois agrupamentos exactos; a freguesia da Lomba, de Gondomar, pertence à ULS Vila Nova de Gaia/Espinho."),
    list("ARS", "As cinco regiões de saúde (Norte, Centro, Lisboa e Vale do Tejo, Alentejo, Algarve) que agrupam as ULS. Não coincidem com as regiões NUTS com o mesmo nome."),
    list("NUTS", "A nomenclatura estatística das regiões. A aplicação usa dois níveis: NUTS I (Continente, Açores, Madeira) e NUTS II (as regiões)."),
    list("Definição das regiões (NUTS 2013 / NUTS 2024)", "As duas versões da nomenclatura que a aplicação oferece, no controlo do topo da página. Agrupam os mesmos 308 municípios de formas diferentes: seis nomes existem nas duas e significam coisas diferentes em cada uma."),
    list("Área Metropolitana de Lisboa", "A região de Lisboa em NUTS 2013. Em NUTS 2024 está dividida em Grande Lisboa e Península de Setúbal."),
    list("Oeste e Vale do Tejo", "Região criada em NUTS 2024 com o Oeste, o Médio Tejo e a Lezíria do Tejo. Apesar do nome, não inclui Lisboa."),
    list("Óbitos das regiões", "Por predefinição, os óbitos de uma região vêm das linhas regionais que o INE publica, que têm a repartição por idade completa. A alternativa é a soma dos seus municípios, a que o INE publica com idades em falta em alguns anos e causas. A população é sempre a soma dos municípios, com a mesma lista aplicada a todos os anos, o que mantém a série contínua apesar da revisão NUTS de 2024."),
    list("Triénio", "Três anos seguidos somados (por exemplo 2022-2024), como nos mapas e nos Indicadores de Planeamento: os óbitos e as pessoas-ano dos três anos, para estabilizar áreas pequenas."),
    list("Revisão da população (2021)", "O INE publica duas estimativas de população que não coincidem; a aplicação usa a revista a partir de 2021. Todas as taxas desde esse ano mudaram, e há um degrau de cerca de 1,7% em 2020/2021.")
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
    list("Fonte de dados", "A escolha entre ler os Ficheiros RDS (rápido) ou consultar o INE em directo (mais lento, para dados não incluídos)."),
    list("Dados corrigidos / originais do INE", "O controlo Dados, no topo da página. Corrigidos: a aplicação corrige erros de publicação do INE, distribui pelas idades os óbitos sem idade e reconstrói valores confidenciais e o último ano. Originais do INE: tal como publicados. A lista completa está em DADOS_INE_LIMITACOES.md."),
    list("Suavização (Bayes empírico)", "Nos mapas e no funil, aproxima o valor de cada área pequena da média de todas, tanto mais quanto menos óbitos tem. Reduz as cores extremas que só se devem ao acaso; o valor observado continua no rato e no CSV.")
  )

  tabPanel(
    "Glossário",
    fluidRow(
      column(
        width = 10, offset = 1,
        h2("Glossário"),
        p("Explicações simples dos termos usados na aplicação. Não é preciso conhecê-los todos para começar."),
        glossary_section("Conceitos de mortalidade", mortality_terms),
        glossary_section("Geografia", geography_terms),
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
        selectInput("cause", "Causa de Morte:", choices = cause_choices),
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
        region_aggregation_note(),
        year_range_slider(
          years = rate_years,
          "years_import",
          "Anos a importar:"
        ),
        data_source_input("data_source"),
        actionButton("go_rates", "Carregar dados"),
        br(), br(),
        actionButton("cancel_rates", "Interromper carregamento"),
        map_link_button("map_from_observed")
      ),
      mainPanel(
        uiOutput("observedNotes"),
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
        h4("Tendência"),
        fluidRow(
          column(6, selectInput(
            "observed_trend_span", "Período da tendência:",
            choices = c("Últimos 10 anos" = "10", "Últimos 5 anos" = "5", "Últimos 15 anos" = "15", "Todos os anos carregados" = "all"),
            selected = "10"
          )),
          column(6, selectInput(
            "observed_trend_joinpoints", "Mudanças de tendência:",
            choices = c("Nenhuma (uma só recta)" = "0", "Até 1" = "1", "Até 2" = "2"),
            selected = "0"
          ))
        ),
        tableOutput("rateTrendTable"),
        uiOutput("rateTrendNote"),
        trend_help(),
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
        uiOutput("beginnerNotes"),
        conditionalPanel(
          "input.beginner_areas_mode == 'separate'",
          h4("Previsão por local"),
          plotly::plotlyOutput("beginnerMultiPlot", height = "560px"),
          uiOutput("beginnerMultiNote"),
          div(style = "overflow-x:auto;", tableOutput("beginnerMultiTable")),
          downloadButton("downloadBeginnerMultiCSV", "Descarregar tabela (CSV)")
        ),
        conditionalPanel(
          "input.beginner_areas_mode != 'separate'",
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
        forecast_table_note(),
        tableOutput("beginnerForecastTable"),
        downloadButton("downloadBeginnerForecastCSV", "Descarregar tabela (CSV)")
        )
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
          choices = cause_choices,
          selected = if ("Todas as causas de morte" %in% diseases) "Todas as causas de morte" else utils::head(diseases, 1),
          multiple = TRUE
        ),
        selectInput("annual_sex", "Sexo:", choices = sex_levels, selected = "HM"),
        selectInput(
          "annual_area",
          "Local:",
          choices = local_area,
          multiple = TRUE,
          selected = character(0)
        ),
        selectizeInput(
          "annual_comparators",
          "Comparar com:",
          choices = annual_comparator_choices(),
          multiple = TRUE,
          selected = ANNUAL_DEFAULT_COMPARATORS,
          options = list(maxItems = 4)
        ),
        helpText(
          "Por omissão, Portugal e a região NUTS II do local (o Norte para Matosinhos, o Alentejo para Beja). ",
          "Pode juntar ou trocar por qualquer região, ULS ou município; cada um é uma coluna."
        ),
        textInput("annual_area_label", "Nome da selecção (opcional):", placeholder = "Ex.: ACES / ULS"),
        selectInput(
          "annual_metric",
          "Métrica:",
          choices = annual_metric_choices,
          selected = "dsr"
        ),
        annual_metric_help(),
        region_aggregation_note(),
        conditionalPanel(
          "input.annual_metric == 'smr' || input.annual_metric == 'isr'",
          selectInput(
            "annual_smr_reference",
            "Referência da padronização indirecta:",
            choices = smr_reference_choices,
            selected = "Portugal"
          ),
          smr_reference_help()
        ),
        selectInput(
          "annual_pooling",
          "Agregação plurianual:",
          choices = POOLING_WINDOW_CHOICES,
          selected = "1"
        ),
        pooling_help(),
        data_source_input("annual_data_source"),
        actionButton("go_annual_metrics", "Carregar métricas"),
        br(), br(),
        actionButton("cancel_annual_metrics", "Interromper carregamento"),
        map_link_button("map_from_annual")
      ),
      mainPanel(
        uiOutput("annualNotes"),
        tableOutput("annualMetricsTable"),
        uiOutput("annualMetricsFootnote"),
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

avoidable_tab_ui <- function() {
  tabPanel(
    "Mortalidade Evitável",
    sidebarLayout(
      sidebarPanel(
        selectInput(
          "avoidable_year",
          "Ano:",
          # By cause: not the latest year, rebuilt for all causes only.
          choices = sort(setdiff(rate_years, calendar_death_years)),
          selected = max(setdiff(rate_years, calendar_death_years))
        ),
        selectInput(
          "avoidable_area",
          "Local de residência:",
          choices = local_area,
          multiple = TRUE,
          selected = "Portugal"
        ),
        textInput("avoidable_area_label", "Nome da selecção (opcional):", placeholder = "Ex.: ULS"),
        selectInput("avoidable_sex", "Sexo:", choices = sex_levels, selected = "HM"),
        selectInput(
          "avoidable_pooling",
          "Agregação plurianual:",
          choices = POOLING_WINDOW_CHOICES,
          selected = "1"
        ),
        pooling_help(),
        region_aggregation_note(),
        data_source_input("avoidable_data_source"),
        actionButton("go_avoidable", "Carregar"),
        br(), br(),
        actionButton("cancel_avoidable", "Interromper carregamento")
      ),
      mainPanel(
        uiOutput("avoidableNotes"),
        helpText(avoidable_scope_note()),
        h4("Repartição dos óbitos com menos de 75 anos"),
        tableOutput("avoidableTable"),
        uiOutput("avoidableNote"),
        downloadButton("downloadAvoidableCSV", "Descarregar tabela (CSV)"),
        br(), br(),
        plotOutput("avoidablePlot", height = "360px"),
        downloadButton("downloadAvoidablePlot", "Descarregar gráfico (PNG)"),
        br(), br(),
        h4("Que causas entram em cada grupo"),
        tableOutput("avoidableCausesTable"),
        downloadButton("downloadAvoidableCausesCSV", "Descarregar lista de causas (CSV)")
      )
    )
  )
}

# Planning indicators, built around one location: each indicator as a chart and
# a table, set against the areas that contain the location where the indicator
# is comparable across sizes, plus a summary of every indicator, the ranking of
# ULS, the pyramid and proportional mortality. Excel downloads close the tab.
planning_tab_ui <- function() {
  years <- tryCatch(
    sort(unique(unlist(lapply(PLANNING_INDICATORS$id, planning_indicator_years)))),
    error = function(e) integer(0)
  )
  if (length(years) == 0) years <- c(1991L, 2025L)
  last <- max(years)

  tabPanel(
    "Indicadores de Planeamento",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        selectInput("planning_area", "Local:", choices = local_area, selected = "Portugal"),
        selectInput("planning_indicator", "Indicador:", choices = planning_indicator_choices(), selected = "ageing_index"),
        uiOutput("planningComparatorsUI"),
        radioButtons(
          "planning_split_mode", "ULS que partilham um município:",
          choices = PLANNING_SPLIT_MODES, selected = PLANNING_DEFAULT_SPLIT_MODE
        ),
        helpText(
          "Lisboa, Loures e Porto estão repartidos por ULS ao nível da freguesia. ",
          "«Município inteiro»: cada ULS leva o município todo — nada é estimado, mas ",
          "as seis ULS sobrepõem-se e a sua soma conta esses municípios mais do que uma vez. ",
          "«Ponderação por freguesias»: cada município é repartido pela população das suas ",
          "freguesias nos Censos de 2021, por grupo etário — as partes somam o país, ",
          "assumindo que as quotas se mantêm. Os agrupamentos exactos de ULS não dependem da escolha. ",
          "Gondomar é a excepção: a freguesia da Lomba (0,8%) é da ULS Vila Nova de Gaia/Espinho, mas ",
          "com «Município inteiro» Gondomar conta todo na ULS Santo António; ver as notas do indicador."
        ),
        radioButtons(
          "planning_portugal", "Portugal:",
          choices = c("Total publicado pelo INE" = "published", "Soma dos 308 municípios" = "municipal"),
          selected = "published"
        ),
        helpText(
          "O total do INE inclui os acontecimentos de residência desconhecida (0,3-0,9% dos ",
          "óbitos), que nenhuma região, ULS ou município contém. A soma dos municípios ",
          "compara igual com igual. Vale para comparadores, significância e funil."
        ),
        selectInput(
          "planning_education_age", "Escolaridade [I24], população:",
          choices = planning_education_age_choices(), selected = 0L
        ),
        sliderInput(
          "planning_years", "Anos:",
          min = min(years), max = last, value = c(max(min(years), last - 19L), last), step = 1, sep = ""
        ),
        helpText(
          "Cada local é a soma dos seus municípios, e cada indicador a razão dessas ",
          "somas. Portugal e o Continente usam as linhas publicadas pelo INE."
        )
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          id = "planning_view",
          tabPanel(
            "Indicador",
            br(),
            uiOutput("planningIndicatorHeader"),
            plotly::plotlyOutput("planningIndicatorPlot", height = "440px"),
            uiOutput("planningIndicatorNotes"),
            h5("Valores"),
            div(style = "overflow-x:auto;", tableOutput("planningIndicatorTable")),
            h5("Tendência nos anos escolhidos"),
            div(style = "overflow-x:auto;", tableOutput("planningTrendTable")),
            uiOutput("planningTrendNote")
          ),
          tabPanel(
            "Resumo do local",
            br(),
            helpText(
              "Todos os indicadores no último período disponível para o local, até ao ",
              "último ano escolhido. As contagens absolutas não têm comparadores."
            ),
            div(style = "overflow-x:auto;", tableOutput("planningSummaryTable")),
            uiOutput("planningSummaryNotes")
          ),
          tabPanel(
            "Comparação entre ULS",
            br(),
            uiOutput("planningRankingNote"),
            plotly::plotlyOutput("planningRankingPlot", height = "820px")
          ),
          tabPanel(
            "Mapa",
            br(),
            fluidRow(
              column(6, selectizeInput(
                "planning_map_indicators", "Indicadores (até 4):",
                choices = planning_map_indicator_choices(), multiple = TRUE,
                selected = c("smr_all", "dsr_premature", "ageing_index", "purchasing_power"),
                options = list(maxItems = PLANNING_MAP_MAX_INDICATORS)
              )),
              column(3, radioButtons("planning_map_level", "Unidades:", choices = PLANNING_MAP_LEVELS, selected = "municipio")),
              column(3, radioButtons(
                "planning_map_scale", "Cores:",
                choices = c("Face a Portugal (= 100): vermelho pior, azul melhor" = "index",
                            "Significância face a Portugal (\u25b2 = \u25bc)" = "significance",
                            "Quintis de cada indicador" = "quintile"),
                selected = "index"
              ))
            ),
            checkboxInput("planning_map_smoothed", "Suavizar as taxas de áreas pequenas (Bayes empírico)", value = FALSE),
            uiOutput("planningMapNote"),
            uiOutput("planningMapWidgets"),
            uiOutput("planningMapClicked"),
            downloadButton("downloadPlanningMapCSV", "Descarregar valores (CSV)"),
            downloadButton("downloadPlanningMapPNG", "Descarregar mapas (PNG)")
          ),
          tabPanel(
            "Funil",
            br(),
            radioButtons(
              "planning_funnel_units", "Unidades:",
              choices = c("ULS" = "ULS", "Municípios" = "Município", "NUTS III" = "NUTS III"),
              selected = "ULS", inline = TRUE
            ),
            radioButtons(
              "planning_funnel_values", "Valores:",
              choices = c("Observados (funil)" = "observed", "Suavizados (Bayes empírico)" = "smoothed"),
              selected = "observed", inline = TRUE
            ),
            uiOutput("planningFunnelNote"),
            plotly::plotlyOutput("planningFunnelPlot", height = "560px")
          ),
          tabPanel(
            "Pirâmide etária",
            br(),
            helpText(
              "Percentagem da população em cada grupo etário e sexo, no último ano ",
              "escolhido. O contorno escuro é o primeiro comparador seleccionado ",
              "(ou Portugal), para comparar a estrutura de áreas de tamanhos diferentes."
            ),
            plotly::plotlyOutput("planningPyramidPlot", height = "540px")
          ),
          tabPanel(
            "Mortalidade proporcional",
            br(),
            radioButtons(
              "planning_proportional_ages", "Idades:",
              choices = c("Todas as idades [I45]" = "all", "Menos de 75 anos [I46]" = "under75"),
              selected = "all", inline = TRUE
            ),
            uiOutput("planningProportionalNote"),
            plotly::plotlyOutput("planningProportionalPlot", height = "520px"),
            div(style = "overflow-x:auto;", tableOutput("planningProportionalTable"))
          ),
          tabPanel(
            "Mortalidade por causa",
            br(),
            radioButtons(
              "planning_cause_sex", "Sexo:",
              choices = c("Ambos" = "HM", "Homens" = "H", "Mulheres" = "M"), selected = "HM", inline = TRUE
            ),
            uiOutput("planningCauseNote"),
            plotly::plotlyOutput("planningCausePlot", height = "560px"),
            div(style = "overflow-x:auto;", tableOutput("planningCauseTable"))
          ),
          tabPanel(
            "Cuidados de saúde primários",
            br(),
            selectInput("planning_sns_indicator", "Indicador (Portal da Transparência do SNS):", choices = planning_sns_choices(), width = "100%"),
            uiOutput("planningSnsNote"),
            plotly::plotlyOutput("planningSnsPlot", height = "420px"),
            div(style = "overflow-x:auto;", tableOutput("planningSnsTable")),
            h5("Todas as ULS no último período completo"),
            plotly::plotlyOutput("planningSnsRanking", height = "820px"),
            h5("Mapa: ACES até 2023, ULS desde 2024"),
            fluidRow(
              column(4, selectInput("planning_sns_map_year", "Ano:", choices = NULL)),
              column(8, radioButtons(
                "planning_sns_map_scale", "Cores:",
                choices = c("Face ao Continente (= 100)" = "index", "Significância (\u25b2 = \u25bc)" = "significance", "Quintis" = "quintile"),
                selected = "index", inline = TRUE
              ))
            ),
            uiOutput("planningSnsMapNote"),
            plotly::plotlyOutput("planningSnsMap", height = "680px")
          ),
          tabPanel(
            "Mortalidade semanal",
            br(),
            radioButtons("planning_weekly_age", "Idades:", choices = WEEKLY_AGE_GROUPS, selected = "all", inline = TRUE),
            uiOutput("planningWeeklyNote"),
            plotly::plotlyOutput("planningWeeklyPlot", height = "420px"),
            div(style = "overflow-x:auto;", tableOutput("planningWeeklyTable")),
            h5("Excesso acumulado ao longo do ano"),
            plotly::plotlyOutput("planningWeeklyCumulative", height = "380px")
          ),
          tabPanel("Notas", br(), planning_method_notes())
        ),
        hr(),
        wellPanel(
          h4("Descarregar em Excel"),
          fluidRow(
            column(
              6,
              downloadButton("downloadPlanningSelectionXLSX", "Local e comparadores"),
              helpText(
                "O local escolhido e os comparadores seleccionados, com todos os ",
                "indicadores nos anos escolhidos: resumo, uma folha por indicador, ",
                "dados com intervalos de confiança, pirâmide etária e mortalidade proporcional."
              )
            ),
            column(
              6,
              downloadButton("downloadPlanningFullXLSX", "Todas as áreas"),
              helpText(
                "Portugal, NUTS I, II e III, ARS, ULS e os 308 municípios, todos os ",
                "indicadores e anos, uma folha por indicador. Demora cerca de dois ",
                "minutos a gerar da primeira vez; depois é reutilizado."
              )
            )
          ),
          hr(),
          h4("Perfil do local"),
          downloadButton("downloadPlanningProfile", "Perfil do local (Word)"),
          helpText(
            "Documento Word editável com o retrato do local escolhido face aos comparadores: ",
            "indicadores-chave com significância face a Portugal, pirâmide etária, evolução, ",
            "posição entre as ULS, mortalidade proporcional e notas de método."
          ),
          uiOutput("planningDataDate")
        )
      )
    )
  )
}

planning_method_notes <- function() {
  tagList(
    h4("Comparadores"),
    tags$ul(
      tags$li("Para cada local, a aplicação propõe as áreas que o contêm, uma por nível: ULS, ARS, NUTS III, NUTS II, NUTS I e Portugal. Uma área com exactamente os mesmos municípios do local (por exemplo a ULS Matosinhos para o município de Matosinhos) é omitida, porque repetiria os mesmos valores."),
      tags$li("Só têm comparadores os indicadores que não dependem do tamanho da área: taxas, proporções, índices e valores por habitante. As contagens (população, nados-vivos, óbitos, beneficiários, pensionistas) mostram-se apenas para o local."),
      tags$li("Cada nível tem sempre a mesma cor, e o intervalo de confiança do local aparece como faixa. Os pontos vazios têm uma marca (* ou \u2020).")
    ),
    h4("Portugal, significância e funil"),
    tags$ul(
      tags$li("ULS que partilham um município: Lisboa, Loures e Porto estão repartidos por ULS ao nível da freguesia (Decreto-Lei n.º 102/2023). «Município inteiro» dá a cada ULS o município todo, sem estimar nada, mas as seis ULS sobrepõem-se e a sua soma conta esses municípios mais do que uma vez (em 2024, mais 18.186 óbitos do que o Continente). «Ponderação por freguesias» reparte cada município pela população das suas freguesias nos Censos de 2021, por grupo etário: as partes somam o município e o país, assumindo que as quotas se mantêm. Os agrupamentos exactos de ULS não dependem da escolha e não assumem nada."),
      tags$li("A freguesia da Lomba (Gondomar, 0,8% dos residentes) pertence à ULS Vila Nova de Gaia/Espinho, mas só a ponderação por freguesias a separa: com «Município inteiro» Gondomar conta todo na ULS Santo António, porque dar o município inteiro também à de Gaia/Espinho somaria 164 mil residentes a esta ULS. Pela mesma razão, o agrupamento exacto ULS Santo António + São João inclui a Lomba."),
      tags$li("Portugal pode ser o total publicado pelo INE, que inclui os acontecimentos de residência desconhecida (0,3-0,9% dos óbitos), ou a soma dos 308 municípios, que compara igual com igual. A escolha vale para comparadores, significância, classificação das ULS, funil e perfil."),
      tags$li("\u25b2 / \u25bc / = : o intervalo de confiança de 95% fica inteiramente acima, abaixo ou inclui o valor de Portugal no mesmo período (critério do PHE Fingertips). Só para indicadores comparáveis com intervalo; não diz se a diferença é boa ou má. Na classificação das ULS: laranja acima, azul abaixo, cinzento sem diferença."),
      tags$li("Funil: cada unidade contra o tamanho do denominador, com os limites do que o acaso produziria à volta de Portugal (95% e 99,8%), calculados pelos quantis exactos da contagem (Poisson ou binomial, Spiegelhalter 2005). Só para taxas de acontecimentos e proporções de nascimentos."),
      tags$li("Escolaridade com idade mínima: Censos de 2011 e 2021 por grupo etário (o INE não a publica por idade e município em 1991 e 2001). Com toda a população, as crianças contam como sem nível completo."),
      tags$li("* junto ao nome de um indicador: nota de método. Ganho médio e sectores [I27, I12] contam no local de trabalho, sem Administração Pública nem trabalhadores por conta própria; a esperança de vida [I10] reproduz o Eurostat e fica 0,8-0,9 anos acima do INE.")
    ),
    h4("Mortalidade padronizada, prematura e evitável"),
    tags$ul(
      tags$li("SMR: óbitos observados sobre os esperados com as taxas por idade de Portugal (a opção escolhida) no mesmo triénio, vezes 100; intervalo de Poisson exacto."),
      tags$li("Taxas padronizadas: População Padrão Europeia de 2013, por 100.000 habitantes; intervalo de Dobson. Prematura: antes dos 75 anos."),
      tags$li("Evitável por prevenção e por cuidados de saúde: listas Eurostat/OCDE de 2019 adaptadas à lista sucinta do INE, antes dos 75 anos. É um limite inferior: seis causas (cerca de 18% dos óbitos antes dos 75) ficam de fora."),
      tags$li("Anos potenciais de vida perdidos: anos que faltavam até aos 70 em cada óbito, por 100.000 residentes com menos de 70 anos."),
      tags$li("Os óbitos por idade e causa de cada município estão incompletos no INE; são completados até ao total de todas as idades. Os que faltam são colocados onde a população do próprio município os torna prováveis (população por idade às taxas nacionais) e depois ajustados para que a soma dos municípios reproduza também as idades em falta no país. \u2021 quando mais de 2% foram redistribuídos; acima de 25%, a esperança de vida e as taxas padronizadas não são apresentadas. As áreas cujo território o INE publica como linha própria - as regiões NUTS e as ULS que coincidem com uma delas - usam essas linhas por idade em vez da soma dos municípios, pelo que não têm óbitos redistribuídos.")
    ),
    h4("Cuidados de saúde primários e mortalidade semanal"),
    tags$ul(
      tags$li("Cuidados de saúde primários: Portal da Transparência do SNS, por ULS desde Janeiro de 2024; base: utentes inscritos. Os rastreios e o exame dos pés acumulam ao longo do ano, a tensão arterial e a HbA1c ao longo do semestre: só os fins de ciclo são comparáveis. As cinco ULS de Lisboa e do Porto somam-se nos dois agrupamentos exactos."),
      tags$li("Mortalidade semanal: INE, óbitos semanais por NUTS III e idade. Esperados: taxas por idade dos anos de base (desde 2023, depois da COVID-19 e na série de população revista) aplicadas à população do ano, com intervalo de previsão de 95%."),
      tags$li("Os dados que mudam com frequência são actualizados semanalmente por uma tarefa agendada; cada actualização fica no histórico dos dados.")
    ),
    h4("Como são calculados"),
    tags$ul(
      tags$li("Índice de envelhecimento: população com 65 e mais anos por 100 com 0-14 anos."),
      tags$li("Índices de dependência: jovens (0-14) ou idosos (65+) por 100 pessoas com 15-64 anos."),
      tags$li("Taxas brutas de natalidade e mortalidade: nados-vivos ou óbitos por 1.000 habitantes, sobre a população média do ano (a média das estimativas a 31 de Dezembro do ano anterior e do próprio ano), como o INE. Reproduzem as taxas municipais do INE em todos os 308 municípios."),
      tags$li("Taxa de mortalidade infantil: óbitos com menos de 1 ano por 1.000 nados-vivos, somando três anos de cada."),
      tags$li("Mortalidade proporcional: óbitos de cada grande grupo de causas sobre o total, no triénio. Para todas as idades [I45] lê os totais municipais por causa, completos. Para as idades abaixo de 75 [I46] precisa da repartição por idade: usa a linha regional do INE onde existe (Portugal, Continente, regiões NUTS e as ULS que coincidem com uma NUTS III) e, nas restantes áreas, a soma dos municípios, que reproduz as linhas do INE com um desvio até 0,3 pontos percentuais - excepto nos triénios que incluem 2014, assinalados com \u00a7."),
      tags$li("Esperança de vida à nascença: tábua de mortalidade abreviada (Chiang II) por triénio, com grupos quinquenais até 85 e mais anos, óbitos de todas as causas e população a meio do ano. Os óbitos sem idade publicada num município são distribuídos pelas idades na proporção dos restantes (marca \u2021 quando excedem 2%); as regiões usam as linhas por idade do próprio INE. Reproduz os valores do Eurostat para Portugal, mas fica cerca de 0,8-0,9 anos acima dos publicados pelo INE, que usa outra metodologia; a ordenação das regiões coincide (correlação 0,97)."),
      tags$li("Mortalidade neonatal, neonatal precoce e pós-neonatal: óbitos com menos de 28 dias, menos de 7 dias e de 28 a 364 dias por 1.000 nados-vivos, no triénio. A neonatal e a pós-neonatal somam a infantil."),
      tags$li("Índice sintético de fecundidade: soma das taxas de fecundidade por grupo quinquenal de idade da mãe (15-49 anos), vezes 5. Os nascimentos de mães com menos de 15 anos contam no grupo 15-19 e os de 50 e mais no grupo 45-49."),
      tags$li("Nascimentos em mães com menos de 20 anos, ou com 35 e mais: proporção do total de nados-vivos no triénio. Pré-termo: menos de 37 semanas de gestação, sobre os nascimentos com duração conhecida."),
      tags$li("Beneficiários do RSI e pensionistas por 1.000 habitantes com 15 e mais anos, como no ficheiro de apoio e no INE. Os beneficiários do RSI, contados ao longo do ano, dividem-se pela população média; os pensionistas, contados a 31 de Dezembro, pela população nessa data."),
      tags$li("Valor médio das pensões: soma do valor das pensões (pensionistas vezes o valor médio de cada município) sobre o total de pensionistas."),
      tags$li("Poder de compra per capita (Portugal = 100): a quota de cada município no poder de compra nacional, dividida pela sua quota implícita de população. Para um agrupamento é a soma das quotas sobre a soma das quotas de população, e não a média dos índices municipais."),
      tags$li("Resíduos urbanos por habitante: toneladas recolhidas vezes 1.000, sobre a população média do ano. Loures e Odivelas têm um serviço conjunto (SIMAR) que o INE regista todo em Loures: o indicador só existe para áreas que contenham os dois municípios.")
    ),
    h4("O último ano (provisório, p)"),
    tags$ul(
      tags$li("O INE publica os óbitos por causa e idade cerca de 15 meses depois do fim do ano. Até lá, o último ano é lido de duas fontes do INE que chegam antes: os óbitos por freguesia de residência (0012542), que dão os óbitos de todas as causas de cada município, e os óbitos semanais por região e idade (0012100), que dão a repartição por idade. Em 2020-2024 a soma das freguesias coincide exactamente com os totais do INE; os óbitos por idade reconstruídos ficaram a 0,13% dos publicados, e a taxa padronizada a 0,03%."),
      tags$li("Com isso, a taxa bruta de mortalidade e os óbitos existem para todas as áreas; a esperança de vida, as taxas padronizadas de todas as causas e a SMR, só para Portugal, o Continente, as regiões NUTS II e as áreas feitas de regiões inteiras (Lisboa e Vale do Tejo), e só para ambos os sexos. Por causa (evitável, prematura por causa, mortalidade por causa) não há valores até o INE publicar."),
      tags$li("Os valores cujo período inclui esse ano têm a marca p. A classificação das ULS, o funil e a mortalidade por causa usam o último ano publicado por causa.")
    ),
    h4("Mudanças de série"),
    tags$ul(
      tags$li("Pensões: em 2017 a Série 2017 da segurança social substitui a Série 1990-2023, com cerca de 5,5% menos pensionistas. Os anos antes e depois não são directamente comparáveis; o gráfico de evolução marca a mudança."),
      tags$li("População: a série revista do INE a partir de 2021 (ver o separador Métricas Anuais e o manual)."),
      tags$li("Óbitos com menos de 1 ano por município: completos desde 2011; em 1995-2001 a soma municipal fica abaixo do total nacional e o valor é marcado com \u2020.")
    ),
    h4("Lacunas nos dados municipais do INE"),
    tags$ul(
      tags$li("Uma célula em branco no INE não é um zero. Nos dados publicados para todos os municípios (resíduos, pensões, trabalhadores, ganho médio), um município sem valor deixa sem valor todas as áreas que o contêm, em vez de desaparecer da soma."),
      tags$li("Trabalhadores por conta de outrem e ganho médio: o total do Continente e o de Portugal publicados pelo INE incluem o que ficou oculto nos municípios, pelo que o valor oculto é o total menos os municípios publicados (assinalado com \u2248). Um só município oculto recupera o seu valor (São João da Madeira em 2015); vários partilham-no, e só as áreas que os contêm a todos têm valor (São João da Madeira e Vale de Cambra em 2013: Norte, Área Metropolitana do Porto e ULS Entre Douro e Vouga; os 19 municípios dos Açores em 2013-2014: a Região Autónoma). Os trabalhadores por sector desses municípios continuam em falta. Nas pensões e nos resíduos não é possível, porque aí os totais do INE também deixam de fora os municípios ocultos."),
      tags$li("Odivelas, Trofa e Vizela foram criados em 1998, a partir de Loures, Santo Tirso e Guimarães. A população está estimada com os limites actuais, mas os nascimentos, os óbitos e o Censo de 1991 até 1998 estão registados no município de origem. Até 1998 (e nos triénios que incluem esses anos) os indicadores só existem para áreas que contenham os dois municípios de cada par."),
      tags$li("Óbitos de todas as causas: em alguns municípios e anos (Vimioso em 2024; Alfândega da Fé e Miranda do Douro em 2015) o INE deixa em branco o total de todas as idades, embora publique os grupos etários. O total é então a soma desses grupos."),
      tags$li("Trabalhadores por sector: o INE oculta dois dos três sectores quando um deles revelaria uma empresa (6 municípios em 2024, 41-49 em 2013-2016). O remanescente é repartido pelos sectores ocultos na proporção do resto da NUTS III nesse ano; as quotas com mais de 1% de trabalhadores estimados são marcadas com \u2248.")
    ),
    h4("De onde vêm os dados"),
    tags$ul(
      tags$li("População: estimativas anuais do INE; a série revista (0012918) a partir de 2021."),
      tags$li("Óbitos: o total de todas as idades por município e causa (0008206 até 2021, 0013166 desde 2022). Não se usa a repartição por idade, que o INE publica incompleta ao nível municipal."),
      tags$li("Nados-vivos e óbitos com menos de 1 ano: os mesmos ficheiros da métrica de mortalidade infantil."),
      tags$li("ULS e ARS: somadas a partir dos municípios que as compõem, com a composição actual aplicada a todos os anos.")
    ),
    h4("Porque podem diferir do ficheiro de apoio aos PLS"),
    tags$ul(
      tags$li("A população de 2021 em diante é a série revista pelo INE; valores calculados com a estimativa anterior ficam desactualizados (por exemplo, índices de envelhecimento mais altos)."),
      tags$li("As ULS que partilham um município ao nível da freguesia (Lisboa, Loures, Porto e, numa freguesia, Gondomar) não podem ser calculadas separadamente sem dados por freguesia. Atribuir o município inteiro a cada uma conta a mesma população duas vezes; a aplicação mostra antes os dois agrupamentos exactos."),
      tags$li("As regiões são somadas a partir dos municípios, pelo que excluem os acontecimentos de residência desconhecida, que só entram nos totais de Portugal e do Continente."),
      tags$li("Poder de compra e mortalidade neonatal: nalgumas ULS os valores do ficheiro não coincidem com os publicados pelo INE para os seus municípios (por exemplo, o poder de compra de Matosinhos em 2021 é 118,1 no INE e 130,6 no ficheiro), e a mortalidade neonatal e pós-neonatal do ficheiro não soma a infantil. A aplicação calcula a partir dos valores municipais do INE.")
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
    uiOutput("forecastNotes"),
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
    forecast_table_note(),
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
