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
# wherever a regional figure is read.
nuts_vintage_control <- function() {
  tags$div(
    class = "nuts-vintage-header",
    style = paste(
      "display:flex; align-items:center; gap:0.6rem; flex-wrap:wrap;",
      "padding:0.35rem 0.9rem; margin:0 0 0.4rem 0;",
      "border-bottom:1px solid rgba(128,128,128,0.25);"
    ),
    tags$label(
      "Definição das regiões:",
      `for` = "nuts_vintage",
      style = "margin:0; font-weight:600; white-space:nowrap;"
    ),
    tags$div(
      style = "min-width:16rem;",
      selectInput(
        "nuts_vintage",
        label = NULL,
        choices = nuts_vintage_choices,
        selected = default_nuts_vintage(),
        width = "100%"
      )
    ),
    tags$span(
      style = "font-size:0.85em; opacity:0.75;",
      "Agrupa os mesmos 308 municípios de outra forma; não altera os dados lidos."
    ),
    # How a region's deaths are built. Sits beside the vintage because both
    # change what a regional figure means, and both apply to every tab.
    tags$label(
      "Óbitos das regiões:",
      `for` = "region_source",
      style = "margin:0 0 0 0.8rem; font-weight:600; white-space:nowrap;"
    ),
    tags$div(
      style = "min-width:17rem;",
      selectInput(
        "region_source",
        label = NULL,
        choices = region_source_choices,
        selected = default_region_source(),
        width = "100%"
      )
    ),
    data_version_badge()
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
    region_aggregation_note(),
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
    list("ULS (Unidade Local de Saúde)", "A unidade de organização do SNS a que corresponde a população de um conjunto de municípios. Na aplicação é somada a partir dos seus municípios. Cinco ULS partilham municípios ao nível da freguesia e aparecem em dois agrupamentos exactos."),
    list("ARS", "As cinco regiões de saúde (Norte, Centro, Lisboa e Vale do Tejo, Alentejo, Algarve) que agrupam as ULS. Não coincidem com as regiões NUTS com o mesmo nome."),
    list("NUTS", "A nomenclatura estatística das regiões. A aplicação usa dois níveis: NUTS I (Continente, Açores, Madeira) e NUTS II (as regiões)."),
    list("Definição das regiões (NUTS 2013 / NUTS 2024)", "As duas versões da nomenclatura que a aplicação oferece, no controlo do topo da página. Agrupam os mesmos 308 municípios de formas diferentes: seis nomes existem nas duas e significam coisas diferentes em cada uma."),
    list("Área Metropolitana de Lisboa", "A região de Lisboa em NUTS 2013. Em NUTS 2024 está dividida em Grande Lisboa e Península de Setúbal."),
    list("Oeste e Vale do Tejo", "Região criada em NUTS 2024 com o Oeste, o Médio Tejo e a Lezíria do Tejo. Apesar do nome, não inclui Lisboa."),
    list("Agregação por municípios", "As regiões são sempre somadas a partir dos seus municípios, com a mesma lista aplicada a todos os anos. Mantém a série contínua apesar da revisão NUTS de 2024, mas os totais não coincidem exactamente com os publicados pelo INE."),
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
        region_aggregation_note(),
        year_range_slider(
          years = rate_years,
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
          # Portugal and Norte are always the first two columns.
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
        actionButton("cancel_annual_metrics", "Interromper carregamento")
      ),
      mainPanel(
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
          choices = sort(rate_years),
          selected = max(rate_years)
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
          "assumindo que as quotas se mantêm. Os agrupamentos exactos de ULS não dependem da escolha."
        ),
        radioButtons(
          "planning_split_mode", "ULS que partilham um município:",
          choices = PLANNING_SPLIT_MODES, selected = PLANNING_DEFAULT_SPLIT_MODE
        ),
        helpText(
          "Lisboa, Loures e Porto estão repartidos por ULS ao nível da freguesia. ",
          "«Município inteiro»: cada ULS leva o município todo — nada é estimado, mas as seis ",
          "ULS sobrepõem-se e a sua soma conta esses municípios mais do que uma vez. ",
          "«Ponderação por freguesias»: cada município é repartido pela população das suas freguesias ",
          "nos Censos de 2021, por grupo etário — as partes somam o país, assumindo que as quotas se mantêm. ",
          "Os agrupamentos exactos de ULS não dependem desta escolha."
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
            div(style = "overflow-x:auto;", tableOutput("planningIndicatorTable"))
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
            "Funil",
            br(),
            radioButtons(
              "planning_funnel_units", "Unidades:",
              choices = c("ULS" = "ULS", "Municípios" = "Município", "NUTS III" = "NUTS III"),
              selected = "ULS", inline = TRUE
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
            plotly::plotlyOutput("planningSnsRanking", height = "820px")
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
                "indicadores e anos, uma folha por indicador. Pode demorar cerca de ",
                "um minuto e meio a gerar da primeira vez."
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
      tags$li("Os óbitos por idade e causa de cada município estão incompletos no INE; são completados até ao total de todas as idades, e os que faltam distribuídos pelas idades com o perfil dos que faltam no país, de modo que a soma dos municípios reproduz Portugal por idade. \u2021 quando mais de 2% foram redistribuídos.")
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
      tags$li("Esperança de vida à nascença: tábua de mortalidade abreviada (Chiang II) por triénio, com grupos quinquenais até 85 e mais anos, óbitos de todas as causas e população a meio do ano. Os óbitos sem idade publicada num município são distribuídos pelas idades na proporção dos restantes (marca \u2021 quando excedem 2%). Reproduz os valores do Eurostat para Portugal, mas fica cerca de 0,8-0,9 anos acima dos publicados pelo INE, que usa outra metodologia; a ordenação das regiões coincide (correlação 0,97)."),
      tags$li("Mortalidade neonatal, neonatal precoce e pós-neonatal: óbitos com menos de 28 dias, menos de 7 dias e de 28 a 364 dias por 1.000 nados-vivos, no triénio. A neonatal e a pós-neonatal somam a infantil."),
      tags$li("Índice sintético de fecundidade: soma das taxas de fecundidade por grupo quinquenal de idade da mãe (15-49 anos), vezes 5. Os nascimentos de mães com menos de 15 anos contam no grupo 15-19 e os de 50 e mais no grupo 45-49."),
      tags$li("Nascimentos em mães com menos de 20 anos, ou com 35 e mais: proporção do total de nados-vivos no triénio. Pré-termo: menos de 37 semanas de gestação, sobre os nascimentos com duração conhecida."),
      tags$li("Beneficiários do RSI e pensionistas por 1.000 habitantes com 15 e mais anos, como no ficheiro de apoio e no INE. Os beneficiários do RSI, contados ao longo do ano, dividem-se pela população média; os pensionistas, contados a 31 de Dezembro, pela população nessa data."),
      tags$li("Valor médio das pensões: soma do valor das pensões (pensionistas vezes o valor médio de cada município) sobre o total de pensionistas."),
      tags$li("Poder de compra per capita (Portugal = 100): a quota de cada município no poder de compra nacional, dividida pela sua quota implícita de população. Para um agrupamento é a soma das quotas sobre a soma das quotas de população, e não a média dos índices municipais."),
      tags$li("Resíduos urbanos por habitante: toneladas recolhidas vezes 1.000, sobre a população média do ano. Loures e Odivelas têm um serviço conjunto (SIMAR) que o INE regista todo em Loures: o indicador só existe para áreas que contenham os dois municípios.")
    ),
    h4("Mudanças de série"),
    tags$ul(
      tags$li("Pensões: em 2017 a Série 2017 da segurança social substitui a Série 1990-2023, com cerca de 5,5% menos pensionistas. Os anos antes e depois não são directamente comparáveis; o gráfico de evolução marca a mudança."),
      tags$li("População: a série revista do INE a partir de 2021 (ver o separador Métricas Anuais e o manual)."),
      tags$li("Óbitos com menos de 1 ano por município: completos desde 2011; em 1995-2001 a soma municipal fica abaixo do total nacional e o valor é marcado com \u2020.")
    ),
    h4("Lacunas nos dados municipais do INE"),
    tags$ul(
      tags$li("Uma célula em branco no INE não é um zero. Nos dados publicados para todos os municípios (resíduos, pensões, trabalhadores, ganho médio), um município sem valor deixa sem valor todas as áreas que o contêm, em vez de desaparecer da soma: por exemplo, os Açores não têm trabalhadores por conta de outrem em 2013-2014, nem ganho médio em 2011-2014."),
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
      tags$li("As ULS que partilham um município ao nível da freguesia (Lisboa, Loures, Porto) não podem ser calculadas separadamente sem dados por freguesia. Atribuir o município inteiro a cada uma conta a mesma população duas vezes; a aplicação mostra antes os dois agrupamentos exactos."),
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
