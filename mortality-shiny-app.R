# =========================================================
# Packages
# =========================================================
if (!requireNamespace("pacman", quietly = TRUE)) {
  stop(
    "Package 'pacman' is required. Install project dependencies first (for example with renv::restore()).",
    call. = FALSE
  )
}
pacman::p_load(
  glue,
  PHEindicatormethods,
  tidyverse,
  shiny,
  forecast,
  strucchange,   # for breakpoints
  memoise        # for caching INE queries
)

pacman::p_load_gh("c-matos/ineptR")

# =========================================================
# Parameters
# =========================================================

year_of_interest <- 1991:2023

local_area <- c(
  "Portugal",
  "Abrantes", "Águeda", "Aguiar da Beira", "Alandroal", "Albergaria-a-Velha",
  "Albufeira", "Alcácer do Sal", "Alcanena", "Alcobaça", "Alcochete",
  "Alcoutim", "Alenquer", "Alfândega da Fé", "Alijó", "Aljezur", "Aljustrel",
  "Almada", "Almeida", "Almeirim", "Almodôvar", "Alpiarça", "Alter do Chão",
  "Alvaiázere", "Alvito", "Amadora", "Amarante", "Amares", "Anadia",
  "Angra do Heroísmo", "Ansião", "Arcos de Valdevez", "Arganil", "Armamar",
  "Arouca", "Arraiolos", "Arronches", "Arruda dos Vinhos", "Aveiro", "Avis",
  "Azambuja", "Baião", "Barcelos", "Barrancos", "Barreiro", "Batalha", "Beja",
  "Belmonte", "Benavente", "Bombarral", "Borba", "Boticas", "Braga",
  "Bragança", "Cabeceiras de Basto", "Cadaval", "Caldas da Rainha",
  "Calheta", "Câmara de Lobos", "Caminha", "Campo Maior", "Cantanhede",
  "Carrazeda de Ansiães", "Carregal do Sal", "Cartaxo", "Cascais",
  "Castanheira de Pêra", "Castelo Branco", "Castelo de Paiva", "Castelo de Vide",
  "Castro Daire", "Castro Marim", "Castro Verde", "Celorico da Beira",
  "Celorico de Basto", "Chamusca", "Chaves", "Cinfães", "Coimbra",
  "Condeixa-a-Nova", "Constância", "Coruche", "Corvo", "Covilhã", "Crato",
  "Cuba", "Elvas", "Entroncamento", "Espinho", "Esposende", "Estarreja",
  "Estremoz", "Évora", "Fafe", "Faro", "Felgueiras", "Ferreira do Alentejo",
  "Ferreira do Zêzere", "Figueira da Foz", "Figueira de Castelo Rodrigo",
  "Figueiró dos Vinhos", "Fornos de Algodres", "Freixo de Espada à Cinta",
  "Fronteira", "Funchal", "Fundão", "Gavião", "Góis", "Golegã", "Gondomar",
  "Gouveia", "Grândola", "Guarda", "Guimarães", "Horta", "Idanha-a-Nova",
  "Ílhavo", "Lagoa", "Lagos", "Lajes das Flores", "Lajes do Pico", "Lamego",
  "Leiria", "Lisboa", "Loulé", "Loures", "Lourinhã", "Lousã", "Lousada",
  "Mação", "Macedo de Cavaleiros", "Machico", "Madalena", "Mafra", "Maia",
  "Mangualde", "Manteigas", "Marco de Canaveses", "Marinha Grande", "Marvão",
  "Matosinhos", "Mealhada", "Mêda", "Melgaço", "Mértola", "Mesão Frio",
  "Mira", "Miranda do Corvo", "Miranda do Douro", "Mirandela", "Mogadouro",
  "Moimenta da Beira", "Moita", "Monção", "Monchique", "Mondim de Basto",
  "Monforte", "Montalegre", "Montemor-o-Novo", "Montemor-o-Velho", "Montijo",
  "Mora", "Mortágua", "Moura", "Mourão", "Murça", "Murtosa", "Nazaré", "Nelas",
  "Nisa", "Nordeste", "Óbidos", "Odemira", "Odivelas", "Oeiras", "Oleiros",
  "Olhão", "Oliveira de Azeméis", "Oliveira de Frades", "Oliveira do Bairro",
  "Oliveira do Hospital", "Ourém", "Ourique", "Ovar", "Paços de Ferreira",
  "Palmela", "Pampilhosa da Serra", "Paredes", "Paredes de Coura",
  "Pedrógão Grande", "Penacova", "Penafiel", "Penalva do Castelo", "Penamacor",
  "Penedono", "Penela", "Peniche", "Peso da Régua", "Pinhel", "Pombal",
  "Ponta Delgada", "Ponta do Sol", "Ponte da Barca", "Ponte de Lima",
  "Ponte de Sor", "Portalegre", "Portel", "Portimão", "Porto", "Porto de Mós",
  "Porto Moniz", "Porto Santo", "Póvoa de Lanhoso", "Póvoa de Varzim",
  "Povoação", "Proença-a-Nova", "Redondo", "Reguengos de Monsaraz", "Resende",
  "Ribeira Brava", "Ribeira de Pena", "Ribeira Grande", "Rio Maior",
  "Sabrosa", "Sabugal", "Salvaterra de Magos", "Santa Comba Dão", "Santa Cruz",
  "Santa Cruz da Graciosa", "Santa Cruz das Flores", "Santa Maria da Feira",
  "Santa Marta de Penaguião", "Santana", "Santarém", "Santiago do Cacém",
  "Santo Tirso", "São Brás de Alportel", "São João da Madeira",
  "São João da Pesqueira", "São Pedro do Sul", "São Roque do Pico",
  "São Vicente", "Sardoal", "Sátão", "Seia", "Seixal", "Sernancelhe", "Serpa",
  "Sertã", "Sesimbra", "Setúbal", "Sever do Vouga", "Silves", "Sines",
  "Sintra", "Sobral de Monte Agraço", "Soure", "Sousel", "Tábua", "Tabuaço",
  "Tarouca", "Tavira", "Terras de Bouro", "Tomar", "Tondela", "Torre de Moncorvo",
  "Torres Novas", "Torres Vedras", "Trancoso", "Trofa", "Vagos",
  "Vale de Cambra", "Valença", "Valongo", "Valpaços", "Velas", "Vendas Novas",
  "Viana do Alentejo", "Viana do Castelo", "Vidigueira", "Vieira do Minho",
  "Vila da Praia da Vitória", "Vila de Rei", "Vila do Bispo", "Vila do Conde",
  "Vila do Porto", "Vila Flor", "Vila Franca de Xira", "Vila Franca do Campo",
  "Vila Nova da Barquinha", "Vila Nova de Cerveira", "Vila Nova de Famalicão",
  "Vila Nova de Foz Côa", "Vila Nova de Gaia", "Vila Nova de Paiva",
  "Vila Nova de Poiares", "Vila Pouca de Aguiar", "Vila Real",
  "Vila Real de Santo António", "Vila Velha de Ródão", "Vila Verde",
  "Vila Viçosa", "Vimioso", "Vinhais", "Viseu", "Vizela", "Vouzela"
)

diseases <- c(
  "Todas as causas de morte",
  "Doenças do aparelho circulatório",
  "Doenças cérebro-vasculares",
  "Doenças isquémicas do coração",
  "Tumores (neoplasmas) malignos",
  "Tumor (neoplasma) maligno da laringe, da traqueia, dos brônquios e dos pulmões",
  "Tumor (neoplasma) maligno da mama",
  "Doenças do aparelho respiratório",
  "Diabetes mellitus"
)

sex_levels   <- c("HM", "H", "M")
exclude_bands <- c("75 - 79 anos", "80 - 84 anos", "85 e mais anos")
forecast_model_choices <- c(
  "ARIMA" = "arima",
  "ETS" = "ets",
  "Passeio Aleatório com Drift" = "rwf",
  "Naive" = "naive",
  "Theta" = "theta",
  "TBATS" = "tbats",
  "Holt" = "holt",
  "Holt (amortecido)" = "holt_damped"
)

age_levels <- c(
  "0 - 4 anos","5 - 9 anos","10 - 14 anos","15 - 19 anos",
  "20 - 24 anos","25 - 29 anos","30 - 34 anos","35 - 39 anos",
  "40 - 44 anos","45 - 49 anos","50 - 54 anos","55 - 59 anos",
  "60 - 64 anos","65 - 69 anos","70 - 74 anos","75 - 79 anos",
  "80 - 84 anos","85 e mais anos"
)

esp2013_df <- tibble(
  age_band = factor(age_levels, levels = age_levels, ordered = TRUE),
  stdpop   = c(
    5000, 5500, 5500, 5500, 6000, 6000, 6500, 7000,
    7000, 7000, 7000, 6500, 6000, 5500, 5000, 4000, 2500, 2500
  )
)

# =========================================================
# Helpers for INE dimensions
# =========================================================

get_cat_id <- function(value, dimension_values, dim_name = "unknown") {
  requested <- as.character(value)

  ids <- dimension_values %>%
    dplyr::mutate(categ_dsg_chr = as.character(categ_dsg)) %>%
    dplyr::filter(categ_dsg_chr %in% requested) %>%
    dplyr::pull(cat_id) %>%
    unique()

  if (length(ids) == 0) {
    stop(
      glue::glue("Não foi possível mapear categorias para {dim_name}: {paste(requested, collapse = ', ')}"),
      call. = FALSE
    )
  }

  ids
}

# Bounded cache for INE dimension metadata (per indicator)
dim_values_cache <- cachem::cache_mem(
  max_size = 30 * 1024^2,
  max_age  = 24 * 60 * 60
)

get_dim_values_cached <- memoise::memoise(
  ineptR::get_dim_values,
  cache = dim_values_cache
)

# General downloader: can include or exclude cause (dim5)
# ---- replace the whole previous download_data() with this ----
download_data <- function(indicator, dims, has_cause = FALSE) {
  dv   <- get_dim_values_cached(indicator)
  cats <- purrr::imap(dims, ~ get_cat_id(.x, dv, dim_name = .y))
  names(cats) <- names(dims)
  
  # Build argument list without dim5 by default
  args <- list(
    indicator = indicator,
    dim1      = cats$dim1,
    dim2      = cats$dim2
  )
  # Only add dim5 if we actually have it
  if ("dim5" %in% names(cats)) {
    args$dim5 <- cats$dim5
  }
  
  # Call get_ine_data with only the arguments we need
  raw <- do.call(ineptR::get_ine_data, args)
  
  out <- raw %>%
    dplyr::transmute(
      year     = dim_1,
      area     = geodsg,
      sex_raw  = dim_3_t,
      sex      = dplyr::recode(
        sex_raw,
        "Total"    = "HM",
        "Homens"   = "H",
        "Mulheres" = "M",
        .default   = sex_raw
      ),
      age_band = dim_4_t,
      value    = as.numeric(valor)
    )
  
  if (has_cause) {
    out <- out %>%
      dplyr::mutate(cause = raw$dim_5_t)
  }
  
  out
}


# =========================================================
# Compute metrics (crude rates + DSR with proper Poisson CIs)
# =========================================================

compute_metrics <- function(df) {
  # Crude rates and Poisson CIs
  crude <- df %>%
    dplyr::group_by(year, sex, cause) %>%
    dplyr::summarise(
      deaths_total = sum(deaths),
      pop_total    = sum(pop),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      crude_rate = dplyr::if_else(pop_total > 0, deaths_total / pop_total * 1e5, NA_real_),
      ci = purrr::map2(
        deaths_total, pop_total,
        ~ if (is.na(.y) || .y <= 0) c(NA_real_, NA_real_) else stats::poisson.test(.x)$conf.int / .y * 1e5
      ),
      crude_lower = purrr::map_dbl(ci, 1),
      crude_upper = purrr::map_dbl(ci, 2)
    ) %>%
    dplyr::select(-ci)
  
  # Age-standardised (DSR)
  dsr <- df %>%
    dplyr::left_join(esp2013_df, by = "age_band") %>%
    dplyr::group_by(year, sex, cause) %>%
    PHEindicatormethods::calculate_dsr(
      x          = deaths,
      n          = pop,
      stdpop     = stdpop,
      multiplier = 1e5,
      confidence = 0.95
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      dsr       = value,
      dsr_lower = lowercl,
      dsr_upper = uppercl
    )
  
  dplyr::left_join(crude, dsr, by = c("year", "sex", "cause"))
}


# =========================================================
# Lazy INE loader: get data only for selected area + cause
# =========================================================

get_data_for <- function(area, cause) {
  # Population (no cause)
  df_pop1 <- download_data(
    "0008273",
    dims      = list(dim1 = year_of_interest, dim2 = area),
    has_cause = FALSE
  ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = 1L)
  
  df_pop2 <- download_data(
    "0003182",
    dims      = list(dim1 = year_of_interest, dim2 = area),
    has_cause = FALSE
  ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = 2L)
  
  df_pop <- dplyr::bind_rows(df_pop1, df_pop2) %>%
    dplyr::group_by(year, area, sex, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)
  
  # Deaths (with cause)
  df_death1 <- download_data(
    "0008206",
    dims      = list(dim1 = year_of_interest, dim2 = area, dim5 = cause),
    has_cause = TRUE
  ) %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(
      age_band = dplyr::case_when(
        age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ age_band
      )
    ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = 1L)
  
  df_death2 <- download_data(
    "0013166",
    dims      = list(dim1 = year_of_interest, dim2 = area, dim5 = cause),
    has_cause = TRUE
  ) %>%
    dplyr::rename(deaths = value) %>%
    dplyr::mutate(
      age_band = dplyr::case_when(
        age_band %in% c("Menos de 1 ano", "1 - 4 anos") ~ "0 - 4 anos",
        TRUE ~ age_band
      )
    ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(source_priority = 2L)
  
  df_death <- dplyr::bind_rows(df_death1, df_death2) %>%
    dplyr::group_by(year, area, sex, cause, age_band) %>%
    dplyr::arrange(dplyr::desc(source_priority), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-source_priority)
  
  # Combine
  df_full <- df_pop %>%
    dplyr::right_join(df_death,
                      by = c("year", "area", "sex", "age_band")) %>%
    tidyr::replace_na(list(deaths = 0)) %>%
    dplyr::mutate(
      age_band = factor(age_band, levels = age_levels, ordered = TRUE)
    )
  
  df_trunc <- df_full %>%
    dplyr::filter(!age_band %in% exclude_bands)
  
  list(full = df_full, trunc = df_trunc)
}

# Bounded cache to avoid repeated downloads for same (area, cause)
data_query_cache <- cachem::cache_mem(
  max_size = 300 * 1024^2,
  max_age  = 6 * 60 * 60
)

get_data_for_cached <- memoise::memoise(get_data_for, cache = data_query_cache)

get_selection_label <- function(selected_areas, custom_label = NULL) {
  label <- if (is.null(custom_label)) "" else trimws(custom_label)
  if (nzchar(label)) {
    return(label)
  }
  if (length(selected_areas) == 1) {
    return(selected_areas[[1]])
  }
  "Soma de locais"
}

safe_filename_token <- function(x) {
  gsub("[^[:alnum:]_]+", "_", iconv(x, to = "ASCII//TRANSLIT"))
}

get_model_labels <- function(model_ids) {
  labels <- names(forecast_model_choices)[match(model_ids, forecast_model_choices)]
  labels[is.na(labels)] <- model_ids[is.na(labels)]
  labels
}

# =========================================================
# UI Helpers
# =========================================================

# Shared control panels ----------------------------------------------------

forecast_controls_panel <- function() {
  tagList(
    selectInput("area2", "Local de residência:", choices = local_area, multiple = TRUE, selected = "Portugal"),
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
    sliderInput(
      "years_fit",
      "Janela de ajuste:",
      min   = min(year_of_interest),
      max   = max(year_of_interest),
      value = range(year_of_interest),
      step  = 1,
      sep   = ""
    ),
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
    selectInput(
      "beginner_training_window",
      "Janela de ajuste:",
      choices = c(
        "Histórico completo" = "full",
        "Últimos 10 anos" = "10",
        "Últimos 15 anos" = "15",
        "Últimos 20 anos" = "20"
      ),
      selected = "full"
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

observed_mortality_tab_ui <- function() {
  tabPanel(
    "Mortalidade Observada",
    sidebarLayout(
      sidebarPanel(
        selectInput("area", "Local de residência:", choices = local_area, multiple = TRUE, selected = "Portugal"),
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
        actionButton("go_rates", "Carregar dados"),
        br(), br(),
        actionButton("cancel_rates", "Interromper carregamento")
      ),
      mainPanel(
        h4("Resumo"),
        tableOutput("rateSummaryTable"),
        br(),
        plotOutput("ratePlot", height = "400px"),
        br(),
        h4("Série anual observada"),
        tableOutput("rateTable")
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
        uiOutput("beginnerForecastSummary"),
        br(),
        uiOutput("beginnerForecastReliability")
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
        tableOutput("forecastSpecTable")
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
    br(),
    fluidRow(
      column(
        width = 6,
        plotOutput("diagnosticAcfPlot", height = "260px")
      ),
      column(
        width = 6,
        plotOutput("diagnosticPacfPlot", height = "260px")
      )
    ),
    br(),
    h4("Teste de Ljung-Box"),
    tableOutput("diagnosticLjungBoxTable"),
    br(),
    h4("Resumo do Modelo"),
    verbatimTextOutput("diagnosticModelSummary")
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
    br(),
    h4("Valores das Métricas"),
    tableOutput("accuracyTable"),
    br(),
    plotOutput("comparisonPlot", height = "360px")
  )
}

advanced_breaks_tab_ui <- function() {
  tabPanel(
    "Quebras e Estrutura",
    forecast_selection_note_ui(),
    uiOutput("breakInterpretation"),
    br(),
    plotOutput("breakPlot", height = "400px"),
    br(),
    tableOutput("breakTable")
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

# =========================================================
# UI
# =========================================================

ui <- navbarPage(
  title = "Mortalidades e Projecções",

  observed_mortality_tab_ui(),
  beginner_forecasting_tab_ui(),
  advanced_forecasting_tab_ui()
)

# =========================================================
# Server
# =========================================================

server <- function(input, output, session) {
  cancel_seq <- reactiveValues(rates = 0L, forecast = 0L)

  observeEvent(input$cancel_rates, {
    cancel_seq$rates <- cancel_seq$rates + 1L
    showNotification("Pedido de interrupção recebido (Taxas).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  observeEvent(input$cancel_forecast, {
    cancel_seq$forecast <- cancel_seq$forecast + 1L
    showNotification("Pedido de interrupção recebido (Projecções).", type = "warning", duration = 3)
  }, ignoreInit = TRUE)

  abort_if_cancelled <- function(kind, token) {
    if (!identical(cancel_seq[[kind]], token)) {
      validate(need(FALSE, "Operação interrompida. A cache foi preservada."))
    }
  }

  get_rate_mapping <- function(rate_type) {
    if (identical(rate_type, "crude")) {
      list(
        value_col = "crude_rate",
        lower_col = "crude_lower",
        upper_col = "crude_upper",
        y_label = "Taxa Bruta por 100.000",
        rate_label = "Bruta"
      )
    } else {
      list(
        value_col = "dsr",
        lower_col = "dsr_lower",
        upper_col = "dsr_upper",
        y_label = "Taxa Padronizada por 100.000",
        rate_label = "Padronizada"
      )
    }
  }

  # Shared historical-series pipeline:
  # 1) Build a stable query spec with the inputs that determine which raw data
  #    need to be downloaded and aggregated.
  # 2) Load one metric bundle (both population scopes, both rate outputs) for
  #    that query.
  # 3) Apply the final population/rate/year filter to get a reusable series.
  make_query_spec <- function(area, area_label, cause, sex) {
    area_key <- sort(unique(area))

    list(
      area_key = area_key,
      area_label = get_selection_label(area_key, area_label),
      cause = cause,
      sex = sex
    )
  }

  make_series_spec <- function(query_spec, population, rate_type, year_range = range(year_of_interest)) {
    c(
      query_spec,
      list(
        population = population,
        rate_type = rate_type,
        year_range = range(year_range)
      )
    )
  }

  load_metric_bundle <- function(query_spec, kind, token) {
    validate(need(length(query_spec$area_key) > 0, "Selecione pelo menos um local de residência."))

    abort_if_cancelled(kind, token)
    incProgress(0.1)

    dat <- get_data_for_cached(query_spec$area_key, query_spec$cause)

    abort_if_cancelled(kind, token)
    incProgress(0.5)

    df_full <- dat$full %>%
      dplyr::filter(sex == query_spec$sex)
    df_trunc <- dat$trunc %>%
      dplyr::filter(sex == query_spec$sex)

    metrics <- dplyr::bind_rows(
      compute_metrics(df_full) %>% dplyr::mutate(População = "Total"),
      compute_metrics(df_trunc) %>% dplyr::mutate(População = "Menos de 75 anos")
    ) %>%
      dplyr::arrange(year)

    abort_if_cancelled(kind, token)
    incProgress(0.4)

    list(
      query_spec = query_spec,
      metrics = metrics
    )
  }

  build_historical_series <- function(metric_bundle, series_spec) {
    rate_map <- get_rate_mapping(series_spec$rate_type)

    series <- metric_bundle$metrics %>%
      dplyr::filter(
        População == series_spec$population,
        year >= series_spec$year_range[1],
        year <= series_spec$year_range[2]
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

  build_accuracy_table <- function(fits) {
    metrics <- c("ME", "RMSE", "MAE", "MAPE", "MASE")

    if (length(fits) == 0) {
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

    acc_list <- lapply(fits, function(fit) {
      acc <- forecast::accuracy(fit)
      acc[1, metrics, drop = FALSE]
    })

    dplyr::bind_rows(lapply(names(acc_list), function(m) {
      df <- as.data.frame(acc_list[[m]])
      df$Model <- m
      df
    }), .id = NULL) %>%
      dplyr::select(Model, dplyr::everything())
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
      geom_line(color = "#0b2e4f", size = 0.8) +
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

  compute_forecast_error_metrics <- function(actual, predicted, scale_denom = NA_real_) {
    validate(need(length(actual) == length(predicted), "As séries observada e prevista têm de ter o mesmo comprimento."))

    err <- predicted - actual
    nonzero_actual <- abs(actual) > 1e-8

    tibble(
      ME = mean(err, na.rm = TRUE),
      RMSE = sqrt(mean(err^2, na.rm = TRUE)),
      MAE = mean(abs(err), na.rm = TRUE),
      MAPE = if (any(nonzero_actual)) mean(abs(err[nonzero_actual] / actual[nonzero_actual]), na.rm = TRUE) * 100 else NA_real_,
      MASE = if (is.finite(scale_denom) && scale_denom > 0) mean(abs(err), na.rm = TRUE) / scale_denom else NA_real_
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

  build_fitted_values_df <- function(obs, fits) {
    if (length(fits) == 0) {
      return(tibble(year = numeric(0), fitted = numeric(0), model = character(0)))
    }

    dplyr::bind_rows(lapply(names(fits), function(model_id) {
      fitted_vals <- as.numeric(stats::fitted(fits[[model_id]]))
      year_count <- min(length(fitted_vals), nrow(obs))

      tibble(
        year = obs$year[seq_len(year_count)],
        fitted = fitted_vals[seq_len(year_count)],
        model = model_id
      )
    }))
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
            size = 0.8
          ) +
          geom_line(
            data = comparison_dat$holdout_actual,
            aes(x = year, y = value),
            color = "black",
            size = 1
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
            size = 1
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

    fitted_df <- build_fitted_values_df(
      obs = comparison_dat$obs,
      fits = comparison_dat$fits
    )

    validate(need(nrow(fitted_df) > 0, "Não existem valores ajustados disponíveis para a comparação de modelos."))

    ggplot() +
      geom_line(
        data = comparison_dat$obs,
        aes(x = year, y = value),
        color = "black",
        size = 1
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
    bp_obj <- tryCatch(
      strucchange::breakpoints(ts_y ~ 1),
      error = function(e) NULL
    )

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
      `Ano da Quebra` = c(break_years, NA_integer_)
    )

    list(
      history = history,
      series = df,
      breakpoints = bp_obj,
      break_index = break_index,
      break_years = break_years,
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
      geom_line(color = "#0b2e4f", size = 1) +
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

    if (break_n == 0) {
      return("A rotina actual de breakpoint não detectou qualquer potencial mudança estrutural, pelo que o nível médio parece relativamente estável no histórico seleccionado.")
    }

    if (break_n == 1) {
      return(glue::glue(
        "Foi detectada uma potencial mudança estrutural em torno de {break_info$break_years[[1]]}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que o nível médio de mortalidade poderá ter mudado nessa altura."
      ))
    }

    glue::glue(
      "Foram detectadas {break_n} potenciais mudanças estruturais em torno de {paste(break_info$break_years, collapse = ', ')}, dividindo o histórico seleccionado em {segment_n} segmentos. Isto sugere que o nível médio poderá não ser estável ao longo de todo o período observado."
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
    value_or_default(dat$fitted_models, character(0))
  }

  get_default_forecast_output_model <- function(dat) {
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
          geom_line(data = obs, aes(x = year, y = value), size = 1) +
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
            size = 1
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
      geom_line(data = obs, aes(x = year, y = value), size = 1) +
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
        size = 1
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
    failures <- dat$failures
    successful_models <- get_successful_model_ids(dat)

    if (nrow(failures) == 0) {
      return(NULL)
    }

    wellPanel(
      h4("Avisos dos modelos"),
      p(
        if (length(successful_models) == 0) {
          "Nenhum dos modelos pedidos pôde ser estimado para esta especificação."
        } else {
          "Alguns dos modelos pedidos não puderam ser estimados. Os resultados abaixo utilizam os modelos que foram ajustados com sucesso."
        }
      ),
      tags$ul(
        lapply(seq_len(nrow(failures)), function(i) {
          tags$li(
            paste0(get_model_labels(failures$Model[[i]]), ": ", failures$Message[[i]])
          )
        })
      )
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
        forward = function(x) as.numeric(x),
        inverse = function(x) as.numeric(x)
      ))
    }

    min_positive <- suppressWarnings(min(values[values > 0], na.rm = TRUE))
    log_offset <- if (is.finite(min_positive)) min_positive / 2 else 1e-6

    if (any(values < 0, na.rm = TRUE)) {
      validate(
        need(FALSE, "A série contém valores negativos e não pode usar a transformação log.")
      )
    }

    list(
      method = "log_offset",
      label = "Transformação log com offset",
      forward = function(x) log(as.numeric(x) + log_offset),
      inverse = function(x) pmax(exp(as.numeric(x)) - log_offset, 0)
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
    model_specs = list()
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

      tibble(
        year  = seq(max(df_vals$year) + 1, by = 1, length.out = horizon),
        model = m,
        mean  = back_transform(fc$mean),
        lower = back_transform(
          if (is.null(dim(fc$lower))) as.numeric(fc$lower) else as.numeric(fc$lower[, ncol(fc$lower)])
        ),
        upper = back_transform(
          if (is.null(dim(fc$upper))) as.numeric(fc$upper) else as.numeric(fc$upper[, ncol(fc$upper)])
        )
      )
    }))

    failure_tbl <- tibble(
      Model = failed_ids,
      Message = vapply(fit_results[failed_ids], `[[`, character(1), "message")
    )

    accuracy_tbl <- build_accuracy_table(fits)

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
      transform_label = transform_setup$label,
      model_specs = model_specs,
      area_label = history$area_label
    )
  }

  get_beginner_training_range <- function(years, training_window) {
    year_end <- max(years)

    if (identical(training_window, "full")) {
      return(range(years))
    }

    window_n <- as.integer(training_window)
    c(max(min(years), year_end - window_n + 1L), year_end)
  }

  get_beginner_training_label <- function(training_window, training_history) {
    if (identical(training_window, "full")) {
      "todo o histórico observado"
    } else {
      paste0("os últimos ", nrow(training_history$series), " anos observados")
    }
  }

  build_beginner_forecast_plot <- function(dat) {
    full_obs <- dat$full_history$series
    train_obs <- dat$obs
    recommended_fc <- dat$fc %>%
      dplyr::filter(model == dat$recommended_model)

    p <- ggplot() +
      geom_line(
        data = full_obs,
        aes(x = year, y = value),
        color = "grey75",
        size = 0.8
      ) +
      geom_line(
        data = train_obs,
        aes(x = year, y = value),
        color = "#1f4e79",
        size = 1
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
        size = 1.1
      ) +
      labs(
        title = paste("Previsão Guiada -", dat$area_label),
        subtitle = if (identical(dat$mode, "compare")) {
          "As linhas cinzentas mostram trajectórias alternativas dos modelos."
        } else {
          NULL
        },
        x = "Ano",
        y = dat$history$y_label
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))

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

    wellPanel(
      h4("Fiabilidade"),
      p(tags$strong(paste("Leitura global:", overall_label))),
      p(training_message),
      p(interval_message),
      p(agreement_message),
      if (!is.null(structural_warning)) {
        p(tags$strong(structural_warning))
      }
    )
  }
  
  # -------------------------
  # Observed Mortality
  # -------------------------
  observed_metric_bundle <- eventReactive(input$go_rates, {
    token <- isolate(cancel_seq$rates)
    query_spec <- make_query_spec(input$area, input$area_label, input$cause, input$sex)

    shiny::withProgress(message = "A obter dados do INE...", value = 0, {
      load_metric_bundle(query_spec, "rates", token)
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
      year_range = range(year_of_interest)
    )
  })

  # `observed_history()` is the final shared series object used by the observed
  # plot/table. It is rebuilt from the cached metric bundle, not from raw data.
  observed_history <- reactive({
    req(input$go_rates > 0)
    build_historical_series(observed_metric_bundle(), observed_series_spec())
  })

  observed_summary <- reactive({
    df <- observed_history()$series
    req(nrow(df) > 0)
    first_row <- dplyr::slice(df, 1)
    last_row <- dplyr::slice(df, nrow(df))

    absolute_change <- last_row$value - first_row$value
    percent_change <- if (isTRUE(all.equal(first_row$value, 0))) {
      NA_real_
    } else {
      (absolute_change / first_row$value) * 100
    }

    tibble(
      Métrica = c(
        "Último ano",
        "Última taxa",
        "Variação absoluta no período observado",
        "Variação percentual no período observado"
      ),
      Valor = c(
        as.character(last_row$year),
        sprintf("%.2f", last_row$value),
        sprintf("%+.2f", absolute_change),
        if (is.na(percent_change)) "N/D" else sprintf("%+.2f%%", percent_change)
      )
    )
  })

  output$rateSummaryTable <- renderTable({
    observed_summary()
  }, bordered = TRUE, spacing = "s")
  
  output$ratePlot <- renderPlot({
    dat <- observed_history()
    df <- dat$series
    
    ggplot(df, aes(x = year, group = 1)) +
      geom_ribbon(
        aes(ymin = lower, ymax = upper),
        fill = "grey80", alpha = 0.4
      ) +
      geom_line(aes(y = value), size = 1) +
      geom_point(aes(y = value), size = 2) +
      labs(
        title = paste(dat$area_label, "-", dat$spec$cause, dat$spec$sex, "(", dat$spec$population, ")"),
        x = "Ano",
        y = dat$y_label
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$rateTable <- renderTable({
    df <- observed_history()$series
    
    df %>%
      dplyr::transmute(
        Ano = year,
        `Taxa (Intervalo de Confiança 95%)` = glue::glue(
          "{round(value, 2)} ({round(lower, 2)}; {round(upper, 2)})"
        )
      )
  }, sanitize.text.function = identity)
  
  # -------------------------
  # Beginner Forecasting
  # -------------------------
  # The guided workflow reuses the currently loaded observed series and only asks
  # for a horizon, a training window, and whether the user wants a simple
  # recommendation or a visual comparison of reasonable alternatives.
  beginner_training_history <- eventReactive(input$go_beginner_forecast, {
    validate(need(input$go_rates > 0, "Carregue primeiro uma série em 'Mortalidade Observada'."))

    base_history <- observed_history()
    training_range <- get_beginner_training_range(
      years = base_history$series$year,
      training_window = input$beginner_training_window
    )

    build_historical_series(
      metric_bundle = observed_metric_bundle(),
      series_spec = make_series_spec(
        query_spec = observed_metric_bundle()$query_spec,
        population = base_history$spec$population,
        rate_type = base_history$spec$rate_type,
        year_range = training_range
      )
    )
  }, ignoreNULL = TRUE)

  beginner_forecast <- eventReactive(input$go_beginner_forecast, {
    guided_result <- run_forecast_models(
      history = beginner_training_history(),
      model_ids = unname(forecast_model_choices),
      horizon = input$beginner_horizon
    )

    validate(
      need(length(guided_result$fitted_models) > 0, "Não foi possível estimar qualquer modelo de projecção para esta série.")
    )

    c(
      guided_result,
      list(
        full_history = observed_history(),
        horizon = input$beginner_horizon,
        mode = input$beginner_mode,
        training_label = get_beginner_training_label(
          training_window = input$beginner_training_window,
          training_history = beginner_training_history()
        )
      )
    )
  }, ignoreNULL = TRUE)

  output$beginnerForecastPlot <- renderPlot({
    req(input$go_beginner_forecast > 0)
    build_beginner_forecast_plot(beginner_forecast())
  })

  output$beginnerForecastSummary <- renderUI({
    req(input$go_beginner_forecast > 0)
    build_beginner_summary_ui(beginner_forecast())
  })

  output$beginnerForecastReliability <- renderUI({
    req(input$go_beginner_forecast > 0)
    build_beginner_reliability_ui(beginner_forecast())
  })

  # -------------------------
  # Advanced Forecasting: Shared Frozen Specification
  # -------------------------
  # Freeze the filtered historical series at click time so the advanced tabs all
  # reflect the same technical specification until the user reruns the model.
  forecast_history <- eventReactive(input$go_forecast, {
    token <- isolate(cancel_seq$forecast)
    query_spec <- make_query_spec(input$area2, input$area_label2, input$cause2, input$sex2)

    shiny::withProgress(message = "A obter dados do INE...", value = 0, {
      metric_bundle <- load_metric_bundle(query_spec, "forecast", token)
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
      model_specs = build_advanced_model_specs()
    )
  }, ignoreNULL = TRUE)

  # These thin wrappers keep the downstream renderers simple and make it clear
  # that all advanced sub-tabs are reading from one frozen result object.
  advanced_forecast_result <- reactive({
    req(input$go_forecast > 0)
    forecast_sel()
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
  output$forecastSpecTable <- renderTable({
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
          spec$population,
          dat$history$rate_label,
          glue::glue("{spec$year_range[1]} - {spec$year_range[2]}"),
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
  }, striped = TRUE, bordered = TRUE, spacing = "s")

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

  output$forecastPlot <- renderPlot({
    build_advanced_forecast_plot(
      dat = advanced_forecast_result(),
      view_mode = input$forecast_output_view,
      selected_model = advanced_forecast_focus_model()
    )
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

  output$diagnosticAcfPlot <- renderPlot({
    plot_diagnostic_correlation(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model(),
      partial = FALSE
    )
  })

  output$diagnosticPacfPlot <- renderPlot({
    plot_diagnostic_correlation(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model(),
      partial = TRUE
    )
  })

  output$diagnosticLjungBoxTable <- renderTable({
    build_ljung_box_table(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model()
    )
  }, digits = 4, striped = TRUE, bordered = TRUE, spacing = "s")

  output$diagnosticModelSummary <- renderText({
    build_diagnostic_model_summary(
      dat = advanced_forecast_result(),
      model_id = advanced_diagnostic_model()
    )
  })

  # -------------------------
  # Advanced Forecasting: Backtesting & Comparison
  # -------------------------
  output$comparisonHoldoutControl <- renderUI({
    if (!identical(input$comparison_validation_mode, "holdout")) {
      return(NULL)
    }

    series_n <- nrow(advanced_forecast_history()$series)
    max_k <- min(10, series_n - 3)

    if (max_k < 1) {
      return(
        p("A validação nos últimos anos requer pelo menos 4 anos observados no histórico de ajuste seleccionado.")
      )
    }

    numericInput(
      "comparison_holdout_k",
      "Últimos k anos para validação:",
      value = min(5, max_k),
      min = 1,
      max = max_k,
      step = 1
    )
  })

  advanced_comparison <- reactive({
    base_result <- advanced_forecast_result()
    mode <- input$comparison_validation_mode

    if (!identical(mode, "holdout")) {
      metric_tbl <- base_result$accuracy
      ranking_tbl <- rank_model_metrics(metric_tbl)

      return(list(
        mode = "insample",
        metrics = metric_tbl,
        ranking = ranking_tbl,
        failures = base_result$failures,
        fits = base_result$fits,
        obs = base_result$obs,
        y_label = base_result$history$y_label
      ))
    }

    max_holdout_k <- min(10, nrow(advanced_forecast_history()$series) - 3)
    validate(
      need(
        max_holdout_k >= 1,
        "A validação nos últimos anos requer pelo menos 4 anos observados no histórico de ajuste seleccionado."
      )
    )

    holdout_k <- value_or_default(input$comparison_holdout_k, min(5, max_holdout_k))
    split <- build_holdout_split(
      history = base_result$history,
      holdout_k = holdout_k
    )

    holdout_result <- run_forecast_models(
      history = split$training_history,
      model_ids = base_result$selected_model_ids,
      horizon = nrow(split$holdout_actual),
      conf_level = base_result$conf_level,
      transform_method = base_result$transform_method,
      model_specs = base_result$model_specs
    )

    metric_tbl <- build_holdout_metric_table(
      training_result = holdout_result,
      holdout_actual = split$holdout_actual
    )

    list(
      mode = "holdout",
      metrics = metric_tbl,
      ranking = rank_model_metrics(metric_tbl),
      failures = holdout_result$failures,
      training_obs = split$training_history$series,
      holdout_actual = split$holdout_actual,
      forecast_df = holdout_result$fc,
      holdout_k = nrow(split$holdout_actual),
      y_label = base_result$history$y_label
    )
  })

  output$comparisonWarnings <- renderUI({
    build_forecast_warning_ui(advanced_comparison())
  })

  output$comparisonRankingTable <- renderTable({
    comparison_dat <- advanced_comparison()
    validate(need(nrow(comparison_dat$ranking) > 0, "Nenhum modelo foi estimado com sucesso para a configuração de comparação seleccionada."))

    comparison_dat$ranking %>%
      dplyr::mutate(Model = get_model_labels(Model)) %>%
      dplyr::rename(Modelo = Model)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$accuracyTable <- renderTable({
    comparison_dat <- advanced_comparison()
    validate(need(nrow(comparison_dat$metrics) > 0, "Nenhum modelo foi estimado com sucesso para a configuração de comparação seleccionada."))

    comparison_dat$metrics %>%
      dplyr::left_join(comparison_dat$ranking %>% dplyr::select(Model, Classificação), by = "Model") %>%
      dplyr::mutate(Model = get_model_labels(Model)) %>%
      dplyr::rename(Modelo = Model) %>%
      dplyr::select(Classificação, Modelo, dplyr::everything()) %>%
      dplyr::arrange(Classificação)
  }, digits = 3)

  output$comparisonPlot <- renderPlot({
    build_comparison_plot(advanced_comparison())
  })
  
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

  output$breakTable <- renderTable({
    advanced_break_analysis()$segments
  }, striped = TRUE, bordered = TRUE, spacing = "s", digits = 2)
}

# =========================================================
# Run app
# =========================================================

shinyApp(ui, server)
