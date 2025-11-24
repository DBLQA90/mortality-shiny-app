# =========================================================
# Packages
# =========================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readxl,
  glue,
  expss,
  PHEindicatormethods,
  tidyverse,
  shiny,
  forecast,
  ineptR,        # for INE data
  strucchange,   # for breakpoints
  memoise        # for caching INE queries
)

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

get_cat_id <- function(value, dimension_values) {
  ids <- dimension_values %>%
    dplyr::filter(categ_dsg %in% value) %>%
    dplyr::pull(cat_id)
  if (length(ids) == 0) "" else ids
}

# General downloader: can include or exclude cause (dim5)
# ---- replace the whole previous download_data() with this ----
download_data <- function(indicator, dims, has_cause = FALSE) {
  dv   <- ineptR::get_dim_values(indicator)
  cats <- purrr::map(dims, ~ get_cat_id(.x, dv))
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
    dplyr::rowwise() %>%
    dplyr::mutate(
      crude_rate  = deaths_total / pop_total * 1e5,
      crude_lower = stats::poisson.test(deaths_total)$conf.int[1] / pop_total * 1e5,
      crude_upper = stats::poisson.test(deaths_total)$conf.int[2] / pop_total * 1e5
    ) %>%
    dplyr::ungroup()
  
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
    dplyr::rename(pop = value)
  
  df_pop2 <- download_data(
    "0003182",
    dims      = list(dim1 = year_of_interest, dim2 = area),
    has_cause = FALSE
  ) %>%
    dplyr::filter(!age_band %in% c("Idade ignorada", "Total")) %>%
    dplyr::rename(pop = value)
  
  df_pop <- dplyr::bind_rows(df_pop1, df_pop2)
  
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
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  
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
    dplyr::summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop")
  
  df_death <- dplyr::bind_rows(df_death1, df_death2)
  
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

# Cache to avoid repeated downloads for same (area, cause)
get_data_for_cached <- memoise::memoise(get_data_for)

# =========================================================
# UI
# =========================================================

ui <- navbarPage(
  title = "PNS Monitorização não oficial",
  
  tabPanel(
    "Taxas de Mortalidade",
    sidebarLayout(
      sidebarPanel(
        selectInput("area", "Local de residência:", choices = local_area),
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
        actionButton("go_rates", "Carregar dados")
      ),
      mainPanel(
        plotOutput("ratePlot", height = "400px"),
        br(),
        tableOutput("rateTable")
      )
    )
  ),
  
  tabPanel(
    "Projecções",
    sidebarLayout(
      sidebarPanel(
        selectInput("area2", "Local de residência:", choices = local_area),
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
          "Anos para ajustar o modelo:",
          min   = min(year_of_interest),
          max   = max(year_of_interest),
          value = range(year_of_interest),
          step  = 1,
          sep   = ""   # avoids 1,991 formatting
        ),
        checkboxGroupInput(
          "models",
          "Modelos de projecção:",
          choices = c(
            "ARIMA"    = "arima",
            "ETS"      = "ets",
            "RW drift" = "rwf",
            "Naive"    = "naive",
            "Theta"    = "theta",
            "TBATS"    = "tbats"
          ),
          selected = c("arima", "ets")
        ),
        sliderInput(
          "horizon",
          "Horizonte Temporal (anos):",
          min   = 1,
          max   = 8,
          value = 7
        ),
        actionButton("go_forecast", "Carregar projecções")
      ),
      mainPanel(
        plotOutput('forecastPlot',height='400px'),
        br(),
        downloadButton("downloadForecastPlot", "Descarregar gráfico (PNG)"),
        br(), br(),
        tableOutput('forecastTable'),
        br(),
        downloadButton("downloadForecastCSV", "Descarregar tabela (CSV)"),
        br(), br(),
        tableOutput('accuracyTable'),
        br(),
        plotOutput('residPlot',height='300px'),
        br(),
        tableOutput('residTable')
      )
    )
  ),
  
  tabPanel(
    "Análise de Quebras",
    sidebarLayout(
      sidebarPanel(
        selectInput("area3", "Local de residência:", choices = local_area),
        selectInput("cause3", "Causa de Morte:", choices = diseases),
        selectInput("sex3", "Sexo:", choices = sex_levels, selected = "HM"),
        radioButtons(
          "population3",
          "População:",
          choices = list("Total" = "Total", "Menos de 75 anos" = "Menos de 75 anos")
        ),
        radioButtons(
          "rate_type3",
          "Taxa:",
          choices = c("Bruta" = "crude", "Padronizada" = "dsr")
        ),
        actionButton("go_breaks", "Carregar análise")
      ),
      mainPanel(
        plotOutput("breakPlot", height = "400px"),
        br(),
        tableOutput("breakTable")
      )
    )
  )
)

# =========================================================
# Server
# =========================================================

server <- function(input, output, session) {
  
  # -------------------------
  # Rates Explorer
  # -------------------------
  rates_sel <- eventReactive(input$go_rates, {
    shiny::withProgress(message = "A obter dados do INE...", value = 0, {
      incProgress(0.1)
      dat <- get_data_for_cached(input$area, input$cause)
      incProgress(0.5)
      
      df1 <- dat$full  %>% dplyr::filter(sex == input$sex)
      df2 <- dat$trunc %>% dplyr::filter(sex == input$sex)
      
      res <- dplyr::bind_rows(
        compute_metrics(df1) %>% dplyr::mutate(População = "Total"),
        compute_metrics(df2) %>% dplyr::mutate(População = "Menos de 75 anos")
      )
      
      incProgress(0.4)
      res
    })
  })
  
  output$ratePlot <- renderPlot({
    req(input$go_rates > 0)
    df <- rates_sel() %>%
      dplyr::filter(População == input$population)
    
    aes_y  <- if (input$rate_type == "crude") "crude_rate"  else "dsr"
    aes_lo <- if (input$rate_type == "crude") "crude_lower" else "dsr_lower"
    aes_hi <- if (input$rate_type == "crude") "crude_upper" else "dsr_upper"
    ylbl   <- if (input$rate_type == "crude")
      "Taxa Bruta por 100.000"
    else
      "Taxa Padronizada por 100.000"
    
    ggplot(df, aes(x = year, group = 1)) +
      geom_ribbon(aes_string(ymin = aes_lo, ymax = aes_hi),
                  fill = "grey80", alpha = 0.4) +
      geom_line(aes_string(y = aes_y), size = 1) +
      geom_point(aes_string(y = aes_y), size = 2) +
      labs(
        title = paste(input$cause, input$sex, "(", input$population, ")"),
        x = "Ano",
        y = ylbl
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$rateTable <- renderTable({
    req(input$go_rates > 0)
    df <- rates_sel() %>%
      dplyr::filter(População == input$population)
    
    vals   <- if (input$rate_type == "crude") df$crude_rate  else df$dsr
    lowers <- if (input$rate_type == "crude") df$crude_lower else df$dsr_lower
    uppers <- if (input$rate_type == "crude") df$crude_upper else df$dsr_upper
    
    df %>%
      dplyr::transmute(
        Ano = year,
        `Taxa (Intervalo de Confiança 95%)` = glue::glue(
          "{round(vals, 2)} ({round(lowers, 2)}; {round(uppers, 2)})"
        )
      )
  }, sanitize.text.function = identity)
  
  # -------------------------
  # Forecast Explorer
  # -------------------------
  forecast_sel <- eventReactive(input$go_forecast, {
    shiny::withProgress(message = "A obter dados do INE...", value = 0, {
      incProgress(0.1)
      dat <- get_data_for_cached(input$area2, input$cause2)
      incProgress(0.4)
      
      df1 <- dat$full  %>% dplyr::filter(sex == input$sex2)
      df2 <- dat$trunc %>% dplyr::filter(sex == input$sex2)
      
      rates_fc <- dplyr::bind_rows(
        compute_metrics(df1) %>% dplyr::mutate(População = "Total"),
        compute_metrics(df2) %>% dplyr::mutate(População = "Menos de 75 anos")
      ) %>%
        dplyr::filter(População == input$population2)
      
      rate_col <- if (input$rate_type2 == "crude") "crude_rate" else "dsr"
      df_vals <- rates_fc %>%
        dplyr::filter(
          year >= input$years_fit[1],
          year <= input$years_fit[2]
        ) %>%
        dplyr::arrange(year) %>%
        dplyr::select(year, value = !!rlang::sym(rate_col)) %>%
        dplyr::mutate(year = as.integer(year))
      
      if (nrow(df_vals) < 3) {
        validate(
          need(FALSE, "Selecione pelo menos 3 anos para ajustar o modelo de projecção.")
        )
      }
      
      ts_data <- stats::ts(df_vals$value,
                           start = min(df_vals$year),
                           frequency = 1)
      
      h <- input$horizon
      
      fits <- lapply(input$models, function(m) switch(
        m,
        arima = forecast::auto.arima(ts_data),
        ets   = forecast::ets(ts_data),
        rwf   = forecast::rwf(ts_data, drift = TRUE),
        naive = forecast::naive(ts_data),
        theta = forecast::thetaf(ts_data),
        tbats = forecast::tbats(ts_data)
      ))
      names(fits) <- input$models
      
      fc_list <- lapply(fits, forecast::forecast, h = h)
      
      obs_df <- df_vals
      fc_df <- dplyr::bind_rows(lapply(names(fc_list), function(m) {
        fc <- fc_list[[m]]
        tibble(
          year  = seq(max(df_vals$year) + 1, by = 1, length.out = h),
          model = m,
          mean  = as.numeric(fc$mean),
          lower = as.numeric(fc$lower[, 2]),
          upper = as.numeric(fc$upper[, 2])
        )
      }))
      
      incProgress(0.5)
      list(obs = obs_df, fc = fc_df, fits = fits)
    })
  }, ignoreNULL = TRUE)
  
  output$forecastPlot <- renderPlot({
    req(input$go_forecast > 0)
    dat <- forecast_sel()
    obs <- dat$obs
    fc  <- dat$fc
    
    ggplot() +
      geom_line(data = obs, aes(x = year, y = value), size = 1) +
      geom_point(data = obs, aes(x = year, y = value), size = 2) +
      geom_ribbon(
        data = fc,
        aes(x = year, ymin = lower, ymax = upper, fill = model),
        alpha = 0.2
      ) +
      geom_line(
        data = fc,
        aes(x = year, y = mean, color = model),
        linetype = "dashed",
        size = 1
      ) +
      labs(
        title = "Horizonte Temporal (anos)",
        x = "Ano",
        y = if (input$rate_type2 == "crude") {
          "Taxa Bruta por 100.000"
        } else {
          "Taxa Padronizada por 100.000"
        }
      ) +
      scale_color_brewer(palette = "Set1", name = "Modelo") +
      scale_fill_brewer(palette = "Set1", name = "Modelo") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  output$downloadForecastPlot <- downloadHandler(
    filename = function() {
      paste0(
        "forecast_plot_",
        gsub(" ", "_", input$area2),
        "_",
        gsub(" ", "_", input$cause2),
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      req(input$go_forecast > 0)
      dat <- forecast_sel()
      obs <- dat$obs
      fc  <- dat$fc
      
      grDevices::png(file, width = 1200, height = 800, res = 150)
      print(
        ggplot() +
          geom_line(data = obs, aes(x = year, y = value), size = 1) +
          geom_point(data = obs, aes(x = year, y = value), size = 2) +
          geom_ribbon(
            data = fc,
            aes(x = year, ymin = lower, ymax = upper, fill = model),
            alpha = 0.2
          ) +
          geom_line(
            data = fc,
            aes(x = year, y = mean, color = model),
            linetype = "dashed",
            size = 1
          ) +
          labs(
            title = "Horizonte Temporal (anos)",
            x = "Ano",
            y = if (input$rate_type2 == "crude")
              "Taxa Bruta por 100.000"
            else
              "Taxa Padronizada por 100.000"
          ) +
          scale_color_brewer(palette = "Set1", name = "Modelo") +
          scale_fill_brewer(palette = "Set1", name = "Modelo") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5))
      )
      grDevices::dev.off()
    }
  )
  
  
  output$forecastTable <- renderTable({
    req(input$go_forecast > 0)
    dat <- forecast_sel()
    obs <- dat$obs
    fc  <- dat$fc
    
    # ensure year is integer
    obs <- obs %>% dplyr::mutate(year = as.integer(year))
    fc  <- fc  %>% dplyr::mutate(year = as.integer(year))
    
    obs_fmt <- obs %>%
      dplyr::transmute(
        Ano       = year,
        Observado = glue::glue("{round(value, 2)}")
      )
    
    fc_fmt <- fc %>%
      dplyr::mutate(
        texto = glue::glue(
          "{round(mean, 2)} ({round(lower, 2)}; {round(upper, 2)})"
        )
      ) %>%
      dplyr::select(year, model, texto)
    
    fc_wide <- fc_fmt %>%
      tidyr::pivot_wider(
        id_cols     = year,
        names_from  = model,
        values_from = texto
      )
    
    dplyr::full_join(obs_fmt, fc_wide, by = c("Ano" = "year")) %>%
      dplyr::arrange(Ano)
  }, sanitize.text.function = identity)
  
  output$downloadForecastCSV <- downloadHandler(
    filename = function() {
      paste0(
        "forecast_",
        gsub(" ", "_", input$area2),
        "_",
        gsub(" ", "_", input$cause2),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(file) {
      req(input$go_forecast > 0)
      dat <- forecast_sel()
      obs <- dat$obs
      fc  <- dat$fc
      
      obs <- obs %>% dplyr::mutate(year = as.integer(year))
      fc  <- fc  %>% dplyr::mutate(year = as.integer(year))
      
      obs_fmt <- obs %>%
        dplyr::transmute(
          Ano       = year,
          Observado = round(value, 2)
        )
      
      fc_fmt <- fc %>%
        dplyr::mutate(
          texto = glue::glue(
            "{round(mean, 2)} ({round(lower, 2)}; {round(upper, 2)})"
          )
        ) %>%
        dplyr::select(year, model, texto)
      
      fc_wide <- fc_fmt %>%
        tidyr::pivot_wider(
          id_cols     = year,
          names_from  = model,
          values_from = texto
        )
      
      out <- dplyr::full_join(obs_fmt, fc_wide, by = c("Ano" = "year")) %>%
        dplyr::arrange(Ano)
      
      utils::write.csv(out, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )
  
  
  output$accuracyTable <- renderTable({
    req(input$go_forecast > 0)
    fits <- forecast_sel()$fits
    metrics <- c("ME", "RMSE", "MAE", "MAPE", "MASE")
    
    acc_list <- lapply(fits, function(fit) {
      acc <- forecast::accuracy(fit)
      acc[1, metrics, drop = FALSE]
    })
    
    tbl <- dplyr::bind_rows(lapply(names(acc_list), function(m) {
      df <- as.data.frame(acc_list[[m]])
      df$Model <- m
      df
    }), .id = NULL)
    
    tbl %>%
      dplyr::select(Model, dplyr::everything())
  }, digits = 2)
  
  output$residPlot <- renderPlot({
    req(input$go_forecast > 0)
    fits <- forecast_sel()$fits
    res_df <- dplyr::bind_rows(lapply(names(fits), function(m) {
      res <- residuals(fits[[m]])
      tibble(
        year  = time(res),
        resid = as.numeric(res),
        model = m
      )
    }))
    
    ggplot(res_df, aes(x = year, y = resid, color = model)) +
      geom_line() +
      facet_wrap(~model, scales = "free_y") +
      labs(
        title = "Resíduos por Modelo",
        x = "Ano",
        y = "Resíduo"
      ) +
      theme_minimal()
  })
  
  output$residTable <- renderTable({
    req(input$go_forecast > 0)
    fits <- forecast_sel()$fits
    tbl <- dplyr::bind_rows(lapply(names(fits), function(m) {
      fit <- fits[[m]]
      res <- residuals(fit)
      lb  <- stats::Box.test(res, type = "Ljung")
      tibble(
        Model       = m,
        Mean        = round(mean(res), 2),
        SD          = round(sd(res), 2),
        LjungBox_p  = round(lb$p.value, 3)
      )
    }))
    tbl
  }, digits = 3)
  
  # -------------------------
  # ITS / Breakpoints
  # -------------------------
  its_sel <- eventReactive(input$go_breaks, {
    shiny::withProgress(message = "A obter dados do INE...", value = 0, {
      incProgress(0.1)
      dat <- get_data_for_cached(input$area3, input$cause3)
      incProgress(0.4)
      
      df1 <- dat$full  %>% dplyr::filter(sex == input$sex3)
      df2 <- dat$trunc %>% dplyr::filter(sex == input$sex3)
      
      df <- dplyr::bind_rows(
        compute_metrics(df1) %>% dplyr::mutate(População = "Total"),
        compute_metrics(df2) %>% dplyr::mutate(População = "Menos de 75 anos")
      ) %>%
        dplyr::filter(População == input$population3) %>%
        dplyr::arrange(year)
      
      incProgress(0.5)
      df
    })
  }, ignoreNULL = TRUE)
  
  output$breakPlot <- renderPlot({
    req(input$go_breaks > 0)
    df <- its_sel()
    y <- if (input$rate_type3 == "crude") df$crude_rate else df$dsr
    ylab <- if (input$rate_type3 == "crude") {
      "Taxa Bruta por 100.000"
    } else {
      "Taxa Padronizada por 100.000"
    }
    
    ts_y <- stats::ts(y, start = min(df$year), frequency = 1)
    bp_obj <- tryCatch(
      strucchange::breakpoints(ts_y ~ 1),
      error = function(e) NULL
    )
    
    plot(ts_y, ylab = ylab, xlab = "Ano", main = "Deteção de Quebras")
    lines(ts_y)
    if (!is.null(bp_obj) && length(bp_obj$breakpoints) > 0) {
      bp_years <- df$year[bp_obj$breakpoints]
      abline(v = bp_years, col = "red", lty = 2)
    }
  })
  
  output$breakTable <- renderTable({
    req(input$go_breaks > 0)
    df <- its_sel()
    y <- if (input$rate_type3 == "crude") df$crude_rate else df$dsr
    ts_y <- stats::ts(y, start = min(df$year), frequency = 1)
    
    bp_obj <- tryCatch(
      strucchange::breakpoints(ts_y ~ 1),
      error = function(e) NULL
    )
    
    if (is.null(bp_obj) || length(bp_obj$breakpoints) == 0) {
      tibble(Message = "Não foram detectadas quebras de série")
    } else {
      yrs <- df$year[bp_obj$breakpoints]
      tibble(Quebra = seq_along(yrs), Ano = yrs)
    }
  }, digits = 0)
}

# =========================================================
# Run app
# =========================================================

shinyApp(ui, server)
