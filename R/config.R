# =========================================================
# Parameters
# =========================================================

population_indicator_current <- "0008273"
population_indicator_legacy  <- "0003182"
death_indicator_legacy       <- "0008206"
death_indicator_current      <- "0013166"

population_indicators <- c(population_indicator_current, population_indicator_legacy)
death_indicators <- c(death_indicator_legacy, death_indicator_current)

fallback_year_of_interest <- 1991:2023

fallback_local_area <- c(
  "Portugal",
  "Norte",
  "Alentejo",
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

fallback_diseases <- c(
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

annual_metric_choices <- c(
  "Óbitos" = "deaths",
  "Mortalidade Bruta" = "crude",
  "Mortalidade Padronizada" = "dsr",
  "Mortalidade Proporcional" = "proportional",
  "AVPP" = "ypll"
)

data_source_choices <- c(
  "INE em directo" = "ine",
  "Ficheiros RDS" = "snapshot"
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
