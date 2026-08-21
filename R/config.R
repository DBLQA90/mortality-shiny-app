# =========================================================
# Parameters
# =========================================================

# Population comes from three indicators, newest first.
#
# 0012918 is not simply a continuation of 0008273. It carries a revised
# estimate that runs progressively higher on the same years - +1.71% for 2021,
# +3.93% for 2022, +5.31% for 2023 on the national total - and both indicators
# were updated within two months of each other, so this is two published series
# rather than stale data on one side.
#
# It is therefore used for the whole of its range, not only for the years the
# older one lacks. Taking 2024 from the revised series while leaving 2023 on
# the old one would put a 5.3% step at the seam, which reads as a real fall in
# mortality and would be inherited by any forecast fitted across it. Using the
# revised series from 2021 moves the seam to 2020/2021, where the two differ by
# 1.7%. See METHODOLOGY.md.
population_indicator_revised <- "0012918"   # NUTS-2024, 2021-2025
population_indicator_current <- "0008273"   # NUTS-2013, 2011-2023
population_indicator_legacy  <- "0003182"   # NUTS-2002, 1991-2013
death_indicator_legacy       <- "0008206"
death_indicator_current      <- "0013166"

population_indicators <- c(
  population_indicator_revised,
  population_indicator_current,
  population_indicator_legacy
)
death_indicators <- c(death_indicator_legacy, death_indicator_current)

# Where the population series changes basis.
#
# The archive splices three indicators, and consecutive ones disagree about the
# years they both publish. Every such handover puts a step into any rate that
# crosses it, so each is declared here and warned about rather than left for the
# reader to discover as an unexplained break.
#
# `year` is the first year served by the new source; the step falls between
# `year - 1` and `year`.
POPULATION_SEAMS <- list(
  list(
    year = 2014L,
    # 0003182 (NUTS-2002, frozen 2014) hands over to 0008273 (NUTS-2013).
    # Measured on the years both publish, Portugal: the totals agree to 0.16%,
    # but the standardised rates differ by 2.49% (2011), 2.74% (2012) and 2.95%
    # (2013), because the two distribute the population across age bands
    # differently and standardisation is sensitive to exactly that.
    #
    # The consequence is a break the eye reads as real. Portugal's standardised
    # rate falls 6.70% from 2013 to 2014 as the archive shows it, against 3.86%
    # with one indicator on both sides - which is in line with the 3.00% of the
    # year before. strucchange puts a breakpoint at 2013|2014 in every fit
    # window from 1991 to 2002, and roughly half of it is this handover.
    effect = "sobretudo nas taxas padronizadas, cerca de 3%; as taxas brutas quase não mudam (0,2%)"
  ),
  list(
    year = 2021L,
    # 0008273 (NUTS-2013) hands over to 0012918 (NUTS-2024), which is a revised
    # estimate rather than a continuation: +1.71% for 2021, +3.93% for 2022,
    # +5.31% for 2023 on the national total.
    effect = "as taxas descem cerca de 1,7%, e a divergência cresce nos anos seguintes"
  )
)

# Which indicator serves a year the sources disagree about.
#
# get_source_year_plan() resolves overlaps by DESCENDING priority: the largest
# number wins, not the first listed. Stated here as a named vector so the
# ordering is visible and testable - written inline it reads as a list of
# indicators with 1, 2, 3 beside them and invites exactly the wrong reading.
#
#   0012918  2021-2025, revised. Wins every year it covers.
#   0003182  1991-2013, NUTS-2002. Keeps the 2011-2013 overlap, as it always
#            has: it is the contiguous source for 1991-2013 and the one the
#            Lisboa label repair was applied to.
#   0008273  2011-2023, NUTS-2013. Serves 2014-2020, between the other two.
population_source_priorities <- c(3L, 2L, 1L)
names(population_source_priorities) <- c(
  population_indicator_revised,
  population_indicator_legacy,
  population_indicator_current
)

# Seams the selected years cross. Only rates are affected: a count does not
# divide by this denominator.
population_seams_crossed <- function(years, metric_id = NULL) {
  if (!is.null(metric_id) && !(metric_id %in% metrics_requiring_population)) {
    return(list())
  }

  years <- as.integer(years)
  years <- years[is.finite(years)]

  if (length(years) == 0) {
    return(list())
  }

  Filter(
    function(seam) any(years < seam$year) && any(years >= seam$year),
    POPULATION_SEAMS
  )
}

# Warning for a rate series that crosses one or more population seams.
#
# The point is not that the data is wrong but that part of the movement at the
# seam is not mortality. A break detected there - by the structural-break
# analysis or by eye - is at least partly the denominator changing basis.
population_revision_warning <- function(years, metric_id = NULL) {
  seams <- population_seams_crossed(years, metric_id)

  if (length(seams) == 0) {
    return(NULL)
  }

  described <- vapply(
    seams,
    function(seam) {
      as.character(glue::glue("{seam$year - 1}/{seam$year} ({seam$effect})"))
    },
    character(1)
  )

  as.character(glue::glue(
    "Atenção: a série atravessa {length(seams)} ",
    "{if (length(seams) == 1) 'mudança' else 'mudanças'} de base da população ",
    "do INE, em {paste(described, collapse = ' e em ')}. Nesses pontos a taxa ",
    "muda em parte por mudança de denominador e não de mortalidade, pelo que ",
    "uma quebra detectada aí não deve ser lida como um acontecimento ",
    "epidemiológico."
  ))
}

# The two sides no longer end together, and now it is population that runs
# ahead: 0012918 publishes 2025 while cause-specific deaths stop at 2024.
# (2025 deaths exist, in 0013331 and 0013332, but without a cause dimension, so
# they cannot serve this app.) The analysable range is the union rather than the
# intersection, with rate metrics guarded for any year missing either side.
fallback_year_of_interest <- 1991:2024
fallback_population_years <- 1991:2025

# Metrics that divide by population, and so cannot be produced for a year that
# has deaths but no population estimate.
metrics_requiring_population <- c("crude", "dsr", "smr", "isr")

# Fallback region list, used only when no NUTS lookup has been built - the live
# vocabulary is derived from the lookup of the selected vintage. `Continente` is
# a NUTS I unit and the two autonomous regions carry the same name at NUTS I and
# NUTS II, so this list is the whole regional vocabulary of NUTS-2024, not only
# its second level. See R/regions.R.
fallback_nuts2_areas <- c(
  "Continente",
  "Norte",
  "Centro",
  "Oeste e Vale do Tejo",
  "Grande Lisboa",
  "Península de Setúbal",
  "Alentejo",
  "Algarve",
  "Região Autónoma dos Açores",
  "Região Autónoma da Madeira"
)

fallback_local_area <- c(
  "Portugal",
  fallback_nuts2_areas,
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
  "Calheta (R.A.A.)", "Calheta (R.A.M.)", "Câmara de Lobos", "Caminha", "Campo Maior", "Cantanhede",
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
  "Ílhavo", "Lagoa", "Lagoa (R.A.A.)", "Lagos", "Lajes das Flores", "Lajes do Pico", "Lamego",
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
  "Mortalidade Padronizada (directa, ESP 2013)" = "dsr",
  "SMR (padronização indirecta, referência = 100)" = "smr",
  "Taxa Padronizada Indirecta (por 100.000)" = "isr",
  "Mortalidade Proporcional" = "proportional",
  "AVPP" = "ypll",
  "Óbitos infantis (< 1 ano)" = "infant_deaths",
  "Mortalidade Infantil (por 1.000 nados-vivos)" = "infant"
)

# Infant mortality needs live births and under-1 deaths, neither of which comes
# from the main pipeline. It is therefore excluded from the population guard -
# it does not use population at all - and has its own coverage check.
#
# The count is offered alongside the rate because at municipal scale the rate is
# often not worth reading: a place with a handful of births produces a figure
# that swings by hundreds per 1,000 on a single death. The count says what
# actually happened and cannot mislead that way, so both are available and the
# reader chooses. The count is also available over a wider span: it needs only
# the under-1 death archive, so it covers 1991-1994 as well, where the rate
# cannot be computed because published births start in 1995.
infant_metric_id <- "infant"
infant_count_metric_id <- "infant_deaths"
infant_metric_ids <- c(infant_metric_id, infant_count_metric_id)

# Metrics reported as whole counts rather than as rates, so they are rounded to
# integers.
count_metric_ids <- c("deaths", "ypll", infant_count_metric_id)

# Counts that are divided by the length of a pooled window, so a 3-year window
# reads as a yearly average instead of a tripled total.
#
# Infant deaths are deliberately left out. Annualising them would put back the
# very problem the count exists to avoid: a municipality with two infant deaths
# in three years would read as "1", a rounded fraction, when what happened is
# two deaths. Accumulating events is the whole reason to pool a count this
# rare. Left as a window total it is also exactly the numerator of the pooled
# rate shown beside it, so the two metrics agree.
annualised_metric_ids <- c("deaths", "ypll")

# Metrics that compare an area against a reference schedule rather than against
# an external standard population. They need the reference area to be loaded
# alongside the selected area.
indirect_metric_ids <- c("smr", "isr")

# Reference geography for indirect standardisation. Portugal is the default so
# an SMR of 100 always means "national average"; the region is offered for
# users who want to know whether a municipality is unusual *within* its region.
smr_reference_choices <- c(
  "Portugal" = "Portugal",
  "Norte" = "Norte"
)

# How regional aggregates (Norte, Alentejo, ...) should be assembled.
#
# INE changed the NUTS boundaries in 2024: Lezíria do Tejo moved out of
# Alentejo, so the regional rows of the current death indicator (0013166,
# NUTS-2024) are not comparable with the historical indicator (0008206,
# NUTS-2013) or with the population indicators (NUTS-2013 / NUTS-2002).
# Municipality boundaries are unaffected, so the two honest options are to
# rebuild the region from its municipalities or to keep INE's own rows and show
# the break where it falls.
region_mode_choices <- c(
  "Aproximação por municípios (série contínua)" = "municipal",
  "Dados originais INE (quebra em 2022)" = "original"
)

# Regions are always built by summing their municipalities, and this is not a
# user-facing choice.
#
# Municipality boundaries have not moved: the archive holds the same 308
# municipalities in every year from 1991 to 2023. Regions are unions of whole
# municipalities, so one fixed membership list reconstructs any region for any
# year. Regions whose definition never changed are unaffected by this; regions
# that INE redrew - Alentejo losing Lezíria do Tejo in NUTS-2024 - come out as
# one continuous series instead of two incompatible halves.
#
# Presenting it as a toggle asked users to arbitrate a question about INE's
# geography vintages that they have no way to judge, and let a wrong answer
# through: reading INE's own rows across 2022 makes Alentejo look like it has
# below-average mortality. The cost is that regional totals no longer match
# INE's published regional rows, and fall short of the national total by the
# deaths INE cannot assign to a municipality (about 0.3-1.0%, see
# METHODOLOGY.md). That is documented rather than offered as an option.
#
# MORTALITY_REGION_MODE=original restores INE's own rows for anyone who needs
# to reproduce a published regional figure.
default_region_mode <- function() {
  requested <- Sys.getenv("MORTALITY_REGION_MODE", unset = "municipal")
  if (requested %in% unname(region_mode_choices)) requested else "municipal"
}

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
