# Extracted from test-regions.R:452

# prequel ----------------------------------------------------------------------
test_lookup <- function() {
  tibble(
    municipality = c("Beja", "Serpa", "Santarém", "Almeirim", "Braga", "Fafe", "Calheta (R.A.M.)"),
    municipality_code = c("1C1401", "1C1408", "1D1418", "1D1402", "111303", "111305", "3003"),
    nuts1 = c(rep("Continente", 6), "Região Autónoma da Madeira"),
    nuts2 = c("Alentejo", "Alentejo", "Oeste e Vale do Tejo", "Oeste e Vale do Tejo",
              "Norte", "Norte", "Região Autónoma da Madeira"),
    nuts3 = c("Baixo Alentejo", "Baixo Alentejo", "Lezíria do Tejo", "Lezíria do Tejo",
              "Cávado", "Ave", "Região Autónoma da Madeira")
  )
}
lookup_file <- function(vintage) sprintf("../../data/nuts_lookup_%s.rds", vintage)

# test -------------------------------------------------------------------------
skip_if_not(file.exists("../../data/uls_lookup.rds"), "ULS lookup not built")
old <- readRDS("../../data/nuts_lookup_2013.rds")
new <- readRDS("../../data/nuts_lookup_2024.rds")
members <- function(lk, level, units) sort(unique(lk$municipality[lk[[level]] %in% units]))
pairs <- list(
    "ULS Alto Minho" = "Alto Minho", "ULS Viseu Dão-Lafões" = "Viseu Dão Lafões",
    "ULS Litoral Alentejano" = "Alentejo Litoral", "ULS Baixo Alentejo" = "Baixo Alentejo",
    "ULS Alto Alentejo" = "Alto Alentejo", "ULS Alentejo Central" = "Alentejo Central"
  )
for (unit in names(pairs)) {
    expect_equal(sort(region_municipalities(unit, new)), members(new, "nuts3", pairs[[unit]]), info = unit)
    expect_equal(sort(region_municipalities(unit, new)), members(old, "nuts3", pairs[[unit]]), info = unit)
    expect_equal(regional_row_plan(unit, "2024", 1991:2024)$year, 1991:2024, info = unit)
  }
