# Extracted from test-regions.R:189

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
skip_if_not(file.exists(lookup_file("2013")), "2013 lookup not built")
old <- readRDS(lookup_file("2013"))
new <- readRDS(lookup_file("2024"))
expect_length(region_choices_for("2013", old), 8)
expect_length(region_choices_for("2024", new), 10)
expect_equal(region_choices_for("2024", new)[[1]], "Continente")
expect_equal(region_choices_for("2024", new)[[4]], "Norte")
areas_new <- area_choices_for("2024", new)
expect_equal(areas_new[[1]], "Portugal")
expect_equal(length(areas_new), 1 + 10 + 308)
