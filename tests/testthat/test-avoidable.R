# Tests for R/avoidable.R (preventable and treatable mortality).

av_frame <- function(years = 2023L, areas = "A") {
  bands <- c("0 - 4 anos", "30 - 34 anos", "60 - 64 anos", "80 - 84 anos")
  causes <- c(
    "Todas as causas de morte",
    "Acidentes de transporte",                 # preventable
    "Diabetes mellitus",                       # treatable
    "Doenças isquémicas do coração",           # unresolved
    "Tumor (neoplasma) maligno do pâncreas"    # in no list
  )
  tidyr::expand_grid(year = years, area = areas, cause = causes, age_band = bands) %>%
    dplyr::mutate(
      sex = "HM",
      pop = 10000,
      deaths = dplyr::case_when(
        cause == "Todas as causas de morte" ~ 100,
        cause == "Acidentes de transporte" ~ 10,
        cause == "Diabetes mellitus" ~ 5,
        cause == "Doenças isquémicas do coração" ~ 20,
        TRUE ~ 3
      )
    )
}

test_that("the two lists are disjoint and leaf-level", {
  expect_true(avoidable_lists_are_disjoint())

  # Nothing that contains another shortlist rubric may appear, or a total would
  # count the same deaths twice.
  aggregates <- c(
    "Todas as causas de morte", "Tumores (neoplasmas)", "Tumores (neoplasmas) malignos",
    "Doenças do aparelho circulatório", "Doenças do aparelho respiratório",
    "Doenças do aparelho digestivo", "Causas externas de lesão e envenenamento",
    "Acidentes", "Algumas doenças infeciosas e parasitárias",
    "Transtornos mentais e comportamentais",
    "Malformações congénitas do aparelho circulatório",
    "Malformações congénitas do sistema nervoso"
  )
  expect_length(intersect(avoidable_group_causes("avoidable"), aggregates), 0)

  # The unresolved set is excluded from both lists, not silently folded in.
  expect_length(intersect(avoidable_group_causes("avoidable"), AVOIDABLE_UNRESOLVED$cause), 0)

  # The total group is exactly the union.
  expect_setequal(
    avoidable_group_causes("avoidable"),
    c(AVOIDABLE_PREVENTABLE, AVOIDABLE_TREATABLE)
  )
})

test_that("everything the tab needs is loaded in one pass", {
  required <- avoidable_required_causes()
  expect_true(all(avoidable_group_causes("avoidable") %in% required))
  expect_true(all(AVOIDABLE_UNRESOLVED$cause %in% required))
  expect_true("Todas as causas de morte" %in% required)
  expect_equal(anyDuplicated(required), 0L)
})

test_that("deaths above the cutoff never enter the calculation", {
  df <- av_frame()
  # Three bands under 75, one above: 3 x 10, not 4 x 10.
  expect_equal(avoidable_deaths(df, "Acidentes de transporte"), 30)
  expect_equal(avoidable_deaths(df, "Todas as causas de morte"), 300)
  expect_equal(avoidable_deaths(df, character(0)), 0)
})

test_that("the breakdown reconciles to all under-75 deaths", {
  b <- avoidable_breakdown(av_frame())
  get <- function(g) b$deaths[b$group == g]

  expect_equal(get("preventable"), 30)
  expect_equal(get("treatable"), 15)
  expect_equal(get("avoidable"), 45)
  expect_equal(get("unresolved"), 60)
  expect_equal(get("total"), 300)

  # The parts must add to the whole, or the tab would be showing a hole it
  # never accounts for.
  expect_equal(get("preventable") + get("treatable") + get("unresolved") + get("other"), get("total"))
  expect_equal(b$share[b$group == "total"], 100)
})

test_that("the standardised rate deduplicates population over causes", {
  df <- av_frame()

  # Population repeats once per cause in the loaded frame. If it were summed
  # blindly the denominator would be five times too large.
  one_cause <- avoidable_group_dsr(df, "Todas as causas de morte", 2023L)
  expect_true(is.finite(one_cause))

  # Same data, one extra cause loaded: the denominator must not move.
  wider <- dplyr::bind_rows(df, df %>% dplyr::mutate(cause = "Pneumonia", deaths = 0))
  expect_equal(avoidable_group_dsr(wider, "Todas as causas de morte", 2023L), one_cause)
})

test_that("a pooled window divides by person-years, not one year", {
  one <- av_frame(2023L)
  five <- av_frame(2019:2023)

  # Deaths and population both grow five-fold, so the rate must stay put.
  # Taking the population of a single year instead made this five times too
  # high for every pooled selection.
  expect_equal(
    avoidable_group_dsr(five, "Todas as causas de morte", 2021L),
    avoidable_group_dsr(one, "Todas as causas de morte", 2023L)
  )

  # Two areas of identical size behave the same way: one combined geography.
  two_areas <- av_frame(2023L, c("A", "B"))
  expect_equal(
    avoidable_group_dsr(two_areas, "Todas as causas de morte", 2023L),
    avoidable_group_dsr(one, "Todas as causas de morte", 2023L)
  )
})

test_that("disjoint groups add up on the standardised scale", {
  df <- av_frame()
  # A DSR is a weighted sum of band rates, so disjoint cause sets are additive.
  expect_equal(
    avoidable_group_dsr(df, avoidable_group_causes("avoidable"), 2023L),
    avoidable_group_dsr(df, AVOIDABLE_PREVENTABLE, 2023L) +
      avoidable_group_dsr(df, AVOIDABLE_TREATABLE, 2023L)
  )
})

test_that("the exclusions are stated, not hidden", {
  note <- avoidable_unresolved_note()
  expect_match(note, "^\\*")
  # Every excluded cause is named.
  for (cause in AVOIDABLE_UNRESOLVED$cause) {
    expect_true(grepl(substr(cause, 1, 25), note, fixed = TRUE))
  }
  expect_match(note, "limite inferior")
  expect_match(avoidable_scope_note(), "75 anos")
})
