# =========================================================
# Avoidable mortality: preventable and treatable
# =========================================================
# Deaths before age 75 that could have been avoided, split the way Eurostat and
# the OECD split them since their 2019 joint revision:
#
#   Prevenível  avoidable through public health and primary prevention -
#               tobacco, alcohol, road safety, suicide prevention.
#   Tratável    avoidable through timely and effective health care -
#               screening, diagnosis, treatment.
#
# The two answer different planning questions, which is the point of separating
# them: one is a public-health brief, the other a health-services brief.
#
# Under 75 is part of the definition, not an option. The app already computes
# under-75 rates, so the same age scope is reused.
#
# ---------------------------------------------------------------------------
# What this can and cannot claim
# ---------------------------------------------------------------------------
# Eurostat defines both lists at ICD-10 level. This archive holds INE's European
# shortlist ("lista sucinta europeia"), which is coarser, so the mapping is a
# correspondence between a fine list and a coarse one. Three consequences, all
# of which the tab states rather than hides:
#
#   1. Only causes whose shortlist rubric maps onto one Eurostat group without
#      a judgement call are included. Six that do not are listed in
#      AVOIDABLE_UNRESOLVED and deliberately left out - together they are about
#      18% of under-75 deaths, so leaving them out is not a rounding decision
#      and has to be reported alongside the totals.
#
#   2. The shortlist is nested: `Tumores (neoplasmas) malignos` contains the
#      individual sites, `Causas externas` contains transport, falls and
#      poisoning. Every list here is leaf-level only. Summing an aggregate
#      beside its own children would double-count on a large scale, so
#      avoidable_group_causes() is the only supported way to build a total.
#
#   3. Coverage is incomplete by construction. The 16 tumour sites the
#      shortlist names sum to less than its own malignant-neoplasm total,
#      because some ICD groups in the Eurostat lists have no shortlist rubric
#      at all. The figures here are therefore a lower bound on avoidable
#      mortality, not an estimate of it.
#
# The result is useful for comparing places and years on a consistent basis. It
# is not a reproduction of Eurostat's published figures for Portugal and must
# not be presented as one.

# Public health and primary prevention.
AVOIDABLE_PREVENTABLE <- c(
  "Doenças pelo vírus da imunodeficiência humana [HIV]",
  "Tumor (neoplasma) maligno do lábio, cavidade bucal e faringe",
  "Tumor (neoplasma) maligno do esófago",
  "Tumor (neoplasma) maligno do fígado e das vias biliares intra-hepáticas",
  "Tumor (neoplasma) maligno da laringe, da traqueia, dos brônquios e dos pulmões",
  "Melanoma maligno da pele",
  "Tumor (neoplasma) maligno da bexiga",
  "Transtornos mentais e comportamentais devidos ao uso de álcool",
  "Dependência de drogas (toxicomania)",
  "Doenças crónicas das vias aéreas inferiores",
  "Doenças crónicas do fígado",
  "Acidentes de transporte",
  "Quedas acidentais",
  "Envenenamento (intoxicação) acidental por drogas, medicamentos e substâncias biológicas",
  "Suicídios e lesões autoprovocadas voluntariamente",
  "Homicídios e lesões provocadas intencionalmente por outras pessoas"
)

# Timely and effective health care.
AVOIDABLE_TREATABLE <- c(
  "Tuberculose",
  "Infecção meningocócica",
  "Meningites (excepto 03 - Infecção meningocócica)",
  "Tumor (neoplasma) maligno do cólon",
  "Tumor (neoplasma) maligno (neoplasma) da junção retossigmoideia, reto, ânus e canal anal",
  "Tumor (neoplasma) maligno da mama",
  "Tumor (neoplasma) maligno do colo do útero",
  "Tumor (neoplasma) maligno do rim, excepto pelve renal",
  "Diabetes mellitus",
  "Pneumonia",
  "Gripe (Influenza)",
  "Asma e estado de mal asmático",
  "Úlcera gástrica, duodenal, péptica de localização não especificada e gastrojejunal",
  "Doenças do rim e ureter",
  "Complicações da gravidez, parto e puerpério",
  "Algumas afecções originadas no período perinatal",
  "Malformações congénitas, deformidades e anomalias cromossómicas"
)

# Left out on purpose, with the reason each cannot be assigned from the
# shortlist alone. Reported on the tab so the exclusion is visible: together
# these are roughly 18% of under-75 deaths, and the two cardiovascular rubrics
# are two thirds of that.
AVOIDABLE_UNRESOLVED <- tibble::tribble(
  ~cause, ~reason,
  "Doenças isquémicas do coração",
  "O Eurostat reparte esta causa em 50% prevenível e 50% tratável ao nível CID-10; a lista sucinta não permite separá-la.",
  "Doenças cérebro-vasculares",
  "Classificada como tratável no conjunto, mas parte é atribuível a prevenção primária (hipertensão, tabaco).",
  "Tumor (neoplasma) maligno do estômago",
  "Consta de algumas revisões da lista e está ausente noutras.",
  "Tumor (neoplasma) maligno do tecido linfático e hematopoético e tecidos relacionados",
  "Apenas o linfoma de Hodgkin e a leucemia infantil constam como tratáveis; a rubrica agrega todos os tumores do tecido linfático.",
  "Tumor (neoplasma) maligno do ovário",
  "Consta de algumas revisões da lista e está ausente noutras.",
  "Hepatite viral",
  "A hepatite B é prevenível por vacinação e a hepatite C é tratável; a rubrica não as separa."
)

# Age bands above the definition's cutoff. Avoidable mortality is an under-75
# measure, so these never enter the calculation.
avoidable_excluded_bands <- function() exclude_bands

avoidable_group_choices <- c(
  "Evitável (total)" = "avoidable",
  "Prevenível" = "preventable",
  "Tratável" = "treatable"
)

# Leaf causes of a group. `avoidable` is the union of the other two, never a
# shortlist aggregate.
avoidable_group_causes <- function(group = "avoidable") {
  switch(
    group,
    preventable = AVOIDABLE_PREVENTABLE,
    treatable = AVOIDABLE_TREATABLE,
    avoidable = unique(c(AVOIDABLE_PREVENTABLE, AVOIDABLE_TREATABLE)),
    character(0)
  )
}

# Every cause the tab needs to load: the two groups plus the unresolved set and
# the all-cause denominator, so one load serves the whole reconciliation.
avoidable_required_causes <- function() {
  unique(c(
    avoidable_group_causes("avoidable"),
    AVOIDABLE_UNRESOLVED$cause,
    "Todas as causas de morte"
  ))
}

# The two lists must stay disjoint: a cause in both would be counted twice in
# the total and would make the parts contradict the whole.
avoidable_lists_are_disjoint <- function() {
  length(intersect(AVOIDABLE_PREVENTABLE, AVOIDABLE_TREATABLE)) == 0
}

# Deaths under 75 for a set of causes, from an already-loaded frame.
avoidable_deaths <- function(df, causes) {
  if (nrow(df) == 0 || length(causes) == 0) {
    return(0)
  }

  rows <- df[df$cause %in% causes & !(as.character(df$age_band) %in% avoidable_excluded_bands()), , drop = FALSE]
  sum(rows$deaths, na.rm = TRUE)
}

# The reconciliation the tab shows: the parts of under-75 mortality, adding to
# the whole. `unresolved` is the deliberately excluded set and `other` is
# everything the lists do not reach - together they are the distance between
# what is reported and all under-75 deaths, which is the honest way to present
# a lower bound.
avoidable_breakdown <- function(df) {
  total <- avoidable_deaths(df, "Todas as causas de morte")
  preventable <- avoidable_deaths(df, AVOIDABLE_PREVENTABLE)
  treatable <- avoidable_deaths(df, AVOIDABLE_TREATABLE)
  unresolved <- avoidable_deaths(df, AVOIDABLE_UNRESOLVED$cause)

  tibble::tibble(
    group = c("preventable", "treatable", "avoidable", "unresolved", "other", "total"),
    label = c(
      "Prevenível",
      "Tratável",
      "Evitável (total)",
      "Não classificada (por resolver) *",
      "Não evitável / fora das listas",
      "Todas as causas, < 75 anos"
    ),
    deaths = c(
      preventable,
      treatable,
      preventable + treatable,
      unresolved,
      max(total - preventable - treatable - unresolved, 0),
      total
    )
  ) %>%
    # `total` is one number, not a column, so if_else() would try to recycle
    # the whole deaths vector into a length-1 condition.
    dplyr::mutate(
      share = if (total > 0) .data$deaths / total * 100 else NA_real_
    )
}

# Under-75 standardised rate for one group, from an already-loaded frame.
#
# The denominator needs care in two directions at once. Population repeats
# across the causes loaded - the same band appears once per cause - so summing
# it blindly would multiply the denominator by the number of causes. But it must
# still be summed across years and areas, because a pooled window divides by
# person-years and a multi-area selection is one combined geography.
#
# Taking it once per (year, area, band) and then summing does both: deduplicated
# over causes, accumulated over everything else. Using first() alone was wrong
# for any pooled window - Beja over five years came out five times too high.
avoidable_group_dsr <- function(df, causes, year = NA_integer_) {
  if (nrow(df) == 0 || length(causes) == 0) {
    return(NA_real_)
  }

  bands <- df[!(as.character(df$age_band) %in% avoidable_excluded_bands()), , drop = FALSE]

  pop_by_band <- bands %>%
    dplyr::distinct(.data$year, .data$area, .data$age_band, .keep_all = TRUE) %>%
    dplyr::group_by(.data$age_band) %>%
    dplyr::summarise(pop = sum(.data$pop, na.rm = TRUE), .groups = "drop")

  deaths_by_band <- bands %>%
    dplyr::filter(.data$cause %in% causes) %>%
    dplyr::group_by(.data$age_band) %>%
    dplyr::summarise(deaths = sum(.data$deaths, na.rm = TRUE), .groups = "drop")

  collapsed <- dplyr::left_join(pop_by_band, deaths_by_band, by = "age_band") %>%
    dplyr::mutate(deaths = dplyr::coalesce(.data$deaths, 0))

  if (nrow(collapsed) == 0 || sum(collapsed$pop, na.rm = TRUE) <= 0) {
    return(NA_real_)
  }

  metric <- compute_metrics(
    collapsed %>%
      dplyr::mutate(year = as.integer(year), sex = "HM", cause = "grupo")
  )

  as.numeric(metric$dsr[[1]])
}

avoidable_unresolved_note <- function() {
  as.character(glue::glue(
    "* Seis causas não podem ser atribuídas a um dos grupos a partir da lista ",
    "sucinta europeia e ficam de fora dos totais: ",
    "{paste(AVOIDABLE_UNRESOLVED$cause, collapse = '; ')}. ",
    "Representam cerca de 18% dos óbitos com menos de 75 anos, pelo que os ",
    "valores apresentados são um limite inferior da mortalidade evitável e não ",
    "uma estimativa dela."
  ))
}

avoidable_scope_note <- function() {
  as.character(glue::glue(
    "Mortalidade evitável é definida para óbitos antes dos 75 anos. As listas ",
    "seguem a revisão conjunta Eurostat/OCDE de 2019, adaptadas à lista sucinta ",
    "europeia que o INE publica, que é menos detalhada: só entram as causas ",
    "cuja correspondência não exige juízo clínico. Os valores servem para ",
    "comparar locais e anos numa base consistente, não para reproduzir os ",
    "números publicados pelo Eurostat."
  ))
}
