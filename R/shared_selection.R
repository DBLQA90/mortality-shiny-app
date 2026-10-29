# =========================================================
# Shared selection across tabs
# =========================================================

# Each tab keeps its own input ids - the server and the smoke tests read them -
# but the user makes one selection: changing the place, cause or sex in any tab
# changes it in every tab linked to the same field. Which tabs share what is
# declared here and nowhere else, so a tab can later leave a field, or a field
# be split into two groups, by editing this list.
#
# `adapt` turns the shared value into one the linked input can hold: the
# planning tab takes a single place, and the tabs that take one cause or one sex
# keep the first of the observed tab's several. A `required` field never
# propagates an empty value, so clearing a field in one tab does not empty the
# others.
shared_link <- function(id, type = c("select", "radio", "slider", "text"), adapt = identity) {
  list(id = id, type = match.arg(type), adapt = adapt)
}

first_only <- function(value) utils::head(value, 1)

SHARED_SELECTION_FIELDS <- list(
  area = list(
    required = TRUE,
    links = list(
      shared_link("area"),
      shared_link("area2"),
      shared_link("beginner_area"),
      shared_link("avoidable_area"),
      shared_link("availability_area"),
      shared_link("planning_area", adapt = first_only)
    )
  ),
  area_label = list(
    required = FALSE,
    links = list(
      shared_link("area_label", "text"),
      shared_link("area_label2", "text"),
      shared_link("beginner_area_label", "text"),
      shared_link("avoidable_area_label", "text")
    )
  ),
  cause = list(
    required = TRUE,
    links = list(
      shared_link("cause"),
      shared_link("cause2", adapt = first_only),
      shared_link("beginner_cause", adapt = first_only),
      shared_link("availability_cause")
    )
  ),
  sex = list(
    required = TRUE,
    links = list(
      shared_link("sex"),
      shared_link("sex2", adapt = first_only),
      shared_link("beginner_sex", adapt = first_only),
      shared_link("avoidable_sex", adapt = first_only),
      shared_link("planning_cause_sex", "radio", adapt = first_only)
    )
  ),
  population = list(
    required = TRUE,
    links = list(
      shared_link("population", "radio"),
      shared_link("population2", "radio"),
      shared_link("beginner_population", "radio")
    )
  ),
  rate_type = list(
    required = TRUE,
    links = list(
      shared_link("rate_type2", "radio"),
      shared_link("beginner_rate_type", "radio")
    )
  ),
  # Only the sliders over the same range of rate years. The planning and
  # avoidable tabs choose years over other ranges, or a single year.
  years = list(
    required = TRUE,
    links = list(
      shared_link("years_import", "slider"),
      shared_link("years_fit", "slider"),
      shared_link("beginner_years_fit", "slider")
    )
  )
)

same_selection <- function(a, b) {
  identical(as.character(if (is.null(a)) character(0) else a), as.character(if (is.null(b)) character(0) else b))
}

update_linked_input <- function(session, link, value) {
  switch(
    link$type,
    select = updateSelectInput(session, link$id, selected = if (is.null(value)) character(0) else value),
    radio  = updateRadioButtons(session, link$id, selected = value),
    slider = updateSliderInput(session, link$id, value = value),
    text   = updateTextInput(session, link$id, value = if (is.null(value)) "" else value)
  )
}

# One field: the first link seeds the shared value, and every tab is aligned to
# it at start. An update sent to a tab comes back as an input change; that echo
# is recognised and dropped, so a quick second change by the user is not
# overwritten by the echo of the first.
sync_shared_field <- function(input, session, field) {
  links <- field$links
  shared <- reactiveVal(isolate(input[[links[[1]]$id]]))
  pending <- new.env(parent = emptyenv())

  push <- function(value, from = NULL) {
    for (link in links) {
      if (identical(link$id, from)) next
      target <- link$adapt(value)
      if (same_selection(isolate(input[[link$id]]), target)) next
      assign(link$id, target, envir = pending)
      update_linked_input(session, link, target)
    }
  }

  if (length(isolate(shared())) > 0) push(isolate(shared()))

  for (link in links) local({
    link <- link
    observeEvent(input[[link$id]], {
      value <- input[[link$id]]
      if (exists(link$id, envir = pending, inherits = FALSE)) {
        expected <- get(link$id, envir = pending)
        rm(list = link$id, envir = pending)
        if (same_selection(value, expected)) return()
      }
      if (isTRUE(field$required) && length(value) == 0) return()
      if (same_selection(value, link$adapt(shared()))) return()
      shared(value)
      push(value, from = link$id)
    }, ignoreInit = TRUE, ignoreNULL = FALSE)
  })

  invisible(shared)
}

register_shared_selection <- function(input, session, fields = SHARED_SELECTION_FIELDS) {
  invisible(lapply(fields, function(field) sync_shared_field(input, session, field)))
}

# The app-wide settings live in a side panel; this one line, always on screen,
# says which are in force, because six region names mean different things under
# the two vintages and a regional figure must not be read without it.
app_settings_summary <- function(vintage, region_source, population_basis, data_source) {
  region_source_short <- c(ine_rows = "óbitos regionais: linhas do INE", municipal_sum = "óbitos regionais: soma dos municípios")
  population_basis_short <- c(mean = "população média do ano", end = "população a 31 de dezembro")
  paste(
    c(
      paste0("Regiões NUTS ", vintage),
      unname(region_source_short[region_source]),
      unname(population_basis_short[population_basis]),
      names(data_source_choices)[match(data_source, data_source_choices)]
    ),
    collapse = " · "
  )
}
