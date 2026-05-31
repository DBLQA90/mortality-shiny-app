get_entry_dir <- function() {
  ofiles <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  ofiles <- ofiles[!is.na(ofiles)]

  if (length(ofiles) > 0) {
    return(dirname(normalizePath(utils::tail(ofiles, 1), mustWork = FALSE)))
  }

  getwd()
}

app_dir <- get_entry_dir()
source(file.path(app_dir, "mortality-shiny-app.R"), local = TRUE, chdir = TRUE)

shiny::shinyApp(ui = ui, server = server)
