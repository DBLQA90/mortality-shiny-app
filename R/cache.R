# =========================================================
# Persistent cache
# =========================================================

persistent_cache_root <- Sys.getenv(
  "MORTALITY_APP_CACHE_DIR",
  file.path(get_app_dir(), ".mortality-shiny-cache")
)
persistent_metadata_cache_max_age <- as.numeric(Sys.getenv(
  "MORTALITY_METADATA_CACHE_MAX_AGE",
  24 * 60 * 60
))
persistent_data_cache_max_age <- as.numeric(Sys.getenv(
  "MORTALITY_DATA_CACHE_MAX_AGE",
  7 * 24 * 60 * 60
))

cache_file_token <- function(key) {
  key_file <- tempfile(fileext = ".rds")
  on.exit(unlink(key_file), add = TRUE)
  saveRDS(key, key_file, version = 2)
  unname(tools::md5sum(key_file))
}

ensure_cache_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

read_persistent_cache <- function(path, max_age = NULL) {
  if (!file.exists(path)) {
    return(NULL)
  }

  if (!is.null(max_age) && is.finite(max_age)) {
    age <- as.numeric(difftime(Sys.time(), file.info(path)$mtime, units = "secs"))
    if (is.na(age) || age > max_age) {
      return(NULL)
    }
  }

  tryCatch(readRDS(path), error = function(e) NULL)
}

write_persistent_cache <- function(value, path) {
  ensure_cache_dir(dirname(path))
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(value, tmp, version = 2)
  file.rename(tmp, path)
  value
}

with_persistent_cache <- function(path, max_age, label, expr) {
  cached <- read_persistent_cache(path, max_age = max_age)
  if (!is.null(cached)) {
    return(cached)
  }

  tryCatch(
    write_persistent_cache(force(expr), path),
    error = function(e) {
      stale <- read_persistent_cache(path, max_age = NULL)
      if (!is.null(stale)) {
        warning(
          glue::glue("Using stale cached {label} because INE could not be reached: {conditionMessage(e)}"),
          call. = FALSE
        )
        return(stale)
      }
      stop(e)
    }
  )
}

metadata_cache_file <- function(indicator) {
  file.path(
    persistent_cache_root,
    "metadata",
    paste0(cache_file_token(list(lang = ine_client$lang, indicator = indicator)), ".rds")
  )
}

data_cache_file <- function(indicator, cats, has_cause) {
  file.path(
    persistent_cache_root,
    "data",
    indicator,
    paste0(cache_file_token(list(
      lang = ine_client$lang,
      indicator = indicator,
      cats = cats,
      has_cause = has_cause
    )), ".rds")
  )
}
