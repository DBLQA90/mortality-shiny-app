# Tests for R/data_versions.R (import log, archive, reconstruction).

with_data_root <- function(code) {
  base <- tempfile("versions-")
  root <- file.path(base, "data")
  dir.create(file.path(root, "snapshots", "births"), recursive = TRUE)
  old <- Sys.getenv("DATA_RUN_ID", unset = NA)
  on.exit({
    if (is.na(old)) Sys.unsetenv("DATA_RUN_ID") else Sys.setenv(DATA_RUN_ID = old)
    unlink(base, recursive = TRUE)
  })
  code(root)
}

births_frame <- function(lisboa) {
  tibble::tibble(year = 2001L, area = c("Lisboa", "Porto", "Portugal"), births = c(lisboa, 2413, 112774), source_indicator = "0000003")
}

test_that("a new file is added and logged with its import date and total", {
  with_data_root(function(root) {
    Sys.setenv(DATA_RUN_ID = "2026-09-17T100000")
    path <- file.path(root, "snapshots", "births", "year_2001.rds")
    expect_equal(versioned_save_rds(births_frame(37208), path, tool = "fetch_births.R"), "added")

    log <- read_import_log(root)
    expect_equal(nrow(log), 1)
    expect_equal(log$relpath, "snapshots/births/year_2001.rds")
    expect_equal(log$dataset, "births")
    expect_equal(log$year, 2001L)
    expect_equal(log$action, "added")
    # The comparable total is Portugal's row.
    expect_equal(log$total_new, 112774)
    expect_equal(log$run_id, "2026-09-17T100000")
  })
})

test_that("identical content leaves the file, the log and the archive untouched", {
  with_data_root(function(root) {
    Sys.setenv(DATA_RUN_ID = "run1")
    path <- file.path(root, "snapshots", "births", "year_2001.rds")
    versioned_save_rds(births_frame(5604), path)
    mtime <- file.info(path)$mtime

    Sys.setenv(DATA_RUN_ID = "run2")
    # Same rows in another order.
    expect_equal(versioned_save_rds(births_frame(5604)[3:1, ], path), "unchanged")
    expect_equal(nrow(read_import_log(root)), 1)
    expect_false(dir.exists(file.path(root, "archive")))
    expect_equal(file.info(path)$mtime, mtime)
  })
})

test_that("a revised file archives the previous version and logs both totals", {
  with_data_root(function(root) {
    path <- file.path(root, "snapshots", "births", "year_2001.rds")
    Sys.setenv(DATA_RUN_ID = "run1")
    versioned_save_rds(births_frame(37208), path)
    Sys.setenv(DATA_RUN_ID = "run2")
    expect_equal(versioned_save_rds(births_frame(5604), path, note = "Lisboa by code"), "replaced")

    archived <- file.path(root, "archive", "run2", "snapshots", "births", "year_2001.rds")
    expect_true(file.exists(archived))
    expect_equal(readRDS(archived)$births[[1]], 37208)
    expect_equal(readRDS(path)$births[[1]], 5604)

    log <- read_import_log(root)
    replaced <- log[log$action == "replaced", ]
    expect_equal(replaced$archived_as, "archive/run2/snapshots/births/year_2001.rds")
    # Portugal is unchanged; the revision shows as one changed row.
    expect_equal(replaced$total_old, replaced$total_new)
    expect_equal(replaced$rows_changed, 1)

    changes <- run_changes(log, "run2")
    expect_equal(changes$dataset, "births")
    expect_equal(changes$change_pct, 0)
    expect_equal(changes$rows_changed, 1)

    runs <- import_runs_summary(log)
    expect_equal(runs$replaced[runs$run_id == "run2"], 1)
    expect_equal(runs$note[runs$run_id == "run2"], "Lisboa by code")
  })
})

test_that("the data can be rebuilt as it stood on an earlier date", {
  with_data_root(function(root) {
    path <- file.path(root, "snapshots", "births", "year_2001.rds")
    later <- file.path(root, "snapshots", "births", "year_2002.rds")
    versioned_save_rds(births_frame(37208), path)
    versioned_save_rds(births_frame(5604), path)
    versioned_save_rds(births_frame(5785), later)

    # Backdate the first import, as if it had happened years ago.
    log_path <- import_log_path(root)
    log <- utils::read.csv(log_path, colClasses = "character")
    log$imported_at[1] <- "2020-01-01T10:00:00+0000"
    log$imported_at[2:3] <- "2026-09-17T10:00:00+0000"
    utils::write.csv(log, log_path, row.names = FALSE)

    out <- tempfile("as-of-")
    result <- data_as_of("2021-06-30", root, out)
    expect_equal(result$restored, 1L)
    expect_equal(result$not_yet_imported, 1L)
    expect_equal(readRDS(file.path(out, "snapshots", "births", "year_2001.rds"))$births[[1]], 37208)
    expect_false(file.exists(file.path(out, "snapshots", "births", "year_2002.rds")))
    expect_false(file.exists(file.path(out, "archive")))

    # Today: everything as it is.
    now <- data_as_of(Sys.Date() + 1, root, tempfile("as-of-"))
    expect_equal(now$restored, 0L)
    expect_equal(readRDS(file.path(now$out_dir, "snapshots", "births", "year_2001.rds"))$births[[1]], 5604)
  })
})

test_that("writes outside a data directory are plain", {
  scratch <- tempfile("scratch-")
  path <- file.path(scratch, "x.rds")
  expect_equal(versioned_save_rds(1:3, path), "unversioned")
  expect_equal(readRDS(path), 1:3)
})

test_that("relative paths map to datasets and years", {
  expect_equal(relpath_dataset("snapshots/deaths/0013166/year_2024/cause_x.rds"), "deaths")
  expect_equal(relpath_dataset("snapshots/planning_extra/pensioners/year_2020.rds"), "planning_extra/pensioners")
  expect_equal(relpath_dataset("uls_lookup.rds"), "uls_lookup")
  expect_equal(relpath_year(c("snapshots/deaths/0013166/year_2024/cause_x.rds", "uls_lookup.rds")), c(2024L, NA))
})

test_that("content comparison ignores row order but not values", {
  a <- births_frame(5604)
  expect_true(snapshot_content_equal(a, a[c(2, 3, 1), c(4, 1, 3, 2)]))
  expect_false(snapshot_content_equal(a, births_frame(5605)))
  expect_false(snapshot_content_equal(a, a[1:2, ]))
})

test_that("totals compare Portugal, both sexes, all causes across versions", {
  pop <- tidyr::expand_grid(area = c("Portugal", "Norte", "Lisboa"), sex = c("H", "M", "HM"), age_band = c("0 - 4 anos", "5 - 9 anos")) %>%
    dplyr::mutate(pop = dplyr::case_when(area == "Portugal" ~ 100, TRUE ~ 40) * ifelse(sex == "HM", 2, 1))
  expect_equal(value_total(pop, "pop"), 400)
  # A version without regional rows or sex breakdown still gives the same total.
  expect_equal(value_total(dplyr::filter(pop, area == "Portugal", sex == "HM"), "pop"), 400)
  # No Portugal row: sum what there is.
  expect_equal(value_total(tibble::tibble(area = c("A", "B"), deaths = c(3, 4)), "deaths"), 7)
  extra <- tibble::tibble(area = "Portugal", category = c("Total", "15 - 19 anos"), value = c(100, 5))
  expect_equal(value_total(extra, "value"), 100)
  expect_equal(relpath_dataset("snapshots/ambiguous_areas.rds"), "ambiguous_areas")
})

test_that("changed rows are matched on keys, ignoring the source edition", {
  old <- tibble::tibble(area = c("A", "B", "C"), sex = "HM", deaths = c(1, 2, 3), source_indicator = "0008206")
  new <- tibble::tibble(area = c("A", "B", "D"), sex = "HM", deaths = c(1, 5, 4), source_indicator = "0013166")
  # B changed, C dropped, D added.
  expect_equal(rows_changed_count(old, new, "deaths"), 3)
  expect_equal(rows_changed_count(old, dplyr::mutate(old, source_indicator = "x"), "deaths"), 0)
})

test_that("a log written with fewer columns is upgraded before appending", {
  with_data_root(function(root) {
    old_columns <- setdiff(import_log_columns, "rows_changed")
    legacy <- as.data.frame(stats::setNames(as.list(rep("x", length(old_columns))), old_columns))
    legacy$year <- "2020"
    utils::write.csv(legacy, import_log_path(root), row.names = FALSE)
    versioned_save_rds(births_frame(1), file.path(root, "snapshots", "births", "year_2001.rds"))
    log <- read_import_log(root)
    expect_equal(nrow(log), 2)
    expect_equal(log$relpath[[2]], "snapshots/births/year_2001.rds")
    expect_equal(log$action[[2]], "added")
  })
})

test_that("exports carry the data version that produced them", {
  stamp <- export_stamp("2026-09-17", "2024")
  expect_match(stamp, "Dados do INE importados até 2026-09-17")
  expect_match(stamp, "regiões NUTS 2024")
  expect_match(stamp, format(Sys.Date(), "%Y-%m-%d"))
  expect_match(export_stamp(NA_character_), "desconhecida")

  # The CSV keeps its shape: the stamp is a comment after the data.
  path <- tempfile(fileext = ".csv")
  helpers_write_csv_utf8(tibble::tibble(area = c("A", "B"), value = c(1.5, 2.5)), path, stamp = stamp)
  lines <- readLines(path)
  expect_equal(length(lines), 4)
  expect_match(lines[[4]], "^# Dados do INE importados até")
  back <- utils::read.csv(path, comment.char = "#")
  expect_equal(nrow(back), 2)
  expect_equal(back$value, c(1.5, 2.5))

  # Without a stamp the file is unchanged.
  plain <- tempfile(fileext = ".csv")
  helpers_write_csv_utf8(tibble::tibble(a = 1), plain)
  expect_equal(length(readLines(plain)), 2)
})

test_that("a chart keeps its own caption and gains the stamp", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(tibble::tibble(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_equal(stamp_ggplot(p, "Dados de 2026-09-17")$labels$caption, "Dados de 2026-09-17")

  with_caption <- p + ggplot2::labs(caption = "Fonte: INE")
  stamped <- stamp_ggplot(with_caption, "Dados de 2026-09-17")
  expect_equal(stamped$labels$caption, "Fonte: INE\nDados de 2026-09-17")

  # Not a ggplot, or no stamp: returned untouched.
  expect_equal(stamp_ggplot(p, NULL)$labels$caption, NULL)
})
