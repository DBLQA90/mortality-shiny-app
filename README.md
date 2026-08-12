# Mortality Shiny App

Unofficial Shiny app for exploring Portuguese mortality indicators from INE.

The app supports observed mortality analysis, guided forecasting, advanced model comparison, diagnostics, and structural break exploration. Results are intended for exploration, decision support, and research workflows, and should be interpreted with appropriate epidemiological and statistical caution.

For a practical tab-by-tab guide to using and interpreting the app, see [USER_MANUAL.md](USER_MANUAL.md).

For calculation details, assumptions, and forecasting notes, see [METHODOLOGY.md](METHODOLOGY.md).

## Current Version Highlights

- Replaced `ineptR` with the CRAN package `ineptr2`.
- Uses INE indicators `0008273` and `0003182` for population.
- Uses INE indicators `0008206` and `0013166` for deaths by cause.
- Detects available years and causes from INE metadata instead of hard-coding them.
- Keeps the app geography list based on the original manual `local_area` vector, with `Norte` also available.
- Lets users select the year range to import from the years available in the source indicators.
- Uses year-range sliders for the observed, guided forecast, and advanced forecast windows.
- Adds a one-year annual metrics comparison tab for Portugal, Norte, and one selected local aggregation, with multiple causes sorted by the selected local value.
- Adds CSV and PNG exports for the app's tables and plots.
- Requests only the years needed from each source indicator.
- Adds persistent local caching for INE metadata and data queries.
- Downloads data in small year/area/cause slices so interrupted or failed runs can reuse data already fetched.
- Prioritises year loading based on the latest slider movement, so leftward changes load recent-to-older and rightward changes load older-to-newer.
- Uses a large INE client timeout for long indicator calls.
- Adds an optional RDS snapshot data source so users can load prebuilt data files instead of querying INE live.
- Splits the app internals into smaller `R/` files for configuration, INE access, snapshot access, metrics, data assembly, and UI helpers.
- Adds a snapshot inventory manifest so the app can discover available chunked RDS files before loading data.
- Adds a `Disponibilidade de Dados` tab to inspect RDS coverage by year, area, cause, and source indicator.
- Shows source indicators used in observed and annual analyses.
- Adds 95% uncertainty intervals to annual metric tables and plots where estimable.
- Adds a base-R dependency installer and first-run package bootstrap.
- Selects the recommended forecast model from out-of-sample accuracy (rolling validation or a single train/test split, with a user-set test-set percentage), instead of in-sample fit.
- Renders the forecast and observed-rate charts as interactive `plotly` widgets: hover a point to read the year and value, plus zoom and pan.
- Adds indirect standardisation: `SMR` (reference = 100) and an indirectly standardised rate per 100,000, with Byar/exact-Poisson intervals and a significance flag against the reference. This is the metric to use for small municipalities, where direct standardisation is unstable or unestimable.
- Adds selectable 3- and 5-year pooling for the annual comparison, summing deaths and person-years so sparse local series become readable without averaging annual rates.
- Builds every region by summing its municipalities, so the NUTS-2024 boundary change does not break the series. Applies to the observed, forecast and annual tabs alike. Regional figures therefore differ slightly from INE's published regional rows; see [METHODOLOGY.md](METHODOLOGY.md).

## Running The App

From a fresh R installation, install/check the runtime packages once:

```r
source("install_dependencies.R")
```

Then run the app from this repository:

```r
shiny::runApp(".")
```

In RStudio, opening this folder and pressing **Run App** works through the repository's `app.R` launcher.

The app startup also checks for missing runtime packages and installs them from CRAN on first run when needed. Set `MORTALITY_INSTALL_MISSING_PACKAGES=false` before launching if you prefer the app to stop and report missing packages instead of installing them automatically.

The required runtime packages are:

- `glue`
- `PHEindicatormethods`
- `tidyverse`
- `shiny`
- `plotly`
- `forecast`
- `ineptr2`
- `strucchange`
- `memoise`
- `cachem`
- `later`

## Reproducible Environment (Optional)

A pinned set of package versions is recorded in `renv.lock` (a known-good snapshot the app and tests were verified against). This is optional: the app still runs against your system library via the bootstrap installer above, and the lockfile does not auto-activate `renv`.

To reproduce the pinned environment instead, install [`renv`](https://rstudio.github.io/renv/) and restore into a project-local library:

```r
install.packages("renv")
renv::restore()
```

`renv.lock` records the recursive dependency closure from CRAN, so restoring installs the exact pinned versions. Regenerate it after changing dependencies with:

```r
renv::snapshot(packages = required_packages)
```

The recorded R version is the one the snapshot was taken with; `renv::restore()` warns but proceeds on a different R version.

## Code Layout

The app entry point is `mortality-shiny-app.R`. Most helper logic is split into smaller files under `R/`:

- `R/dependencies.R`: runtime package list and first-run installation helper.
- `R/config.R`: indicator IDs, default choices, age groups, standard population weights.
- `R/cache.R`: persistent metadata/data cache helpers.
- `R/ine_client.R`: INE metadata and live data download helpers.
- `R/snapshots.R`: flat/chunked RDS readers, source priority handling, and snapshot inventory helpers.
- `R/metrics.R`: mortality-rate, direct-standardisation, and AVPP calculations.
- `R/standardisation.R`: indirect standardisation (SMR/ISR) and multi-year pooling (unit tested).
- `R/regions.R`: NUTS vintage handling and municipal rebuilds of regional aggregates (unit tested).
- `R/forecast_helpers.R`: pure forecast-metric and out-of-sample validation helpers (unit tested).
- `R/data_access.R`: shared data assembly for snapshot and live INE sources.
- `R/ui_helpers.R`: reusable Shiny UI panels and tabs.

## Running Tests

Unit tests for the calculation and forecast-helper modules live in `tests/testthat/`. Run them with:

```sh
Rscript run_tests.R
```

The runner loads only the network-free calculation modules (`R/config.R`, `R/metrics.R`, `R/forecast_helpers.R`), so the tests do not contact INE. They require `testthat`, `tidyverse`, `PHEindicatormethods`, and `forecast`.

## Optional RDS Snapshots

The app can load data from prebuilt RDS files instead of querying INE live. In the app controls, choose `Ficheiros RDS` under `Fonte de dados`; choose `INE em directo` to query INE instead. This selector is available in the observed mortality, guided forecast, annual metrics, and advanced model specification loading controls.

By default, the app first supports flat snapshot files:

- `data/snapshots/population.rds`
- `data/snapshots/deaths.rds`

For larger datasets, the repository uses chunked files:

- `data/snapshots/population/year_<year>.rds`
- `data/snapshots/deaths/<indicator>/year_<year>/cause_<cause-token>.rds`

By default, the app first uses local snapshot files under `data/snapshots` when both population and death snapshots are present. This means a full repository download works from its own local RDS files without reading the same files back from GitHub. If local snapshot files are absent or only partially present, the app falls back to the configured GitHub raw snapshot directory. The manifest lists available chunks, areas, years, causes, row counts, and source priorities, allowing the app to choose relevant RDS files before reading the data itself. After adding or changing local snapshot chunks, refresh the manifest with:

```sh
Rscript tools/update_snapshot_inventory.R
```

You can also point to another location with environment variables:

- `MORTALITY_SNAPSHOT_DIR`
- `MORTALITY_USE_LOCAL_SNAPSHOTS=false` to skip local `data/snapshots` and use the remote snapshot directory when no explicit snapshot directory is set
- `MORTALITY_POPULATION_SNAPSHOT_RDS`
- `MORTALITY_DEATHS_SNAPSHOT_RDS`
- `MORTALITY_SNAPSHOT_RDS` for one combined RDS list containing `population` and `deaths`
- `MORTALITY_SNAPSHOT_INVENTORY_RDS` for a custom inventory manifest path
- `MORTALITY_DEFAULT_DATA_SOURCE=snapshot` if you want the app to open with `Ficheiros RDS` selected by default

To build separate snapshot files from INE:

```sh
Rscript tools/build_ine_snapshot.R out=data/snapshots years=2022:2023 areas=Portugal\|Norte causes="Todas as causas de morte|Diabetes mellitus"
```

Use `ALL` for years, areas, or causes when you intentionally want a broad snapshot. Large snapshots can take a long time to build and may be too large for normal GitHub commits, so consider GitHub Releases or another file host for production-size files.

For the slow historical deaths indicator `0008206`, the repository includes a resumable portal exporter. It uses INE's own web table export, then reshapes the CSV into the same chunked RDS layout used by the app:

```sh
Rscript tools/build_0008206_snapshot_from_portal.R out=data/snapshots years=ALL areas=ALL causes=ALL max_batches=26
```

This writes small files under `data/snapshots/deaths/0008206/`. It is intended for manual maintenance if the local snapshot archive needs to be rebuilt. The repository no longer runs scheduled GitHub Actions jobs for this backfill.

Population and API-backed death indicators can also be built directly into the same chunked layout:

```sh
Rscript tools/build_population_snapshot_chunks.R years=2019:2023 areas=Portugal\|Norte out=data/snapshots
Rscript tools/build_death_snapshot_chunks.R indicator=0013166 years=2022:2023 areas=Portugal\|Norte causes=ALL out=data/snapshots
```

Chunked death files are stored per indicator:

```text
data/snapshots/deaths/<indicator>/year_<year>/cause_<cause-token>.rds
```

When overlapping years exist, the app prefers the current death indicator `0013166` over the historical `0008206` snapshot for the same year, cause, area, sex, and age band. If a higher-priority indicator lacks a requested area, lower-priority chunks can still fill those rows.

The portal exporter fetches a table from the INE portal, requests CSV, parses the returned file, combines `Menos de 1 ano` and `1 - 4 anos` into `0 - 4 anos`, and writes one RDS file per year and cause. Defaults are conservative: latest available `0008206` year, `Portugal|Norte`, and all causes. Use `areas=ALL`, `years=2019:2022`, `area_batch_size=12`, or `max_batches=1` to control how much work is done per run.

The repository currently includes complete population chunks for the configured app locations, complete `0013166` chunks for 2022-2023 where INE returns location data, and complete `0008206` chunks for 1991-2022 across the configured locations.

## App Modules

### Observed Mortality

Explore historical mortality rates by geography, cause of death, sex, and population scope.

Outputs include:

- mortality rates per 100,000
- Poisson 95% confidence intervals for crude rates
- directly standardised rates using ESP 2013
- time-series plots
- summary and annual data tables
- source indicators used in the loaded data

### Guided Forecasting

Provides a guided forecasting workflow with simpler controls and reasonable defaults.

The user can choose:

- residence location and optional selection name
- cause of death
- sex
- population scope
- rate
- data source
- forecast horizon
- training window
- how the recommended model is chosen: rolling validation (default) or a single train/test split, with a test-set size in percent of the selected years (falls back to in-sample fit when the series is too short)
- recommended model or model comparison mode

Guided and advanced forecasts can project up to 30 years beyond the last observed year.

### Advanced Forecasting

Provides a fuller modelling workflow, including:

- ARIMA, ETS, random walk with drift, naive, Theta, TBATS, Holt, and damped Holt models
- custom training windows
- confidence interval controls
- optional log transform
- forecast tables and downloadable outputs
- residual diagnostics
- backtesting and model comparison
- structural break analysis

### Annual Metrics

Compares one selected year across Portugal, Norte, and a selected local aggregation. The tab shows one metric at a time for one or more causes of death, ordered from highest to lowest by the selected local aggregation.

Available metrics:

- deaths
- crude mortality
- directly standardised mortality (ESP 2013)
- SMR, indirectly standardised against a selectable reference (`Portugal` or `Norte`), expressed with the reference as 100
- indirectly standardised rate per 100,000
- proportional mortality, using all-cause deaths as the denominator for each location
- years of potential life lost before age 70

Each metric can be computed for a single year or pooled over a rolling 3- or
5-year window. The region control chooses whether `Norte` and `Alentejo` are
rebuilt from their municipalities or read from INE's own rows.

Annual tables show point estimates with 95% intervals where the interval can be estimated. A separate source table reports the population and death indicators used for each location/cause.

### Data Availability

Inspect the local RDS snapshot inventory before loading an analysis. The tab summarises available population and death chunks, then checks selected years, areas, and causes against `data/snapshots/snapshot_inventory.rds`.

Coverage states:

- `Disponível`: all selected areas are present for that year/cause.
- `Parcial`: at least one selected area is present, but others are missing.
- `Indisponível`: none of the selected areas are present.

The coverage table can be exported as CSV.

## Data Sources

Data can be read from the prebuilt RDS snapshots in this repository or fetched live from INE through `ineptr2`.

Population indicators:

- `0008273`
- `0003182`

Deaths by cause indicators:

- `0008206`
- `0013166`

The app harmonises age bands, recodes infant mortality into the `0-4` age group, excludes total or ignored age categories where needed, and can compute rates for the full population or the population under 75 years.

## Methods Summary

The app sums deaths and population over the selected geography before calculating rates, so multiple selected areas are interpreted as one combined area.

Main calculations:

- crude mortality: deaths divided by population, multiplied by 100,000
- crude confidence intervals: exact Poisson intervals scaled to the selected population
- standardised mortality: direct standardisation with European Standard Population 2013 weights; the under-75 scope is the conventional premature-mortality rate standardised to the ESP 0-74 sub-population and is labelled as such
- proportional mortality: selected-cause deaths divided by all-cause deaths for the same year, sex, and geography
- AVPP: years of potential life lost before age 70, approximated from age-band midpoints, with Dobson intervals for sparse counts

Forecasts are exploratory extrapolations of annual mortality-rate series using models from the `forecast` package. Model comparison uses common forecast accuracy metrics such as RMSE, MAE, MAPE, and MASE.

## Caching And Performance

The first request for a new area, cause, and year range can still take time, especially when it needs historical deaths from indicator `0008206`.

To reduce repeated delays, the app uses:

- in-memory caching during the Shiny session
- persistent RDS files on disk
- separate metadata and data cache expiry windows
- granular data slices so partial downloads survive interruptions
- pending Shiny event servicing between slices, so cancellation is checked before the next INE request

By default, cache files are written to `.mortality-shiny-cache` next to the app file.

Optional environment variables:

- `MORTALITY_APP_CACHE_DIR`: custom cache directory
- `MORTALITY_METADATA_CACHE_MAX_AGE`: metadata cache age in seconds, default 24 hours
- `MORTALITY_DATA_CACHE_MAX_AGE`: data cache age in seconds, default 7 days

If an INE request fails but a stale cached file exists, the app will use the stale file and show a warning.

## Snapshot Maintenance

`.github/workflows/refresh-snapshots.yml` runs the snapshot maintenance tasks on
a GitHub runner, weekly and on demand. The work is almost entirely spent waiting
on INE, and the driver is resumable: each run makes what progress fits in its
time budget, commits it, and the next run continues.

```sh
Rscript tools/refresh_snapshots.R task=all minutes=300
```

| Task | What it does |
|---|---|
| `fixareas` | Repairs geographies stored under an ambiguous INE label (`Lisboa`, `Calheta`, `Lagoa`) by refetching them by category code |
| `deaths2024` | Asks INE which death years exist and fetches those the archive lacks |
| `nuts2` | Backfills regional rows, for the `Dados originais INE` region mode only |
| `ambiguous` | Reports municipalities INE labels ambiguously, without guessing |
| `inventory` | Rebuilds the snapshot manifest |

Note that INE publishes no municipality-level population by age and sex beyond
2023, so 2024 supports counts, proportional mortality and AVPP but not rates.

## Known Issues

- **INE blocks GitHub's hosted runners.** Requests from Azure IP ranges time out
  after ~135 s, while the same call answers a normal connection in under two
  seconds. The workflow checks this and fails fast; run
  `tools/refresh_snapshots.R` locally, or register a self-hosted runner.
- **2024 has deaths but no population.** INE publishes no municipality-level
  population by age and sex beyond 2023, so 2024 supports deaths, proportional
  mortality and AVPP, but crude, standardised and SMR values are refused for
  that year with an explanatory message.
- Regional totals are municipal sums and do not match INE's published regional
  figures, nor do they add up to the national total: the national row includes
  deaths INE cannot assign to a municipality, about 0.3-1.0% depending on year.
  `MORTALITY_REGION_MODE=original` restores INE's own rows, but historical
  chunks carry them only for `Norte` and `Alentejo`.
- Selecting overlapping areas (for example a region and one of its own
  municipalities) sums them and double-counts. The app warns but does not block.

The `Lisboa`, `Calheta` and `Lagoa` label defects described in
[METHODOLOGY.md](METHODOLOGY.md) were repaired in the committed archive: the
2013/2014 seam now shows no area moving more than 20%, and the municipal sum
for 2000 matches the national row exactly.

## Limitations

- INE API calls can be slow or temporarily unavailable.
- Indicator `0008206` is particularly slow for some historical mortality requests.
- Small municipalities and rare causes can produce sparse counts and unstable rates.
- Forecasts are sensitive to short time series, low counts, and structural breaks.
- Historical INE data revisions are not versioned inside the app.

This is a non-official tool.
