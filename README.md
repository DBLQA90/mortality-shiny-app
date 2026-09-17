# Mortality Shiny App

Unofficial Shiny app for exploring Portuguese mortality indicators from INE.

The app supports observed mortality analysis, guided forecasting, advanced model comparison, diagnostics, and structural break exploration. Results are intended for exploration, decision support, and research workflows, and should be interpreted with appropriate epidemiological and statistical caution.

For a practical tab-by-tab guide to using and interpreting the app, see [USER_MANUAL.md](USER_MANUAL.md).

For calculation details, assumptions, and forecasting notes, see [METHODOLOGY.md](METHODOLOGY.md).

## Current Version Highlights

- Replaced `ineptR` with the CRAN package `ineptr2`.
- Uses INE indicators `0012918`, `0008273` and `0003182` for population. `0012918` is the NUTS-2024 series, 2021-2025, and is a **revision** rather than a continuation: it runs +1.71% above the older series for 2021, rising to +5.31% for 2023. It is used across its whole range so the single seam falls at 2020/2021 rather than a 5.3% cliff at 2023/2024, and the app warns when a rate series crosses it. Every rate from 2021 onward changed as a result — Portugal's 2023 crude mortality reads 1,056 per 100,000 against 1,112 before.
- Uses INE indicators `0008206` and `0013166` for deaths by cause.
- Detects available years and causes from INE metadata instead of hard-coding them.
- Derives the geography list from the NUTS lookup of the selected vintage, so the places offered are exactly the places the app can aggregate; the manual `local_area` vector remains only as a fallback.
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
- Adds infant mortality (deaths under 1 year per 1,000 live births), 1995-2024, from two datasets the main pipeline cannot supply: under-1 deaths, which the standard ingest folds into `0-4`, and live births, since no population indicator has an under-1 band. Reconciles with INE's published national series throughout. Offered both as a rate and as a plain death count, since at municipal scale the rate often is not one; rates on fewer than 1,000 births are marked with an asterisk rather than hidden.
- Corrects the AVPP weight of the `0 - 4 anos` band using those under-1 counts. Most deaths in the band are infants — 254 of 286 nationally in 2024 — and the band midpoint credited each with 67.5 lost years instead of 69.5.
- Adds a `Mortalidade Evitável` tab: deaths under 75 split into **preventable** (public health, primary prevention) and **treatable** (timely, effective care), following the 2019 Eurostat/OECD split. Only causes whose shortlist rubric maps without a clinical judgement are included — the six that cannot (ischaemic heart disease, cerebrovascular, stomach, lymphatic/haematopoietic, ovary, viral hepatitis) are shown as their own line marked `*`, not folded in or dropped, because they are ~18% of under-75 deaths. The figures are a lower bound on avoidable mortality, consistent for comparing places and years, and are not a reproduction of Eurostat's published figures.
- Builds every region by summing its municipalities, so the NUTS-2024 boundary change does not break the series. Applies to the observed, forecast and annual tabs alike.
- **Adds ULS and ARS as a geography in every tab**: 34 ULS, 5 ARS and two exact groups for the five ULS that share Lisboa, Loures or Porto at parish level. Built from the PNS2030 planning workbook's municipality mapping (`tools/build_uls_lookup.R`). ULS do not nest in NUTS, so they are selectable alongside NUTS regions rather than derived from them. Seven ULS and two ARS coincide with an INE subregion or region and use its complete rows; the rest are municipal sums.
- **Takes regional deaths from INE's own regional rows by default.** INE publishes cause-specific municipal deaths with incomplete age breakdowns, so summing municipalities under-counts - every year, worst in small regions (lung cancer 2013: Açores −18%, Madeira −12%, Alentejo −9%) and catastrophically in 2014 (−31% to −84%). Regional rows are complete. Redrawn regions are composed from NUTS III subregion rows in the years their own row does not exist, so every region is on INE rows in every year except Grande Lisboa and Península de Setúbal before 2022. Each composition equals the directly published row exactly in 2022, the year both indicators overlap. `Soma dos municípios` remains available as the second option in the page header; population is the municipal sum under both.
- Adds the NUTS I level: `Continente`, plus `R.A. dos Açores` and `R.A. da Madeira`, which the islands already carried at NUTS II. Built like every other region, from municipalities, under both vintages. The three partition the country exactly — 119,589 + 2,366 + 2,875 = 124,830 for 2021, which is INE's national row to the death.
- Repaired `0013166` for 2022 and 2023, which were missing three municipalities entirely (`Calheta (R.A.A.)`, `Calheta (R.A.M.)`, `Lagoa (R.A.A.)`) because they were fetched area-by-area against a stale area list. Açores read 7.0% low and Madeira 6.5% low for those years. Both now reconcile to the national row exactly.
- Makes the NUTS vintage a selection: **NUTS 2013** (7 regions, `Área Metropolitana de Lisboa`, Lezíria inside `Alentejo`) or **NUTS 2024** (9 regions, `Grande Lisboa` + `Península de Setúbal`, `Oeste e Vale do Tejo`). One selector in the page header, applying app-wide. Both cover the same 308 municipalities, so switching regroups the data rather than changing what is read, and either vintage gives a series continuous across 2022. Under NUTS 2013 the municipal sums reproduce INE's published regional rows exactly (`Norte` 37,121 and `Alentejo` 11,742 for 2021, seven regions summing to the national 124,830).

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

The runner loads only the network-free calculation modules, so the tests do not contact INE. They require `testthat`, `tidyverse`, `PHEindicatormethods`, and `forecast`.

Those cover the calculations in isolation. To exercise the Shiny server itself -
inputs set, event observers fired, outputs read - against the committed
snapshots, run the end-to-end smoke test:

```sh
Rscript tests/smoke_app.R
```

It needs no browser and makes no network calls. Wiring bugs live here rather
than in the calculation modules: it caught a pooled window reaching into a year
with no population estimate, and rate metrics being offered for years that have
none.

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
- infant mortality per 1,000 live births (1995-2024), asterisked where the period has fewer than 1,000 births
- infant deaths under 1 year, as a count (1991-2024), which is what to read at municipal scale
- directly standardised mortality (ESP 2013)
- SMR, indirectly standardised against a selectable reference (`Portugal` or `Norte`), expressed with the reference as 100
- indirectly standardised rate per 100,000
- proportional mortality, using all-cause deaths as the denominator for each location
- years of potential life lost before age 70

Each metric can be computed for a single year or pooled over a rolling 3- or
5-year window. Regions are always rebuilt from their municipalities; the header
control chooses which NUTS vintage groups them.

Annual tables show point estimates with 95% intervals where the interval can be estimated. A separate source table reports the population and death indicators used for each location/cause.

### Planning Indicators

Demographic and mortality indicators in the form the local health plans (PLS)
report them, one column per selected area (Portugal, Continente, NUTS regions,
ARS, ULS, municipalities): population structure, ageing and dependency indices,
births and crude birth rate, deaths and crude death rate, triennial infant
mortality, and proportional mortality by the 13 large cause groups. Each carries
its reference in the DRS/PNS2030 support workbook (I1, I3-I9, I13-I17, I28,
I32, I33, I35, I37-I42, I45, I64, I65). Socio-economic, birth and neonatal
components are fetched by `tools/fetch_planning_extra.R` into
`data/snapshots/planning_extra`.

The tab is built around one location: each indicator as a chart and a table,
with the areas containing the location (its ULS, ARS, NUTS III/II/I, Portugal)
as comparators for rates, shares and indices but not for counts; a summary of
every indicator; the ULS ranking; the pyramid; and proportional mortality. Two
Excel downloads close the tab: the location and its comparators, and every area
the app can build (the automated counterpart of the PLS workbook). Areas are ratios of sums over their municipalities;
Portugal and Continente use INE's published rows. Deaths come from the
municipal all-ages totals (`data/snapshots/death_totals`), which are complete,
rather than from the age breakdown, which is not. The engine is
`R/planning_indicators.R`.

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
- AVPP: years of potential life lost before age 70, approximated from age-band midpoints, with Dobson intervals for sparse counts. `0 - 4 anos` is split into `< 1 ano` (midpoint 0.5) and `1 - 4 anos` (midpoint 3) using the under-1 death counts, because most of the band's deaths are infants and the band midpoint understates their lost years

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
Rscript tools/refresh_snapshots.R task=all minutes=300 recent=2 note="INE 2026 release"
```

Tasks run one at a time: INE answers parallel requests with `429 Too Many
Requests` and then refuses connections for hours.

| Task | What it does |
|---|---|
| `deaths` | Deaths by cause and age (`0013166`): missing years, and the last `recent` years re-checked |
| `population` | Resident population (`0012918`): missing years and re-check |
| `deathtotals` | Municipal death totals by cause, all ages |
| `regional` | INE's regional death rows used for regions and some ULS |
| `infant` | Live births, under-1 deaths by cause, and complete under-1 counts |
| `planning` | RSI, pensions, purchasing power, waste, births by mother's age and gestation, under-1 deaths by age |
| `ambiguous` | Reports municipalities INE labels ambiguously, without guessing |
| `inventory` | Rebuilds the snapshot manifest |
| `fixareas` | Explicit only: re-runs the completed Lisboa/Calheta/Lagoa repair of old death chunks |
| `nuts2` | Explicit only. **Do not run**: five region names denote two or three NUTS levels at once, so those rows are double- or triple-counted; see [METHODOLOGY.md](METHODOLOGY.md) |

### Data versions

INE revises what it has published, so the same analysis can give different
figures a year apart. Every data file is written through `R/data_versions.R`:

- **`data/import_log.csv`** records each file written, with its import date, the
  tool and run that wrote it, whether it was new or replaced a different
  version, how many rows changed value, and Portugal's total before and after.
- **`data/archive/<run id>/`** keeps the previous version of every file a run
  replaced, under its original path. Re-fetching a year INE has not revised
  writes nothing.
- **`data/snapshots/REFRESH_STATUS.md`** ends with a table of what the run
  revised.
- Revisions made before the log existed (the Lisboa repair, the population
  revision, the births fix) were recorded from git by
  `tools/backfill_import_log.R`; their previous versions are restored from git.

To repeat an analysis on the data as it stood on a date:

```sh
Rscript tools/data_as_of.R date=2026-09-01
MORTALITY_SNAPSHOT_DIR=.mortality-shiny-cache/data_as_of/2026-09-01/snapshots Rscript -e 'shiny::runApp()'
```

The app shows the latest import date in its header, and the import history -
per dataset, per run, and the values each run revised - in the Data
Availability tab.

Note that population now runs *ahead* of deaths: `0012918` publishes 2025 while
cause-specific deaths stop at 2024. 2025 is therefore selectable for infant
mortality only, and refused for everything else with an explanatory message.

## Known Issues

- **INE blocks GitHub's hosted runners.** Requests from Azure IP ranges time out
  after ~135 s, while the same call answers a normal connection in under two
  seconds. The workflow checks this and fails fast; run
  `tools/refresh_snapshots.R` locally, or register a self-hosted runner.
- **2025 has population but no deaths by cause.** `0013166` ends at 2024, so
  2025 supports infant mortality only. 2025 deaths exist in `0013331` and
  `0013332` but carry no cause dimension. All rates work for 2024.
- **The population series changes basis at 2021**, by about 1.7%. See
  [METHODOLOGY.md](METHODOLOGY.md); the app warns when a series crosses it.
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
