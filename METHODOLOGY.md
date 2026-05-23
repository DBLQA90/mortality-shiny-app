# Methodology

This document describes how the app prepares INE data, computes mortality metrics, and fits exploratory forecasts. The app is non-official and intended for analysis support, not as a substitute for validated statistical production workflows.

## Data Sources

The app fetches data from INE through `ineptr2`.

Population indicators:

- `0008273`
- `0003182`

Deaths by cause indicators:

- `0008206`
- `0013166`

Available years, areas, and causes of death are read from INE metadata where possible. The app keeps the original manual area list for user-facing location choices, with `Norte` added explicitly.

When a requested year can be found in more than one indicator, the app builds a non-overlapping source-year plan and de-duplicates loaded rows by year, area, sex, cause, and age band. Data are downloaded in small year, area, and cause slices so partial downloads can be cached and reused.

## Geography

Users can select one or more local areas. When more than one local area is selected, the app treats the selection as one combined geography by summing deaths and population before calculating rates. A custom label can be supplied for this aggregate.

`Portugal` and `Norte` are used as fixed comparator geographies in the annual metrics tab.

## Age Groups

The app harmonises INE age bands into these groups:

- `0 - 4 anos`
- `5 - 9 anos`
- `10 - 14 anos`
- `15 - 19 anos`
- `20 - 24 anos`
- `25 - 29 anos`
- `30 - 34 anos`
- `35 - 39 anos`
- `40 - 44 anos`
- `45 - 49 anos`
- `50 - 54 anos`
- `55 - 59 anos`
- `60 - 64 anos`
- `65 - 69 anos`
- `70 - 74 anos`
- `75 - 79 anos`
- `80 - 84 anos`
- `85 e mais anos`

Deaths recorded as `Menos de 1 ano` and `1 - 4 anos` are recoded into `0 - 4 anos`. `Total` and `Idade ignorada` are excluded from age-specific calculations.

For the `Menos de 75 anos` population option, the app excludes `75 - 79 anos`, `80 - 84 anos`, and `85 e mais anos` before computing rates.

## Mortality Metrics

All rates are calculated after filtering by selected year, area, cause of death, sex, and population scope.

### Deaths

`Óbitos` is the sum of deaths over the selected age bands, areas, and sex.

### Crude Mortality

Crude mortality is expressed per 100,000 inhabitants:

```text
crude rate = deaths / population * 100,000
```

Exact Poisson 95% confidence intervals are calculated for crude rates using `poisson.test()` and then scaled by the population denominator.

### Age-Standardised Mortality

Directly standardised mortality is calculated with `PHEindicatormethods::calculate_dsr()` using the European Standard Population 2013 weights embedded in the app.

The standard population weights used by the app are:

```text
0-4: 5000
5-9: 5500
10-14: 5500
15-19: 5500
20-24: 6000
25-29: 6000
30-34: 6500
35-39: 7000
40-44: 7000
45-49: 7000
50-54: 7000
55-59: 6500
60-64: 6000
65-69: 5500
70-74: 5000
75-79: 4000
80-84: 2500
85+: 2500
```

The multiplier is 100,000 and confidence intervals are requested at 95%.

### Proportional Mortality

Proportional mortality is calculated for a selected cause as:

```text
proportional mortality = deaths for selected cause / deaths from all causes * 100
```

The denominator is loaded for `Todas as causas de morte` for the same year, sex, and geography. In the annual metrics tab this denominator is loaded whenever `Mortalidade Proporcional` is selected, even if `Todas as causas de morte` is not one of the selected causes.

### Years Of Potential Life Lost

`AVPP` uses 70 years as the cutoff. Because INE data are grouped by age band, the app approximates age at death using each age band's midpoint:

```text
AVPP = sum(deaths in age band * max(70 - age midpoint, 0))
```

For `0 - 4 anos`, the midpoint is 2.5. For five-year age bands, the midpoint is the average of the lower and upper bound. Age groups with midpoints at or above 70 contribute zero years lost.

## Annual Metrics Tab

The annual metrics tab compares one selected metric for one selected year across:

- `Portugal`
- `Norte`
- the selected local area or aggregate of local areas

Users can select multiple causes of death. The table and plot are ordered from highest to lowest according to the value in the selected local area or aggregate. This ordering is intended to help identify which causes contribute most in the local geography, while keeping national and regional comparators visible.

## Forecasting

Forecasting is exploratory. It uses annual mortality-rate series from the observed mortality pipeline, after the selected geography, cause, sex, population scope, rate type, and fitting window are applied.

The app supports these model families through the `forecast` package:

- ARIMA
- ETS
- random walk with drift
- naive forecast
- Theta
- TBATS
- Holt
- damped Holt

The guided forecast tab uses simpler controls and recommends among available models using in-sample accuracy. The advanced forecast tab exposes model families, training windows, confidence interval level, optional transformation, diagnostics, backtesting, and structural-break exploration.

### Transformations

The model runner can fit models on transformed values and back-transform forecasts for display. The default workflow uses a log offset transformation where configured by the app controls. This can improve stability for positive rates but does not remove the need to inspect fit quality.

### Model Comparison

Model accuracy is summarised with:

- ME
- RMSE
- MAE
- MAPE
- MASE

The current recommendation logic prioritises lower RMSE, then MAE, then MASE, then MAPE, using the first available metric in that order.

Backtesting can evaluate forecasts against a holdout period from the end of the observed series. Holdout errors are computed against the observed values for overlapping years.

### Diagnostics And Structural Breaks

Diagnostics include residual plots, ACF, PACF, Ljung-Box tables, and model summaries where available.

Structural breaks are explored with `strucchange::breakpoints()` on the selected annual rate series. This is a screening tool: detected breakpoints should be interpreted with epidemiological context, data revisions, coding changes, and small-number instability in mind.

## Caching And Interruption

The app uses both in-memory and persistent RDS caching. Metadata and data have separate cache expiry settings.

Data requests are intentionally granular. If a long request is interrupted or an INE call fails, completed slices remain cached and can be reused in later runs. If a stale cached slice exists and a live request fails, the app may use the stale slice and show a warning.

## RDS Snapshot Source

The app can use prebuilt RDS files as an alternative to live INE requests. This is intended for faster app use when the relevant INE data have already been downloaded and normalised.

The snapshot source expects either separate files:

- `population.rds`
- `deaths.rds`

or one combined RDS file containing a list with `population` and `deaths` elements.

The population table must contain:

- `year`
- `area`
- `sex`
- `age_band`
- `pop`

The deaths table must contain:

- `year`
- `area`
- `sex`
- `cause`
- `age_band`
- `deaths`

When `Ficheiros RDS` is selected in the app, the same downstream metric calculations are used. The only difference is that rows are filtered from the snapshot files instead of requested from INE. `INE em directo` keeps the live INE path. The selector is available in the loading controls for observed mortality, annual metrics, and the advanced model specification; guided forecasts reuse the observed series that was already loaded.

The helper script `tools/build_ine_snapshot.R` can build these files from INE. Snapshot files should be rebuilt when INE updates the underlying indicators or when the app needs years, areas, or causes not present in the existing snapshot.

For large or slow indicators, the app can also read chunked files:

- `data/snapshots/population/year_<year>.rds`
- `data/snapshots/deaths/<indicator>/year_<year>/cause_<cause-token>.rds`

When `data/snapshots/snapshot_inventory.rds` is present, the app reads that manifest before loading chunk files. The inventory records the dataset, indicator, year, cause, relative path, row count, available areas, sexes, age bands, and source priority. This lets the app avoid unnecessary chunk discovery and makes missing data easier to identify. The manifest is rebuilt with:

```sh
Rscript tools/update_snapshot_inventory.R
```

`tools/build_0008206_snapshot_from_portal.R` is the preferred route for historical `0008206` death chunks. Instead of calling the INE API for each cause slice, it uses the INE web portal's BDDXplorer CSV export, then normalises that CSV into the same columns used by the app. The accompanying GitHub Actions workflow uses this portal exporter on a schedule, builds a limited number of missing area batches per run, and commits new chunks back to the repository. This avoids requiring one long INE download to complete successfully.

`tools/build_0008206_snapshot_chunks.R` remains available as an API-based fallback for the same chunked layout. It is usually slower for `0008206`, but can still be useful if the portal export route is temporarily unavailable.

`tools/build_population_snapshot_chunks.R` creates yearly population chunks. `tools/build_death_snapshot_chunks.R` creates the same per-year/per-cause death chunks for API-backed indicators such as `0013166`. When multiple death indicators contain the same year and cause, the snapshot reader resolves priority at row level by year, area, sex, cause, and age band. This means `0013166` is used ahead of `0008206` where it exists, while `0008206` can still fill areas not present in `0013166`.

## Interpretation Notes

- Small local areas and rare causes can produce unstable rates and wide confidence intervals.
- Direct standardisation reduces age-structure confounding but does not correct for all comparability problems.
- AVPP is approximate because age at death is inferred from grouped age bands.
- Forecasts extrapolate past rate patterns and should not be read as targets or official projections.
- INE data can be revised; the app cache is not a formal versioned data archive.
