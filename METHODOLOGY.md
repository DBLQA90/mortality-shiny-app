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

### NUTS Vintages And Regional Aggregates

The four source indicators do not share one geography version:

| Indicator | Role | Geography | Years |
|---|---|---|---|
| `0003182` | population | NUTS-2002 | 1991-2013 |
| `0008273` | population | NUTS-2013 | 2011-2023 |
| `0008206` | deaths | NUTS-2013 | 1980-2022 |
| `0013166` | deaths | NUTS-2024 | 2022-2024 |

Municipality boundaries are identical across all three vintages, but the
regional groupings are not: NUTS-2024 moved Lezíria do Tejo out of `Alentejo`
into the new `Oeste e Vale do Tejo`. Reading INE's own regional rows across the
2022 seam therefore compares two different Alentejos. For 2022 all-cause deaths
`0008206` reports 11,327 and `0013166` reports 7,898, and because the app
prefers the newer indicator while population remains NUTS-2013, the regional
rate is understated by roughly 30%. `Portugal` and `Norte` are unaffected -
they are identical under both vintages.

The `Regiões` control offers the two honest responses:

- `Aproximação por municípios` rebuilds the region by summing its
  municipalities, using one fixed NUTS-2024 membership list for every year.
  Regions are unions of whole municipalities, so no parish-level data is
  needed. The series is continuous and means "this region as defined today".
  Where a municipality cannot be matched by name the shortfall is reported, not
  absorbed silently.
- `Dados originais INE` keeps INE's own regional rows and warns when the
  selected years span 2022, leaving the change of definition visible.

The difference is not cosmetic. For Alentejo in 2022, against Portugal, the
original rows give an SMR of 77.7 - implying below-average mortality in the
country's oldest region - while the municipal rebuild gives 113.3.

Municipality membership is read from `data/nuts_lookup.rds`, derived from INE's
own hierarchical geography codes and rebuilt with:

```sh
Rscript tools/build_nuts_lookup.R
```

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

`calculate_dsr()` normalises by the sum of the supplied standard weights. For the all-age scope this is the full ESP-2013 (weights summing to 100,000). For the `Menos de 75 anos` scope only the 0-74 age bands are supplied, so the routine returns the conventional premature-mortality rate standardised to the ESP-2013 0-74 sub-population. This under-75 rate is a valid rate per 100,000, but it uses a different standard age structure from the all-age rate, so the two are not directly comparable. The app makes this explicit by labelling the under-75 standardised rate as `padrão ESP 0-74` rather than rescaling it onto the all-age standard.

If the direct-standardisation routine cannot estimate a valid interval for very sparse selected data, the app reports the value as unavailable rather than silently substituting another method.

### Indirect Standardisation (SMR and ISR)

Direct standardisation needs stable age-specific rates in the area being
standardised. In a small municipality most age bands contain one or two deaths,
so the directly standardised rate becomes unstable, its interval very wide, and
for the sparsest selections `calculate_dsr()` cannot estimate an interval at
all. Indirect standardisation is the conventional alternative for small-area
comparison and is offered as the `SMR` and `Taxa Padronizada Indirecta`
metrics.

Expected deaths apply the reference area's age-specific rates to the local age
structure:

```text
expected = sum(local population in band i * reference deaths in band i / reference population in band i)
SMR      = observed / expected * 100
```

An SMR of 100 means the area experienced exactly the number of deaths expected
if it had the reference area's age-specific rates. The reference is selectable
(`Portugal` by default, `Norte` also available) and is loaded for the same
period, sex and cause as the area being compared, so both sides of the ratio
rest on identical data.

Intervals come from `PHEindicatormethods::calculate_ISRatio()`, which uses
Byar's method and exact Poisson limits below 10 observed deaths. The app also
reports whether the interval excludes the reference value, as
`Acima da referência`, `Abaixo da referência` or
`Sem diferença significativa`.

`Taxa Padronizada Indirecta` is the same comparison expressed as a rate per
100,000: the SMR multiplied by the reference crude rate. It carries no extra
information beyond the SMR, but puts the result on a scale that can be read
next to the crude and directly standardised rates.

Age bands present locally but absent from the reference cannot produce an
expected count. They are excluded and named in the result rather than being
silently treated as zero risk.

### Multi-Year Pooling

Any metric can be computed over a rolling 3- or 5-year window instead of a
single year. Deaths and population are both summed over the window, so the
denominator becomes person-years and the pooled value stays a valid rate:

```text
pooled rate = sum(deaths over window) / sum(population over window) * 100,000
```

This is not the mean of the annual rates, which would weight small years
equally with large ones. Windows are centred on the target year. At the ends of
the series the window is truncated rather than dropped, so the most recent year
stays visible; the reported period label names the real span and `n_years`
records how many years actually contributed.

Pooling trades resolution for stability. For Barrancos, all causes, 2022, the
single-year SMR is 213.9 (153.5-290.2); pooled over 2018-2022 it is 145.2
(123.1-170.1) - an interval roughly a third as wide, showing that 2022 was an
unusually bad year rather than the local level. Genuine year-to-year signal is
smoothed away in the same operation, so pooled and annual values answer
different questions.

Pooling applies to the annual comparison tab. Forecasting always uses unpooled
annual series, because a moving-average filter induces autocorrelation that
would invalidate the forecast intervals.

### Proportional Mortality

Proportional mortality is calculated for a selected cause as:

```text
proportional mortality = deaths for selected cause / deaths from all causes * 100
```

The denominator is loaded for `Todas as causas de morte` for the same year, sex, and geography. In the annual metrics tab this denominator is loaded whenever `Mortalidade Proporcional` is selected, even if `Todas as causas de morte` is not one of the selected causes.

Annual proportional mortality intervals use an exact binomial interval for selected-cause deaths over all-cause deaths.

### Years Of Potential Life Lost

`AVPP` uses 70 years as the cutoff. Because INE data are grouped by age band, the app approximates age at death using each age band's midpoint:

```text
AVPP = sum(deaths in age band * max(70 - age midpoint, 0))
```

For `0 - 4 anos`, the midpoint is 2.5. For five-year age bands, the midpoint is the average of the lower and upper bound. Age groups with midpoints at or above 70 contribute zero years lost.

Annual AVPP intervals use the Dobson et al. (1991) method for a weighted sum of Poisson counts. With `estimate = sum(deaths * years_lost)`, `variance = sum(deaths * years_lost^2)`, and `O` the total number of premature deaths (before the cutoff), the exact Poisson confidence limits of `O` are scaled by `sqrt(variance / O)` and centred on the estimate. This yields asymmetric limits that behave better for sparse local counts than a plain normal approximation; deaths at or after the cutoff contribute no years lost and are excluded from `O`. When there are no premature deaths the interval is reported as zero. These intervals remain approximate because age at death is inferred from grouped age-band midpoints.

### Source Transparency

Loaded rows keep their source indicator where this can be identified. Live INE loads retain the indicator code used for each row. Chunked death snapshots infer the death indicator from the snapshot path, so overlapping years can show `0013166`, `0008206`, or both depending on row-level fallback. Existing population RDS chunks do not contain their original indicator code, so the app labels them as `RDS population` unless future chunks include a `source_indicator` column.

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

The guided forecast tab uses simpler controls and recommends a model using out-of-sample forecast accuracy (see below). The advanced forecast tab exposes model families, training windows, confidence interval level, optional transformation, diagnostics, backtesting, and structural-break exploration.

Both forecasting tabs allow horizons up to 30 years beyond the last observed year. Longer horizons widen the gap between a statistical extrapolation and interpretable epidemiological expectation, so the final years should be read with extra caution, especially for sparse local series or unstable causes.

Each requested model is estimated independently. If one model fails, the app keeps the successful models and shows a model warning with the technical error message returned by the estimator. If all requested models fail, the forecast is treated as an error condition: the app shows `Erro detectado na previsão` and does not present forecast values as valid results for that selection.

Missing or incomplete source data can also invalidate a forecast. When RDS snapshots are selected, the app checks the snapshot inventory for the requested years, areas, and causes before fitting. If coverage is partial or unavailable, it shows a warning so the user can distinguish a modelling failure from a data-availability problem.

### Transformations

The model runner can fit models on transformed values and back-transform forecasts for display. The default workflow uses a log offset transformation where configured by the app controls. This can improve stability for positive rates but does not remove the need to inspect fit quality.

The offset is a data-dependent pseudo-count equal to half the smallest positive rate in the fitting series (or `1e-6` when every value is zero). Its value is shown in the transformation label (for example in the advanced model specification table), and when the series contains zeros the app flags that the offset materially affects the back-transformed forecast and intervals, since this is the case where a small additive constant has the largest effect.

### Model Comparison

Model accuracy is summarised with:

- ME
- RMSE
- MAE
- MAPE
- MASE

The current recommendation logic prioritises lower RMSE, then MAE, then MASE, then MAPE, using the first available metric in that order.

### Model Selection (out-of-sample)

The recommended model is chosen from out-of-sample forecast accuracy rather than in-sample fit, so short annual series do not simply reward the most flexible model. The most recent portion of the selected series is used as an evaluation region; its size is set by the user as a percentage of the available years.

Two schemes are offered in both the guided and advanced tabs:

- **Rolling validation (default):** for each year in the evaluation region the model is re-fitted on all earlier years and scored on a one-step-ahead forecast; the errors are pooled across origins. This uses the limited data efficiently and is less sensitive to any single split.
- **Single split:** the model is fitted once on the earlier years and scored on a single multi-step forecast over the whole evaluation region.

When the selected series is too short to leave at least three training years and one test year, selection falls back to the in-sample accuracy table and a note is shown. Out-of-sample errors are computed on the original rate scale (per 100,000), matching the forecasts and the holdout metrics. The evaluation refits models on the transformed modelling scale but scores back-transformed predictions.

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

By default, the app first reads local snapshot files under `data/snapshots` when both population and death snapshots are present. If the local snapshot files are absent or only partially present, the app reads the manifest and chunk files from the configured GitHub raw snapshot directory. `MORTALITY_SNAPSHOT_DIR` can force a specific local or remote folder, and `MORTALITY_USE_LOCAL_SNAPSHOTS=false` can skip the default local folder. A stray local inventory without snapshot data is not enough to redirect the app away from GitHub. The inventory records the dataset, indicator, year, cause, relative path, row count, available areas, sexes, age bands, and source priority. This lets the app avoid unnecessary chunk discovery and makes missing data easier to identify. The manifest is rebuilt with:

```sh
Rscript tools/update_snapshot_inventory.R
```

The `Disponibilidade de Dados` tab uses the same inventory to classify selected RDS coverage as available, partial, or unavailable. This is an inventory-level check: it tells whether the necessary chunks and requested areas are present before the app reads the data rows for an analysis.

When a user starts an analysis with `Ficheiros RDS`, the app repeats this inventory-level check for the active selection. Partial or unavailable coverage is shown as a warning before the detailed rows are loaded. If the requested rows are truly missing, the load still stops with an explicit message rather than silently producing a partial result.

`tools/build_0008206_snapshot_from_portal.R` is the preferred route for manually rebuilding historical `0008206` death chunks. Instead of calling the INE API for each cause slice, it uses the INE web portal's BDDXplorer CSV export, then normalises that CSV into the same columns used by the app. The repository no longer runs scheduled GitHub Actions jobs for this backfill because the snapshot archive is now committed.

`tools/build_population_snapshot_chunks.R` creates yearly population chunks. `tools/build_death_snapshot_chunks.R` creates the same per-year/per-cause death chunks for API-backed indicators such as `0013166`. When multiple death indicators contain the same year and cause, the snapshot reader resolves priority at row level by year, area, sex, cause, and age band. This means `0013166` is used ahead of `0008206` where it exists, while `0008206` can still fill areas not present in `0013166`.

## Interpretation Notes

- Small local areas and rare causes can produce unstable rates and wide confidence intervals.
- Direct standardisation reduces age-structure confounding but does not correct for all comparability problems.
- AVPP is approximate because age at death is inferred from grouped age bands.
- Forecasts extrapolate past rate patterns and should not be read as targets or official projections.
- INE data can be revised; the app cache is not a formal versioned data archive.
