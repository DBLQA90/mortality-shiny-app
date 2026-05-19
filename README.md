# Mortality Shiny App

Unofficial Shiny app for exploring Portuguese mortality indicators from INE.

The app supports observed mortality analysis, guided forecasting, advanced model comparison, diagnostics, and structural break exploration. Results are intended for exploration, decision support, and research workflows, and should be interpreted with appropriate epidemiological and statistical caution.

## What changed in v5

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

## Running The App

Install `pacman` once if needed:

```r
install.packages("pacman")
```

Then run the app from this repository:

```r
shiny::runApp("mortality-shiny-app.R")
```

`pacman::p_load()` will load or install the required runtime packages, including:

- `glue`
- `PHEindicatormethods`
- `tidyverse`
- `shiny`
- `forecast`
- `ineptr2`
- `strucchange`
- `memoise`
- `later`

## App Modules

### Observed Mortality

Explore historical mortality rates by geography, cause of death, sex, and population scope.

Outputs include:

- mortality rates per 100,000
- Poisson 95% confidence intervals for crude rates
- directly standardised rates using ESP 2013
- time-series plots
- summary and annual data tables

### Annual Metrics

Compares one selected year across Portugal, Norte, and a selected local aggregation. The tab shows one metric at a time for one or more causes of death, ordered from highest to lowest by the selected local aggregation.

Available metrics:

- deaths
- crude mortality
- standardised mortality
- proportional mortality, using all-cause deaths as the denominator for each location
- years of potential life lost before age 70

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

### Beginner Forecasting

Provides a guided forecasting workflow with simpler controls and reasonable defaults.

The user can choose:

- forecast horizon
- training window
- recommended model or model comparison mode

## Data Sources

All data is fetched from INE through `ineptr2`.

Population indicators:

- `0008273`
- `0003182`

Deaths by cause indicators:

- `0008206`
- `0013166`

The app harmonises age bands, recodes infant mortality into the `0-4` age group, excludes total or ignored age categories where needed, and can compute rates for the full population or the population under 75 years.

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

## Limitations

- INE API calls can be slow or temporarily unavailable.
- Indicator `0008206` is particularly slow for some historical mortality requests.
- Small municipalities and rare causes can produce sparse counts and unstable rates.
- Forecasts are sensitive to short time series, low counts, and structural breaks.
- Historical INE data revisions are not versioned inside the app.

This is a non-official tool.
