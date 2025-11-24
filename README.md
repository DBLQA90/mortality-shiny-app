PNS Monitorização — Mortality Monitoring Shiny App

This Shiny application provides an interactive way to explore mortality rates in Portugal using data obtained directly from INE (Instituto Nacional de Estatística).
It was designed to support public-health surveillance, epidemiological analysis, and forecasting of mortality trends at the national and municipal level.

Main Features
1. Mortality Rates (Taxas de Mortalidade)
Download population and death counts dynamically from INE
Compute:
-Crude mortality rates
-Age-standardised rates (DSR, using ESP 2013)
-Poisson 95% confidence intervals
Visualise trends across years
Export tables as needed

2. Forecasting (Projecções)
Fit time series models using selected years:
-ARIMA
-ETS
-Random Walk with Drift
-Naive
-Theta
-TBATS
Combine observed values with projected values
Export:
-Forecast table as CSV
-Forecast graph as PNG

3. Breakpoint Detection (Análise de Quebras)
Detect structural changes in mortality trends
Works with both crude and age-standardised mortality rates
Visualise breakpoints and export a clean summary table

Performance Considerations
Access to INE’s API is relatively slow and sometimes rate-limited.
To avoid long delays every time a user adjusts a setting, the app uses a controlled workflow:
Data does NOT load automatically
Data is only downloaded after clicking:
-“Carregar dados”
-“Carregar projecções”
-“Carregar análise”
The results are cached with memoise to avoid repeated downloads

Data Sources
All data is fetched live from INE:
Population
-Indicators: 0008273, 0003182
Deaths by cause
-Indicators: 0008206, 0013166
Adjustments applied:
-Harmonisation of age bands
-Infant mortality recoded as 0–4 anos
-Exclusion of "Total" and "Idade ignorada" categories

Limitations
INE API is slow and may temporarily fail
Some municipalities have sparse mortality data, especially for specific causes
Breakpoint detection may be unstable for series with low counts
Historical INE revisions are not version-controlled by the app
