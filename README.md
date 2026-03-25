Overview

The application is organised into three main components:

Observed Mortality
Beginner Forecasting (Guided)
Advanced Forecasting (Full Modelling Framework)

All computations are performed dynamically from INE data, with caching mechanisms to improve performance.

1. Observed Mortality

This module allows users to explore historical mortality data.

Features
Selection of:
Geographic area (single or multiple municipalities)
Cause of death
Sex
Population scope:
Total population
Population under 75 years
Rate type:
Crude mortality rate
Age-standardised rate (ESP 2013)
Automatic computation of:
Mortality rates (per 100,000)
Poisson 95% confidence intervals (crude rates)
Directly standardised rates (DSR) using ESP 2013
Outputs:
Summary table
Time series plot
Annual data table

2. Beginner Forecasting (Guided)

This module provides a simplified, user-friendly forecasting workflow.

Features
Uses the series already loaded in Observed Mortality
Minimal configuration required:
Forecast horizon
Training window:
Full history
Last 10 / 15 / 20 years
Mode:
Recommended model
Model comparison
Outputs:
Forecast plot
Summary interpretation
Reliability feedback (based on data characteristics)
Purpose

Designed for users who want:

Quick projections
Reasonable defaults
Minimal statistical configuration

3. Advanced Forecasting

This is a full modelling framework for time series analysis and evaluation.

It is divided into multiple sub-modules:

3.1. Model Specification

Users define the full forecasting setup:

Model families:
ARIMA (auto or manual)
ETS (auto or manual)
Random Walk with Drift
Naive
Theta
TBATS
Holt / Holt (damped)
Additional controls:
Training window (custom years)
Forecast horizon
Confidence level
Data transformation:
Log with offset
None
Manual parameter specification available for:
ARIMA (p, d, q, seasonal terms)
ETS (error, trend, seasonality)

3.2. Forecast Results

Visualisation of:
Observed vs forecasted values
Prediction intervals
Two modes:
Single model
Model comparison
Outputs:
Forecast table
Summary metrics
Download options:
CSV (data)
PNG (plot)

3.3. Diagnostics

Model validation tools:

Residual analysis:
Time series of residuals
ACF / PACF plots
Statistical testing:
Ljung–Box test
Model summaries:
Full statistical output

3.4. Backtesting and Model Comparison

Two validation approaches:

In-sample performance
Holdout validation (last k years)
Metrics
ME
RMSE
MAE
MAPE
MASE
Outputs
Model ranking
Accuracy tables
Comparative plots

3.5. Structural Break Analysis

Detection of structural changes in mortality trends:

Uses breakpoint detection methods
Identifies:
Break years
Segmented periods
Changes in mean level
Outputs
Breakpoint plot
Segment summary table
Interpretation text
Data Sources

All data is fetched dynamically from INE using the ineptR package:

Population
Indicators:
0008273
0003182
Deaths by cause
Indicators:
0008206
0013166
Adjustments applied
Harmonisation of age bands
Infant mortality recoded as 0–4 years
Exclusion of:
“Total”
“Idade ignorada”
Optional exclusion of older age bands (≥75 years)
Performance and Architecture
Controlled Data Loading

Data is only fetched when explicitly requested:

“Carregar dados”
“Carregar projecções”
Caching
INE queries cached using memoise
Separate caches for:
Dimension metadata
Data queries
Benefits
Avoids repeated API calls
Improves responsiveness
Reduces INE rate-limit issues
Limitations
INE API:
Can be slow
May temporarily fail or rate-limit requests
Small-area data:
Some municipalities have sparse counts
Leads to instability in rates and forecasts
Forecasting:
Sensitive to:
Short time series
Structural breaks
Low counts
Breakpoint detection:
Requires sufficient data length
May produce unstable results in noisy series
Data revisions:
Historical INE updates are not version-controlled within the app
Notes
This is an unofficial tool (“não oficial”)
Intended for:
Exploration
Decision support
Research workflows
Results should be interpreted with appropriate epidemiological and statistical caution
