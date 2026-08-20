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

Available years and causes of death are read from INE metadata where possible. The user-facing location list is derived from the NUTS lookup of the selected vintage, so the places offered and the places the app can actually aggregate are the same list; a hand-written list is kept only as a fallback for when no lookup has been built.

When a requested year can be found in more than one indicator, the app builds a non-overlapping source-year plan and de-duplicates loaded rows by year, area, sex, cause, and age band. Data are downloaded in small year, area, and cause slices so partial downloads can be cached and reused.

## Geography

Users can select one or more local areas. When more than one local area is selected, the app treats the selection as one combined geography by summing deaths and population before calculating rates. A custom label can be supplied for this aggregate.

`Portugal` and `Norte` are used as fixed comparator geographies in the annual metrics tab.

### Ambiguous INE Labels

The snapshot builders resolve a geography through the INE category *label*.
Several labels are not unique, and when a label matches more than one category
the download returns all of them and they are summed into one row. Three cases
affect this archive:

| Indicator | Label | Matches | Consequence |
|---|---|---|---|
| `0003182` | `Lisboa` | region `17` + município `1711106` | 1991-2013 population is region + município |
| `0008206`, `0008273` | `Calheta` | `2004501` (Açores) + `3003101` (Madeira) | two municipalities added together |
| `0008206`, `0008273` | `Lagoa` | `1500806` (Algarve) + `2004201` (Açores) | two municipalities added together |

`Lisboa` is the damaging one. Deaths are always the município, so dividing them
by a denominator that also contains the region understates Lisboa's mortality
roughly six-fold for 1991-2013, and produces a spurious six-fold jump at the
2013/2014 source seam (crude rate ~211 per 100,000 in 2013 against ~1,231 in
2015). Any trend or forecast covering those years is affected.

A scan of all 308 areas across the 2013/2014 seam found exactly these
discontinuities, with a median population ratio of 1.0139 elsewhere, so the
rest of the archive is unaffected.

`0008273` and `0013166` name the Lisbon region `Área Metropolitana de Lisboa`,
so `Lisboa` is unambiguous there and population from 2014 onward is correct.

The repair requests each affected geography by its unique category code and
stores it under an unambiguous label, splitting the conflated municipalities
into `Calheta (R.A.A.)`, `Calheta (R.A.M.)`, `Lagoa` (Algarve) and
`Lagoa (R.A.A.)`:

```sh
Rscript tools/fix_ambiguous_areas.R            # add dry_run=true to preview
```

The repair has been applied to the committed archive. Afterwards, Lisboa's crude
rate runs 1474.6 (2000), 1367.5 (2013), 1240.7 (2014), 1280.7 (2022) - continuous
across the source seam that previously produced the six-fold jump - and a rescan
of all areas across 2013/2014 finds none moving more than 20%. Re-run the tool
after rebuilding any chunk from scratch, since the builders still resolve
geographies by label.

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

**The app therefore always builds a region by summing its municipalities**,
using one fixed membership list for every year. That much is not a user-facing
option. *Which* membership list is - see *Choosing the NUTS vintage* below.

Three facts make it sound:

- Regions are unions of whole municipalities, so no parish-level data is needed.
- Municipality boundaries have not moved: the archive holds the same 308
  municipalities in every year from 1991 to 2023.
- Regions whose definition never changed are unaffected. Only regions INE
  redrew differ, and there the aggregate is the point.

The alternative - reading INE's own regional rows - is not merely less tidy, it
is wrong across the 2022 seam. For Alentejo in 2022, against Portugal, INE's
rows give an SMR of 77.7, implying below-average mortality in the country's
oldest region; the municipal aggregate gives 113.3. As a series, the
standardised rate reads 1072 (2021), 722 (2022), 676 (2023) from INE's rows
against 1093, 1060, 997 when aggregated. Presenting that as a toggle asked
users to arbitrate a question about NUTS vintages they have no way to judge,
and let the wrong answer through.

What this costs:

- Regions inherit only deaths that INE could assign to a municipality. INE
  publishes `Ignorado` and `Estrangeiro` geographies at every level - deaths
  with unknown or foreign municipality of residence - and the archive excludes
  them, because they are not places and must not join a regional aggregate.
  For most years this is nothing at all: 1991-2023 the municipal rows sum to
  INE's national row exactly. **2024 is the exception**, with 545 deaths (0.46%)
  unattributed - 355 on the mainland and 190 on the islands. Re-fetching 2024
  returned byte-identical data, so this is INE's published position for a
  recent year rather than an incomplete download.
- A region means "this territory as defined by the selected vintage", applied to
  every year. That is a deliberate choice, not a reconstruction of what was
  published at the time.

Setting `MORTALITY_REGION_MODE=original` restores INE's own rows for anyone who
needs to reproduce a published regional figure. Note that those rows are not
trustworthy for five regions: see *A defect in INE's own regional rows* below.

### NUTS I: Continente, Açores And Madeira

INE's geography has a NUTS I level above the regions: `Continente` and the two
autonomous regions. The islands carry the same name at NUTS I, II and III -
they are the same territory - so `Continente` is the only label this level adds.
It is built like every other region, by summing its 278 mainland municipalities,
and is available under both vintages (the 2024 reform did not move anything
between the mainland and the islands, so the two vintages give identical
figures for all three).

The three NUTS I units partition the country, and the partition closes exactly
against INE's national row for the years whose municipal rows are complete:

| All-cause deaths | Continente | Açores | Madeira | Sum | INE national row |
|---|---:|---:|---:|---:|---:|
| 2021 | 119,589 | 2,366 | 2,875 | 124,830 | 124,830 |
| 2022 | 118,517 | 2,712 | 3,104 | 124,333 | 124,333 |
| 2023 | 113,164 | 2,369 | 2,791 | 118,324 | 118,324 |

NUTS III is deliberately not offered in the selectors - it would add 21 entries
to a list already 316 long, and every one sits inside an offered NUTS II region
- but `region_municipalities()` resolves a NUTS III name if one is supplied.

### Overlapping Selections

Selected areas are summed into one geography, so an overlapping selection has to
be flagged. The two cases behave differently, and the warning says which:

- **`Portugal` with anything else** is genuinely double-counted. Portugal is not
  a region label, so it is never expanded into municipalities: its own row is
  loaded and added to whatever else is selected. `Portugal + Beja` for 2021
  gives 125,382 against a true 124,830.
- **A region with something inside it** is *absorbed*, not double-counted. The
  region is expanded into a sorted unique set of municipalities, so
  `Alentejo + Beja` yields Alentejo (8,248, not 8,800) and
  `Continente + Centro` yields Continente. Nothing is counted twice, but the
  user does not get the combination they asked for.

Disjoint selections sum as expected: `Centro + Alentejo` gives 31,464, which is
23,216 + 8,248.

Containment is detected by comparing municipality *sets*. An earlier version
intersected a region's municipalities with the other selected area *names*,
which caught `Alentejo + Beja` but never `Continente + Norte`, because `Norte`
is a region label and so appears in nobody's municipality list. That was
harmless only while every offered region sat at one NUTS level and was therefore
disjoint from the others; NUTS I regions contain NUTS II ones.

### Choosing The NUTS Vintage

Which grouping is used is the user's choice. The app ships one lookup per
vintage and a single app-wide selector in the page header:

| | Regions | Lisbon | Lezíria do Tejo, Oeste, Médio Tejo |
|---|---:|---|---|
| NUTS 2013 | 7 | one region, `Área Metropolitana de Lisboa` | inside `Alentejo` and `Centro` |
| NUTS 2024 | 9 | `Grande Lisboa` + `Península de Setúbal` | form `Oeste e Vale do Tejo` |

Both lookups cover the same 308 municipalities under the same labels, verified
by test. Switching vintage therefore never changes which data is read, only how
it is grouped - and every year remains available under either, so both series
are continuous across the 2022 seam.

Six region names exist in both vintages and mean different things in each:
NUTS-2013 `Centro` has 100 municipalities, NUTS-2024 `Centro` has 77. The
selector is in the header rather than in a tab so the active definition is on
screen wherever a regional figure is read, and a selection naming a region the
new vintage does not have is dropped and reported rather than silently
misread.

Under NUTS 2013 the municipal aggregates reproduce INE's published regional
rows **exactly**, where those rows exist in the archive. For 2021 all-cause
deaths: `Norte` 37,121 and `Alentejo` 11,742 from the municipal sum against
37,121 and 11,742 from INE's own rows, and the seven regions sum to 124,830,
which is the national total to the death. That is the strongest available check
that the lookup is right: it is the same arithmetic INE did, from the same
parts.

`Grande Lisboa` (2024) plus `Península de Setúbal` (2024) reproduce
`Área Metropolitana de Lisboa` (2013) exactly - 33,288 deaths in 2021 - which
is the reform's own definition, recovered from the data rather than asserted.

### A Defect In INE's Own Regional Rows

Five region names denote both a NUTS II and a NUTS III unit, and the autonomous
regions are a NUTS I level as well. The archive stores rows by label, so those
levels collapse into one row and the deaths are counted two or three times:

| Region, 2024 | INE rows as stored | Municipal sum | Ratio |
|---|---:|---:|---:|
| Grande Lisboa | 41,430 | 20,715 | 2.00 |
| Península de Setúbal | 18,240 | 9,120 | 2.00 |
| Algarve | 11,268 | 5,616 | 2.01 |
| Região Autónoma da Madeira | 7,722 | 2,531 | 3.05 |
| Região Autónoma dos Açores | 7,362 | 2,308 | 3.19 |

The four regions with no name collision are all within 0.3%. The app never
reads these rows, because regions are always summed from municipalities and the
area selectors only ever offer NUTS II names and municipalities. The check that
the municipal figures are the correct ones: `Continente` (113,357) plus the two
islands' single-counted figures reconstructs the 2024 national total of 118,386.

This is a further reason not to use `MORTALITY_REGION_MODE=original`, and it
means `tools/refresh_snapshots.R task=nuts2` would deepen the problem rather
than fix it.

### Rebuilding The Lookups

Municipality membership is read from `data/nuts_lookup_2013.rds` and
`data/nuts_lookup_2024.rds`, derived from INE's own hierarchical geography
codes - a municipality's code is prefixed by its NUTS III and NUTS II parents -
and rebuilt with:

```sh
Rscript tools/build_nuts_lookup.R indicator=0013166 out=data/nuts_lookup_2024.rds
Rscript tools/build_nuts_lookup.R indicator=0008206 out=data/nuts_lookup_2013.rds \
  canonical=data/nuts_lookup_2024.rds
```

`canonical=` relabels municipalities to the archive's vocabulary, because
0008206 publishes two municipalities called `Calheta` while 0013166
distinguishes `Calheta (R.A.A.)` from `Calheta (R.A.M.)`. Matching is by
geography code first and by name second: the code encodes the hierarchy, so it
is stable only for municipalities the reform did not move (132 of 308), and the
ambiguous island labels are all in that group, so no name is ever matched
against two candidates. Anything resolving by neither route is an error rather
than a guess.

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

### Infant Mortality

`Mortalidade Infantil` is deaths under one year of age per 1,000 live births.
Neither part comes from the main pipeline, and neither could:

- **Numerator.** Both death indicators publish `Menos de 1 ano` as its own age
  band, but the ingest recodes it into `0 - 4 anos`, so the main death archive
  cannot separate an infant death from a death at age four. For Portugal 2024 it
  holds 286 deaths in `0 - 4 anos`, being 254 infant deaths plus 32 at ages one
  to four. `tools/fetch_infant_deaths.R` writes a parallel dataset holding only
  the under-1 band, leaving every existing rate untouched.
- **Denominator.** Live births, not population. Neither population indicator has
  an under-1 age band, so infant deaths per under-1 population is not computable
  at all. `tools/fetch_births.R` writes live births, assembled from three INE
  vintages.

```text
infant mortality rate = deaths under 1 year / live births * 1,000
```

The interval is an exact Poisson interval on the death count, scaled by births.
Births are treated as a fixed denominator, the usual convention: the sampling
variation that matters is in the small number of deaths.

Coverage is 1995-2024, bounded by births. The three source vintages are used
only for the years the others do not cover, and they agree where they overlap -
`0008084` and `0012434` both report 83,671 live births for Portugal in 2022 -
which is what confirms the dimension handling is right in each. The assembled
series runs continuously across both seams: 2.94 (2013), 2.87 (2014), 2.44
(2020), 2.43 (2021) per 1,000.

The national series reconciles with INE's published figures throughout: 7.43 in
1995, 5.52 in 2000, 3.51 in 2005, 2.53 in 2010, 3.00 in 2024.

At municipality level the metric is extremely sparse - a small municipality may
record fewer than ten births in a year, so a single death moves the rate by
hundreds per 1,000 and the interval is correspondingly enormous. Barrancos in
2024 had 9 live births and no infant deaths, giving 0.0 with an upper limit of
409.9. Multi-year pooling helps but cannot manufacture events that did not
happen.

The app handles this in two ways, neither of which suppresses a value.

**The count is offered as its own metric.** `Óbitos infantis (< 1 ano)` reports
the number of deaths under one year with an exact Poisson interval and no
denominator at all. It is the honest answer at municipal scale: it says what
happened and cannot be misread as a comparable rate. Because it needs only the
numerator it also covers 1991-1994, where the rate cannot be computed. Unlike
the other counts it is *not* annualised when a window is pooled - a
municipality with two infant deaths in three years would read as `1`, a rounded
fraction, when what happened is two deaths. Left as a window total it is also
exactly the numerator of the pooled rate beside it.

**Rates on thin denominators are marked.** A rate computed on fewer than 1,000
live births in the period is shown with an asterisk, in the table, in the CSV
and on the chart. The threshold is not a significance test but a statement about
resolution: below 1,000 births, one additional death moves the rate by more than
one whole unit per 1,000 - larger than the entire national rate of about 3.
Ranking such places, or reading a change between years, is reading noise. Most
Portuguese municipalities fall below the threshold, which is the point: the mark
describes the ordinary case rather than singling out a few outliers. The value
is still shown, still exact, and its interval already states the uncertainty;
the mark only stops a reader skimming the table from treating it as comparable.

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

For five-year age bands, the midpoint is the average of the lower and upper bound. Age groups with midpoints at or above 70 contribute zero years lost.

`0 - 4 anos` is a special case, because most of its deaths are not spread across it. Its nominal midpoint of 2.5 credits every death in the band with 67.5 lost years, but infant deaths sit at the very bottom of it and lose close to the whole 70 — and they are the majority of the band, not a fringe of it: of the 286 deaths Portugal recorded in `0 - 4 anos` in 2024, 254 were infants. Applying the band midpoint to all of them understates AVPP.

The under-1 counts fetched for infant mortality make the correction possible, so the app applies it. `split_infant_age_band()` divides the band into two before the weights are applied:

| Band | Midpoint | Years lost against a cutoff of 70 |
|---|---:|---:|
| `< 1 ano` | 0.5 | 69.5 |
| `1 - 4 anos` | 3.0 | 67.0 |

The `1 - 4 anos` midpoint is 3, not 2.5: the band spans exact ages 1 up to 5. Deaths are unchanged by the split; only their weights move. The under-1 archive covers 1980-2024, so the correction applies to every year the app offers, but it is applied only when the counts cover the whole pooled window — with partial coverage the band is left undivided and AVPP falls back to its previous behaviour, rather than silently assigning the uncovered years' infant deaths to ages one to four.

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

Back-transforming a forecast of `log(rate)` with `exp()` returns the **median**,
not the mean: for a lognormal, `E[Y] = exp(mu + sigma^2 / 2)`. Reported directly,
the point forecast would therefore sit below the expected value, and the gap
widens with the forecast variance - worst at the long horizons the app permits.
The app recovers the standard deviation from the interval the model reported and
adds the variance term, so the point forecast is the expected value. Interval
limits are quantiles and map through a monotone transform unchanged, so they are
back-transformed directly and are unaffected. The correction can be switched off
in the advanced tab, which returns the median. In a 20-year projection the two
differ by roughly 5%; at one year ahead, by about 0.3%.

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

Structural breaks are explored with `strucchange::breakpoints()` on the selected
annual rate series, fitted as a **segmented trend** (`rate ~ time`): each segment
gets its own intercept and slope, so a break is reported where the level or the
rate of change shifts.

A mean-only model (`rate ~ 1`) is the wrong tool for these series. Portuguese
mortality falls steadily - the national standardised rate is down about 42% since
1991 - and with no slope to fit, a mean-shift model explains a smooth decline by
cutting it into a staircase of level shifts. On the national all-cause series it
reports five breaks (1994, 1999, 2005, 2009, 2013), none of which is an event;
the segmented trend reports one, in 2013. Where a series is too short to support
per-segment trends the app falls back to the mean-only model and says so.

This remains a screening tool: detected breakpoints should be interpreted with
epidemiological context, data revisions, coding changes, and small-number
instability in mind.

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
