# Methodology

This document describes how the app prepares INE data, computes mortality metrics, and fits exploratory forecasts. The app is non-official and intended for analysis support, not as a substitute for validated statistical production workflows.

## Data Sources

The app fetches data from INE through `ineptr2`.

Population indicators:

- `0012918` (NUTS-2024, 2021-2025 — a revision of the series below, not a continuation; see *The 2021 population revision*)
- `0008273` (NUTS-2013, 2011-2023)
- `0003182` (NUTS-2002, 1991-2013)

Deaths by cause indicators:

- `0008206` (NUTS-2013, 1980-2022)
- `0013166` (NUTS-2024, 2022-2024)

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
discontinuities, with a median population ratio of 1.0139 elsewhere. That scan
used a 20% threshold, chosen to catch this six-fold error, and it did. It is too
coarse to see the systematic 3-7% step the same seam also carries for other
reasons - see *The 2013/2014 source seam* below.

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
| `0012918` | population | NUTS-2024 | 2021-2025 |
| `0008206` | deaths | NUTS-2013 | 1980-2022 |
| `0013166` | deaths | NUTS-2024 | 2022-2024 |

Municipality boundaries are identical across every vintage, but the
regional groupings are not: NUTS-2024 moved Lezíria do Tejo out of `Alentejo`
into the new `Oeste e Vale do Tejo`. Reading INE's own regional rows across the
2022 seam therefore compares two different Alentejos. For 2022 all-cause deaths
`0008206` reports 11,327 and `0013166` reports 7,898, and because the app
prefers the newer indicator while population remains NUTS-2013, the regional
rate is understated by roughly 30%. `Portugal` and `Norte` are unaffected -
they are identical under both vintages.

**The app therefore always defines a region by its municipalities**, using one
fixed membership list for every year - the population is always their sum. *Which*
membership list is a user choice (see *Choosing the NUTS vintage*), and so is
where the deaths come from (see *Regional deaths: INE's rows or the municipal
sum*).

Three facts make it sound:

- Regions are unions of whole municipalities, so no parish-level data is needed.
- Municipality boundaries have not moved: the archive holds the same 308
  municipalities in every year from 1991 to 2023.
- Regions whose definition never changed are unaffected. Only regions INE
  redrew differ, and there the aggregate is the point.

Reading INE's regional rows *without regard to vintage* is wrong across the 2022
seam. For Alentejo in 2022, against Portugal, NUTS-2024 rows divided by NUTS-2013
population give an SMR of 77.7, implying below-average mortality in the country's
oldest region; the municipal aggregate gives 113.3. As a series, the standardised
rate reads 1072 (2021), 722 (2022), 676 (2023) from those rows against 1093, 1060,
997 when aggregated. The fix is not to avoid INE's rows but to use each only for
the territory it describes - which is what the regional-rows source does, taking a
row only where its vintage matches the selected definition.

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

### Regional Deaths: INE's Rows Or The Municipal Sum

**Summing municipalities under-counts cause-specific deaths, in every year.**
INE publishes complete municipal totals but incomplete municipal age breakdowns,
and the missing detail is concentrated where counts are small. Because every
rate in the app is built from age bands, a region rebuilt from its
municipalities' age bands comes out short. Lung cancer, share of each region's
age-banded deaths lost by summing its municipalities, against INE's own row:

| Year | Norte | Centro | Alentejo | Açores | Madeira |
|---|---:|---:|---:|---:|---:|
| 2002 | −1.4% | −2.3% | −5.5% | −10.1% | −2.3% |
| 2010 | −1.4% | −1.8% | −8.9% | −14.8% | −11.0% |
| 2013 | −1.2% | −2.1% | −8.9% | −18.4% | −11.8% |
| **2014** | **−30.7%** | **−50.7%** | **−64.6%** | **−71.4%** | **−83.5%** |
| 2018 | −0.5% | −1.3% | −7.6% | −18.1% | −4.0% |
| 2021 | 0.0% | −0.2% | −3.4% | −11.6% | −6.0% |

2014 is the extreme of a persistent bias, not an isolated defect. All-cause
deaths are barely affected - the cells are large - which is why every all-cause
reconciliation in this document closes exactly. The municipal *totals* are
complete; it is the breakdown by age that is not.

INE's regional rows do not have the problem: their age bands sum to their
totals exactly, in every year fetched, 2014 included. So the app offers two ways
to build a region's deaths, chosen once in the page header:

| | Deaths come from | Trade-off |
|---|---|---|
| **Linhas regionais do INE** (default) | INE's row for the territory wherever one exists; the municipal sum otherwise | Accurate; redrawn regions step where the source switches |
| **Soma dos municípios** | Always the municipal sum | One consistent source and no seams; cause-specific figures biased low |

Population is the sum of municipalities under both. It is complete; the defect is
only in the age breakdown of deaths.

**Every region is built from INE rows in every year** except the two Lisbon
regions before 2022:

- **Continente, Norte, Algarve, Açores, Madeira** are the same territory under
  both vintages, so their own rows serve every year - `0008206` to 2021 and
  `0013166` from 2022, the death archive's own precedence.
- **The redrawn regions are composed from subregion (NUTS III) rows** in the
  years their own row does not exist. Oeste, Médio Tejo and Lezíria do Tejo were
  subregions under both definitions, which is what makes this possible:

| Vintage | Years | Region | Built from |
|---|---|---|---|
| NUTS 2013 | 2023-2024 | Centro | NUTS-2024 Centro + Oeste + Médio Tejo |
| NUTS 2013 | 2023-2024 | Alentejo | NUTS-2024 Alentejo + Lezíria do Tejo |
| NUTS 2013 | 2023-2024 | Área Metropolitana de Lisboa | Grande Lisboa + Península de Setúbal |
| NUTS 2024 | 1991-2021 | Alentejo | its four NUTS-2013 subregions |
| NUTS 2024 | 1991-2021 | Centro | six NUTS-2013 subregions **+ Sertã + Vila de Rei** |
| NUTS 2024 | 1991-2021 | Oeste e Vale do Tejo | NUTS-2013 Oeste + Médio Tejo + Lezíria **− Sertã − Vila de Rei** |

Sertã and Vila de Rei moved from Médio Tejo to Beira Baixa in 2024, which is the
only reason two compositions are not pure subregion sums; their own municipal rows
make the correction, about 1% of Centro. Every composition was checked two ways:
against the municipality lookups, where each reproduces the region's membership
exactly (a unit test), and against **2022, the one year both indicators publish**,
where each composed value equals the directly published row to the death -
all-cause and lung cancer, all six compositions.

- **Grande Lisboa and Península de Setúbal before 2022** cannot be composed:
  under NUTS 2013 the Lisbon metropolitan area was a single subregion. They fall
  back to the municipal sum and the app warns about the change of source. In
  practice the effect is negligible - their municipalities are large enough to
  keep their age detail - but it is not guaranteed for every cause.

**A municipality selected on its own cannot be repaired this way**: there is no
finer row to fall back on. The app warns instead - for any cause-specific
municipal figure, and more strongly when 2014 is included.

The rows are in `data/snapshots/regional_deaths/<indicator>/year_<y>.rds`, fetched
by `tools/fetch_regional_deaths.R` with one request per year. They are keyed by
**geography code**, never by label: Algarve, the Lisbon metropolitan area and both
autonomous regions carry the same name at two or three NUTS levels, and that is
exactly why the label-keyed regional rows in the municipal archive are
multi-counted. The fetcher checks that each row's age bands sum to its total, and
records the count of cells where they do not; it has been zero in every year.

The substitution is `substitute_regional_deaths()`, applied immediately after
every death load - observed rates, both forecast tabs, annual metrics including
the SMR reference and the proportional-mortality denominator, and avoidable
mortality - so all tabs build a region the same way. A nested selection covers
each municipality once: `Continente` with `Norte` substitutes Continente only. A
composition is applied only when every one of its subregion rows is present;
otherwise that year keeps the municipal sum rather than dropping part of the
territory.

The effect on a series is not subtle. Lung cancer, crude rate per 100,000,
NUTS-2024 Alentejo:

| | 2012 | 2013 | 2014 | 2015 | 2021 |
|---|---:|---:|---:|---:|---:|
| sum of municipalities | 44 | 39 | **16** | 40 | 45 |
| INE rows | 48 | 44 | **48** | 46 | 47 |

### Health-System Geography: ULS And ARS

Health planning in Portugal is organised by **Unidade Local de Saúde** (ULS),
grouped into the five former ARS regions. These do not nest inside NUTS: five
ULS straddle a NUTS II boundary (Entre Douro e Vouga, Guarda, Estuário do Tejo,
Médio Tejo, Região de Leiria), and ARS Norte and NUTS Norte are different sets of
municipalities. So ULS and ARS are a second, independent geography, built the
same way as NUTS regions - as unions of whole municipalities - and offered in
every area selector.

Membership comes from the `Var` sheet of the PNS2030 workbook *Indicadores de
Apoio ao Planeamento Local em Saúde*, keyed by INE's 2024 municipality codes, via
`tools/build_uls_lookup.R` into `data/uls_lookup.rds`. The workbook itself is not
committed. ARS follow from the first digit of the ULS code. A ULS has no NUTS
vintage: the same municipalities apply under either definition.

**Five ULS cannot be built individually.** Lisboa, Loures and Porto are each
divided between two ULS at parish level, and nothing below municipality exists in
the app's data. They are offered as the two smallest unions that contain only
whole municipalities, which are exact:

| Offered as | Municipalities |
|---|---|
| ULS Santo António + São João | Gondomar, Maia, Porto, Valongo |
| ULS Loures/Odivelas + Santa Maria + São José | Lisboa, Loures, Mafra, Odivelas |

The workbook handles these differently, and wrongly for counts: it assigns the
*whole* shared municipality to each ULS that touches it. ULS Santo António is
Gondomar plus all of Porto; ULS São João is Maia, Valongo and all of Porto. Porto's
3,000 deaths in 2019 are counted twice, so the workbook's 14 Norte ULS sum to 38,271
against its own ARS Norte total of 35,278. Rates for those five ULS mix in the whole
of the shared municipality in the same way.

The result is: **34 individual ULS, 2 groups and 5 ARS**, covering the 278 mainland
municipalities with each municipality counted once at the ULS level and once at the
ARS level (unit-tested). Checked against the workbook's own death counts (I37),
which the population revision does not affect: in 2023, 34 of 39 ULS and ARS match
exactly and the rest differ by one death. The five ARS sum to Continente exactly.

**Regional deaths for ULS.** Most ULS are finer than any NUTS III subregion - they
are frequently two or three municipalities - so INE publishes no row for them.
Seven coincide exactly with a subregion in both vintages (Alto Minho, Viseu
Dão-Lafões, Alentejo Litoral, Baixo Alentejo, Alto Alentejo, Alentejo Central,
Algarve), and ARS Alentejo and ARS Algarve with a NUTS II region; those take INE's
complete rows like any region. The rest use the municipal sum, with the
cause-specific warning. The split falls fortunately: the uncovered ULS are mostly
urban, whose large municipalities keep their age detail, while the rural ULS where
the loss is worst are largely the covered ones. Not all - **ULS Guarda, with 13
small municipalities and no INE row, shows zero lung-cancer deaths in 2014**, and
figures like that are why the warning names the areas it applies to.

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

### The 2013/2014 Source Seam

The population archive splices three indicators, and consecutive ones disagree
about the years they both publish. `0003182` (NUTS-2002, frozen in 2014) serves
1991-2013 and `0008273` (NUTS-2013) serves 2014-2020, so there is a handover
between 2013 and 2014.

Measured on the years both publish, for Portugal:

| | 2011 | 2012 | 2013 |
|---|---:|---:|---:|
| crude rate, `0008273` against `0003182` | −0.16% | −0.16% | −0.16% |
| **standardised rate** | **−2.49%** | **−2.74%** | **−2.95%** |

The totals agree to a sixth of a per cent. The *standardised* rates differ by
about 3%, because the two indicators distribute the population across age bands
differently and direct standardisation is sensitive to exactly that. The 65+
share barely moves (19.85% against 19.99% in 2013), so it is redistribution
within the bands, not a gross ageing shift.

The consequence is a break that reads as real:

| Portugal, standardised rate | |
|---|---:|
| 2012 → 2013, one consistent indicator | −3.00% |
| **2013 → 2014, as the archive shows it** | **−6.70%** |
| 2013 → 2014, `0008273` on both sides | −3.86% |

Roughly half the apparent fall is the handover. `strucchange` places a
breakpoint at `2013|2014` in Portugal's standardised series in every fit window
from 1991 through 2002. On a consistent basis the transition is ordinary.

At municipality level it is larger. Year-over-year population change across the
305 municipalities present in both indicators:

| Transition | IQR | moving more than 2% |
|---|---:|---:|
| 2012-13 | 1.09 | 19 |
| **2013-14** | **2.19** | **126** |
| 2014-15 | 0.94 | 9 |

Corvo −13.8%, Tabuaço −9.3%, Aguiar da Beira +8.4%, Faro +7.2%, Lisboa +6.3%,
Porto +5.1%, against ±1-3% either side.

This was missed by the earlier check of the same seam, which reported that no
area moved more than 20%. That is true, and the threshold was chosen to catch
the six-fold `Lisboa` conflation, which it did. It cannot see a systematic 3-7%
step.

**The seam is not removable.** `0008273` starts in 2011, so the handover can sit
at 2010/2011, 2011/2012, 2012/2013 or 2013/2014, and the discrepancy is of the
same character wherever it is placed - moving it to 2010/2011 would reduce the
step from about 2.95% to about 2.49%. There is no municipal population series
covering 1991-2010 on the newer basis. The archive therefore keeps `0003182` for
1991-2013, which is also the range the `Lisboa` repair was applied to, and warns
about the seam instead.

Note that for 2011-2013 this means using a superseded estimate: INE's current
national long series agrees with `0008273` for those three years, not with
`0003182`. The difference is 0.16% on the total.

### The 2021 Population Revision

INE publishes two overlapping population estimates, and they do not agree.
`0008273` (NUTS-2013) runs to 2023; `0012918` (NUTS-2024) covers 2021-2025 and
carries a revised estimate that is progressively higher on the same years:

| Portugal, resident population | `0012918` | `0008273` | |
|---|---:|---:|---:|
| 2021 | 10,599,117 | 10,421,117 | +1.71% |
| 2022 | 10,929,704 | 10,516,621 | +3.93% |
| 2023 | 11,204,347 | 10,639,726 | +5.31% |
| 2024 | 11,387,222 | not published | |
| 2025 | 11,424,031 | not published | |

Both were updated within two months of each other, so this is two published
series rather than stale data on one side. The divergence growing year on year
is the signature of revised migration estimates.

**The archive uses `0012918` for the whole of its range**, not only for the
years `0008273` lacks. Taking 2024 from the revised series while leaving 2023 on
the old one would put a 5.3% step at the seam: mortality rates would fall about
5% from 2023 to 2024 for reasons that have nothing to do with mortality, and any
forecast fitted across it would inherit the step. Using the revised series from
2021 moves the single seam to 2020/2021, where the two differ by 1.7%.

### Warning About Both Seams

Both handovers are declared in `POPULATION_SEAMS` in `R/config.R` and warned
about the same way. A rate series or a pooled window that crosses one raises a
warning naming the year and the size and direction of the effect;
`population_revision_warning()` is skipped for counts, AVPP and infant
mortality, none of which divide by this denominator.

The structural-break analysis carries a second, sharper warning. When a detected
breakpoint falls on a seam - at `year - 1` or `year`, since a breakpoint is the
last year of a segment - the narrative says so explicitly, because the breaks
tab is where a reader meets the year and would otherwise have no way to
distinguish an artefact from an event. The same note is appended to the
reliability warning shown with a guided forecast.

The practical consequence is that **every rate from 2021 onward changed** when
this was adopted. Portugal's crude mortality for 2023 reads 1,056 per 100,000
against 1,112 on the previous basis.

Before 2021 nothing was revised. INE's national long series `0001223`
(1970-2025, updated 22 June 2026) matches this archive exactly for 1995, 2000,
2005, 2010, 2014, 2016, 2018, 2019 and 2020 - to the person, +0.00% in every
case. The revision applies only from 2021, which is why the seam sits there and
why there is no revised municipal series to splice onto the earlier years.

One small exception, unrelated to the revision: **2011-2013 run 0.16% below
INE's current national figure**, because the app serves those years from
`0003182` (NUTS-2002, last updated June 2014) rather than `0008273`, which
publishes them too. That preference is deliberate - `0003182` is the contiguous
source for 1991-2013 and the one the `Lisboa` label repair was applied to, so
switching would move the source seam from 2013/2014 to 2010/2011 rather than
remove it. At a sixth of one per cent it is an order of magnitude below the
2021 step the app already warns about.

Which indicator serves which year is set by `population_source_priorities` in
`R/config.R`. Note that `get_source_year_plan()` resolves overlaps by
*descending* priority - the largest number wins, not the first listed.

### Which Years Support What

The two sides of a rate no longer end together, and it is now population that
runs ahead:

| | Range | Notes |
|---|---|---|
| Population | 1991-2025 | `0012918` publishes 2025 |
| Deaths by cause | 1980-2024 | `0013166` ends at 2024 |
| Live births | 1995-2025 | |
| Under-1 deaths | 1980-2025 | 2025 is all-cause, both sexes only |

So 2025 is selectable for infant mortality and refused for everything else, and
the year selector stops there rather than offering 2026. The observed and
forecast tabs, which work only in rates, bound their sliders by `rate_years` -
the intersection of population and deaths - rather than by either side alone.

2025 deaths do exist in `0013331` and `0013332`, but without a cause-of-death
dimension, so they cannot serve this app.

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

Pooling applies to the annual comparison tab and to avoidable mortality.
Forecasting always uses unpooled annual series, because a moving-average filter
induces autocorrelation that would invalidate the forecast intervals.

### The Pooled Denominator

A pooled rate divides by person-years, not by one year's population. Verified
against the archive by hand for Lisboa 2020-2022: 22,189 deaths over 1,792,724
person-years, which is 2.961 times the 2021 population, giving a crude rate of
1237.73 per 100,000 - exactly what the app reports, and deliberately not the
simple mean of the three annual rates (1243.07).

How the population is accumulated depends on the shape of the frame, and
getting it wrong produces a plausible number rather than an error:

- The annual tab loads **one cause at a time**, so each (year, area, band)
  appears once and the population is simply summed.
  `collapse_annual_cause_data()` now refuses a frame holding more than one
  cause, because summing there would count the same population once per cause.
- The avoidable tab loads **forty causes at once**, so the population repeats
  once per cause and must be deduplicated - taken once per (year, area, band)
  and then summed. Taking it once per band alone was wrong for every pooled
  window, and made a five-year selection for Beja read 2,173.7 against a true
  417.9.

Every pooled metric was audited against its unpooled value across several
geographies. The ratio of a pooled rate to a single-year rate sits near 1 in
every case; a denominator error of this class would show up as a ratio near the
window length. `Óbitos infantis` is the one deliberate exception, being a window
total rather than an annual average.

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

### Infant Mortality: Sources Keyed By Code

Two defects in the infant datasets were corrected in September 2026.

**Births were matched by label.** In `0000003` (NUTS-2002, 1995-2013) the label
`Lisboa` names both the NUTS II region and the municipality, and the fetcher
summed by label: the municipality received 37,208 births in 2001 against a true
5,604, deflating the infant rate of Lisboa and of every region or ULS containing
it. The same summing merged the two `Calheta` and the two `Lagoa`
municipalities. `tools/fetch_births.R` now maps every row by geography code: the
last four digits of a municipal code are its DICO, identical in every NUTS
vintage. A year is rejected unless all 308 municipalities map. After the fix the
municipal sum equals Portugal to within 20 births in every year.

**Under-1 deaths were incomplete at municipal level.** `data/snapshots/infant_deaths`
takes the `Menos de 1 ano` band of the cause-of-death indicators, whose municipal
age breakdown is incomplete: in 2014 the municipalities add up to 112 against a
national 236. INE also publishes under-1 deaths as a subject of their own
(`0008181`, NUTS-2013, from 2011; `0012541`, NUTS-2024, from 2022), by
municipality, sex and age in days and months, and those close exactly against
the national total in every year. `tools/fetch_infant_totals.R` writes them to
`data/snapshots/infant_totals`, keyed by code; all-cause requests read them
wherever a year has them. Cause-specific requests, which they cannot answer,
still read the band-derived dataset, and so does the AVPP correction, which
splits the `0 - 4 anos` band of that same breakdown.

Before 2011 no complete municipal count exists. The band-derived municipal sum
reaches 99-100% of Portugal in 2002-2010 but only about 85% in 1995-2001.
`infant_undercount_years()` identifies years below 97%, the annual tab warns
when a non-national selection touches them, and the planning tab marks the
value with `†`. Portugal reads its own published row and is unaffected.

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

## Planning Indicators Tab

The `Indicadores de Planeamento` tab reproduces the demographic and mortality
indicators of the DRS/PNS2030 support workbook for local health plans
(`Indicadores_Apoio_PLS`), computed from the app's own snapshots for any
location. The engine is `R/planning_indicators.R`.

| Indicator | Workbook | Definition |
|---|---|---|
| Resident population | I1 | INE annual estimate |
| Share aged 0-14, 65+, 75+ | I1 | group / total x 100 |
| Ageing index | I4 | pop 65+ / pop 0-14 x 100 |
| Youth dependency | I5 | pop 0-14 / pop 15-64 x 100 |
| Old-age dependency | I6 | pop 65+ / pop 15-64 x 100 |
| Live births | I7 | count |
| Crude birth rate | I8 | births / mean population x 1,000 |
| Deaths | I37 | count, all causes and ages |
| Crude death rate | I38 | deaths / mean population x 1,000 |
| Life expectancy at birth and at 65, total and by sex | I10 | abridged life table per triennium (below) |
| Total fertility rate | I9 | sum over mother's age 15-49 of births / women x 5 |
| Births to mothers under 20, 35+ | I32, I33 | share of live births, three years |
| Preterm births | I35 | under 37 weeks / births of known duration x 100, three years |
| Low birth weight | I36 | under 2,500 g / births of known weight x 100, three years |
| Census population and change | I2 | census population; change against the previous census |
| Education level | I24 | share of the census population per completed level |
| Illiteracy rate | I26 | illiterate aged 10+ / population aged 10+ x 100 |
| RSI beneficiaries | I13, I14 | count; per 1,000 residents aged 15+ (mean population) |
| Social security pensioners | I15, I16 | count; per 1,000 residents aged 15+ (31 December) |
| Mean pension | I17 | sum(pensioners x mean) / sum(pensioners) |
| Purchasing power per capita | I28 | sum(share) / sum(share / index) x 100 |
| Urban waste per inhabitant | I64, I65 | tonnes x 1,000 / mean population |
| Average monthly earnings | I27 | sum(earnings x employees) / sum(employees) |
| Employees and sector shares | I12 | count; share of employees per sector |
| Infant mortality | I39 | under-1 deaths / live births x 1,000, pooled over three years |
| Neonatal, early neonatal, post-neonatal | I40-I42 | deaths <28 d, <7 d, 28-364 d / live births x 1,000, three years |
| Late fetal, perinatal mortality | I43, I44 | stillbirths 28+ weeks, and those plus deaths <7 d, / (live births + stillbirths) x 1,000, three years |
| Proportional mortality | I45, I46 | deaths per large cause group / all deaths x 100, three years, all ages and under 75 |
| Population pyramid | I3 | share by five-year band and sex |

### Life expectancy at birth (I10)

`R/life_expectancy.R` builds an abridged period life table per area, sex and
triennium, with the formulas and variance of
`PHEindicatormethods::phe_life_expectancy()` (Chiang II; Silcocks' variance for
the open interval; suppression when person-years are 5,000 or fewer or the 95%
interval exceeds 20 years). A test reproduces that function exactly on its own
age structure. Two adaptations to the app's data:

- **First band 0-4.** The death archive has no under-1 band with a matching
  population, so 0-4 is one interval whose `a` (fraction lived by those dying)
  combines 0.1 years for infant deaths and 2.5 years for deaths at 1-4, weighted
  by the complete under-1 counts.
- **Unrecorded ages.** INE's municipal breakdown of all-cause deaths by age
  misses some deaths (4,130 in 2014, about 0.5% in 2013, 2015 and 2024). Each
  municipality's missing deaths, the difference to its complete total, are spread
  over its ages in proportion to its recorded ones (the national profile if none
  are recorded). A value is marked `‡` when this exceeds 2% of the triennium's
  deaths; that happens only at municipal level, chiefly for triennia containing
  2014.

Deaths are pooled over the three years; person-years are the sum of mid-year
populations, each the mean of consecutive end-of-year estimates.

Both ages come from one pass of the table: `abridged_life_table()` returns life
expectancy and its standard error for every band, and each indicator reads the
band it needs, with PHE's suppression rule applied at that age.

**Validation.** For Portugal the tables agree with Eurostat's (`demo_mlexpec`):
2017-2019 gives 81.9 at birth and 20.7 at 65 (Eurostat 2019: 82.0 and 20.6);
men 78.8 / 18.6 (79.0 / 18.7); women 84.9 / 22.4 (84.8 / 22.3). INE's published
tables (Metodologia 2007: `0001724` for Portugal, `0013473` and `0008459` for
NUTS III) are systematically lower. With INE's pre-revision population restored
for 2021-2023, the 26 NUTS III of the current edition differ by +0.80 (2020-2022)
and +0.87 (2021-2023) years, standard deviation 0.3, correlation 0.97, rank
correlation 0.94-0.95. The offset is not the open interval (a Gompertz extension
beyond 85 would raise, not lower, the values) nor the population revision; it is
INE's method. App values are comparable with each other, not with INE's.

### Proportional mortality under 75 (I46)

I45 reads the complete all-ages totals; I46 needs the age breakdown, so each
area takes the better of two sources (`R/planning_under75.R`): INE's own
regional row where one exists - Portugal, Continente, every NUTS region, and the
ULS that coincide with a NUTS III unit, matched by identical membership so a
territory published under one of its two names still counts - and otherwise the
sum of its municipalities. An area is built entirely from one source: if any of
the year's cause files is missing from INE's rows, the whole area falls back to
municipal sums, so its shares always add up.

How good the municipal sum is, measured against INE's rows: for 2020-2022 it
reproduces the under-75 deaths of Alto Minho and Algarve exactly and their
shares to within 0.08 and 0.29 points; Portugal's shares are within 0.21 points.
The exception is 2014, where INE published an age for only 79.7% of municipal
deaths (52.7% for the worst cause group); Alto Minho's 2012-2014 shares are then
off by up to 2.2 points. Areas not on INE rows are marked in the three triennia
containing 2014.

Deaths are deliberately not rescaled to the complete totals. Measured against
INE's rows, rescaling makes the shares worse: deaths with no published age sit
mostly at older ages, so spreading them proportionally moves too many below 75
(Portugal 2020-2022, malignant tumours: 38.69% exact, 38.48% from municipal
sums, 39.20% rescaled). The same bias applies, mildly, to the life-expectancy
repair, which affects only flagged municipalities.

All 381 areas and 32 triennia take about 30 seconds, and the all-areas Excel
file about 85 seconds in total.

### Socio-economic, birth and neonatal sources

`tools/fetch_planning_extra.R` writes every measure as additive municipal
components to `data/snapshots/planning_extra/<measure>/year_<year>.rds`. Each
measure is published in up to three editions, one per NUTS vintage; for every
year the newest edition covering it wins, and rows are mapped by DICO.

| Measure | Editions (oldest first) | Years |
|---|---|---|
| RSI beneficiaries | `0004299`, `0008251`, `0013417` | 2007-2025 |
| Pensioners | `0004294`, `0010271`, `0013395`, `0014534` (Série 2017) | 2004-2025 |
| Mean pension | `0004149`, `0010266`, `0013398`, `0014532` (Série 2017) | 2004-2025 |
| Urban waste collected (t, by collection type) | `0000482`, `0009612`, `0012769` | 1995-2024 |
| Purchasing power per capita / share | `0001354`+`0001355`, `0008614`+`0008615`, `0014580`+`0014581` | biennial, 1993-2023 |
| Births by mother's age | `0005952`, `0008092`, `0012441` | 1995-2025 |
| Births by gestation | `0005950`, `0008084`, `0012434` | 1995-2025 |
| Under-1 deaths by age | `0008181`, `0012541` | 2011-2025 |
| Births by weight | `0005611`, `0008088`, `0012438` | 1995-2025 |
| Perinatal deaths | `0003527`, `0008173`, `0012549` | 1995-2025 |
| Census population | `0014353` | 1991, 2001, 2011, 2021 |
| Census population by age | `0014164` | 1991, 2001, 2011, 2021 |
| Census education level | `0014380` | 1991, 2001, 2011, 2021 |
| Census illiteracy rate | `0014375` | 1991, 2001, 2011, 2021 |
| Average monthly earnings | `0009047`, `0012656` | 2011-2024 |
| Employees by sector | `0010378`, `0012648` | 2013-2024 |

Every dimension other than area and the measure's own is pinned to its total
category; a response whose other dimension has no total is refused rather than
summed. Birth indicators pad their labels with Unicode spaces, which are trimmed
before matching.

Non-additive published figures are rebuilt from additive parts. A mean pension
is kept as pensioners and pensioners x mean. Purchasing power per capita is an
index with Portugal = 100: a municipality's share of national purchasing power
divided by its index is its implied share of population, so an area's index is
the sum of shares over the sum of implied population shares, times 100 -
reproducing INE's published Continente (100.63) and Norte (92.90) exactly.

The mother's-age dimension carries overlapping categories together: single
years, five-year groups, a 15-49 group, and `50 - 54`, `50 e mais` and
`55 e mais` side by side. Only the five-year groups below the lowest open group
and that open group are read, so no birth is counted twice. For the fertility
index, births to mothers under 15 are counted in 15-19 and those of 50 and over
in 45-49; births of unknown mother's age are left out. Women are counted at
mid-year, as the mean of the end-of-year estimates of the previous year and the
current one. With that denominator the index reproduces INE's published series
for Portugal (`0001293`) to two decimals in 2018-2020 and 2022-2025; 2021 reads
1.32 against 1.30 because its mid-year mean straddles the population revision.
The end-of-year estimate alone reads up to 0.03 low.

Earnings and employees come from Quadros de Pessoal (MTSSS/GEP): employees
only, counted at their workplace. The workbook's I12 uses census employment by
residence, a different universe, and no census series of employment by sector
exists per municipality for 2011. An area's mean earnings weight the municipal
means by employees; for a single-municipality ULS the result equals the
workbook's exactly (ULS Matosinhos 2013-2018). For a ULS of several
municipalities the workbook repeats the first municipality alphabetically - its
ULS Alto Minho series is Arcos de Valdevez's, year for year - which the app does
not reproduce.

The census series code municipalities as one character plus the DICO (five
characters) rather than the usual seven, which the fetcher handles. INE
publishes the illiteracy rate per municipality but not the count, so an area's
rate is its municipalities' rates weighted by population aged 10 and over -
equal to the ratio of the implied counts. The education series counts only
people with a completed level, so "no level completed" is the census population
less that total. Validated against workbook v26: Continente census population
exact for 2011 and 2021, illiteracy 5.19% and 3.04% (workbook 5.187, 3.047),
education shares within 0.05 points. Regional rows differ (Norte 2011:
3,689,682 here, the published Census figure, against 3,737,768 in the workbook),
which points to an older region definition there.

INE publishes perinatal deaths per municipality but not stillbirths on their
own, so stillbirths of 28 or more weeks are the perinatal deaths less the deaths
under 7 days of the same municipality, which are known from 2011. I43 and I44
therefore start with the 2011-2013 triennium, and both divide by live births
plus those stillbirths, as INE defines the rates.

The workbook's RSI and pensioner rates divide by residents aged 15 and over (its
2021 denominators match the 15+ population of the superseded estimate to within
0.02%), and so do INE's own indicators (`0013420`, `0014599`), which the app
matches to about 1% per municipality in 2024. An earlier version of this note
said INE used 15-64; that was wrong.

Pensions change series in 2017 (Série 1990-2023 to Série 2017, about 5.5% fewer
pensioners); `PLANNING_SERIES_BREAKS` records it and the evolution chart marks
it.

### Location, comparators and export

The tab is built around one location. `planning_comparators()` proposes, for
that location, one containing area per level above it - ULS, ARS, NUTS III,
NUTS II, NUTS I, Portugal - by testing which units' municipalities include all
of the location's. The two geographies interleave (a ULS can lie inside a NUTS
III, a NUTS III inside an ARS), so levels are ranked Município < ULS < NUTS III
< ARS = NUTS II < NUTS I < Portugal. A containing unit with exactly the
location's municipalities, or a nearer comparator's, is dropped as it would
repeat the same values; Portugal is always kept, since it reads INE's national
row. Comparators are offered only for indicators flagged `comparable` in
`PLANNING_INDICATORS` (rates, shares, indices, per-capita values); counts are
shown for the location alone.

Components are summed with a membership matrix (areas x municipalities), one
matrix product per year, and indicators are computed for all areas and years in
one vectorised pass, with exact Poisson and Clopper-Pearson intervals from their
closed forms (identical to `poisson.test()` and `binom.test()`). All 381 areas
x 35 years x 27 indicators take about 11 seconds, which is what allows the
all-areas Excel file (`R/planning_export.R`) to be built on demand; it is cached
per NUTS vintage and data import date.

Chart forms (`R/planning_charts.R`): lines for comparable indicators, with the
location's 95% interval as a band; bars with error bars for counts; sorted bars
for the ULS ranking; a share-based pyramid with the comparator as an outline;
bars and dots for cause groups. Each level has a fixed colour from the reference
categorical palette, in its validated order; direct labels are drawn at line
ends only up to four series, beyond which the legend carries identity.

### ULS that share a municipality: whole or parish weights

Lisboa, Loures and Porto are divided between ULS at parish level, so six ULS
cannot be built from whole municipalities. `data-raw/uls_parish.csv` records the
assignment (Decreto-Lei n.º 102/2023, which defines each ULS by the ACES it
integrates, plus those ACES' parish lists); `tools/build_uls_parish.R` matches
it to INE's parish names and fails unless each municipality is covered exactly
once. The decree puts three Lisboa parishes (Ajuda, Alcântara, Belém) in ULS
Lisboa Ocidental, which the support workbook does not record;
`tools/build_uls_lookup.R` now adds that pair before deriving the exact groups,
so the Lisboa group is `ULS Lisboa Ocidental + Loures/Odivelas + Santa Maria +
São José` (six municipalities). The lookup keeps those ULS individually as
`ULS (partilhada)`, beside the groups.

`R/planning_parish.R` provides the two readings (`PLANNING_SPLIT_MODES`):
`whole`, where each ULS takes every municipality it serves (nothing estimated,
but the six overlap: in 2024 they sum to 18,186 deaths and 1.8 million people
more than the Continente), and `parish`, where each shared municipality is
divided by the parishes' census share. Births and deaths are not estimated: INE publishes them by parish every year
(`tools/fetch_parish_vitals.R`; the current parishes match from 2014, and the
parishes of a municipality sum exactly to its total), so each ULS takes its
parishes' own counts, year by year - in 2024 ULS São José holds 50.2% of
Lisboa's deaths and 61.3% of Loures'. Population has no annual parish figure,
so it uses the 2021 census share by age group (`tools/fetch_census_parish.R`,
`0012364`), as does everything derived from population, and as do births and
deaths before 2014 (women 15-49 and an expected-deaths share respectively).
For quantities carried by age band, the census gives the shape over ages and
the registers the level: `planning_band_product()` scales each municipality's
band weights so they sum to that year's registered share (raking), which makes
the six ULS reproduce the parish registers exactly. `planning_membership_matrix()` takes a mode
and a basis and returns fractional weights; `planning_components()` multiplies
each column with the matrix of its own basis, and `planning_band_product()`
does the same per age band for life expectancy and the standardised module. In
the parish reading the 39 ULS sum exactly to the Continente in every component.

### Portugal benchmark, significance and funnel

`PLANNING_PORTUGAL_MUNICIPAL` ("Portugal (soma dos municípios)") is a pseudo-area
with every municipality as members and no published row, so every path
(components, life table, I45, I46) builds it as a municipal sum; the I46 alias
lookup is told never to swap it for INE's national row. The tab's `Portugal`
option swaps it in for `Portugal` in the comparators, the significance
benchmark, the ranking, the funnel and the profile. The all-areas export
carries both rows.

Significance (`planning_add_significance()`) follows PHE Fingertips: higher or
lower when the whole 95% interval lies above or below the benchmark's value in
the same period, similar otherwise. It is computed only for `comparable`
indicators with an interval; counts are never compared. Colours in the ranking
are a diverging pair (orange above, blue below, grey similar) validated for
colour-vision deficiency (worst adjacent OKLab distance 16.8); orange rather
than red because the desirable direction depends on the indicator.

The funnel (`planning_funnel_limits()`) uses exact quantiles of the count under
the benchmark rate - Poisson for event rates, binomial for the birth
proportions - interpolated between integers (Spiegelhalter, Stat Med 2005), at
95% and 99.8%. A first version used the exact interval around the expected
count, which flagged 136 of 308 municipalities as significantly low on infant
mortality because a zero count in a small unit fell below it; the quantile
method puts 275 inside the limits. Only indicators in `PLANNING_FUNNEL_MODELS`
have a funnel.

### Education by age

`census_education_by_age` (`tools/fetch_census_education_age.R`) holds the 2011
(`0006350`) and 2021 (`0012364`) census population by five-year age group and
highest completed level, per municipality; post-secondary counts as secondary,
as in the historical series. The 2021 table codes municipalities by the bare
DICO (four characters) and parishes by six. With `education_min_age > 0` each
level is divided by everyone of that age and over; 1991 and 2001 have no age
breakdown and stay empty. With 0 (the default, as in the workbook) the
historical series `0014380` is used unchanged.

### Standardised, premature and avoidable mortality

`R/planning_standardised.R` adds, per triennium and area: the SMR against the
benchmark (indirect, Portugal = 100, exact Poisson interval on the observed
count), the directly standardised rate (ESP-2013, Dobson interval with exact
Poisson limits - PHEindicatormethods uses Byar's approximation; they agree to
0.002% at 120 deaths and exactly at regional counts), the same under 75
(premature), the preventable and treatable rates under 75 (the lists of
`R/avoidable.R`, a lower bound), premature and avoidable death counts, and years
of potential life lost before 70 per 100,000 residents under 70 (band
midpoints; infant deaths weighted at 69.5 years). `planning_cause_standardised()`
gives observed, expected, SMR and both standardised rates for the 13 I45 groups
and all causes, by sex.

Deaths by age and cause per municipality are incomplete at INE (97.9% of
circulatory deaths in the municipal bands in 2023; 65% of suicides in 2014),
while the national row and the municipal all-ages totals are complete. Each
municipality and cause is completed to its total (death_totals), and the
missing deaths are spread with the national gap profile - Portugal's row by age
less the municipal sum - rather than the municipality's own profile. The gap is
concentrated at young ages (INE suppresses small cells), and the own-profile
version read 2012-2014 premature mortality of the municipal sum at 338.8
against Portugal's 350.5; with the gap profile the two agree (350.5), and Norte
2022-2024 deaths under 75 reproduce INE's regional row exactly (33,155).
Alentejo 2014-2016 is 2% off INE's composed row. Values with more than 2% of
their deaths spread carry `‡`. The regional rows INE publishes by age could
replace the municipal sums for NUTS areas; not done yet.

### Primary care from the SNS Transparency portal

`tools/fetch_sns.R` exports five datasets of transparencia.sns.gov.pt
(Opendatasoft Explore API v2.1, `/exports/json`, no key) whole, in long form,
through the versioning layer. `R/planning_sns.R` keeps January 2024 onwards (ULS
units; the ACES before do not map onto ULS), maps the portal's labels onto the
app's ULS (`SNS_UNIT_RENAMES`: spelling variants, and the five Lisboa/Porto ULS
into the two exact groups), and rebuilds every proportion as numerator and
denominator (the denominator from count / proportion when only those are
published) so ARS and Continente are exact sums. Clopper-Pearson intervals;
significance against the Continente.

Several indicators accumulate and reset. From the published series: blood
pressure and HbA1c control rise from about 10% in January to 50-65% in June and
restart in July; the foot exam and the three screenings rise from January to
December. So `SNS_INDICATORS$cycle` marks each as month, semester or year, and
only ends of cycle are compared; the chart breaks the line at each reset. A
field the portal renames leaves its indicator out with a warning. An area is
shown only when it is a union of whole ULS (or a municipality, which shows its
ULS; Portugal shows the Continente).

### Weekly deaths and excess mortality

`tools/fetch_weekly_deaths.R` fetches `0012100` (NUTS 2024, 2021-) and `0010112`
(NUTS 2013, 2018-2024), both sexes, all weeks in one request each. Older-edition
years are kept only for regions whose annual counts agree between the editions
within 0.5% in 2021-2024 (not Centro, Alentejo, Médio Tejo, Beira Baixa).
`planning_weekly_excess()` applies the mean age-specific weekly rate of the
baseline years (under 65, 65-74, 75-84, 85+; population from the annual
estimates, the latest available for the current year) to the year's
population; the 95% prediction band uses the between-year variance of the
baseline rates, at least Poisson, times (1 + 1/n). Baseline: up to five years,
from 2023 only. A first version used 2018-2019 as well; over the superseded
population estimate their rates were 4-7% higher, and Portugal showed a
spurious deficit of 4-7% in every year from 2023. Week 53 borrows week 52.
Deaths of unknown age (0.014%) are left out on both sides.

### Location profile

`write_planning_profile()` (`R/planning_profile.R`) builds a .docx with
officer/flextable: the indicators significantly above and below the benchmark,
the full table on landscape pages, the pyramid, six trends as small multiples,
the location's ULS ranked on nine event rates (the smallest ULS containing the
location, even when it is not offered as a comparator for having the same
municipalities), I45 and the notes.

### Aggregation

Every area is a ratio of sums. The components (population by broad age group,
births, deaths, infant deaths) are summed over the area's municipalities, and
the indicator is computed from those sums, never as an average of municipal
indicators. Portugal and Continente use INE's published rows where the dataset
has them, because those include events whose municipality of residence is
unknown; any other area is its municipal sum and so excludes them (0.3-0.9% of
deaths). ULS and ARS use the membership in `data/uls_lookup.rds`, applied to
every year.

For each year the engine builds one compact table of components per area label
present in the files (about 310 rows) and caches it for the session, so a
ranking of all ULS or a 35-year series is a sum over that table.

### Mean or end-of-year population

Events counted over a year - births, deaths, RSI beneficiaries, waste - are
divided by the year's mean population, the mean of the end-of-year estimates of
the previous year and of this one. That is INE's convention and the workbook's:
with it the app reproduces INE's 2024 crude birth and death rates and waste per
inhabitant in all 308 municipalities (to the published decimal), and the
workbook's Continente I8, I38 and I14 of 2015 and 2019 exactly. The end-of-year
estimate alone read about 0.5-0.8% low in 2024, while the population grew.
Stocks counted on 31 December (pensioners, the age structure) use that day's
estimate. The first year of the series has no previous estimate and uses its
own. Life expectancy and the fertility index already used mid-year population.

### Gaps in INE's municipal data

A blank cell at INE is not a zero, and older fetches stored blanks as 0. The
audit of 2026-09-21 found four kinds, each now handled explicitly:

- **Complete datasets.** Waste, pensions, employees and earnings are published
  for every municipality (`PLANNING_COMPLETE_BLOCKS`). A municipality that has a
  population but no value there makes every area containing it missing, rather
  than dropping out of the sum: the Azores have no employees in 2013-2014 and
  no earnings in 2011-2014; 1995-2006 waste lacks many island municipalities;
  four municipalities lack pensions in 2005. Before the fix, those regions read
  as zero employees and their NUTS II totals were understated.
- **Joint reporting** (`PLANNING_JOINT_REPORTING`). Odivelas, Trofa and Vizela
  were created in 1998 out of Loures, Santo Tirso and Guimarães. INE's
  population is back-cast to today's boundaries, but births (to 1997, 1998
  partly), deaths (to 1998) and the 1991 census sit with the parent. Until 1998
  every area holding one of a pair without the other is missing, for every
  dataset but the population - including life expectancy and I45/I46 in
  triennia that include those years. Before the fix Odivelas had a life
  expectancy of 88.7 and Trofa 92.4 in 1997-1999, and Loures a birth rate over a
  population without Odivelas. The Loures-Odivelas waste service (SIMAR) is
  recorded wholly under Loures, with Odivelas blank, in every year: waste per
  inhabitant exists only for areas holding both (INE shows Loures 402 kg, the
  joint figure, and Odivelas blank; the app had shown 708 and 0).
- **Blank all-ages death totals.** INE leaves the all-ages cell blank for a few
  municipality-years while filling in the age bands: Vimioso 2024,
  Alfândega da Fé and Miranda do Douro 2015, and single municipalities in
  1993-1999. `repair_death_totals()` takes the larger of the published total and
  the sum of its bands for all causes and the I45 groups, and marks the cell
  `repaired`. Vimioso's 2024 death rate goes from 0 to 23.2, INE's value; ULS
  Nordeste gains 87 deaths (4%), which had inflated its proportional mortality.
  `tools/fetch_death_totals.R` now keeps blanks as `NA`.
- **Suppressed employment sectors.** INE hides two of the three sectors when one
  would reveal a single employer (6 municipalities in 2024, 41-49 a year in
  2013-2016: Vizela, Marinha Grande, Boticas...). Read as zeros, all of such a
  municipality's employees fell in its one published sector and its shares
  summed to as little as 29%. `planning_estimate_suppressed_sectors()` splits
  the hidden remainder between the hidden sectors in the proportions of the rest
  of the NUTS III that year; shares where more than 1% of the area's employees
  were estimated carry the flag `≈`.

### Why deaths come from the all-ages totals

INE's municipal breakdown of cause-specific deaths by age is incomplete, worst
in 2014. The all-ages total per municipality and cause is complete. Counts,
crude rates and all-ages proportions need no age, so this tab reads
`data/snapshots/death_totals` (`tools/fetch_death_totals.R`): the age `Total`
row of `0008206` (1991-2021) and `0013166` (2022-2024), mapped to municipalities
by geography code. Municipal sums close against the national row to within the
deaths of unknown residence in every year.

### Uncertainty

Counts and event rates carry exact Poisson intervals on the count, the
denominator treated as fixed; proportional mortality carries exact binomial
intervals. Population-structure indices carry none: population estimates are
not a sample of events and INE publishes no error for them. A triennial infant
rate on fewer than 1,000 births is marked `*`, as elsewhere in the app.

### Proportional mortality groups

The 13 groups of the workbook's I45, in INE's shortlist wording: infectious and
parasitic diseases, malignant neoplasms, blood, endocrine, nervous system,
circulatory, respiratory, digestive, musculoskeletal, genitourinary, perinatal,
ill-defined, external causes. Each is a chapter-level rubric, so none contains
another. What they leave out (mental disorders, skin, pregnancy, congenital
malformations) is shown as `Restantes causas`, so the column sums to 100%.

### Differences from the workbook

- **Population revision.** From 2021 the app uses the revised series
  `0012918`. Indicators built on the superseded estimate are out of date: the
  ageing index of Alto Minho in 2024 is 270.3 on the old estimate and 240.2 on
  the revised one.
- **ULS sharing a municipality.** Lisboa, Loures and Porto are split between ULS
  at parish level. The workbook assigns each whole municipality to every ULS
  that serves part of it, so the ULS of ARS Norte add up to about 3,000 more
  deaths than ARS Norte itself. The app offers the two exact groups instead.
- **Population denominators of 2024.** Every population-based indicator
  differs from the workbook by one factor per area (Continente 1.067 for birth,
  death, waste and RSI rates alike; Algarve 1.18), which is the population
  revision. Indicators without a population denominator agree: mean pension
  2024 Continente 7,697 vs 7,696.62; pensioners and RSI beneficiaries 2021 exact.
- **Fertility index.** The workbook's 2024 values (Continente 1.41) predate
  the population revision; INE's own revised figure for Portugal is 1.27, which
  the app reproduces.
- **Births to mothers under 20 (I32).** The workbook counts mothers aged 15-19
  only (Continente 2022-2024: 1.82%); the app follows the indicator's title and
  includes mothers aged 10-14 (1.87%).
- **Purchasing power and neonatal mortality of some ULS.** The workbook's ULS
  values disagree with INE's municipal values: Matosinhos (a single municipality)
  2021 purchasing power is 118.06 at INE and 130.57 in the workbook; Alto Minho
  is 55.76, neither the weighted index (82.2) nor the plain mean (77.3). Its
  neonatal and post-neonatal rates do not add up to its infant rate (Alto Minho
  2022-2024: 1.1 + 0.9 vs 2.4), which suggests misaligned rows. In the app the
  two always add up.
- **Edition of 2022.** The app reads 2022 from `0013166` (NUTS-2024). For
  Continente 2020-2022 it gets 356,355 deaths against the workbook's 356,333
  (+22, 0.006%); circulatory (95,307), genitourinary (11,855) and perinatal
  (342) agree exactly. The residual is consistent with a revision between the
  two INE tables for 2022.
- **Naming.** The workbook spells some units inconsistently between sheets
  (`Loures-Odivelas` / `Loures/Odivelas`, `Almada-Seixal` / `Almada/Seixal`),
  which breaks lookups by name. The app keys every unit by municipality code.

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

## Data Versions

INE revises published figures - provisional years become final, population
series are re-estimated, indicators are replaced by new editions - so an
analysis repeated later can legitimately give different numbers. The app keeps
the history needed to tell a revision from an error (`R/data_versions.R`).

**Writing.** Every tool writes data through `versioned_save_rds()`:

- content identical to the stored file (compared independently of row order) is
  not rewritten and leaves no log entry;
- a new file is written and logged as `added`;
- a different file is first copied to `data/archive/<run id>/<relative path>`,
  then replaced, and logged as `replaced`.

**The log** (`data/import_log.csv`) has one row per write: run id, import
timestamp, relative path, dataset, year, action, row counts, value column, the
comparable total before and after, the number of rows whose value changed, the
archived location, the tool and a note. The comparable total is Portugal, both
sexes, all causes and the `Total` category, whichever of these the file carries;
summing every row would count regional and sex rows, which differ between INE
editions even when no value does (the population revision reads +1.7%, +3.9%
and +5.3% for 2021-2023 on this total, and over +100% on a sum of all rows).
Death chunks are one file per cause and causes nest, so a year's death total is
taken from its all-cause file. A correction between
municipalities leaves Portugal's total unchanged, so the log also counts
changed rows, matched on every column except the value and the source
indicator.

**Runs.** `tools/refresh_snapshots.R` sets one run id (`DATA_RUN_ID`) for every
tool it launches, re-fetches the last `recent` calendar years (default 2) as
well as missing ones, and appends to `REFRESH_STATUS.md` a table of what the run
revised.

**Before the log.** `tools/backfill_import_log.R` built the log from git: each
file's import date is the commit that gave it its form (files assembled over
many commits in May 2026 take the last of them), and every change from
2026-08-01 is a `replaced` entry whose previous version is referenced as
`git:<parent commit>:<path>` rather than copied. That covers the
Lisboa/Calheta/Lagoa repair (2,146 death chunks, 2026-08-11), the adoption of
the revised population series (2026-08-20), the 0013166 refetch (2026-08-20),
the regional-row refetch (2026-09-16) and the births fix (2026-09-17).

**Exports.** Every CSV ends with a comment line, and every PNG carries a
caption, naming the data import date, the NUTS vintage in force and the export
date (`export_stamp()` in `R/helpers.R`). The Excel files state the same on
their read-me sheet. The comment line sits after the data and starts with `#`,
so the file still parses.

**Reconstruction.** `data_as_of(date)` (`tools/data_as_of.R`) rebuilds the data
directory as it stood at the end of a day: files unchanged since are
hard-linked, files replaced later are restored from the version their first
later replacement archived (from `data/archive` or git), and files first added
later are omitted. Checked against the history: on 2026-09-01 Lisboa has 37,208
births in 2001 (before the fix) and Portugal's 2022 population is 10,929,704
(revised); on 2026-08-15 births are absent and the 2022 population is 10,516,621
(before the revision).

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
