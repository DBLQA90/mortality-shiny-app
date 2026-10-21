# Snapshot refresh status

Last run: 2026-10-21 10:07:19 WEST
Task: `current`  |  Duration: 5.3 min

```
[10:07:19] Refresh started; task=current, budget=120 min
[10:07:19] == Task: weekly deaths ==
[10:07:19] -> weekly deaths 0012100 / 0010112
   0012100: 296100 rows, 2021-W01 to 2026-W40, 35 regions
   0010112: 337260 rows, 2018-W01 to 2024-W52, 33 regions
   Done. Failed: 0
[10:09:50]    done
[10:09:50] == Task: births and deaths by parish ==
[10:09:50] -> parish births and deaths
     births 0012450: 16980 rows, 2021-2025, 3396 parishes
     births 0008234: 30940 rows, 2014-2023, 3094 parishes
     births 0006239: 93366 rows, 1996-2013, 5187 parishes
   births: 132004 rows, 1996-2025, 5489 parishes
     deaths 0012542: 16980 rows, 2021-2025, 3396 parishes
     deaths 0008235: 30940 rows, 2014-2023, 3094 parishes
     deaths 0005604: 93366 rows, 1996-2013, 5187 parishes
   deaths: 132004 rows, 1996-2025, 5489 parishes
   Done. Failed: 0
[10:11:14]    done
[10:11:14] == Task: calendar-year deaths by age (weekly + parish) ==
[10:11:14] -> calendar deaths by age
     2025: Portugal 121840, Continente 116499 deaths
   Done.
[10:11:15]    done
[10:11:15] == Task: SNS Transparency portal ==
[10:11:15] -> SNS primary care
   utentes-inscritos-em-cuidados-de-saude-primarios: 6450 rows, 2016-01 to 2026-06, 96 units
   rastreios-oncologicos: 7184 rows, 2015-01 to 2026-08, 77 units
   diabetes: 7192 rows, 2015-01 to 2026-08, 77 units
   hipertensao: 7184 rows, 2015-01 to 2026-08, 77 units
   saude-da-mulher-e-crianca: 7192 rows, 2015-01 to 2026-08, 77 units
   Done. Failed: 0
[10:11:58]    done
[10:11:58] == Task: deaths by cause and age ==
[10:11:59] Nothing to fetch or re-check for 0013166.
[10:11:59] == Task: population ==
[10:11:59] -> population 0012918 2025
   Fetching 0012918 for 2025
     2025: 340 areas, Portugal = 11 424 031
   Done: 1 year files written.
[10:12:03]    done
[10:12:03] == Task: municipal death totals ==
[10:12:03] -> death totals
   Done: 0 written, 0 failed.
[10:12:05]    done
[10:12:05] == Task: births and under-1 deaths ==
[10:12:05] -> live births
   
   == 0000003: 19 year(s) ==
   
   == 0008084: 7 year(s) ==
   
   == 0012434: 5 year(s) ==
     2025: ok (310 areas, Portugal 87764)
   
   Done: 1 written, 30 already present, 0 failed.
[10:12:12]    done
[10:12:12] -> under-1 deaths by cause
   Fetching under-1 deaths for 46 year(s), one request each
     2025 come from 0012540: all causes, both sexes only.
     2025 (0012540, total only): ok - 340 rows, Portugal all-cause 246
   Done: 1 written, 45 already present, 0 failed.
[10:12:16]    done
[10:12:16] -> complete under-1 death counts
   Years: 2011-2025 (15)
     0012541 2025: 308 municipalities | municipal sum 246, Portugal 246
   Done: 1 written, 0 failed.
[10:12:34]    done
[10:12:34] Finished in 5.3 min.
```

## Changes in this run

Run id: `2026-10-21T100719`  |  Re-checked years from 2025  |  0 file(s) added, 1 replaced

Revised by INE since the previous import (previous versions in `data/archive/2026-10-21T100719/`):

Totals are Portugal, both sexes, all causes where the files carry them; a
correction between municipalities shows in the changed rows, not the total.

| Dataset | Year | Files | Rows changed | Total before | Total after | Change |
|---|---|---|---|---|---|---|
| weekly_deaths/0012100 | NA | 1 | 4 609 | 6 917 130 | 7 000 080 | +1.20% |
