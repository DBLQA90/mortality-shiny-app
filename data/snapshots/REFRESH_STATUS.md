# Snapshot refresh status

Last run: 2026-10-26 09:38:24 WET
Task: `current`  |  Duration: 29.1 min

```
[09:38:24] Refresh started; task=current, budget=120 min
[09:38:24] == Task: weekly deaths ==
[09:38:25] -> weekly deaths 0012100 / 0010112
   0012100: 296100 rows, 2021-W01 to 2026-W40, 35 regions
   0010112: 337260 rows, 2018-W01 to 2024-W52, 33 regions
   Done. Failed: 0
[09:40:53]    done
[09:40:53] == Task: births and deaths by parish ==
[09:40:53] -> parish births and deaths
   Failed to perform HTTP request.
   Caused by error in `curl::curl_fetch_memory()`:
   ! Timeout was reached [www.ine.pt]:
   Operation too slow. Less than 1 bytes/sec transferred the last 600 seconds
     births 0012450: 16980 rows, 2021-2025, 3396 parishes
     births 0008234: 30940 rows, 2014-2023, 3094 parishes
     births 0006239: 93366 rows, 1996-2013, 5187 parishes
   births: 132004 rows, 1996-2025, 5489 parishes
     deaths 0012542: 16980 rows, 2021-2025, 3396 parishes
     deaths 0008235: 30940 rows, 2014-2023, 3094 parishes
     deaths 0005604: 93366 rows, 1996-2013, 5187 parishes
   deaths: 132004 rows, 1996-2025, 5489 parishes
   Done. Failed: 0
[09:53:42]    done
[09:53:42] == Task: calendar-year deaths by age (weekly + parish) ==
[09:53:42] -> calendar deaths by age
     2025: Portugal 121840, Continente 116499 deaths
   Done.
[09:53:43]    done
[09:53:43] == Task: SNS Transparency portal ==
[09:53:43] -> SNS primary care
   utentes-inscritos-em-cuidados-de-saude-primarios: 6450 rows, 2016-01 to 2026-06, 96 units
   rastreios-oncologicos: 7223 rows, 2015-01 to 2026-09, 77 units
   diabetes: 7231 rows, 2015-01 to 2026-09, 77 units
   hipertensao: 7223 rows, 2015-01 to 2026-09, 77 units
   saude-da-mulher-e-crianca: 7231 rows, 2015-01 to 2026-09, 77 units
   Done. Failed: 0
[09:54:12]    done
[09:54:12] == Task: deaths by cause and age ==
[09:54:15] Nothing to fetch or re-check for 0013166.
[09:54:15] == Task: population ==
[09:54:16] -> population 0012918 2025
   Fetching 0012918 for 2025
     2025: 340 areas, Portugal = 11 424 031
   Done: 1 year files written.
[09:54:18]    done
[09:54:18] == Task: municipal death totals ==
[09:54:18] -> death totals
     0008461 2010: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2011: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2012: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2013: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2014: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
   503 - Service Unavailable
     0008461 2015: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2016: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2017: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2018: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2019: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2020: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2021: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
     0008461 2022: 308 municipalities | all-cause HM: municipal sum 0, Portugal 0
   Done: 21 written, 0 failed.
[10:07:01]    done
[10:07:01] == Task: births and under-1 deaths ==
[10:07:01] -> live births
   
   == 0000003: 19 year(s) ==
   
   == 0008084: 7 year(s) ==
   
   == 0012434: 5 year(s) ==
     2025: ok (310 areas, Portugal 87764)
   
   Done: 1 written, 30 already present, 0 failed.
[10:07:08]    done
[10:07:08] -> under-1 deaths by cause
   Fetching under-1 deaths for 46 year(s), one request each
     2025 come from 0012540: all causes, both sexes only.
     2025 (0012540, total only): ok - 340 rows, Portugal all-cause 246
   Done: 1 written, 45 already present, 0 failed.
[10:07:12]    done
[10:07:12] -> complete under-1 death counts
   Years: 2011-2025 (15)
     0012541 2025: 308 municipalities | municipal sum 246, Portugal 246
   Done: 1 written, 0 failed.
[10:07:32]    done
[10:07:32] Finished in 29.1 min.
```

## Changes in this run

Run id: `2026-10-26T093824`  |  Re-checked years from 2025  |  21 file(s) added, 4 replaced

Revised by INE since the previous import (previous versions in `data/archive/2026-10-26T093824/`):

Totals are Portugal, both sexes, all causes where the files carry them; a
correction between municipalities shows in the changed rows, not the total.

| Dataset | Year | Files | Rows changed | Total before | Total after | Change |
|---|---|---|---|---|---|---|
| sns/diabetes | NA | 1 | 156 | 104 685 954 | 105 864 845 | +1.13% |
| sns/hipertensao | NA | 1 |  78 |  38 344 401 |  38 685 692 | +0.89% |
| sns/rastreios-oncologicos | NA | 1 | 234 | 462 121 166 | 466 970 544 | +1.05% |
| sns/saude-da-mulher-e-crianca | NA | 1 | 156 |   7 071 911 |   7 133 896 | +0.88% |
