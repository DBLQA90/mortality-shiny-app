# Snapshot refresh status

Last run: 2026-11-02 09:35:02 WET
Task: `current`  |  Duration: 27.7 min

```
[09:35:02] Refresh started; task=current, budget=120 min
[09:35:02] == Task: weekly deaths ==
[09:35:02] -> weekly deaths 0012100 / 0010112
   0012100: 296100 rows, 2021-W01 to 2026-W40, 35 regions
   0010112: 337260 rows, 2018-W01 to 2024-W52, 33 regions
   Done. Failed: 0
[09:48:07]    done
[09:48:07] == Task: births and deaths by parish ==
[09:48:07] -> parish births and deaths
     births 0012450: 16980 rows, 2021-2025, 3396 parishes
     births 0008234: 30940 rows, 2014-2023, 3094 parishes
     births 0006239: 93366 rows, 1996-2013, 5187 parishes
   births: 132004 rows, 1996-2025, 5489 parishes
     deaths 0012542: 16980 rows, 2021-2025, 3396 parishes
     deaths 0008235: 30940 rows, 2014-2023, 3094 parishes
     deaths 0005604: 93366 rows, 1996-2013, 5187 parishes
   deaths: 132004 rows, 1996-2025, 5489 parishes
   Done. Failed: 0
[09:54:38]    done
[09:54:38] == Task: calendar-year deaths by age (weekly + parish) ==
[09:54:38] -> calendar deaths by age
     2025: Portugal 121840, Continente 116499 deaths
   Done.
[09:54:40]    done
[09:54:40] == Task: SNS Transparency portal ==
[09:54:40] -> SNS primary care
   utentes-inscritos-em-cuidados-de-saude-primarios: 6450 rows, 2016-01 to 2026-06, 96 units
   rastreios-oncologicos: 7223 rows, 2015-01 to 2026-09, 77 units
   diabetes: 7231 rows, 2015-01 to 2026-09, 77 units
   hipertensao: 7223 rows, 2015-01 to 2026-09, 77 units
   saude-da-mulher-e-crianca: 7231 rows, 2015-01 to 2026-09, 77 units
   Done. Failed: 0
[09:55:47]    done
[09:55:47] == Task: deaths by cause and age ==
[09:56:05] Nothing to fetch or re-check for 0013166.
[09:56:05] == Task: population ==
[09:56:15] -> population 0012918 2025
   Fetching 0012918 for 2025
     2025: 340 areas, Portugal = 11 424 031
   Done: 1 year files written.
[09:57:14]    done
[09:57:14] == Task: municipal death totals ==
[09:57:14] -> death totals
   Done: 0 written, 0 failed.
[09:57:25]    done
[09:57:25] == Task: corrections applied on reading ==
[09:57:25] -> source corrections
     0008206 2011 Nisa: published 24 (H 83, M 82) -> 165 (homens e mulheres publicados; ambos os sexos = a sua soma)
     0008206 2011 Gavião: published 34 (H 91, M 42) -> 91 (sexo e idade repartidos como em 2010 e 2012)
     0008206 2011 Fronteira: published 1 (H 60, M 26) -> 60 (sexo e idade repartidos como em 2010 e 2012)
     0008206 2011 Monforte: published 31 (H 49, M 25) -> 49 (sexo e idade repartidos como em 2010 e 2012)
     0008206 2011 Marvão: published 49 (H 59, M 28) -> 59 (sexo e idade repartidos como em 2010 e 2012)
     0008206: 261203 blank totals recovered
     0013166: 16624 blank totals recovered
     0008461: 128934 blank totals recovered
     by method: diferença 123270, repartição do total por sexo 231876, resíduo da NUTS II 51615
   All-cause municipal totals differing from the parish deaths by more than 2, after correction: 1
     year          area deaths parish
   1 2019 Ponte de Lima    469    472
[10:02:10]    done
[10:02:10] == Task: births and under-1 deaths ==
[10:02:10] -> live births
   
   == 0000003: 19 year(s) ==
   
   == 0008084: 7 year(s) ==
   
   == 0012434: 5 year(s) ==
     2025: ok (310 areas, Portugal 87764)
   
   Done: 1 written, 30 already present, 0 failed.
[10:02:19]    done
[10:02:19] -> under-1 deaths by cause
   Fetching under-1 deaths for 46 year(s), one request each
     2025 come from 0012540: all causes, both sexes only.
     2025 (0012540, total only): ok - 340 rows, Portugal all-cause 246
   Done: 1 written, 45 already present, 0 failed.
[10:02:22]    done
[10:02:22] -> complete under-1 death counts
   Years: 2011-2025 (15)
     0012541 2025: 308 municipalities | municipal sum 246, Portugal 246
   Done: 1 written, 0 failed.
[10:02:44]    done
[10:02:44] Finished in 27.7 min.
```

## Changes in this run

Run id: `2026-11-02T093502`  |  Re-checked years from 2025  |  0 file(s) added, 3 replaced

Revised by INE since the previous import (previous versions in `data/archive/2026-11-02T093502/`):

Totals are Portugal, both sexes, all causes where the files carry them; a
correction between municipalities shows in the changed rows, not the total.

| Dataset | Year | Files | Rows changed | Total before | Total after | Change |
|---|---|---|---|---|---|---|
| sns/diabetes | NA | 1 | 0 | 105 864 845 | 105 744 282 | -0.11% |
| sns/saude-da-mulher-e-crianca | NA | 1 | 0 |   7 133 896 |   7 127 853 | -0.08% |
