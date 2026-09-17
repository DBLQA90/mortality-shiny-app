# Snapshot refresh status

Last run: 2026-09-17 16:58:30 WEST
Task: `all`  |  Duration: 8 min

```
[16:58:30] Refresh started; task=all, budget=120 min
[16:58:30] == Task: deaths by cause and age ==
[16:58:31] Missing years: none; re-checking: 2024
[16:58:31] -> deaths 0013166 2024 (re-check)
   Fetching 0013166 2024: 66 causes, all areas (1 request each)
     10 causes written ...
     20 causes written ...
     30 causes written ...
     40 causes written ...
     50 causes written ...
     60 causes written ...
   Done: 66 written, 0 already present, 0 failed.
[17:01:17]    done
[17:01:18] == Task: population ==
[17:01:18] -> population 0012918 2024, 2025
   Fetching 0012918 for 2024, 2025
     2024: 340 areas, Portugal = 11 387 222
     2025: 340 areas, Portugal = 11 424 031
   Done: 2 year files written.
[17:01:23]    done
[17:01:23] == Task: municipal death totals ==
[17:01:23] -> death totals
     0013166 2024: 308 municipalities | all-cause HM: municipal sum 118 302, Portugal 118 396
   Done: 1 written, 0 failed.
[17:01:35]    done
[17:01:35] == Task: regional death rows ==
[17:01:35] -> regional deaths
   == 0008206: 32 years, 22 territories ==
   == 0013166: 3 years, 19 territories ==
     2024: 67716 rows | all-cause HM: Continente 113 357, NUTS II sum 118 385 | incomplete cells 0
   Done: 1 written, 0 failed.
[17:01:43]    done
[17:01:43] == Task: births and under-1 deaths ==
[17:01:43] -> live births
   
   == 0000003: 19 year(s) ==
   
   == 0008084: 7 year(s) ==
   
   == 0012434: 5 year(s) ==
     2024: ok (310 areas, Portugal 84642)
     2025: ok (310 areas, Portugal 87764)
   
   Done: 2 written, 29 already present, 0 failed.
[17:01:52]    done
[17:01:52] -> under-1 deaths by cause
   Fetching under-1 deaths for 46 year(s), one request each
     2025 come from 0012540: all causes, both sexes only.
     2024 (0013166): ok - 67320 rows, Portugal all-cause 254
     2025 (0012540, total only): ok - 340 rows, Portugal all-cause 246
   Done: 2 written, 44 already present, 0 failed.
[17:02:00]    done
[17:02:00] -> complete under-1 death counts
   Years: 2011-2025 (15)
     0012541 2024: 308 municipalities | municipal sum 253, Portugal 254
     0012541 2025: 308 municipalities | municipal sum 246, Portugal 246
   Done: 2 written, 0 failed.
[17:02:24]    done
[17:02:24] == Task: planning-tab components ==
[17:02:24] -> planning extra (RSI, pensions, purchasing power, waste, births, neonatal)
   == purchasing_power_share: 1993-2023 (16 years) ==
   
   == births_by_mother_age: 1995-2025 (31 years) ==
     2024 0012441: 308 municipalities, 59 categories, Portugal total 84 642
     2025 0012441: 308 municipalities, 59 categories, Portugal total 87 764
   
   == births_by_gestation: 1995-2025 (31 years) ==
     2024 0012434: 308 municipalities, 8 categories, Portugal total 84 642
     2025 0012434: 308 municipalities, 8 categories, Portugal total 87 764
   
   == infant_deaths_by_age: 2011-2025 (15 years) ==
     2024 0012541: 308 municipalities, 21 categories, Portugal total 254
     2025 0012541: 308 municipalities, 21 categories, Portugal total 246
   
   Done. Failed: 0
[17:05:53]    done
[17:05:53] == Task: ambiguous municipality labels ==
[17:05:53] 0008206: Lagoa (1500806); Calheta (2004501); Lagoa (2004201); Calheta (3003101)  <-- AMBIGUOUS: Calheta, Lagoa
[17:05:54] 0013166: Lagoa (1500806); Calheta (R.A.A.) (2004501); Lagoa (R.A.A.) (2004201); Calheta (R.A.M.) (3003101)  <-- unambiguous
[17:05:54] 0008273: Lagoa (1500806); Calheta (2004501); Lagoa (2004201); Calheta (3003101)  <-- AMBIGUOUS: Calheta, Lagoa
[17:05:55] 0003182: Lagoa (1500806); Calheta (R.A.A.) (2004501); Lagoa (R.A.A) (2004201); Calheta (R.A.M.) (3003101)  <-- unambiguous
[17:05:55] Wrote data/snapshots/ambiguous_areas.rds
[17:05:55] == Task: rebuild inventory ==
[17:05:55] -> snapshot inventory
   Wrote 2345 snapshot inventory rows to /home/dblqa/Desktop/DGS/Mortality/scripts/mortality-app-smr/data/snapshots/snapshot_inventory.rdsPopulation chunks: 35; death chunks: 2310
[17:06:32]    done
[17:06:32] Finished in 8 min.
```

## Changes in this run

Run id: `2026-09-17T165830`  |  Re-checked years from 2024  |  0 file(s) added, 0 replaced

No stored value was revised.
