# Snapshot refresh status

> **Note (2026-08-20):** the archive has been updated since the run logged below,
> by targeted tools rather than by `tools/refresh_snapshots.R`, so the log that
> follows describes only the 2026-08-11 run. Changes since then:
>
> | What | How |
> |---|---|
> | `0013166` 2022-2023 refetched (three municipalities were missing) | `tools/fetch_death_year.R` |
> | Population 2021-2025 replaced with the revised `0012918` series | `tools/fetch_population_year.R` |
> | Under-1 deaths extended to 2025 | `tools/fetch_infant_deaths.R` |
> | Inventory manifest rebuilt | `tools/update_snapshot_inventory.R` |
>
> Coverage now: population 1991-2025, deaths by cause 1991-2024, live births
> 1995-2025, under-1 deaths 1980-2025.

Last run: 2026-08-11 17:07:24 WEST
Task: `all`  |  Duration: 38.7 min

```
[17:07:24] Refresh started; task=all, budget=150 min
[17:07:24] == Task: repair ambiguous geographies ==
[17:07:24] -> fix ambiguous areas (Lisboa, Calheta, Lagoa)
      2010: ok (66 cause files from 1 request)
      2011: ok (66 cause files from 1 request)
      2012: ok (66 cause files from 1 request)
      2013: ok (66 cause files from 1 request)
      2014: ok (66 cause files from 1 request)
      2015: ok (66 cause files from 1 request)
      2016: ok (66 cause files from 1 request)
      2017: ok (66 cause files from 1 request)
      2018: ok (66 cause files from 1 request)
      2019: ok (66 cause files from 1 request)
      2020: ok (66 cause files from 1 request)
      2021: ok (66 cause files from 1 request)
      2022: ok (66 cause files from 1 request)
   
   Repaired 8580 chunk(s); 0 already correct; 0 failed.
[17:40:56]    done
[17:40:56] == Task: latest death years ==
[17:40:57] Missing years: 2024
[17:40:57] -> deaths 0013166 2024
   Fetching 0013166 2024: 66 causes, all areas (1 request each)
     10 causes written ...
     20 causes written ...
     30 causes written ...
     40 causes written ...
     50 causes written ...
     60 causes written ...
   Done: 66 written, 0 already present, 0 failed.
[17:45:24]    done
[17:45:24] == Task: ambiguous municipality labels ==
[17:45:25] 0008206: Lagoa (1500806); Calheta (2004501); Lagoa (2004201); Calheta (3003101)  <-- AMBIGUOUS: Calheta, Lagoa
[17:45:26] 0013166: Lagoa (1500806); Calheta (R.A.A.) (2004501); Lagoa (R.A.A.) (2004201); Calheta (R.A.M.) (3003101)  <-- unambiguous
[17:45:26] 0008273: Lagoa (1500806); Calheta (2004501); Lagoa (2004201); Calheta (3003101)  <-- AMBIGUOUS: Calheta, Lagoa
[17:45:27] 0003182: Lagoa (1500806); Calheta (R.A.A.) (2004501); Lagoa (R.A.A) (2004201); Calheta (R.A.M.) (3003101)  <-- unambiguous
[17:45:27] Wrote data/snapshots/ambiguous_areas.rds
[17:45:27] == Task: regional (NUTS II) rows ==
[17:45:27] Note: rarely needed. Years fetched by fetch_death_year.R already include INE's regional rows, because the all-areas response carries every NUTS level; and the municipal region mode derives regions by summing municipalities. This task only backfills regions into older chunks that were built area-by-area.
[17:45:27] 0013166 regions: Continente | Norte | Centro | Oeste e Vale do Tejo | Grande Lisboa | Península de Setúbal | Alentejo | Algarve | Região Autónoma dos Açores | Região Autónoma da Madeira
[17:45:27] -> regional deaths 0013166
   ERROR: error in running command
[17:45:27]    done
[17:45:28] 0008206 regions: Continente | Norte | Centro | Área Metropolitana de Lisboa | Alentejo | Algarve | Região Autónoma dos Açores | Região Autónoma da Madeira
[17:45:28] -> regional deaths 0008206
   ERROR: error in running command
[17:45:28]    done
[17:45:28] -> regional population
   ERROR: error in running command
[17:45:28]    done
[17:45:28] == Task: rebuild inventory ==
[17:45:28] -> snapshot inventory
   Wrote 2343 snapshot inventory rows to /home/dblqa/Desktop/DGS/Mortality/scripts/mortality-app-smr/data/snapshots/snapshot_inventory.rdsPopulation chunks: 33; death chunks: 2310
[17:46:06]    done
[17:46:06] Finished in 38.7 min.
```
