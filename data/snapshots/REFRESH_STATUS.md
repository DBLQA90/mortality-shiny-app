# Snapshot refresh status

Last run: 2026-08-11 15:33:42 UTC
Task: `all`  |  Duration: 30.6 min

```
[15:33:42] Refresh started; task=all, budget=300 min
[15:33:42] == Task: repair ambiguous geographies ==
[15:33:42] -> fix ambiguous areas (Lisboa, Calheta, Lagoa)
   ! Timeout was reached [www.ine.pt]:
   Failed to connect to www.ine.pt port 443 after 135149 ms: Couldn't connect to server
     0008206 Calheta (R.A.M.): no overlapping years in archive; skipped.
   Failed to perform HTTP request.
   Caused by error in `curl::curl_fetch_memory()`:
   ! Timeout was reached [www.ine.pt]:
   Failed to connect to www.ine.pt port 443 after 135148 ms: Couldn't connect to server
     0008206 Lagoa: no overlapping years in archive; skipped.
   Failed to perform HTTP request.
   Caused by error in `curl::curl_fetch_memory()`:
   ! Timeout was reached [www.ine.pt]:
   Failed to connect to www.ine.pt port 443 after 135145 ms: Couldn't connect to server
     0008206 Lagoa (R.A.A.): no overlapping years in archive; skipped.
   
   Repaired 0 chunk(s); 0 already correct; 0 failed.
[15:46:14]    done
[15:46:14] == Task: latest death years ==
[15:48:29] Could not read INE years for 0013166: Failed to retrieve metadata for indicator '0013166'.
ℹ The API may be unavailable.
[15:48:29] No missing death years for 0013166.
[15:48:29] == Task: ambiguous municipality labels ==
[15:57:30] == Task: regional (NUTS II) rows ==
[15:57:30] Note: rarely needed. Years fetched by fetch_death_year.R already include INE's regional rows, because the all-areas response carries every NUTS level; and the municipal region mode derives regions by summing municipalities. This task only backfills regions into older chunks that were built area-by-area.
[15:59:45] Could not read INE geographies for 0013166: Failed to retrieve metadata for indicator '0013166'.
ℹ The API may be unavailable.
[16:02:00] Could not read INE geographies for 0008206: Failed to retrieve metadata for indicator '0008206'.
ℹ The API may be unavailable.
[16:04:15] Could not read INE geographies for 0008273: Failed to retrieve metadata for indicator '0008273'.
ℹ The API may be unavailable.
[16:04:15] == Task: rebuild inventory ==
[16:04:15] -> snapshot inventory
   Error in library(tidyverse) : there is no package called ‘tidyverse’
   Calls: suppressPackageStartupMessages -> withCallingHandlers -> library
   Execution halted
[16:04:16]    FAILED (exit 1)
[16:04:16] Finished in 30.6 min.
```
