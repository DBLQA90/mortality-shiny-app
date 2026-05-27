#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
cd "$repo_root"

app_file="${APP_FILE:-mortality-shiny-app.R}"
out_dir="${OUT_DIR:-data/snapshots}"
years="${YEARS:-1991:2005,2008}"
areas="${AREAS:-ALL}"
causes="${CAUSES:-ALL}"
priority_areas="${PRIORITY_AREAS:-Norte|Portugal}"
max_batches="${MAX_BATCHES:-26}"
area_batch_size="${AREA_BATCH_SIZE:-12}"
year_batch_size="${YEAR_BATCH_SIZE:-1}"
timeout_seconds="${TIMEOUT:-240}"
curl_retries="${CURL_RETRIES:-4}"
curl_retry_sleep="${CURL_RETRY_SLEEP:-15}"
export_retries="${EXPORT_RETRIES:-3}"
export_retry_sleep="${EXPORT_RETRY_SLEEP:-30}"
iterations="${ITERATIONS:-0}"
sleep_seconds="${SLEEP_SECONDS:-5}"
keep_raw="${KEEP_RAW:-0}"
auto_commit="${AUTO_COMMIT:-0}"
auto_push="${AUTO_PUSH:-0}"
branch="${BRANCH:-main}"
log_dir="${LOG_DIR:-.snapshot-download-logs}"

mkdir -p "$log_dir"

if [[ ! -f "$app_file" ]]; then
  echo "Cannot find app file: $app_file" >&2
  exit 1
fi

if [[ ! -f tools/build_0008206_snapshot_from_portal.R ]]; then
  echo "Cannot find tools/build_0008206_snapshot_from_portal.R" >&2
  exit 1
fi

if [[ ! -f tools/update_snapshot_inventory.R ]]; then
  echo "Cannot find tools/update_snapshot_inventory.R" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required but was not found in PATH." >&2
  exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required but was not found in PATH." >&2
  exit 1
fi

update_inventory() {
  if [[ ! -d "$out_dir" ]]; then
    return 0
  fi

  echo "Updating snapshot inventory in $out_dir"
  MORTALITY_SNAPSHOT_DIR="$out_dir" Rscript tools/update_snapshot_inventory.R
}

sync_progress() {
  if [[ "$auto_commit" != "1" ]]; then
    return 0
  fi

  if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "AUTO_COMMIT=1 was requested, but this folder is not inside a git worktree." >&2
    return 1
  fi

  git add "$out_dir"

  if git diff --cached --quiet; then
    echo "No new snapshot files to commit."
    return 1
  fi

  git commit -m "Fill missing 0008206 snapshot chunks"

  if [[ "$auto_push" == "1" ]]; then
    git pull --rebase --autostash origin "$branch"
    git push origin "HEAD:$branch"
  fi

  return 0
}

finish_on_interrupt() {
  echo
  echo "Interrupted. Completed RDS chunks remain in $out_dir."
  update_inventory || true
  sync_progress || true
  exit 130
}

trap finish_on_interrupt INT TERM

cat <<EOF
Missing 0008206 downloader
  Output:          $out_dir
  Years:           $years
  Areas:           $areas
  Causes:          $causes
  Priority areas:  $priority_areas
  Area batch size: $area_batch_size
  Max batches/run: $max_batches
  Curl retries:    $curl_retries (sleep base ${curl_retry_sleep}s)
  Export retries:  $export_retries (sleep base ${export_retry_sleep}s)
  Iterations:      $iterations (0 means until complete)
  Auto commit:     $auto_commit
  Auto push:       $auto_push

The current online gap is 0008206 local-area data for 1991-2005, plus
the remaining missing areas in 2008. Defaults therefore target
YEARS=1991:2005,2008, all areas, and all causes.
EOF

batch=0
while true; do
  batch=$((batch + 1))
  timestamp="$(date +%Y%m%d-%H%M%S)"
  log_file="$log_dir/0008206-missing-$timestamp.log"

  echo
  echo "Starting 0008206 missing-data batch $batch"
  echo "Log: $log_file"

  set +e
  Rscript tools/build_0008206_snapshot_from_portal.R \
    "app=$app_file" \
    "out=$out_dir" \
    "years=$years" \
    "areas=$areas" \
    "causes=$causes" \
    "priority_areas=$priority_areas" \
    "max_batches=$max_batches" \
    "area_batch_size=$area_batch_size" \
    "year_batch_size=$year_batch_size" \
    "timeout=$timeout_seconds" \
    "curl_retries=$curl_retries" \
    "curl_retry_sleep=$curl_retry_sleep" \
    "export_retries=$export_retries" \
    "export_retry_sleep=$export_retry_sleep" \
    "keep_raw=$keep_raw" 2>&1 | tee "$log_file"
  status="${PIPESTATUS[0]}"
  set -e

  update_inventory || true
  sync_progress || true

  if [[ "$status" -ne 0 ]]; then
    echo "Batch $batch exited with status $status."
    echo "Completed chunks, if any, were kept in $out_dir."
    exit "$status"
  fi

  if grep -q "Done. Batches built: 0." "$log_file"; then
    echo "All requested missing 0008206 chunks are already present in $out_dir."
    break
  fi

  if [[ "$iterations" -gt 0 && "$batch" -ge "$iterations" ]]; then
    echo "Reached ITERATIONS=$iterations. Stopping."
    break
  fi

  sleep "$sleep_seconds"
done

echo
echo "Done. RDS chunks are in $out_dir."
echo "Upload/commit $out_dir/deaths/0008206 and $out_dir/snapshot_inventory.rds when ready."
