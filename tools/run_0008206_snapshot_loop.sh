#!/usr/bin/env bash
set -Eeuo pipefail

repo_root="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$repo_root"

app_file="${APP_FILE:-mortality-shiny-app.R}"
out_dir="${OUT_DIR:-data/snapshots}"
branch="${BRANCH:-main}"
max_chunks="${MAX_CHUNKS:-50}"
years="${YEARS:-ALL}"
causes="${CAUSES:-ALL}"
areas="${AREAS:-ALL}"
priority_areas="${PRIORITY_AREAS:-Portugal}"
population_area_batch_size="${POPULATION_AREA_BATCH_SIZE:-50}"
death_area_batch_size="${DEATH_AREA_BATCH_SIZE:-50}"
iterations="${ITERATIONS:-0}"
sleep_seconds="${SLEEP_SECONDS:-5}"
auto_commit="${AUTO_COMMIT:-1}"
auto_push="${AUTO_PUSH:-1}"
log_dir="${LOG_DIR:-.snapshot-download-logs}"

mkdir -p "$log_dir"

if [[ ! -f "$app_file" ]]; then
  echo "Cannot find app file: $app_file" >&2
  exit 1
fi

if [[ ! -f tools/build_0008206_snapshot_chunks.R ]]; then
  echo "Cannot find tools/build_0008206_snapshot_chunks.R" >&2
  exit 1
fi

sync_progress() {
  if [[ "$auto_commit" != "1" ]]; then
    return 0
  fi

  git add "$out_dir"

  if git diff --cached --quiet; then
    echo "No new snapshot files to commit."
    return 1
  fi

  git commit -m "Update 0008206 snapshot chunks"

  if [[ "$auto_push" == "1" ]]; then
    git pull --rebase --autostash origin "$branch"
    git push origin "HEAD:$branch"
  fi

  return 0
}

trap 'echo; echo "Interrupted. Syncing completed chunks before exit..."; sync_progress || true; exit 130' INT TERM

batch=0
while true; do
  batch=$((batch + 1))
  timestamp="$(date +%Y%m%d-%H%M%S)"
  log_file="$log_dir/0008206-batch-$timestamp.log"

  echo "Starting 0008206 snapshot batch $batch"
  echo "  years=$years causes=$causes areas=$areas priority_areas=$priority_areas max_chunks=$max_chunks"
  echo "  log=$log_file"

  set +e
  Rscript tools/build_0008206_snapshot_chunks.R \
    "app=$app_file" \
    "out=$out_dir" \
    "max_chunks=$max_chunks" \
    "years=$years" \
    "causes=$causes" \
    "areas=$areas" \
    "priority_areas=$priority_areas" \
    "population_area_batch_size=$population_area_batch_size" \
    "death_area_batch_size=$death_area_batch_size" 2>&1 | tee "$log_file"
  status="${PIPESTATUS[0]}"
  set -e

  sync_progress_status=0
  sync_progress || sync_progress_status=$?

  if [[ "$status" -ne 0 ]]; then
    echo "Batch $batch exited with status $status. Completed chunks, if any, were synced."
    exit "$status"
  fi

  if grep -q "All requested 0008206 chunks already exist" "$log_file"; then
    echo "All requested chunks already exist. Done."
    break
  fi

  if [[ "$sync_progress_status" -ne 0 ]]; then
    echo "No new files were produced by this successful batch. Done."
    break
  fi

  if [[ "$iterations" -gt 0 && "$batch" -ge "$iterations" ]]; then
    echo "Reached ITERATIONS=$iterations. Done."
    break
  fi

  sleep "$sleep_seconds"
done
