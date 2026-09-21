#!/usr/bin/env bash
# Scheduled data refresh, run by the systemd user timer in tools/systemd/.
#
# Fetches what changes often - INE weekly deaths, the SNS Transparency portal,
# and INE's latest deaths, population and under-1 deaths (re-checking the last
# year for revisions) - through tools/refresh_snapshots.R task=current. Every
# file goes through R/data_versions.R: a revised file is archived, and the run
# is recorded in data/import_log.csv with its date.
#
# It must run on a network INE serves: INE refuses GitHub's hosted runners.
#
# REFRESH_TASK overrides the task (e.g. sns for a quick test).
#
# The data files are left changed in the working tree. Set AUTO_COMMIT=1 (in
# ~/.config/mortality-refresh.env) to commit them, and AUTO_PUSH=1 to push that
# commit - only data/ is ever committed, and only on the branch checked out.
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
state="${XDG_STATE_HOME:-$HOME/.local/state}/mortality-refresh"
mkdir -p "$state"
log="$state/refresh-$(date +%Y%m%d-%H%M%S).log"
[ -f "$HOME/.config/mortality-refresh.env" ] && . "$HOME/.config/mortality-refresh.env"

cd "$repo" || exit 1
{
  echo "== $(date -Is) refresh started in $repo"
  Rscript tools/refresh_snapshots.R task="${REFRESH_TASK:-current}" minutes="${REFRESH_MINUTES:-120}" recent=1 note="scheduled refresh $(date +%F)"
  status=$?
  echo "== $(date -Is) refresh finished with status $status"

  if [ "${AUTO_COMMIT:-0}" = "1" ] && [ -n "$(git status --porcelain -- data)" ]; then
    git add -- data
    git commit -q -m "Data: scheduled refresh $(date +%F)" -- data && echo "committed data changes"
    if [ "${AUTO_PUSH:-0}" = "1" ]; then git push -q && echo "pushed"; fi
  fi
} >>"$log" 2>&1

# Keep the last 30 logs.
ls -1t "$state"/refresh-*.log 2>/dev/null | tail -n +31 | xargs -r rm -f
ln -sf "$log" "$state/latest.log"
exit "${status:-1}"
