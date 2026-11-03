#!/usr/bin/env bash
# Rebuild the INLA small-area forecasts (tools/precompute_inla.R) when their
# data change, run by the systemd user timer in tools/systemd/ (daily and after
# boot) and started by tools/scheduled_refresh.sh when a refresh ends.
#
# The forecasts start from the last year of causes of death; the app serves
# them only while that is still the last year (R/inla_forecasts.R), and falls
# back to the mean of six models until they are rebuilt. The inputs - the
# municipal deaths by cause and age, the population, the death totals that
# fill the ages, and the source corrections - are fingerprinted; when the
# fingerprint differs from the one the forecasts were built from, they are
# rebuilt in data/inla_forecasts.build/, and swapped in only when all 96
# (8 causes x 3 sexes x 4 switch combinations) are there. The old ones are
# kept in data/inla_forecasts.previous/.
#
# Resumable: the build keeps its finished fits, so a run cut short by a
# shutdown continues at the next one (the timer fires after boot). A build
# whose inputs changed meanwhile starts again. Takes about 1.5 hours and 6 GB;
# it waits while the data refresh runs.
#
# After the first rebuild from 2025 or later, the 2021 step is backtested
# once (step2021_test below).
#
# INLA_REBUILD=1 rebuilds even if the inputs are unchanged (after a change of
# the model or of the engine in R/).
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
state="${XDG_STATE_HOME:-$HOME/.local/state}/mortality-refresh"
mkdir -p "$state"
log="$state/inla-$(date +%Y%m%d-%H%M%S).log"
status_file="$state/inla-status.txt"
[ -f "$HOME/.config/mortality-refresh.env" ] && . "$HOME/.config/mortality-refresh.env"
cd "$repo" || exit 1

live=data/inla_forecasts
build=data/inla_forecasts.build
expected=96

say() { echo "$(date '+%F %T') $*" | tee -a "$status_file"; }

# The 2021 step (switch B) is left out with the COVID years separated: at
# origin 2024 it cannot be told from them. Once the forecasts start from 2025
# or later, it is backtested once from two years before - C alone, B and C (the
# 2014 step), and B and C with the 2021 step - against the two years observed
# since (tools/backtest_inla.R); about an hour. The report says whether the
# 2021 step can join B (R/inla_forecasts.R, tools/precompute_inla.R).
step2021_test() {
  local bt=$(( $1 - 2 )) report="$state/inla-step2021-$(( $1 - 2 )).txt" dir v
  [ "$1" -ge 2025 ] && [ ! -f "$report" ] || return 0
  say "2021 step backtest from $bt started"
  for v in "C B0C1 0" "BC B1C1 0" "BCs2021 B1C1 1"; do
    set -- $v
    dir="$repo/data/inla_backtests/${bt}_$1"
    MORTALITY_INLA_DIR="$dir" Rscript tools/precompute_inla.R combos="$2" origin="$bt" step2021="$3" horizon=2 || return 1
  done
  Rscript tools/backtest_inla.R dirs="data/inla_backtests/${bt}_C,data/inla_backtests/${bt}_BC,data/inla_backtests/${bt}_BCs2021" \
    report="$report" && say "2021 step backtest from $bt done: $report"
}

(
  echo "== $(date -Is) INLA check in $repo"
  waited=0
  while [ "$(systemctl --user is-active mortality-refresh.service 2>/dev/null)" = "active" ] && [ "$waited" -lt 240 ]; do
    sleep 60; waited=$((waited + 1))
  done

  fp="$(Rscript tools/inla_inputs_fingerprint.R)" || { say "fingerprint FAILED"; exit 2; }
  if [ "${INLA_REBUILD:-0}" != "1" ] && [ -f "$live/inputs.md5" ] && [ "$(cat "$live/inputs.md5")" = "$fp" ] && [ ! -d "$build" ]; then
    echo "up to date ($fp)"
    step2021_test "$(Rscript -e "cat(readRDS('$live/manifest.rds')\$origin)")"
    exit $?
  fi
  if [ -d "$build" ] && [ "$(cat "$build/inputs.md5" 2>/dev/null)" != "$fp" ]; then
    echo "inputs changed since the build started: starting again"
    rm -rf "$build"
  fi
  mkdir -p "$build"
  echo "$fp" > "$build/inputs.md5"
  say "INLA rebuild started (inputs $fp)"

  MORTALITY_INLA_DIR="$repo/$build" Rscript tools/precompute_inla.R
  n=$(ls "$build" | grep -c '__B[01]C[01]\.rds$')
  if [ "$n" -ne "$expected" ] || [ ! -f "$build/manifest.rds" ]; then
    say "INLA rebuild incomplete: $n of $expected forecasts; the next run continues"
    exit 1
  fi
  origin=$(Rscript -e "cat(readRDS('$build/manifest.rds')\$origin)")
  rm -rf data/inla_forecasts.previous
  [ -d "$live" ] && mv "$live" data/inla_forecasts.previous
  mv "$build" "$live"
  say "INLA rebuild done: $n forecasts from $origin swapped in"
  step2021_test "$origin"
) >>"$log" 2>&1
status=$?

ls -1t "$state"/inla-2*.log 2>/dev/null | tail -n +31 | xargs -r rm -f
ln -sf "$log" "$state/latest-inla.log"
exit "$status"
