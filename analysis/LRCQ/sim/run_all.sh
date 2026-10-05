#!/usr/bin/env bash
# Run the LRCQ simulation grid in two parallel lanes. Restartable.
cd "$(dirname "$0")/.."
out=${1:-/home/user/data/sim/out}; reps=${2:-20}
lane() { for s in "$@"; do Rscript sim/run_sim.R "$s" 1 "$reps" "$out"; done; }
lane base q30 q100 ukb sparse dense inf strat &
lane strat_miss rb_local rb_cross het ref500 tag08 tag02 q1000 &
wait
