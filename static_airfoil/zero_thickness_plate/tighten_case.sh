#!/bin/bash
# tighten_case.sh <extraIters> <caseDir> [<caseDir> ...]
# Continue finished steady cases with fvSolution.tight (tolerance 1e-12, relTol 0.01 on every solve).
# The production tolerances (p 1e-7, U 1e-8) stall: once the residual is ~8e-8 most solves do 0
# iterations and the state freezes 2.5e-7 to 3e-6 in C_L away from the fixed point (tight_test/,
# 2026-10-19).  About 2400-3000 tight iterations bring C_L to within ~1e-9.
# Restarts from the latest time, writes once at the end; the old log is kept as log.steady_to<old>.
HERE=$(cd "$(dirname "$0")" && pwd)
N=$1; shift
for C in "$@"; do
  [ -f "$C/system/controlDict" ] || { echo "SKIP $C (no case)"; continue; }
  grep -qa "^End" "$C/log.steady" 2>/dev/null || { echo "SKIP $C (previous run not finished)"; continue; }
  old=$(ls -d $C/[0-9]* | xargs -n1 basename | sort -g | tail -1)
  END=$((old + N))
  sed -i -e 's/^startFrom .*/startFrom       latestTime;/' -e "s/^endTime .*/endTime         $END;/" \
         -e "s/^writeInterval .*/writeInterval   $END;/" "$C/system/controlDict"
  cp "$HERE/fvSolution.tight" "$C/system/fvSolution"
  mv "$C/log.steady" "$C/log.steady_to$old"
  echo "$C: $old -> $END (tight)"
done
