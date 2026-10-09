#!/bin/bash
# run_batch.sh <openfoam bashrc> <ranks> <reDir> [<reDir> ...]
# Plate Re sweep (instructions/Quadratic_hypothesis.md): runs every steady_aoa* in each
# re<R>/ folder, one <ranks>-rank job at a time, lowest angle first; skips cases that already
# ended ("End" in log.steady).  Cases were built on the main machine by plate_case/run_sweep.py
# --build-only (6000 SIMPLE iterations, writePrecision 12); force histories are copied back and
# judged there by run_sweep.py --report-only, which extends any case that has not converged.
cd "$(dirname "$0")" || exit 1
source "$1" >/dev/null 2>&1; NP=$2; shift 2
for D in "$@"; do
  for C in $(ls -d $D/steady_aoa* | sort -t a -k4 -g); do
    [ -d "$C" ] || continue
    grep -qa "^End" "$C/log.steady" 2>/dev/null && continue
    ( cd "$C"
      sed -i "s/^numberOfSubdomains .*/numberOfSubdomains  $NP;/" system/decomposeParDict
      [ -d processor0 ] || decomposePar > log.decomposePar 2>&1
      mpirun -np $NP foamRun -parallel > log.steady 2>&1
      reconstructPar -latestTime > log.reconstructPar 2>&1 )
    ok="ran"; grep -qa "^End" "$C/log.steady" || ok="NO-End"
    printf '%s  %s %s  %s\n' "$(date '+%H:%M:%S')" "$C" "$ok" \
           "$(tail -1 "$C"/postProcessing/forceCoeffs/0/*.dat 2>/dev/null | awk '{printf "iter=%s Cd=%.6f Cl=%.6f", $1, $3, $4}')"
  done
  echo "$(date '+%F %T')  DONE $D"
done
echo "$(date '+%F %T')  BATCH COMPLETE: $*"
