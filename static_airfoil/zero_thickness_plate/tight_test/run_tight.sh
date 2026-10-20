#!/bin/bash
# run_tight.sh <bashrc> <ranks> <case...>
# Step 0 of the a5/a7 plan: continue converged cases with tight tolerances (fvSolution.tight).
# Run on the OpenFOAM version that produced the case (cross-version C_L differs by ~4e-7).
BASHRC=$1; NP=$2; shift 2
source $BASHRC
cd "$(dirname "$0")"
for C in "$@"; do
  ( cd $C && grep -qa "^End" log.steady 2>/dev/null && exit 0
    sed -i "s/^numberOfSubdomains .*/numberOfSubdomains  $NP;/" system/decomposeParDict
    [ -d processor0 ] || decomposePar > log.decomposePar 2>&1
    mpirun -np $NP foamRun -parallel > log.steady 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1 )
  echo "$(date '+%F %T') DONE $C"
done
echo "TIGHT TEST COMPLETE"
