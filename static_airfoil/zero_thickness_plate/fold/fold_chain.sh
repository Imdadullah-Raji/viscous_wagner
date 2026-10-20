#!/bin/bash
# fold_chain.sh <openfoam bashrc> <ranks> <dir> <angle> [<angle> ...]
# Fold test for the zero-thickness plate (2026-10-20): Padé puts the nearest singularity of the steady
# C_L(alpha) on the real axis near 10.3-10.8 deg at Re 500.  A fold there would end the steady branch,
# and would usually show as hysteresis between a chain climbing from below and one descending from above.
# Runs the angles IN THE ORDER GIVEN (run_batch.sh alone sorts ascending), one at a time, each seeded
# (WARM=1) from the nearest finished angle in <dir>: so an ascending list continues from below, a
# descending one from above.  <dir>/QUEUE lists the order for progress.py.
cd "$(dirname "$0")" || exit 1
BASHRC=$1; NP=$2; D=$3; shift 3
echo "$*" > $D/QUEUE
for a in "$@"; do
  ONLY="$a" WARM=1 ./run_batch.sh $BASHRC $NP $D | grep -v "DONE\|BATCH COMPLETE"
done
echo "$(date '+%F %T')  CHAIN COMPLETE: $D ($*)"
