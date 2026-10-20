#!/bin/bash
# run_quarters.sh <bashrc> <ranks> <folder> <tighten_case.sh> <angle...>
# New cases for the a5/a7 data (2026-10-20), already seeded by mapFields from a production neighbour.
# Two phases, because tight tolerances from the start cost ~3.5 s/iteration (every corrector must
# reduce its residual 100x while the fields are still far off):
#   1. production fvSolution, 6000 iterations (cheap; ends at the ~1e-7 tolerance stall);
#   2. tighten_case.sh 3000: fvSolution.tight from the latest time.
cd "$(dirname "$0")" || exit 1
BASHRC=$1; NP=$2; D=$3; TIGHTEN=$4; shift 4
ONLY="$*" ./run_batch.sh $BASHRC $NP $D
for a in "$@"; do $TIGHTEN 3000 $D/steady_aoa$a; done
ONLY="$*" ./run_batch.sh $BASHRC $NP $D
echo "$(date '+%F %T')  QUARTERS COMPLETE: $D"
