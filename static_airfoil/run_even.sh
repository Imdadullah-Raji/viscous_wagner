#!/bin/bash
# run_even.sh -- AoA 22 then 24 at Re=500, in sequence, 6 ranks each.
#
# Warm-started up the angle ladder (20 -> 22 -> 24) so only the first run at this
# Re ever paid a cold transient.  Sequential, not concurrent: this box has 6
# PHYSICAL cores and two 4-rank jobs oversubscribe them at ~1.7x cost per step.
#
# RESTART-SAFE.  Re-running this script after an interruption picks up where it
# stopped: a case is built and seeded only if it has no processor0 yet, and
# run_staged.sh resumes from its last phase-1 chunk boundary.  Do NOT "simplify"
# by always calling make_case.sh -- that begins with rm -rf on the case and would
# throw away both the progress and the warm start.
set -e
cd "$(dirname "$0")"
MESH=meshes/M2_ns160_nr180_b025.msh
NP=6

for pair in "22 20" "24 22"; do
  set -- $pair; AOA=$1; SRC=$2
  echo "=============================================================="
  if [ -d re500_aoa$AOA/processor0 ]; then
    echo "$(date '+%F %T')  AoA $AOA  (already seeded -- resuming)"
  else
    echo "$(date '+%F %T')  AoA $AOA  (warm start from AoA $SRC)"
    ./make_case.sh   re500_aoa$AOA 500 $AOA $MESH
    ./warm_start.sh  re500_aoa$SRC latestTime re500_aoa$AOA $NP
  fi
  ./run_staged.sh  re500_aoa$AOA $NP
  echo "$(date '+%F %T')  AoA $AOA done"
  python3 analysis/limit_cycle.py re500_aoa$AOA
done
echo "$(date '+%F %T')  BOTH ANGLES COMPLETE"
