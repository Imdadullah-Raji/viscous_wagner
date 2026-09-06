#!/bin/bash
# run_chain.sh <AoA>:<srcAoA> [<AoA>:<srcAoA> ...]
#
#   ./run_chain.sh 28:30 26:28      # 28 warm-started from 30, then 26 from 28
#
# Warm-started marches along the angle ladder so only the first run at a Re ever
# pays a cold transient.  Sequential, not concurrent: 6 PHYSICAL cores, and two
# 4-rank jobs oversubscribe them at ~1.7x cost per step.
#
# CHOOSING THE SOURCE.  Seed from a case on the SAME side of the period-doubling
# bifurcation as the target, and from as small an angle step as possible.  Near
# onset the subharmonic's growth rate passes through zero, so seeding a
# period-2 angle from a period-1 state means waiting for a mode that grows with
# an e-folding time of 50 t or worse: measured sigma=+0.021/t at AoA 30, which is
# ~376 t from the noise floor to saturation, against T1_MAX=240.  Marching down
# from an established period-2 case instead leaves the subharmonic only needing
# to adjust its level.
#
# RESTART-SAFE.  Re-running picks up where it stopped: a case is built and seeded
# only if it has no processor0 yet, and run_staged.sh resumes from its last
# restart point (phase-1 chunk boundary, or phase 2 via the .phase2 marker).  Do
# NOT "simplify" by always calling make_case.sh -- that begins with rm -rf on the
# case and would throw away both the progress and the warm start.
set -e
cd "$(dirname "$0")"
MESH=meshes/M2_ns160_nr180_b025.msh
NP=6
[ $# -gt 0 ] || { echo "usage: run_chain.sh <AoA>:<srcAoA> [...]"; exit 1; }

for pair in "$@"; do
  AOA=${pair%%:*}; SRC=${pair##*:}
  echo "=============================================================="
  if [ -d re500_aoa$AOA/processor0 ]; then
    echo "$(date '+%F %T')  AoA $AOA  (already seeded -- resuming)"
  else
    [ -d re500_aoa$SRC ] || { echo "ABORT: source case re500_aoa$SRC does not exist"; exit 1; }
    echo "$(date '+%F %T')  AoA $AOA  (warm start from AoA $SRC)"
    ./make_case.sh   re500_aoa$AOA 500 $AOA $MESH
    ./warm_start.sh  re500_aoa$SRC latestTime re500_aoa$AOA $NP
  fi
  ./run_staged.sh  re500_aoa$AOA $NP
  echo "$(date '+%F %T')  AoA $AOA done"
  python3 analysis/limit_cycle.py re500_aoa$AOA
done
echo "$(date '+%F %T')  CHAIN COMPLETE: $*"
