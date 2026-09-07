#!/bin/bash
# run_onset.sh  --  bracket and locate the SHEDDING onset below AoA 20.
#
#   ./run_onset.sh 18 16
#
# BOTH ANGLES ARE SEEDED FROM AoA 20, NOT FROM EACH OTHER.  This is the one
# thing that must not be "simplified" into the run_chain.sh ladder.  The
# measurement wanted here is how a finite-amplitude shedding oscillation EVOLVES
# at the target angle -- it must therefore be given one.  Chaining 16 off 18
# would seed it from whatever 18 ended at, and if 18 decays to a steady state
# then 16 inherits no oscillation at all, sits at its fixed point, and reports
# "no shedding" while proving nothing about its stability.
#
# WHY A "no shedding" ABORT IS A RESULT, NOT A FAILURE.  run_staged.sh exits 1
# when phase 1 finishes with SHEDDING=no, because for a POD run that is a dead
# end.  Here it is the answer.  The loop below records it and carries on to the
# next angle instead of taking the whole chain down with it.
#
# Sequential at 6 ranks: 6 PHYSICAL cores, and concurrent jobs oversubscribe.
cd "$(dirname "$0")"
MESH=meshes/M2_ns160_nr180_b025.msh
# Seed case, overridable:  SRC=re500_aoa16 ./run_onset.sh 14
# Default is AoA 20 -- the strongest well-characterised period-1 state, right for
# seeding angles a few degrees away.  For an angle expected to be BELOW onset,
# override to the nearest shedding angle instead: the decay only becomes a linear
# eigenvalue once the amplitude is small, so starting from the smallest available
# saturated amplitude (AoA 16, A=0.0694, vs AoA 20's 0.1403) reaches that regime
# sooner, and the nearer angle also means less far-field readjustment to ride out.
SRC=${SRC:-re500_aoa20}
NP=6

# dt IS FROZEN, AND THE SAME FROZEN VALUE IS USED AT EVERY ANGLE.  Two reasons,
# both specific to what this sweep measures.
#
#  (1) analysis/growth_rate.py fits an exponential to the C_L envelope across
#      the whole phase-1 record.  The scheme's numerical damping depends on dt,
#      so letting the Courant controller re-tune between chunks would put a step
#      into the quantity being fitted -- an artefact indistinguishable from the
#      physical decay it is there to measure.
#  (2) Comparing sigma BETWEEN angles only means something if the numerical
#      damping is identical at both.  Sharing one dt makes the difference
#      between sigma(18) and sigma(16) physics rather than discretisation.
#
# The value is AoA 18's own phase-0-derived dt, which ran a full 20-t chunk at
# Co_max=0.417 with 7.8% step reversals before the controller raised it and the
# next chunk diverged.  Lower angles separate less and so run at lower Co still,
# which makes this conservative at 16 and safe to reuse.
export DT_FIXED=0.0005572582892

[ $# -gt 0 ] || { echo "usage: run_onset.sh <AoA> [<AoA> ...]"; exit 1; }
[ -d "$SRC" ] || { echo "ABORT: seed case $SRC does not exist"; exit 1; }

for AOA in "$@"; do
  CASE=re500_aoa$AOA
  echo "=============================================================="
  if [ -d $CASE/processor0 ]; then
    echo "$(date '+%F %T')  AoA $AOA  (already seeded -- resuming)"
  else
    echo "$(date '+%F %T')  AoA $AOA  (warm start from $SRC)"
    ./make_case.sh  $CASE 500 $AOA $MESH   || { echo "make_case FAILED for $AOA"; continue; }
    ./warm_start.sh $SRC latestTime $CASE $NP || { echo "warm_start FAILED for $AOA"; continue; }
  fi

  echo "  dt frozen at $DT_FIXED for the whole record"
  if ./run_staged.sh $CASE $NP; then
    echo "$(date '+%F %T')  AoA $AOA: SHEDDING, phase 2 written"
  else
    echo "$(date '+%F %T')  AoA $AOA: run_staged exited non-zero -- see below"
  fi

  echo "--- AoA $AOA limit-cycle verdict ---"
  python3 analysis/limit_cycle.py $CASE || true
  echo "--- AoA $AOA envelope growth rate (sigma < 0 => below onset) ---"
  python3 analysis/growth_rate.py $CASE || true
done
echo "$(date '+%F %T')  ONSET SWEEP COMPLETE: $*"
