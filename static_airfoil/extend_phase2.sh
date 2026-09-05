#!/bin/bash
# extend_phase2.sh <caseDir> <t_conv> <nOrbits> [nCores]
#
# Extend a finished phase 2 to nOrbits FUNDAMENTAL periods measured from t_conv.
#
# Why this is not just "re-run run_staged.sh": the run is already past its gate,
# so re-entering the driver would re-walk phase 1.  This picks up the existing
# case exactly where it stopped and only moves endTime.
#
# Two things are deliberately NOT recomputed:
#
#   dt              -- changing it costs an Euler step (backward needs two old
#                      time levels, a restart has one) and, worse, would break
#                      the snapshot grid below.
#   writeInterval   -- kept byte-identical to what phase 2 already used, so the
#                      new snapshots continue the SAME uniform time grid as the
#                      existing ones.  A POD set with a seam in its sampling
#                      interval is not a uniformly weighted ensemble.  Re-deriving
#                      it from a fresh period estimate would shift it by a step or
#                      two (measured: 224 -> 225 x dt) for no benefit.
#
# endTime is snapped UP to a whole number of write intervals from the restart
# point, so it lands on both a write and a time step.
set -e

CASE=$(cd "$1" && pwd)
T_CONV=$2
NORB=$3
NP=${4:-6}
ROOT=$(cd "$(dirname "$0")" && pwd)

source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
cd "$CASE"
exec > >(tee -a driver.log) 2>&1
echo "=================================================================="
echo "$(date '+%F %T')  EXTEND phase 2  case=$CASE  np=$NP  target=$NORB orbits"

# OpenFOAM names a time directory from the VALUE, so a write that lands on a
# whole number is "213", not "213.0".  The obvious glob [0-9]*.[0-9]* silently
# skips those -- it undercounted this case by one and made a perfectly uniform
# snapshot series look like it had a 2x gap at t=213.  Match on the value.
count_snaps() {   # t_min -> number of time dirs strictly after t_min
  ls -d [0-9]* 2>/dev/null \
    | awk -v t="$1" '/^[0-9]+(\.[0-9]+)?$/ && $1+0 > t+1e-9' | wc -l
}
DT=$(awk '/^deltaT/{print $2}' system/controlDict | tr -d ';')
WI=$(awk '/^writeInterval/{print $2}' system/controlDict | tr -d ';')
T_NOW=$(ls -d processor0/[0-9]* | sed 's|.*/||' | sort -g | tail -1)

eval $(python3 "$ROOT/analysis/limit_cycle.py" .)
[ "$SHEDDING" = yes ] || { echo "ABORT: SHEDDING=no"; exit 1; }
[ "$QUASIPERIODIC" = no ] || { echo "ABORT: QUASIPERIODIC=yes -- not a closed orbit"; exit 1; }

read T2 NSTEP NNEW HRS <<< $(awk -v p=$PERIOD -v tc=$T_CONV -v n=$NORB -v now=$T_NOW \
                                 -v w=$WI -v d=$DT 'BEGIN{
  want = tc + n*p;                      # target end, in whole orbits from t_conv
  k = int((want-now)/w); if (now + k*w < want) k++;   # snap UP to a whole write
  t2 = now + k*w;
  printf "%.10g %d %d %.2f", t2, k*int(w/d+0.5), k, k*int(w/d+0.5)*0.0775/3600}')

echo "  period-$ORDER orbit: T_orbit=$PERIOD  T_carrier=$SAMPLE_PERIOD  subharm=$SUBHARM_RATIO"
echo "  dt=$DT  writeInterval=$WI (unchanged, $(awk -v w=$WI -v s=$SAMPLE_PERIOD 'BEGIN{printf "%.2f", s/w}') samples/carrier period)"
echo "  restart t=$T_NOW -> endTime $T2   (+$NNEW writes, $NSTEP steps, ~$HRS h)"
echo "  orbits from t_conv=$T_CONV: $(awk -v a=$T_NOW -v t=$T_CONV -v p=$PERIOD 'BEGIN{printf "%.2f", (a-t)/p}') -> $(awk -v a=$T2 -v t=$T_CONV -v p=$PERIOD 'BEGIN{printf "%.2f", (a-t)/p}')"

sed -i -e "s/^startFrom .*/startFrom       latestTime;/" \
       -e "s/^endTime .*/endTime         $T2;/" system/controlDict

mpirun -np $NP foamRun -parallel > log.phase2_ext 2>&1 || true
grep -q "^End" log.phase2_ext || { echo "FAILED"; tail -30 log.phase2_ext; exit 1; }

R=$(bash "$ROOT/analysis/sawtooth_check.sh" "$CASE" "$T_NOW" \
    | grep -o 'step-reversals=[0-9.]*' | cut -d= -f2)
echo "  sawtooth check (phase 2 extension): step-reversals=${R:-n/a}%"
awk -v r="${R:-0}" 'BEGIN{ exit (r+0 > 60) }' \
  || { echo "ABORT: 2*dt sawtooth -- result is garbage."; exit 1; }

echo "--- reconstructPar (new times only) ---"
reconstructPar -newTimes -time "${T_CONV}:" > log.reconstructPar_ext 2>&1 || true
tail -3 log.reconstructPar_ext
echo "snapshots: $(count_snaps 0)   $(du -sh . | cut -f1)"
python3 "$ROOT/analysis/limit_cycle.py" .
echo "$(date '+%F %T')  EXTENSION DONE  case=$CASE"
