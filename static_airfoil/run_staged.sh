#!/bin/bash
# run_staged.sh <caseDir> [nCores]
#
# Two-stage static-airfoil run.
#
#   Phase 0  impulsive start at a deliberately small fixed dt, to measure the
#            real Courant number.  A single large dt from uniform flow is what
#            diverged the earlier attempt; cell-size estimates were off by 2.3x,
#            so dt must be sized from a MEASUREMENT, not a prediction.
#   Phase 1  transient, FIELD-FREE.  writeInterval == endTime means exactly one
#            write per chunk (the restart point) and nothing else.  C_L/C_D
#            still accumulate every step because forceCoeffs appends to a single
#            .dat and never mints a time directory.  Repeat until
#            analysis/limit_cycle.py says the limit cycle is reached.
#   Phase 2  limit cycle, WITH fields, at >=32 snapshots per MEASURED period,
#            for NCYCLES cycles.
#
# adjustTimeStep is 'no' throughout, and dt is never changed after phase 1
# begins: the write interval is snapped to an integer multiple of dt instead.
# That keeps writes landing exactly on time steps without an adjustable-dt
# controller anywhere in the loop.
set -e

CASE=$(cd "$1" && pwd)
NP=${2:-4}
ROOT=$(cd "$(dirname "$0")" && pwd)

DT0=4e-5            # impulsive-start dt
T0=0.5              # end of impulsive start
TARGET_CO=0.35      # deliberately below the 0.6 we could afford at t=0.5:
                    # Co is measured while the flow is still largely attached and
                    # climbs as the separated shear layers develop, so this leaves
                    # ~2x headroom.  The chunk loop below re-tunes it from there.
CO_HI=0.80          # shrink dt if a chunk exceeds this
CO_LO=0.15          # grow dt if a chunk stays under this (only once the wake exists)
DT_MAX=1.0e-3
DT_MIN=2.0e-5
CHUNK=20            # phase-1 chunk length (time units)
T1_MIN=60           # no convergence verdict before this
T1_MAX=240          # give up extending the transient here
NCYCLES=32
PER_CYCLE=32
TOL=0.01

source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
cd "$CASE"
exec > >(tee -a driver.log) 2>&1
echo "=================================================================="
echo "$(date '+%F %T')  case=$CASE  np=$NP"

patch_cd() {   # startFrom deltaT endTime writeControl writeInterval
  sed -i -e "s/^startFrom .*/startFrom       $1;/" \
         -e "s/^deltaT .*/deltaT          $2;/" \
         -e "s/^endTime .*/endTime         $3;/" \
         -e "s/^writeControl .*/writeControl    $4;/" \
         -e "s/^writeInterval .*/writeInterval   $5;/" \
         -e "s/^adjustTimeStep .*/adjustTimeStep  no;/" system/controlDict
}

run() {        # logfile
  mpirun -np $NP foamRun -parallel > "$1" 2>&1 || true
  grep -q "^End" "$1" || { echo "FAILED: $1"; tail -30 "$1"; exit 1; }
}

# The 2*dt sawtooth passes every other health check -- "End" in the log, continuity
# errors ~1e-15, no warning -- so it must be tested for explicitly after every phase.
# A clean run sits near 27% step reversals; a corrupted one pegs at 100%.
sawtooth_gate() {   # t_min  label
  R=$(bash "$ROOT/analysis/sawtooth_check.sh" "$CASE" "$1" \
      | grep -o 'step-reversals=[0-9.]*' | cut -d= -f2)
  echo "  sawtooth check ($2): step-reversals=${R:-n/a}%"
  awk -v r="${R:-0}" 'BEGIN{ exit (r+0 > 60) }' \
    || { echo "ABORT: 2*dt sawtooth detected after $2 -- result is garbage."; exit 1; }
}

# ------------------------------------------------------------- decompose ---
if [ ! -d processor0 ]; then
  decomposePar > log.decomposePar 2>&1
  grep -q "^End" log.decomposePar || { echo "decomposePar FAILED"; tail -20 log.decomposePar; exit 1; }
  echo "decomposed onto $NP subdomains"
fi

# --------------------------------------------------------------- phase 0 ---
echo "--- PHASE 0: impulsive start, dt=$DT0 -> t=$T0 (no fields) ---"
patch_cd startTime $DT0 $T0 runTime $T0
S=$(date +%s); run log.phase0; E=$(date +%s)
sawtooth_gate 0.2 "phase 0"
CO=$(grep -E "^Courant" log.phase0 | tail -1 | awk '{print $6}')
NSTEP=$(echo "$T0/$DT0" | bc -l | cut -d. -f1)
SPS=$(awk -v e=$E -v s=$S -v n=$NSTEP 'BEGIN{printf "%.4f",(e-s)/n}')
echo "phase 0 done in $((E-S))s   Co_max=$CO   ${SPS}s/step"

# dt for the rest of the run, from the MEASURED Courant number at t=0.5 (the
# settled value -- sampling during the impulsive spike, which peaks ~14x higher
# around t=1e-3, gives a dt roughly 7x too small).  Snapped so CHUNK is an exact
# integer number of steps, which keeps every chunk boundary, and therefore
# t_conv, exactly on a time step.
snap_dt() {  # target_dt -> dt that divides CHUNK exactly
  awk -v x=$1 -v hi=$DT_MAX -v lo=$DT_MIN -v ch=$CHUNK 'BEGIN{
    if(x>hi)x=hi; if(x<lo)x=lo;
    n=int(ch/x); if(n*x<ch)n++; printf "%.10g", ch/n}'
}
DT=$(snap_dt $(awk -v d=$DT0 -v c=$CO -v t=$TARGET_CO 'BEGIN{print d*t/c}'))
echo "phase-1/2 dt = $DT  (target Co=$TARGET_CO from measured Co=$CO)"
awk -v s=$SPS -v d=$DT -v tm=$T1_MAX -v nc=$NCYCLES 'BEGIN{
  printf "rough ETA: %.1f h to t=%g, plus phase 2\n", s*(tm/d)/3600, tm}'

# --------------------------------------------------------------- phase 1 ---
echo "--- PHASE 1: transient, field-free, chunks of $CHUNK ---"
T=$T0
CONVERGED=no
while awk -v t=$T -v m=$T1_MAX 'BEGIN{exit !(t<m)}'; do
  T=$(awk -v t=$T -v c=$CHUNK -v t0=$T0 'BEGIN{printf "%.10g", (int((t-t0)/c)+1)*c+t0}')
  patch_cd latestTime $DT $T runTime $T          # interval==endTime -> one write
  echo "  running to t=$T  (dt=$DT) ..."
  run log.phase1_t$T

  # Re-tune dt from what the chunk actually did.  Only ever BETWEEN chunks, so
  # dt stays rigorously fixed within any single foamRun -- the discipline the
  # variable-dt controller broke.  Each change costs one Euler step out of
  # ~100k, which is why this is cheap enough to do at all.
  sawtooth_gate $(awk -v t=$T -v c=$CHUNK 'BEGIN{print t-c}') "chunk to t=$T"

  COC=$(grep -E "^Courant" log.phase1_t$T | awk '{print $6}' | sort -g | tail -1)
  NEWDT=$DT
  if awk -v c=$COC -v h=$CO_HI 'BEGIN{exit !(c>h)}'; then
    NEWDT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$TARGET_CO 'BEGIN{print d*t/c}'))
    echo "    Co_max=$COC > $CO_HI  ->  dt $DT -> $NEWDT"
  elif awk -v c=$COC -v l=$CO_LO -v t=$T -v m=$T1_MIN 'BEGIN{exit !(c<l && t>=m/3)}'; then
    NEWDT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$TARGET_CO 'BEGIN{print d*t/c}'))
    echo "    Co_max=$COC < $CO_LO  ->  dt $DT -> $NEWDT"
  fi
  DT=$NEWDT

  if awk -v t=$T -v m=$T1_MIN 'BEGIN{exit !(t>=m)}'; then
    eval $(python3 "$ROOT/analysis/limit_cycle.py" . --tol $TOL)
    echo "    t=$T  shedding=$SHEDDING  f=$FREQ  T=$PERIOD  Cl=$CL_MEAN+-$CL_P2P" \
         "d(mean,p2p,f)=($D_MEAN,$D_P2P,$D_FREQ)  converged=$CONVERGED"
    [ "$CONVERGED" = yes ] && break
  fi
done

T_CONV=$T
if [ "$CONVERGED" != yes ]; then
  echo "WARNING: no limit cycle by t=$T1_MAX; proceeding with the last estimate."
  eval $(python3 "$ROOT/analysis/limit_cycle.py" . --tol $TOL)
fi
[ "$SHEDDING" = yes ] || { echo "ABORT: no shedding detected -- steady solution?"; exit 1; }

# --------------------------------------------------------------- phase 2 ---
# Snap the write interval DOWN to a whole number of time steps.  Rounding down
# means at least PER_CYCLE samples per cycle, never fewer.
read WI NPC T2 <<< $(awk -v p=$PERIOD -v n=$PER_CYCLE -v d=$DT -v tc=$T_CONV -v nc=$NCYCLES 'BEGIN{
  m=int((p/n)/d); if(m<1)m=1; w=m*d;
  printf "%.10g %.2f %.10g", w, p/w, tc+nc*p}')
echo "--- PHASE 2: limit cycle with fields ---"
echo "  period=$PERIOD  writeInterval=$WI  (${NPC} samples/cycle)  t: $T_CONV -> $T2  ($NCYCLES cycles)"
patch_cd latestTime $DT $T2 runTime $WI
run log.phase2
sawtooth_gate $T_CONV "phase 2"

# ----------------------------------------------------------- reconstruct ---
# Only the phase-2 window: the phase-1 chunk-end writes are restart points, not
# snapshots, and must not end up in the POD set.
echo "--- reconstructPar (t > $T_CONV) ---"
reconstructPar -time "${T_CONV}:" > log.reconstructPar 2>&1 || true
tail -3 log.reconstructPar
echo "snapshots: $(ls -d [0-9]*.[0-9]* 2>/dev/null | wc -l)   $(du -sh . | cut -f1)"
echo "$(date '+%F %T')  ALL DONE  case=$CASE"
