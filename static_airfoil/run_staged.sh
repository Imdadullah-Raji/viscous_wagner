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
# Re-tuning targets. TARGET_CO above is applied to the SETTLED Courant number
# measured at t=0.5 (a typical value). The re-tune below instead sees a chunk's
# MAXIMUM Courant, so it needs its own, higher target -- reusing 0.35 there drove
# the maximum to 0.35 and the median to 0.19, costing 2.1x in wall time for no
# accuracy anyone asked for.
RETUNE_CO=0.75      # re-tune aims to put the chunk MAXIMUM here
CO_HI=0.80          # shrink dt if a chunk's max exceeds this
CO_LO=0.45          # grow dt if a chunk's max stays under this.
                    # Must sit below RETUNE_CO or the controller has a dead zone:
                    # with CO_LO=0.15 a shrink landed max Co at 0.34, too high to
                    # ever trigger growth again, so dt could never recover.
# DT_MAX was 1.0e-3 until AoA 18 diverged at exactly that value (2026-09-07).
# The chunk before it ran clean at 5.573e-4 with Co_max=0.417; the controller
# saw 0.417 < CO_LO, scaled dt by 0.75/0.417 = 1.80, and hit this clamp.  The
# next chunk blew up 191 steps later -- at Co_max=0.60, i.e. BELOW both
# RETUNE_CO and CO_HI.  Co never saw it coming and could not have: mean/max Co
# was 1/166, so the Courant number is set by a handful of trailing-edge cells
# and says nothing about stability in the 64-degree non-orthogonal cells where
# the lagged non-orthogonal correction actually fails.  The ceiling is
# therefore an absolute dt, not a Courant target.  6.0e-4 sits below 6.331e-4,
# the largest dt ever run to convergence on this mesh (AoA 26).
DT_MAX=6.0e-4
DT_GROW_MAX=1.5    # cap on one re-tune's growth ratio, so the controller walks
                   # up to an unstable dt with a recoverable overshoot instead
                   # of leaping past every value it has evidence for
DT_MIN=2.0e-5
CHUNK=20            # phase-1 chunk length (time units)
T1_MIN=60           # no convergence verdict before this
T1_MAX=240          # give up extending the transient here
NCYCLES=20
PER_CYCLE=32
# TOL is overridable:  TOL=0.002 ./run_staged.sh <case> 6
# The default 0.01 is right for a case that settles quickly, but near the
# shedding onset the approach to the limit cycle is slow (critical slowing
# down) and 1% window-to-window agreement can be reached while the amplitude is
# still creeping.  Measured at AoA 12: the gate passed at t=60.5 with
# D_P2P=0.00586 while the envelope was still growing at +3.2e-3 /t -- 300x the
# calibrated saturation floor of ~1e-5 -- which would have put a ~4% amplitude
# drift through the phase-2 recording.  Tighten it rather than override the
# gate by hand; overriding it by eye got AoA 14 wrong in the other direction.
TOL=${TOL:-0.01}

# DT_FIXED -- hold dt at this value and disable the Courant re-tune entirely.
#   DT_FIXED=5.5726e-4 ./run_staged.sh NACA0012/re500_aoa18 6
# For a POD run the re-tune is free money: dt only has to keep Co in band, and
# a change between chunks costs one Euler step.  For a STABILITY measurement it
# is not free.  analysis/growth_rate.py fits an exponential to the C_L envelope
# over the whole phase-1 record, and the numerical damping of the scheme is
# dt-dependent, so a re-tune partway through the record puts a step in the very
# quantity being fitted.  Freezing dt costs some wall time and buys a record
# whose damping is homogeneous end to end.
DT_FIXED=${DT_FIXED:-}

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

# OpenFOAM names a time directory from the VALUE, so a write that lands on a
# whole number is "213", not "213.0".  The obvious glob [0-9]*.[0-9]* silently
# skips those -- it undercounted this case by one and made a perfectly uniform
# snapshot series look like it had a 2x gap at t=213.  Match on the value.
count_snaps() {   # t_min -> number of time dirs strictly after t_min
  ls -d [0-9]* 2>/dev/null \
    | awk -v t="$1" '/^[0-9]+(\.[0-9]+)?$/ && $1+0 > t+1e-9' | wc -l
}
run() {        # logfile
  mpirun -np $NP foamRun -parallel > "$1" 2>&1 || true
  grep -q "^End" "$1" || { echo "FAILED: $1"; tail -30 "$1"; exit 1; }
}

# Same thing, but returns 1 instead of exiting, and QUARANTINES what the failed
# attempt wrote.  A diverged chunk leaves postProcessing/forceCoeffs/<restart>/
# holding the blow-up -- values of order 1e200.  A retry restarts from the same
# time and so APPENDS to that same directory, and analysis/limit_cycle.py
# concatenates every directory it finds.  Left in place it poisons C_L for the
# rest of the campaign; it reported CL_MEAN=2.5e206 for AoA 18.  Move it aside
# before retrying, never after.
try_run() {    # logfile  restartTime
  mpirun -np $NP foamRun -parallel > "$1" 2>&1 || true
  grep -q "^End" "$1" && return 0
  echo "    DIVERGED: $1"
  mkdir -p .diverged
  local stamp=$(date +%H%M%S)
  [ -d "postProcessing/forceCoeffs/$2" ] &&
    mv "postProcessing/forceCoeffs/$2" ".diverged/forceCoeffs_${2}_$stamp"
  mv "$1" ".diverged/$(basename $1).$stamp"
  return 1
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

# Snapped so CHUNK is an exact integer number of steps, which keeps every chunk
# boundary -- and therefore t_conv -- exactly on a time step.
snap_dt() {  # target_dt -> dt that divides CHUNK exactly
  awk -v x=$1 -v hi=$DT_MAX -v lo=$DT_MIN -v ch=$CHUNK 'BEGIN{
    if(x>hi)x=hi; if(x<lo)x=lo;
    n=int(ch/x); if(n*x<ch)n++; printf "%.10g", ch/n}'
}

# ---------------------------------------------------------------- resume ---
# Phase 1 writes one field snapshot per chunk purely as a restart point, so an
# interrupted run (killed terminal, power cut, Ctrl-C) can always be picked up at
# the last chunk boundary instead of starting over. Without this the driver would
# re-run phase 0 from t=0 and silently discard both the progress AND any warm
# start that seeded time 0.
LATEST=$(ls -d processor0/[0-9]* 2>/dev/null | sed 's|.*/||' | sort -g | tail -1)

# Interrupted during phase 2?  Finish it on the ORIGINAL dt, writeInterval and
# endTime, so the snapshots either side of the interruption stay on one uniform
# grid and T_CONV keeps its real value.
if [ -f .phase2 ] && [ -n "$LATEST" ]; then
  . ./.phase2
  if awk -v l="$LATEST" -v t="$T_CONV" 'BEGIN{exit !(l>t+1e-9)}'; then
    echo "--- RESUME PHASE 2: snapshot at t=$LATEST, t_conv=$T_CONV, endTime=$T2 ---"
    echo "    period-$ORDER orbit, dt=$DT, writeInterval=$WI (both held fixed)"
    if awk -v l="$LATEST" -v e="$T2" 'BEGIN{exit !(l<e-1e-9)}'; then
      patch_cd latestTime $DT $T2 runTime $WI
      run log.phase2_resume
      sawtooth_gate "$LATEST" "phase 2 resume"
    else
      echo "    already at endTime -- reconstructing only"
    fi
    echo "--- reconstructPar (new times only, t > $T_CONV) ---"
    reconstructPar -newTimes -time "${T_CONV}:" > log.reconstructPar 2>&1 || true
    tail -3 log.reconstructPar
    echo "snapshots: $(count_snaps $T_CONV)   $(du -sh . | cut -f1)"
    python3 "$ROOT/analysis/limit_cycle.py" .
    echo "$(date '+%F %T')  ALL DONE (resumed phase 2)  case=$CASE"
    exit 0
  fi
fi

if [ -n "$LATEST" ] && awk -v l="$LATEST" -v t=$T0 'BEGIN{exit !(l>t+1e-9)}'; then
  T=$LATEST
  DT=$(awk '/^deltaT/{print $2}' system/controlDict | tr -d ';')
  [ -n "$DT_FIXED" ] && DT=$DT_FIXED
  echo "--- RESUME: restart point at t=$T (skipping phase 0), dt=$DT ---"
  # Apply the same Courant re-tune the chunk loop would have applied, using the
  # last completed chunk's log. On resume after a controller fix this is what
  # lets dt recover instead of inheriting a stale, over-conservative value.
  LASTLOG=$(ls -t log.phase1_t* 2>/dev/null | head -1)
  if [ -n "$DT_FIXED" ]; then
    echo "    DT_FIXED=$DT_FIXED -- re-tune disabled"
  elif [ -n "$LASTLOG" ]; then
    COC=$(grep -E "^Courant" "$LASTLOG" | awk '{print $6}' | sort -g | tail -1)
    if awk -v c="$COC" -v h=$CO_HI 'BEGIN{exit !(c>h)}'; then
      DT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$RETUNE_CO 'BEGIN{print d*t/c}'))
      echo "    last chunk Co_max=$COC > $CO_HI  ->  dt -> $DT"
    elif awk -v c="$COC" -v l=$CO_LO 'BEGIN{exit !(c<l)}'; then
      DT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$RETUNE_CO 'BEGIN{print d*t/c}'))
      echo "    last chunk Co_max=$COC < $CO_LO  ->  dt -> $DT"
    else
      echo "    last chunk Co_max=$COC in band [$CO_LO,$CO_HI], dt unchanged"
    fi
  fi
else
  # ------------------------------------------------------------- phase 0 ---
  echo "--- PHASE 0: impulsive start, dt=$DT0 -> t=$T0 (no fields) ---"
  patch_cd startTime $DT0 $T0 runTime $T0
  S=$(date +%s); run log.phase0; E=$(date +%s)
  sawtooth_gate 0.2 "phase 0"
  CO=$(grep -E "^Courant" log.phase0 | tail -1 | awk '{print $6}')
  NSTEP=$(awk -v t=$T0 -v d=$DT0 'BEGIN{printf "%.0f", t/d}')  # awk, not bc: bc cannot parse "4e-5"
  SPS=$(awk -v e=$E -v s=$S -v n=$NSTEP 'BEGIN{printf "%.4f",(e-s)/n}')
  echo "phase 0 done in $((E-S))s   Co_max=$CO   ${SPS}s/step"

  # dt for the rest of the run, from the MEASURED Courant number at t=0.5 (the
  # settled value -- sampling during the impulsive spike, which peaks ~14x higher
  # around t=1e-3, gives a dt roughly 7x too small).
  DT=$(snap_dt $(awk -v d=$DT0 -v c=$CO -v t=$TARGET_CO 'BEGIN{print d*t/c}'))
  if [ -n "$DT_FIXED" ]; then
    echo "phase-1/2 dt = $DT_FIXED  (DT_FIXED; Co-derived value would have been $DT)"
    DT=$DT_FIXED
  else
    echo "phase-1/2 dt = $DT  (target Co=$TARGET_CO from measured Co=$CO)"
  fi
  awk -v s=$SPS -v d=$DT -v tm=$T1_MAX 'BEGIN{
    printf "rough ETA: %.1f h to t=%g, plus phase 2\n", s*(tm/d)/3600, tm}'
  T=$T0
fi

# --------------------------------------------------------------- phase 1 ---
echo "--- PHASE 1: transient, field-free, chunks of $CHUNK ---"
CONVERGED=no
while awk -v t=$T -v m=$T1_MAX 'BEGIN{exit !(t<m)}'; do
  T=$(awk -v t=$T -v c=$CHUNK -v t0=$T0 'BEGIN{printf "%.10g", (int((t-t0)/c)+1)*c+t0}')
  TPREV=$(awk -v t=$T -v c=$CHUNK 'BEGIN{printf "%.10g", t-c}')
  # Bounded retry.  A chunk that diverges has not lost anything -- the restart
  # point at TPREV is the previous chunk's write and is untouched -- so the
  # cheap response is to halve dt and go again rather than end the run and wait
  # for a human.  Bounded at 2 halvings (dt/4) so a genuinely sick case still
  # stops instead of grinding down to DT_MIN.
  ATTEMPT=0
  while :; do
    patch_cd latestTime $DT $T runTime $T       # interval==endTime -> one write
    echo "  running to t=$T  (dt=$DT) ..."
    try_run log.phase1_t$T $TPREV && break
    ATTEMPT=$((ATTEMPT+1))
    [ $ATTEMPT -le 2 ] || { echo "FAILED: diverged twice at t=$T; giving up"; exit 1; }
    DT=$(snap_dt $(awk -v d=$DT 'BEGIN{print d/2}'))
    echo "    retry $ATTEMPT/2 at dt=$DT"
  done

  # Re-tune dt from what the chunk actually did.  Only ever BETWEEN chunks, so
  # dt stays rigorously fixed within any single foamRun -- the discipline the
  # variable-dt controller broke.  Each change costs one Euler step out of
  # ~100k, which is why this is cheap enough to do at all.
  sawtooth_gate $(awk -v t=$T -v c=$CHUNK 'BEGIN{print t-c}') "chunk to t=$T"

  COC=$(grep -E "^Courant" log.phase1_t$T | awk '{print $6}' | sort -g | tail -1)
  NEWDT=$DT
  if [ -n "$DT_FIXED" ]; then
    echo "    Co_max=$COC  (dt held fixed at $DT)"
  elif awk -v c=$COC -v h=$CO_HI 'BEGIN{exit !(c>h)}'; then
    NEWDT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$RETUNE_CO 'BEGIN{print d*t/c}'))
    echo "    Co_max=$COC > $CO_HI  ->  dt $DT -> $NEWDT  (target max Co=$RETUNE_CO)"
  elif awk -v c=$COC -v l=$CO_LO -v t=$T -v m=$T1_MIN 'BEGIN{exit !(c<l && t>=m/3)}'; then
    NEWDT=$(snap_dt $(awk -v d=$DT -v c=$COC -v t=$RETUNE_CO -v g=$DT_GROW_MAX 'BEGIN{
                        x=d*t/c; if(x>d*g)x=d*g; print x}'))
    echo "    Co_max=$COC < $CO_LO  ->  dt $DT -> $NEWDT  (target max Co=$RETUNE_CO)"
  fi
  DT=$NEWDT

  if awk -v t=$T -v m=$T1_MIN 'BEGIN{exit !(t>=m)}'; then
    eval $(python3 "$ROOT/analysis/limit_cycle.py" . --tol $TOL)
    echo "    t=$T  shedding=$SHEDDING  period-$ORDER  f_orbit=$FREQ" \
         "f_carrier=$CARRIER_FREQ  Cl=$CL_MEAN+-$CL_P2P" \
         "d(mean,p2p,f,sub)=($D_MEAN,$D_P2P,$D_FREQ,$D_SUB)  converged=$CONVERGED"
    [ "$CONVERGED" = yes ] && break
  fi
done

T_CONV=$T
if [ "$CONVERGED" != yes ]; then
  echo "WARNING: no limit cycle by t=$T1_MAX; proceeding with the last estimate."
  eval $(python3 "$ROOT/analysis/limit_cycle.py" . --tol $TOL)
fi
[ "$SHEDDING" = yes ] || { echo "ABORT: no shedding detected -- steady solution?"; exit 1; }

# NO_FIELDS=1 -- stop here: forces only, no phase-2 snapshots.  For a run whose
# only question is "does it shed, and what is the limit cycle" (C_L/C_D are
# complete in postProcessing/forceCoeffs/).  The phase-1 chunk writes remain as
# restart points; extend_phase2.sh is NOT applicable since there is no phase 2.
if [ -n "${NO_FIELDS:-}" ]; then
  echo "$(date '+%F %T')  NO_FIELDS: stopping after phase 1 at t=$T_CONV" \
       "(converged=$CONVERGED, period-$ORDER, f_carrier=$CARRIER_FREQ)"
  exit 0
fi

# --------------------------------------------------------------- phase 2 ---
# Two different periods, and they are NOT interchangeable once the wake has
# period-doubled (ORDER>1, see analysis/limit_cycle.py):
#
#   SAMPLE_PERIOD (the carrier) sets the WRITE INTERVAL.  It is the fastest
#     oscillation in the flow, so it is what the snapshot rate has to resolve.
#     PER_CYCLE samples per carrier period, whatever the order.
#   PERIOD (the fundamental, = ORDER x carrier) sets the RUN LENGTH.  It is the
#     time for the orbit to actually close, so NCYCLES of it is NCYCLES
#     statistically independent realisations.  Using the carrier here would
#     record NCYCLES/ORDER orbits -- half the data at ORDER=2.
#
# Snap the write interval DOWN to a whole number of time steps.  Rounding down
# means at least PER_CYCLE samples per cycle, never fewer.
read WI NPC T2 <<< $(awk -v sp=$SAMPLE_PERIOD -v p=$PERIOD -v n=$PER_CYCLE \
                         -v d=$DT -v tc=$T_CONV -v nc=$NCYCLES 'BEGIN{
  m=int((sp/n)/d); if(m<1)m=1; w=m*d;
  printf "%.10g %.2f %.10g", w, sp/w, tc+nc*p}')
echo "--- PHASE 2: limit cycle with fields ---"
echo "  period-$ORDER orbit: T_orbit=$PERIOD  T_carrier=$SAMPLE_PERIOD"
echo "  writeInterval=$WI  (${NPC} samples/carrier period)"
echo "  t: $T_CONV -> $T2  ($NCYCLES orbits)"
# Record what phase 2 is doing before it starts.  Phase-1 chunk boundaries and
# phase-2 snapshots are both just time directories -- nothing on disk
# distinguishes them -- so a resume that finds a phase-2 snapshot would treat it
# as a chunk boundary and set T_CONV to the interruption time, silently
# recording the wrong number of orbits.  This box loses power; that is not
# hypothetical.
printf 'T_CONV=%s\nDT=%s\nWI=%s\nT2=%s\nORDER=%s\nPERIOD=%s\nSAMPLE_PERIOD=%s\n' \
       "$T_CONV" "$DT" "$WI" "$T2" "$ORDER" "$PERIOD" "$SAMPLE_PERIOD" > .phase2
patch_cd latestTime $DT $T2 runTime $WI
run log.phase2
sawtooth_gate $T_CONV "phase 2"

# ----------------------------------------------------------- reconstruct ---
# Only the phase-2 window: the phase-1 chunk-end writes are restart points, not
# snapshots, and must not end up in the POD set.
echo "--- reconstructPar (t > $T_CONV) ---"
reconstructPar -time "${T_CONV}:" > log.reconstructPar 2>&1 || true
tail -3 log.reconstructPar
echo "snapshots: $(count_snaps 0)   $(du -sh . | cut -f1)"
echo "$(date '+%F %T')  ALL DONE  case=$CASE"
