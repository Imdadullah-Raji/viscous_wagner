#!/bin/bash
# run_bisearch_pert.sh <aoa> [<aoa> ...]
#
# Perturbation runs for the bifurcation search.  For each angle alpha:
#   1. build a TRANSIENT case at alpha
#   2. seed it with the converged steady field at alpha+0.5, so the initial
#      condition is  x_s(alpha) + [x_s(alpha+0.5) - x_s(alpha)]  -- a small,
#      smooth, physical perturbation ABOUT THE ALPHA FIXED POINT
#   3. integrate in chunks, dt frozen, forceCoeffs only, NO field snapshots
#   4. after each chunk print sigma from analysis/eigenvalue.py
#
# dt is frozen at the onset-sweep value so sigma is comparable across angles and
# against re500_aoa14/16/18.  A re-tune mid-record would put a step into the
# envelope being fitted -- the scheme's numerical damping is dt-dependent.
#
# THERE IS A STARTUP STAGE, AND SKIPPING IT DIVERGED THE FIRST ATTEMPT.
# The reasoning that looked right -- "the initial field is a converged steady
# solution, not uniform flow, so no impulsive start is needed" -- is wrong.
# mapFields transfers VOLUME FIELDS ONLY, so phi is not transferred; the solver
# rebuilds it from the mapped U and the first pressure solve has to absorb the
# divergence error.  That adjustment is impulsive whatever the field looks like.
# Run at the measurement dt (5.57e-4, 14x the proven startup dt) C_L went
# 0.4585 -> -0.29 in ONE step and the case diverged by t=0.093 with Co ~ 1e103.
# Every warm start in this repo survives because run_staged.sh always opens with
# phase 0 at dt=4e-5; this now does the same.
set -e
cd "$(dirname "$0")"
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
MESH=meshes/M2_ns160_nr180_b025.msh
NP=6
DT=0.0005572582892
CHUNK=20
TMAX=${TMAX:-60}
D=bi_point_search

for A in "$@"; do
  SRC=$(awk -v a="$A" 'BEGIN{printf "%.1f", a+0.5}')
  C=$D/pert_$A
  echo "=============================================================="
  echo "$(date '+%F %T')  AoA $A   perturbation = x_s($SRC) - x_s($A)"
  [ -d "$D/steady_$SRC" ] || { echo "ABORT: no base flow $D/steady_$SRC"; continue; }

  CLS=$(tail -1 $D/steady_$A/postProcessing/forceCoeffs/0/*.dat | awk '{print $4}')
  CLP=$(tail -1 $D/steady_$SRC/postProcessing/forceCoeffs/0/*.dat | awk '{print $4}')
  echo "  Cl_steady($A)=$CLS   Cl_steady($SRC)=$CLP   perturbation dCl=$(awk -v a=$CLS -v b=$CLP 'BEGIN{printf "%.6f", b-a}')"

  if [ ! -d "$C/processor0" ]; then
    ./make_case.sh "$C" 500 "$A" "$MESH" > "$C.build.log" 2>&1
    ST=$(ls -d $D/steady_$SRC/[0-9]* | xargs -n1 basename | sort -g | tail -1)
    ( cd "$C" && mapFields "../steady_$SRC" -consistent -sourceTime "$ST" > log.mapFields 2>&1 )
    grep -qa "^End" "$C/log.mapFields" || { echo "  mapFields FAILED"; tail -5 "$C/log.mapFields"; continue; }
    echo "  seeded from $D/steady_$SRC @ $ST"
    ( cd "$C" && decomposePar > log.decomposePar 2>&1 )
  else
    echo "  (already seeded -- resuming)"
  fi

  # --- startup: rebuild phi gently at the proven phase-0 dt -------------
  T=$(ls -d $C/processor0/[0-9]* 2>/dev/null | xargs -n1 basename | sort -g | tail -1)
  T=${T:-0}
  if awk -v t="$T" 'BEGIN{exit !(t<0.5-1e-9)}'; then
    ( cd "$C" && sed -i -e 's/^startFrom .*/startFrom       latestTime;/' \
        -e 's/^deltaT .*/deltaT          4e-05;/' -e 's/^endTime .*/endTime         0.5;/' \
        -e 's/^writeControl .*/writeControl    runTime;/' -e 's/^writeInterval   .*/writeInterval   0.5;/' \
        -e 's/^adjustTimeStep .*/adjustTimeStep  no;/' system/controlDict
      mpirun -np $NP foamRun -parallel > log.startup 2>&1 || true )
    grep -qa "^End" "$C/log.startup" || { echo "  STARTUP FAILED"; tail -5 "$C/log.startup"; continue; }
    CO=$(grep -a "^Courant" "$C/log.startup" | awk '{print $6}' | sort -g | tail -1)
    CLNOW=$(tail -1 $C/postProcessing/forceCoeffs/*/*.dat | awk '{print $4}')
    echo "  startup to t=0.5 at dt=4e-5: Co_max=$CO  Cl=$CLNOW"
    awk -v c="$CO" 'BEGIN{exit !(c>1)}' && { echo "  ABORT: startup Co>1"; continue; }
    T=0.5
  fi
  while awk -v t="$T" -v m="$TMAX" 'BEGIN{exit !(t<m-1e-9)}'; do
    T=$(awk -v t="$T" -v c=$CHUNK 'BEGIN{printf "%.10g", (int((t-0.5)/c)+1)*c+0.5}')
    ( cd "$C" && sed -i -e 's/^startFrom .*/startFrom       latestTime;/' \
        -e "s/^deltaT .*/deltaT          $DT;/" -e "s/^endTime .*/endTime         $T;/" \
        -e 's/^writeControl .*/writeControl    runTime;/' -e "s/^writeInterval   .*/writeInterval   $T;/" \
        -e 's/^adjustTimeStep .*/adjustTimeStep  no;/' system/controlDict
      mpirun -np $NP foamRun -parallel > log.chunk_$T 2>&1 || true )
    grep -qa "^End" "$C/log.chunk_$T" || { echo "  FAILED at t=$T"; break; }
    echo "  t=$T  $(python3 analysis/eigenvalue.py "$C" --cl-steady "$CLS" 2>/dev/null \
        | awk -F= '/^SIGMA_RMS|^R2_RMS|^FITSPAN_RMS|^OSC_FRAC|^FREQ/{printf "%s=%s  ", $1, $2}')"
  done
  echo "$(date '+%F %T')  AoA $A done"
  python3 analysis/eigenvalue.py "$C" --cl-steady "$CLS" --table || true
done
echo "$(date '+%F %T')  BISEARCH PERT COMPLETE: $*"
