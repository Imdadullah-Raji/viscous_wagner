#!/bin/bash
# run_batch.sh <openfoam bashrc> <ranks> <dir> [<dir> ...]
# Zero-thickness plate sweep on any machine: runs every steady_aoa* in each <dir>, one
# <ranks>-rank job at a time, lowest angle first; skips cases that already ended ("End" in
# log.steady).  Cases are built (mesh, baffles, checkMesh) on the main machine by
# build_sweep.sh, so only decomposePar / foamRun / reconstructPar run here.  Force histories
# are copied back and judged there with analysis/steady_summary.py.
# ONLY="<a> <a> ..." restricts every listed folder to those angles (load rebalancing).
# WARM=1 seeds each case from the converged case at the nearest angle (steady fixed point unchanged,
# the expensive cold start skipped).
cd "$(dirname "$0")" || exit 1
source "$1" >/dev/null 2>&1; NP=$2; shift 2
for D in "$@"; do
  for C in $(ls -d $D/steady_aoa* | awk -F'steady_aoa' '{print $NF, $0}' | sort -g | cut -d' ' -f2); do
    [ -d "$C" ] || continue
    if [ -n "$ONLY" ]; then case " $ONLY " in *" ${C##*steady_aoa} "*) ;; *) continue ;; esac; fi
    grep -qa "^End" "$C/log.steady" 2>/dev/null && continue
    if [ -n "$WARM" ] && [ ! -d "$C/processor0" ]; then
      # WARM=1: seed from the finished case at the nearest angle in the same folder (same mesh, so
      # mapFields -consistent is a 1:1 copy; target BCs, i.e. the freestream angle, are kept).
      a=${C##*steady_aoa}; best=""; bd=1e9
      for S in $D/steady_aoa*; do
        [ "$S" = "$C" ] && continue
        grep -qa "^End" "$S/log.steady" 2>/dev/null || continue
        t=$(ls -d $S/[0-9]* 2>/dev/null | xargs -n1 basename | sort -g | tail -1); [ "$t" = 0 ] && continue
        d=$(awk -v x=${S##*steady_aoa} -v y=$a 'BEGIN{d=x-y; print (d<0?-d:d)}')
        if awk -v d=$d -v b=$bd 'BEGIN{exit !(d<b)}'; then best=$S; bd=$d; bt=$t; fi
      done
      if [ -n "$best" ]; then
        mapFields "$best" -case "$C" -consistent -sourceTime $bt > "$C/log.mapFields" 2>&1 \
          && echo "$(date '+%H:%M:%S')  $C seeded from $best ($bt)" || echo "$(date '+%H:%M:%S')  $C mapFields FAILED, cold start"
      fi
    fi
    ( cd "$C"
      sed -i "s/^numberOfSubdomains .*/numberOfSubdomains  $NP;/" system/decomposeParDict
      [ -d processor0 ] || decomposePar > log.decomposePar 2>&1
      mpirun -np $NP foamRun -parallel > log.steady 2>&1
      reconstructPar -latestTime > log.reconstructPar 2>&1 )
    ok="ran"; grep -qa "^End" "$C/log.steady" || ok="NO-End"
    printf '%s  %s %s  %s\n' "$(date '+%H:%M:%S')" "$C" "$ok" \
           "$(tail -1 "$C"/postProcessing/forceCoeffs/0/*.dat 2>/dev/null | awk '{printf "iter=%s Cd=%.6f Cl=%.6f", $1, $3, $4}')"
  done
  echo "$(date '+%F %T')  DONE $D"
done
echo "$(date '+%F %T')  BATCH COMPLETE: $*"
