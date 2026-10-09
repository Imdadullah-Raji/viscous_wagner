#!/bin/bash
# Remote half of the plate LSB sweep (oct9_workflow.md), for banana-cfd (OpenFOAM 13, 6 cores).
# Copied there as ~/OpenFOAM/raji/lsb_plate_fine_re500/ without processor*/log.*.
# Same sequence as static_airfoil/run_steady.sh: one 6-rank job at a time.
cd "$(dirname "$0")" || exit 1
source /opt/openfoam13/etc/bashrc >/dev/null 2>&1
NP=6
[ $# -gt 0 ] || set -- 8.1 8.2 8.3
for A in "$@"; do
  C=steady_aoa$A
  ( cd "$C"
    [ -d processor0 ] || decomposePar > log.decomposePar 2>&1
    mpirun -np $NP foamRun -parallel > log.steady 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1
    foamPostProcess -func writeCellCentres -latestTime > log.writeCellCentres 2>&1 )
  ok="ran"; grep -qa "^End" "$C/log.steady" || ok="NO-End"
  printf '%s  AoA %-3s %s  %s\n' "$(date '+%H:%M:%S')" "$A" "$ok" \
         "$(tail -1 "$C"/postProcessing/forceCoeffs/0/*.dat 2>/dev/null \
            | awk '{printf "iter=%s Cd=%.6f Cl=%.6f", $1, $3, $4}')"
done
echo "$(date '+%F %T')  REMOTE SWEEP COMPLETE"
