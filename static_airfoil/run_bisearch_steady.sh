#!/bin/bash
# run_bisearch_steady.sh <aoa> [<aoa> ...]
# Steady (SIMPLE) base flows for the bifurcation search, into NACA0012/bi_point_search/.
# Sequential, 6 MPI ranks each -- six concurrent serial jobs measured 4.2x slower
# per case (they each hold the whole mesh and saturate memory bandwidth).
cd "$(dirname "$0")"
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
MESH=meshes/M2_ns160_nr180_b025.msh
NP=6; NITER=6000
for A in "$@"; do
  C=NACA0012/bi_point_search/steady_$A
  if [ -d "$C/processor0" ]; then echo "$(date '+%H:%M:%S')  steady $A already built, skipping"; continue; fi
  ./make_steady.sh "$C" 500 "$A" "$MESH" $NITER > /dev/null 2>&1 || { echo "build FAILED $A"; continue; }
  ( cd "$C"
    decomposePar > log.decomposePar 2>&1
    mpirun -np $NP foamRun -parallel > log.steady 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1 )
  ok=ran; grep -qa "^End" "$C/log.steady" || ok=NO-End
  printf '%s  steady %-5s %s  %s\n' "$(date '+%H:%M:%S')" "$A" "$ok" \
    "$(tail -1 "$C"/postProcessing/forceCoeffs/0/*.dat 2>/dev/null \
       | awk '{printf "Cl=%.6f Cd=%.6f", $4, $3}')"
done
echo "$(date '+%F %T')  BISEARCH STEADY COMPLETE: $*"
