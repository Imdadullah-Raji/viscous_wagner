#!/bin/bash
# [STEADY_DIR=NACA0012/steady_re1000] run_steady.sh <aoa> [<aoa> ...]
#
# Cases live in $STEADY_DIR/steady_aoa<A> (default NACA0012/steady_re500).
#
# Runs steady (SIMPLE) cases SEQUENTIALLY at 6 MPI ranks each.
#
# NOT six concurrent serial jobs, which was the obvious thing and was 4.2x
# slower PER CASE: measured 0.254 s/iter solo against 1.055 s/iter with six
# running, i.e. 5.7 iter/s aggregate versus ~20 iter/s for a single 6-rank job.
# Six independent processes each hold the whole 63k-cell mesh and compete for
# memory bandwidth; MPI splits it into ~10.5k-cell subdomains that stay in
# cache.  Same cores, same total work, several times the throughput.
#
# NO `set -u`: /opt/openfoam-dev/etc/bashrc references unset variables and would
# exit the script the moment it is sourced under nounset.
#
# Convergence is NOT assumed -- checked afterwards on the C_L trace itself
# (analysis/steady_summary.py), not on solver residuals alone, because a
# residual can stall low while the force is still creeping.
cd "$(dirname "$0")/${STEADY_DIR:-NACA0012/steady_re500}" || exit 1
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
NP=6

for A in "$@"; do
  C=steady_aoa$A
  ( cd "$C"
    [ -d processor0 ] || decomposePar > log.decomposePar 2>&1
    mpirun -np $NP foamRun -parallel > log.steady 2>&1
    reconstructPar -latestTime > log.reconstructPar 2>&1 )
  ok="ran"; grep -qa "^End" "$C/log.steady" || ok="NO-End"
  printf '%s  AoA %-3s %s  %s\n' "$(date '+%H:%M:%S')" "$A" "$ok" \
         "$(tail -1 "$C"/postProcessing/forceCoeffs/0/*.dat 2>/dev/null \
            | awk '{printf "iter=%s Cd=%.6f Cl=%.6f", $1, $3, $4}')"
done
echo "$(date '+%F %T')  STEADY SWEEP COMPLETE: $*"
