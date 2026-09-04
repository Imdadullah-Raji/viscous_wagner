#!/bin/bash
# Two-phase FIXED-dt run. No adjustTimeStep anywhere - that controller is what
# corrupted the first attempt (../FINDINGS.md Part 4).
#
# Phase 1: dt=4e-5 through the impulsive start (measured peak Co = 0.699).
#          A single large fixed dt from uniform flow is what diverged
#          diag/a5_fixeddt v1 (Co -> 1e103), so this phase is deliberately small.
# Phase 2: dt sized from the MEASURED Courant number at the end of phase 1,
#          targeting Co ~ 0.6. Cell-size estimates were off by 2.3x, so this
#          must be measured, not predicted.
set -e
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
C=/home/raji/Research/Thesis/static_airfoil/aoa10
cd $C
DT1=4e-5
TARGET_CO=0.6

patch_cd() {  # startFrom deltaT endTime writeInterval
  sed -i -e "s/^startFrom .*/startFrom       $1;/" \
         -e "s/^deltaT .*/deltaT          $2;/" \
         -e "s/^endTime .*/endTime         $3;/" \
         -e "s/^writeInterval .*/writeInterval   $4;/" \
         -e "s/^adjustTimeStep .*/adjustTimeStep  no;/" system/controlDict
}

echo "=== PHASE 1: dt=$DT1 -> t=0.5 ==="
patch_cd startTime $DT1 0.5 0.5
mpirun -np 6 foamRun -parallel > log.phase1 2>&1
grep -q "^End" log.phase1 || { echo "PHASE 1 FAILED"; tail -20 log.phase1; exit 1; }

CO=$(grep -E "^Courant" log.phase1 | tail -1 | awk '{print $6}')
DT2=$(awk -v d=$DT1 -v c=$CO -v t=$TARGET_CO 'BEGIN{
  x=d*t/c; if(x>2.5e-4)x=2.5e-4; if(x<8e-5)x=8e-5; printf "%.3e", x}')
echo "PHASE 1 done. Co_end=$CO at dt=$DT1  ->  phase-2 dt=$DT2 (target Co=$TARGET_CO)"

# Gate: do not burn hours if the 2dt sawtooth is already present.
R=$(bash ../analysis/sawtooth_check.sh $C 0.2 | grep -o 'step-reversals=[0-9.]*' | cut -d= -f2)
echo "PHASE 1 sawtooth check: step-reversals=${R}%"
awk -v r="$R" 'BEGIN{ if (r+0 > 60) exit 1 }' || { echo "ABORT: sawtooth after phase 1"; exit 1; }

echo "=== PHASE 2: dt=$DT2 -> t=30 ==="
patch_cd latestTime $DT2 30 5
mpirun -np 6 foamRun -parallel > log.phase2 2>&1
grep -q "^End" log.phase2 || { echo "PHASE 2 FAILED"; tail -20 log.phase2; exit 1; }

echo "=== reconstructPar ==="
reconstructPar -latestTime > log.reconstructPar 2>&1
tail -3 log.reconstructPar
echo "ALL DONE"
