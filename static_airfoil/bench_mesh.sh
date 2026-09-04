#!/bin/bash
# bench_mesh.sh <mesh.msh> <label> -- build a throwaway case and measure the two
# numbers that set the wall-clock budget: Courant at a reference dt, and s/step.
set -e
SCRATCH=/tmp/claude-1000/-home-raji-Research-Thesis-static-airfoil/59d8fa6c-b299-42a0-849f-f1e49516eb73/scratchpad
C=$SCRATCH/bench_$2
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
cd /home/raji/Research/Thesis/static_airfoil
./make_case.sh $C 500 20 $1 > $SCRATCH/mk_$2.log 2>&1
CELLS=$(grep -m1 'cells:' $C/log.checkMesh | awk '{print $2}')
NONORTH=$(grep -m1 'non-orthogonality' $C/log.checkMesh | awk '{print $4}')
cd $C
decomposePar > log.dp 2>&1
sed -i -e 's/^endTime .*/endTime         0.008;/' -e 's/^writeInterval .*/writeInterval   1000;/' system/controlDict
S=$(date +%s%N); mpirun -np 4 foamRun -parallel > log.b 2>&1; E=$(date +%s%N)
N=$(grep -c "^Time = " log.b)
CO=$(grep -E "^Courant" log.b | tail -1 | awk '{print $6}')
SPS=$(awk -v e=$E -v s=$S -v n=$N 'BEGIN{printf "%.4f",(e-s)/1e9/n}')
# dt at Co=0.6, and hours to t=180 on 4 cores
awk -v l="$2" -v c="$CELLS" -v no="$NONORTH" -v co="$CO" -v sps="$SPS" 'BEGIN{
  dt=4e-5*0.6/co; h=sps*(180/dt)/3600;
  printf "%-22s cells=%-7s nonOrth=%-6s Co@4e-5=%-9.4f dt@Co0.6=%-9.2e %6.3f s/step  -> %6.1f h to t=180\n", l,c,no,co,dt,sps,h}'
