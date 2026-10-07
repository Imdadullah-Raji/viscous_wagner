#!/bin/bash
# NACA0006 steady (SIMPLE) sweep at Re = 500, integer angles 0-20 deg: thickness test of
# the domain of the theta + theta^3 fit.  Same M2 recipe (h1 = 0.05/sqrt(Re)), same
# make_steady.sh, only the section changes (meshing/mesher.py --naca 0006).
# Strictly sequential, 6 MPI ranks per case.
cd "$(dirname "$0")/.."
MSH=meshes/M2_ns160_nr180_b025_naca0006_re500.msh
ANGLES="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
for A in $ANGLES; do
  ./make_steady.sh steady_naca0006_re500/steady_aoa$A 500 $A $MSH 6000 > /dev/null 2>&1 \
    || { echo "BUILD FAIL AoA=$A"; exit 1; }
done
echo "$(date '+%F %T') built"
STEADY_DIR=steady_naca0006_re500 ./run_steady.sh $ANGLES
echo "$(date '+%F %T') NACA0006 SWEEP DONE"
