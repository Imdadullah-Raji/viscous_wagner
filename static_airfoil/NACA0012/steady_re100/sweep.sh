#!/bin/bash
# Re=100 steady sweep: build every case, then run sequentially at 6 ranks.
cd "$(dirname "$0")/.."
MSH=meshes/M2_ns160_nr180_b025_re100.msh
ANGLES="0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 14 16 18 20"
for A in $ANGLES; do ./make_steady.sh steady_re100/steady_aoa$A 100 $A $MSH 6000 || exit 1; done
STEADY_DIR=steady_re100 ./run_steady.sh $ANGLES
