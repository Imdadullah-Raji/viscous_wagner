#!/bin/bash
# Tip-resolution study for the zero-thickness plate, Re = 500, at 6 deg (the
# angle where the edge's viscous scale is smallest, ~c/(Re*alpha)^2 = 3.7e-4c).
# T1..T4 halve the tip cell (2.24e-3 -> 2.8e-4c), raising nplate/nspokes to keep
# the growth ratios ~1.03 / ~1.06.  One mesh is then chosen for the whole 0-6 deg
# sweep: a3 is a curvature, so the angles must share one mesh.
cd "$(dirname "$0")/.." || exit 1
Z=zero_thickness_plate
while pgrep -x foamRun > /dev/null; do sleep 30; done      # wait for the running sweep
for T in T1 T2 T3 T4; do
  MSH=$(ls $Z/meshes/ztp_${T}_*.msh)
  D=$Z/tip_study_re500/$T; mkdir -p $D
  [ -d $D/steady_aoa6 ] || BAFFLES=$Z/createBafflesDict ./make_steady.sh $D/steady_aoa6 500 6 $MSH 6000 || continue
  STEADY_DIR=$D ./run_steady.sh 6
done
echo "$(date '+%F %T') TIP STUDY DONE"
