#!/bin/bash
# Plate Re=500 shedding-onset test at 9.5 deg (oct8_workflow.md: 8.5 was steady -> test 9.5),
# forces only.  Waits for the mesh-convergence queue so one 6-rank job runs at a time.
# Waits on the MESH-CONV DONE marker, not on pgrep of a command line (that self-matched).
cd /home/raji/Research/Thesis/static_airfoil
until grep -q "MESH-CONV DONE" ../plate_case/mesh_conv_re500/run_steady.log; do sleep 30; done
NO_FIELDS=1 ./run_staged.sh ../plate_case/unsteady_re500/aoa9.5 6
