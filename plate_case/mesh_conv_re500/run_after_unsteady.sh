#!/bin/bash
# Mesh-convergence steady runs (oct8_workflow.md), queued behind the plate 8.5 deg unsteady
# run so only one 6-rank job is ever active.  Fine meshes: h1 = 0.025/sqrt(Re), nr 260
# (BL cells within 5/sqrt(Re): 44 -> 72); tangential counts unchanged.
#   plate    : plate_case/mesh_conv_re500/steady_aoa{0,2,6,8}
#   NACA0006 : static_airfoil/steady_naca0006_re500_fine/steady_aoa{0,2,6,8}
SA=/home/raji/Research/Thesis/static_airfoil
# (The original wait loop here pgrep'd for the driver's command line; the launching
# shells contain that same text, so it never returned.  The 8.5 deg run was stopped by
# hand at t~61 (steady) and this script was then started directly.)
echo "$(date '+%F %T') starting mesh-convergence runs"
STEADY_DIR=../plate_case/mesh_conv_re500 $SA/run_steady.sh 0 2 6 8
STEADY_DIR=steady_naca0006_re500_fine     $SA/run_steady.sh 0 2 6 8
echo "$(date '+%F %T') MESH-CONV DONE"
