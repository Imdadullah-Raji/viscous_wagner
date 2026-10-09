#!/bin/bash
# Launch from a normal local terminal: bash mesh_convergence_M1/run_study.sh
set -e
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/matplotlib-mesh-convergence-M1
exec /home/raji/miniconda3/envs/research/bin/python "$(dirname "$0")/study.py" "$@"
