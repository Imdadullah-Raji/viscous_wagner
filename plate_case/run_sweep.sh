#!/bin/bash
# nohup bash run_sweep.sh --re 1000 [--ranks 6] > re1000/run.log 2>&1 < /dev/null &
set -e
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /home/raji/miniconda3/envs/research/bin/python "$(dirname "$0")/run_sweep.py" "$@"
