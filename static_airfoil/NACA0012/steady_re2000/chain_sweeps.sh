#!/bin/bash
# Re=2000 then Re=100 steady sweeps, strictly sequential (6 MPI ranks total).
cd "$(dirname "$0")/.."
steady_re2000/sweep.sh > steady_re2000/run_steady.log 2>&1
steady_re100/sweep.sh  > steady_re100/run_steady.log  2>&1
echo "$(date '+%F %T') CHAIN DONE" >> steady_re100/run_steady.log
