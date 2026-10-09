#!/bin/bash
# Rebalance (2026-10-09 13:45): banana-cfd stalled ~12 min during 7.97, so 7.99 moved here.
cd "$(dirname "$0")" || exit 1
until grep -q "STEADY SWEEP COMPLETE: 7.91" run_local.log; do sleep 15; done
cd ../.. && STEADY_DIR=lsb/plate_fine_re500 ./run_steady.sh 7.99
