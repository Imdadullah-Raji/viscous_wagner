#!/bin/bash
# Load rebalance (2026-10-09): banana-cfd is ~3x slower per case, so 8.4 and 8.5 moved here.
# Waits for the first local sweep's completion marker, then runs them the same way.
cd "$(dirname "$0")" || exit 1
until grep -q "STEADY SWEEP COMPLETE: 7.5" run_local.log; do sleep 20; done
cd ../.. && STEADY_DIR=lsb/plate_fine_re500 ./run_steady.sh 8.4 8.5
