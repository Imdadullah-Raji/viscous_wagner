#!/bin/bash
# queue_re700_ext.sh -- cfd-test, 2026-10-20.  Re 700 tight 4, 4.25 and 4.5 deg were not converged at 9000
# (4: overshot, still returning; 4.25: ratio 0.93, slow; 4.5: ~7e-8 from the fixed point).  Waits for the quarter-angle queue
# (run_re700q.log) to finish so the 4 ranks are never shared, then adds 3000 tight iterations to all three.
cd "$(dirname "$0")" || exit 1
until grep -q "QUARTERS COMPLETE" run_re700q.log 2>/dev/null; do sleep 300; done
./tighten_case.sh 3000 re700/steady_aoa4 re700/steady_aoa4.25 re700/steady_aoa4.5
ONLY="4 4.25 4.5" ./run_batch.sh /opt/openfoam14/etc/bashrc 4 re700
