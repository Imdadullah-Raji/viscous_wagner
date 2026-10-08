#!/bin/bash
# Onset bisection, plate Re=500: 8.5 and 9.5 steady, 15 sheds -> test 12 deg, forces only.
# Starts when the 15 deg driver has finished (its NO_FIELDS/ABORT line or a MANUAL STOP note).
cd /home/raji/Research/Thesis/static_airfoil
until grep -q "NO_FIELDS: stopping\|ABORT\|MANUAL STOP" ../plate_case/unsteady_re500/run_aoa15.log; do sleep 30; done
while pgrep -x foamRun >/dev/null; do sleep 10; done      # let the last solver exit
NO_FIELDS=1 ./run_staged.sh ../plate_case/unsteady_re500/aoa12 6
