#!/bin/bash
# Plate Re=500 at 15 deg, forces only -- launched ONLY if the 9.5 deg run is judged steady
# (user, 2026-10-08: "crank it up to 15 if 9.5 is not a shedding case").  The trigger is
# the "MANUAL STOP ... steady" line appended to run_aoa9.5.log when that run is stopped;
# if 9.5 sheds, its driver ends with a NO_FIELDS line instead and this script exits.
cd /home/raji/Research/Thesis/static_airfoil
L=../plate_case/unsteady_re500/run_aoa9.5.log
until grep -q "MANUAL STOP.*steady\|NO_FIELDS: stopping\|ABORT" $L; do sleep 30; done
grep -q "MANUAL STOP.*steady\|ABORT: no shedding" $L || { echo "$(date '+%F %T') 9.5 deg was not steady; 15 deg NOT run"; exit 0; }
NO_FIELDS=1 ./run_staged.sh ../plate_case/unsteady_re500/aoa15 6
