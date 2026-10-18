#!/bin/bash
# T4 (finest tip, 154k cells) was still oscillating at 6000 SIMPLE iterations (C_L +-2e-4,
# CREEPING); continue it to 12000 from its 6000 state, after this machine's sweep share.
cd "$(dirname "$0")" || exit 1
until grep -q "BATCH COMPLETE" steady_re500/run_batch_local.log 2>/dev/null; do sleep 30; done
source /opt/openfoam-dev/etc/bashrc > /dev/null 2>&1
cd tip_study_re500/T4/steady_aoa6
sed -i -e 's/^startFrom .*/startFrom       latestTime;/' -e 's/^endTime .*/endTime         12000;/' \
       -e 's/^writeInterval .*/writeInterval   12000;/' system/controlDict
for p in processor*; do cp system/controlDict $p/../system/controlDict 2>/dev/null; done
mpirun -np 6 foamRun -parallel > log.steady_to_12000 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
echo "$(date '+%F %T') T4 EXTENDED: $(grep -c '^End' log.steady_to_12000) End"
