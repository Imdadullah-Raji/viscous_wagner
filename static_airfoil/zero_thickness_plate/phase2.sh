#!/bin/bash
# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
# phase2.sh <log> <bashrc> <ranks> <folder> <angles...>
# Rebalance 2026-10-18: after this machine's current queue logs its BATCH COMPLETE line, run the
# extra angles (moved off machineB) with WARM=1.  Starts nothing while the current queue runs.
cd "$(dirname "$0")" || exit 1
LOG=$1; RC=$2; NP=$3; DIR=$4; shift 4
n0=$(grep -c "BATCH COMPLETE" "$LOG" 2>/dev/null)
until [ "$(grep -c "BATCH COMPLETE" "$LOG" 2>/dev/null)" -gt "$n0" ]; do sleep 30; done
echo "$(date '+%F %T') phase 2: $*" >> "$LOG"
WARM=1 ONLY="$*" exec ./run_batch.sh "$RC" "$NP" "$DIR" >> "$LOG" 2>&1
