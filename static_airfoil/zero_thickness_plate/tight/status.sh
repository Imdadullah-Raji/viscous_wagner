#!/bin/bash
# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
# status.sh -- one-shot status of the a5/a7 tight-convergence runs on all three machines (read-only).
# Runs status_here.sh in each machine's tight folder: per queue log, cases "ran" / "NO-End"
# (failed), the last line, and the case currently running.  Push status_here.sh first if missing.
cd "$(dirname "$0")" || exit 1
U=$CFD_USER
date '+%F %T %Z'
echo "== local";    ./status_here.sh
echo "== machineB";    timeout 40 ssh -o BatchMode=yes $HOST_B $REMOTE_B_DIR/ztp_tight/status_here.sh < /dev/null
echo "== cfd-test"; timeout 90 ssh -o BatchMode=yes -l $U cfd-test sudo /home/raji/cases/ztp_tight/status_here.sh < /dev/null 2>&1 \
  | grep -v -i "numpy\|tcp_upload\|^WARNING\|^$"
