#!/bin/bash
# swap_warm.sh <bashrc> <ranks> <folder> <log> [ONLY angles...]
# Stop the running run_batch LOOP (not the subshell running the current case), wait for that case
# to finish, then continue the queue with WARM=1 (seed from the nearest converged angle).
cd "$(dirname "$0")" || exit 1
RC=$1; NP=$2; DIR=$3; LOG=$4; shift 4; ONLYA="$*"
for p in $(pgrep -f "run_batch.sh /opt"); do
  pp=$(ps -o ppid= -p $p | tr -d ' ')
  ps -o args= -p $pp | grep -q "run_batch.sh /opt" || { kill $p; echo "$(date '+%F %T') stopped loop $p" >> $LOG; }
done
while pgrep -f "run_batch.sh /opt" > /dev/null; do sleep 20; done
echo "$(date '+%F %T') switching to WARM=1" >> $LOG
WARM=1 ONLY="$ONLYA" exec ./run_batch.sh $RC $NP $DIR >> $LOG 2>&1
