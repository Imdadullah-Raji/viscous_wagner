#!/bin/bash
# status_here.sh -- queue status in the folder it lives in (copied to every machine's ztp_tight/).
cd "$(dirname "$0")" || exit 1
for L in run_*.log; do
  printf "  %-22s ran=%-3s NO-End=%-2s | %s\n" $L $(grep -c " ran " $L) $(grep -c NO-End $L) "$(tail -1 $L | cut -c1-70)"
done
for c in $(ls -td re*/steady_aoa*/log.steady 2>/dev/null | head -3); do
  grep -qa "^End" $c || echo "  running: ${c%/log.steady}  $(grep -a '^Time =' $c | tail -1)"
done
uptime | sed 's/.*load/  load/'
