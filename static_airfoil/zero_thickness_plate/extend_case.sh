#!/bin/bash
# extend_case.sh <endTime> <caseDir> [<caseDir> ...]
# Prepare finished steady cases for continuation (2026-10-19: Re 800/900 at 4-6 deg were still
# creeping at 6000 iterations).  Restart from the latest time, run to <endTime>, write once at the end.
# The first run's log is kept as log.steady_to<oldEnd> so run_batch.sh does not skip the case.
# The force history continues in postProcessing/forceCoeffs/<oldEnd>/ (steady_summary.py joins them).
# Decomposed data are reused when present (run on the same rank count).
END=$1; shift
for C in "$@"; do
  [ -f "$C/system/controlDict" ] || { echo "SKIP $C (no case)"; continue; }
  old=$(grep "^endTime" "$C/system/controlDict" | tr -dc 0-9)
  grep -qa "^End" "$C/log.steady" 2>/dev/null || { echo "SKIP $C (first run not finished)"; continue; }
  sed -i -e 's/^startFrom .*/startFrom       latestTime;/' -e "s/^endTime .*/endTime         $END;/" \
         -e "s/^writeInterval .*/writeInterval   $END;/" "$C/system/controlDict"
  mv "$C/log.steady" "$C/log.steady_to$old"
  echo "$C: $old -> $END ($(ls -d $C/processor0 2>/dev/null | wc -l) decomposed)"
done
