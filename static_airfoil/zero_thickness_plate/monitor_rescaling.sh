#!/bin/bash
# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
# Monitor for the zero-thickness plate Re sweep (elliptic grid, 7 Re x 13 angles = 91 cases).
# One line per finished / failed / stalled case, keyed by Re+angle (cases copied between machines are
# reported once); a line when an Re completes; machine-unreachable events.  Exits when all 91 are done.
# Reads only: it never starts anything on any machine.
#   INIT=1 ./monitor_rescaling.sh   record what is already finished, silently, then exit
#   ./monitor_rescaling.sh          watch
cd "$(dirname "$0")" || exit 1
S=.monitor_state; mkdir -p $S
Z=$PWD
F='grep -v -i -e numpy -e ^WARNING -e tcp_upload -e increase.the.perf -e cloud.google'
probe='for C in ell_re_scaling/steady_re*/steady_aoa*; do [ -d $C ] || continue; a=${C##*steady_aoa}; r=${C%/steady_aoa*}; r=${r##*steady_re}; L=$C/log.steady; [ -f $L ] || continue;
  e=$(grep -c "^End" $L); f=$(grep -c -E "FOAM FATAL|Floating point exception \(core|sigFpe::sigHandler|\bnan\b" $L);
  it=$(grep "^Time =" $L | tail -1 | tr -dc 0-9); last=$(tail -1 $C/postProcessing/forceCoeffs/0/*.dat 2>/dev/null | awk "{printf \"Cl=%.6f Cd=%.6f\", \$4, \$3}");
  echo "$r $a $e $f $it $last"; done'
declare -A fails
while true; do
  for M in local machineB cfdtest; do
    case $M in
      local)   out=$(cd $Z && eval "$probe") ;;
      machineB)   out=$(timeout 60 ssh -o BatchMode=yes -o ConnectTimeout=20 $HOST_B "cd $REMOTE_B_DIR/ztp_sweep && $probe" < /dev/null 2>/dev/null) ;;
      cfdtest) out=$(timeout 120 ssh -o BatchMode=yes -l $CFD_USER cfd-test "sudo bash -c '$(echo "cd /home/raji/cases/ztp_sweep && $probe" | sed "s/'/'\\\\''/g")'" < /dev/null 2>/dev/null | $F) ;;
    esac
    if [ -z "$out" ]; then
      fails[$M]=$(( ${fails[$M]:-0} + 1 ))
      [ -z "$INIT" ] && [ ${fails[$M]} -eq 3 ] && [ ! -f $S/$M.unreach ] && { touch $S/$M.unreach; echo "$(date +%H:%M) UNREACHABLE $M (3 polls)"; }
      continue
    fi
    fails[$M]=0; rm -f $S/$M.unreach
    while read -r r a e f it rest; do
      [ -n "$a" ] || continue
      key=$S/$r.$a
      if [ "$e" -ge 1 ] && [ ! -f $key.done ]; then
        touch $key.done; [ -n "$INIT" ] && continue
        echo "$(date +%H:%M) DONE $M Re=$r aoa=$a iter=$it $rest"
        n=$(ls $S/$r.*.done 2>/dev/null | wc -l); [ "$n" -eq 13 ] && echo "$(date +%H:%M) Re=$r COMPLETE (13/13)"
      fi
      [ -n "$INIT" ] && continue
      if [ "$f" -ge 1 ] && [ ! -f $key.fail ]; then touch $key.fail; echo "$(date +%H:%M) FAILED $M Re=$r aoa=$a (fatal/nan) iter=$it"; fi
      if [ "$e" -eq 0 ] && [ -n "$it" ]; then
        prev=$(cat $key.$M.it 2>/dev/null); echo "$it" > $key.$M.it
        if [ "$prev" = "$it" ] && [ ! -f $key.stall ]; then touch $key.stall; echo "$(date +%H:%M) STALLED $M Re=$r aoa=$a at iter $it"; fi
      fi
    done <<< "$out"
  done
  [ -n "$INIT" ] && { echo "initialised: $(ls $S/*.done 2>/dev/null | wc -l) finished cases recorded"; exit 0; }
  n=$(ls $S/*.done 2>/dev/null | wc -l)
  if [ "$n" -ge 91 ]; then echo "$(date +%H:%M) ALL 91 Re-SWEEP CASES DONE"; exit 0; fi
  sleep 240
done
