#!/bin/bash
# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
# 2026-10-20 rebalance: Re 500 4-5 deg (never reached machineB intact) -> machineB 4 ranks;
# Re 700 quarter angles (new, seeded) -> cfd-test 4 ranks, two-phase via run_quarters.sh.
cd "$(dirname "$0")" || exit 1
U=$CFD_USER
R5="4 4.25 4.5 4.75 5"
Q="0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25"
( for a in $R5; do
    rsync -a --delete --exclude 'processor*' re500/steady_aoa$a/ $HOST_B:$REMOTE_B_DIR/ztp_tight/re500/steady_aoa$a/ < /dev/null || echo "PUSH FAIL machineB $a"
  done
  echo "$(date '+%F %T') pushed Re 500 to machineB"
  timeout 40 ssh -o BatchMode=yes $HOST_B "cd $REMOTE_B_DIR/ztp_tight && ONLY='$R5' setsid nohup ./run_batch.sh /opt/openfoam14/etc/bashrc 4 re500 > run_re500_b.log 2>&1 < /dev/null &" < /dev/null
  echo "$(date '+%F %T') launched machineB Re 500 4-5" ) &
( rsync -a --rsync-path="sudo -u raji rsync" -e "ssh -o BatchMode=yes -l $U" run_quarters.sh cfd-test:/home/raji/cases/ztp_tight/ < /dev/null
  for a in $Q; do
    rsync -a --exclude 'processor*' --rsync-path="sudo -u raji rsync" -e "ssh -o BatchMode=yes -l $U" \
      re700/steady_aoa$a cfd-test:/home/raji/cases/ztp_tight/re700/ < /dev/null || echo "PUSH FAIL cfd $a"
  done
  echo "$(date '+%F %T') pushed Re 700 quarters to cfd-test"
  timeout 60 ssh -o BatchMode=yes -l $U cfd-test "sudo -u raji bash -c 'cd /home/raji/cases/ztp_tight && setsid nohup ./run_quarters.sh /opt/openfoam14/etc/bashrc 4 re700 ./tighten_case.sh $Q > run_re700q.log 2>&1 < /dev/null &'" < /dev/null
  echo "$(date '+%F %T') launched cfd-test Re 700 quarters" ) &
wait
echo PUSH-AND-LAUNCH DONE
