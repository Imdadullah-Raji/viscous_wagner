#!/bin/bash
# recover_queue.sh <sweepRoot> <folder> [<angle> ...]      (2026-10-18 incident recovery)
# Listed angles are reset even if their log says End (an orphaned solver can finish into a log
# the duplicate run had already truncated).
# swap_warm.sh killed the run_batch loop AND the subshell of the running case; the solver kept
# running orphaned and a new queue restarted the same case, so two solvers wrote one case.
# This stops every queue/solver whose working directory is under <sweepRoot>/<folder> (no other
# process is touched), and resets every unfinished case there to its unrun state (0/, constant/,
# system/ kept; processor*, postProcessing, logs and time directories > 0 removed).
ROOT=$(cd "$1" && pwd); DIR=$2; T="$ROOT/$DIR"; shift 2; FORCE=" $* "
for p in $(pgrep -f 'run_batch[.]sh|swap_warm[.]sh'); do
  [ "$p" = "$$" ] && continue
  case "$(readlink /proc/$p/cwd 2>/dev/null)" in "$ROOT"*) kill "$p" && echo "killed queue $p";; esac
done
for p in $(pgrep -x 'mpirun|foamRun|orterun|prterun|reconstructPar|decomposePar'); do
  case "$(readlink /proc/$p/cwd 2>/dev/null)" in "$T"/*) kill "$p" && echo "killed solver $p ($(readlink /proc/$p/cwd))";; esac
done
sleep 5
for p in $(pgrep -x 'mpirun|foamRun|orterun|prterun'); do
  case "$(readlink /proc/$p/cwd 2>/dev/null)" in "$T"/*) kill -9 "$p" && echo "kill -9 $p";; esac
done
for C in "$T"/steady_aoa*; do
  [ -d "$C" ] || continue
  case "$FORCE" in *" ${C##*steady_aoa} "*) ;; *) grep -qa "^End" "$C/log.steady" 2>/dev/null && continue;; esac
  [ -e "$C/log.steady" ] || [ -d "$C/processor0" ] || continue
  rm -rf "$C"/processor* "$C"/postProcessing "$C"/log.steady "$C"/log.decomposePar "$C"/log.reconstructPar "$C"/log.mapFields
  for t in "$C"/[0-9]*; do [ "$(basename "$t")" != 0 ] && rm -rf "$t"; done
  echo "reset $C"
done
