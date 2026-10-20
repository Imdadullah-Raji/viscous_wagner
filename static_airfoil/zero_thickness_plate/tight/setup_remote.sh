#!/bin/bash
# setup_remote.sh <sweepRoot> <Re> <angles...>
# Run on a remote machine from its ztp_tight/ folder.  Copies the already-converged cases from the
# machine's own sweep folder into ztp_tight/re<Re>/ (no processor dirs; the latest time is
# reconstructed) and prepares them for 3000 tight iterations with tighten_case.sh.
# The three tight_test cases (re500_aoa3, re900_aoa2, re900_aoa3.5) are copied as they are: they
# already ran tight to 12000 and keep their End, so run_batch.sh skips them.
cd "$(dirname "$0")" || exit 1
ROOT=$1; RE=$2; shift 2
mkdir -p re$RE
for a in "$@"; do
  if [ -d re${RE}_aoa$a ]; then src=re${RE}_aoa$a; tight=0
  elif [ "$RE" = 500 ]; then src=$ROOT/ell_re500/steady_aoa$a; tight=1
  else src=$ROOT/ell_re_scaling/steady_re$RE/steady_aoa$a; tight=1; fi
  rsync -a --delete --exclude 'processor*' $src/ re$RE/steady_aoa$a/
  [ $tight = 1 ] && ./tighten_case.sh 3000 re$RE/steady_aoa$a
done
echo "SETUP DONE re$RE"
