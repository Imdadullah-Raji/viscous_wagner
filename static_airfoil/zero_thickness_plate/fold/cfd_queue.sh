#!/bin/bash
# cfd_queue.sh -- cfd-test, 2026-10-20: one 4-rank queue so the two jobs never share cores.
#   1. Re 700 tight extensions at 4, 4.25, 4.5 deg (ztp_tight/queue_re700_ext.sh);
#   2. the descending fold chain at Re 500: cold start at 12 deg, then 11 ... 9, each seeded from above.
cd "$(dirname "$0")/.." || exit 1
ztp_tight/queue_re700_ext.sh
ztp_fold/fold_chain.sh /opt/openfoam14/etc/bashrc 4 re500_desc 12 11 10.5 10 9.5 9
