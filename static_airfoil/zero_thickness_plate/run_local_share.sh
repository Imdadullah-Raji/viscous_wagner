#!/bin/bash
# This machine's share of the Re 500 sweep, after the tip study frees the cores.
cd "$(dirname "$0")" || exit 1
until grep -q "TIP STUDY DONE" tip_study.log; do sleep 30; done
ONLY="1.25 2.5 3.5 4.5" exec ./run_batch.sh /opt/openfoam-dev/etc/bashrc 6 steady_re500
