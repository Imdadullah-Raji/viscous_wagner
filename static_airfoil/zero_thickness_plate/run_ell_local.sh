#!/bin/bash
cd "$(dirname "$0")" || exit 1
ONLY="0.5 1.25 2 2.75 3.5 4.25 5 5.75" exec ./run_batch.sh /opt/openfoam-dev/etc/bashrc 6 ell_re500
