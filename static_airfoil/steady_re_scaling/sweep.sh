#!/bin/bash
# Re-scaling check of the steady cubic term: SIMPLE sweeps at
# Re = 200 300 400 600 700 800 900 on 0-8 deg (0.5 deg grid, same as the
# Re 100/500/1000/2000 fits), M2 recipe with h1 = 0.05/sqrt(Re).
# Cases: steady_re_scaling/steady_re<RE>/steady_aoa<A>.
# Strictly sequential, 6 MPI ranks per case.  Whole degrees for every Re run
# first, so a coarse-grid answer exists after ~half the wall time; the
# half-degree angles fill in afterwards.
cd "$(dirname "$0")/.."
RES="200 300 400 600 700 800 900"
PASS1="0 1 2 3 4 5 6 7 8"
PASS2="0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5"
for RE in $RES; do
  MSH=meshes/M2_ns160_nr180_b025_re$RE.msh
  mkdir -p steady_re_scaling/steady_re$RE      # make_steady.sh writes <case>.build.log here
  for A in $PASS1 $PASS2; do
    ./make_steady.sh steady_re_scaling/steady_re$RE/steady_aoa$A $RE $A $MSH 6000 \
      > /dev/null 2>&1 || { echo "BUILD FAIL Re=$RE AoA=$A"; exit 1; }
  done
  echo "$(date '+%F %T') built Re=$RE"
done
for PASS in "$PASS1" "$PASS2"; do
  for RE in $RES; do
    echo "== Re=$RE  angles: $PASS"
    STEADY_DIR=steady_re_scaling/steady_re$RE ./run_steady.sh $PASS
  done
done
echo "$(date '+%F %T') RE-SCALING SWEEP DONE"
