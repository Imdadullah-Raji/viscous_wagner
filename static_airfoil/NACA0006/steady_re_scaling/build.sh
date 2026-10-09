#!/bin/bash
# NACA0006 Re sweep (instructions/Quadratic_hypothesis.md): builds steady SIMPLE cases on
# 0-6 deg in 0.25-deg steps.
#   Re 200..900 -> NACA0006/steady_re_scaling/steady_re<RE>/steady_aoa<A>
#   Re 500      -> NACA0006/steady_re500/steady_aoa<A>, non-integer angles only (the integer
#                  angles 0-20 already exist there)
# Mesh: the M2 NACA0006 recipe, Re-specific (h1 = 0.05/sqrt(Re)):
#   mesher.py --naca 0006 --re RE -r 30 --nr 180 --ns 160 --nte 32 --bump 0.25
# 6000 iterations, writePrecision 12 (as the plate sweep).  Running is separate: the cases are
# shared between machines with plate_case-style run_batch.sh.
cd "$(dirname "$0")/../.." || exit 1
ANGLES=$(python3 -c "print(' '.join(f'{i*0.25:g}' for i in range(25)))")
build() {  # build <caseDir> <Re> <A>
  [ -f "$1/constant/polyMesh/boundary" ] && return
  ./make_steady.sh "$1" "$2" "$3" meshes/M2_ns160_nr180_b025_naca0006_re$2.msh 6000 > /dev/null 2>&1 \
    || { echo "BUILD FAIL $1"; return; }
  grep -q "Mesh OK" "$1/log.checkMesh" || echo "CHECKMESH FAIL $1"
  sed -i 's/^writePrecision .*/writePrecision  12;/' "$1/system/controlDict"
}
for RE in 200 300 400 600 700 800 900; do
  D=NACA0006/steady_re_scaling/steady_re$RE; mkdir -p $D
  for A in $ANGLES; do build $D/steady_aoa$A $RE $A; done
  echo "$(date '+%F %T') built Re=$RE"
done
for A in $ANGLES; do
  case $A in *.*) build NACA0006/steady_re500/steady_aoa$A 500 $A ;; esac
done
echo "$(date '+%F %T') built Re=500 infill; BUILD DONE"
