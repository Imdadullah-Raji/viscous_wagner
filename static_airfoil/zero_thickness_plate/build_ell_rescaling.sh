#!/bin/bash
# Zero-thickness plate Re sweep on the elliptic grid, one mesh per Re: h1 = 0.05/sqrt(Re) and
# n_nu = 800*sqrt(Re/500) (even), so the tip cells keep Re 500's E800 shape (skewness 3.78):
#   Re 200 300 400 600 700 800 900, 0-6 deg at 0.5 deg (13 angles)
#   -> zero_thickness_plate/ell_re_scaling/steady_re<RE>/steady_aoa<A>
# Re 500 is ell_re500/ (0.25-deg grid).  0 deg of each Re starts cold (10000 iterations); the
# others run 6000, warm-started from the nearest converged angle (run_batch.sh WARM=1).
cd "$(dirname "$0")/.." || exit 1
Z=zero_thickness_plate
for R in 200 300 400 600 700 800 900; do
  D=$Z/ell_re_scaling/steady_re$R; mkdir -p $D
  for i in $(seq 0 12); do
    A=$(python3 -c "print(f'{$i*0.5:g}')")
    echo "$R $A"
  done
done | xargs -P 3 -n 2 sh -c '
  R=$0; A=$1; Z=zero_thickness_plate; C=$Z/ell_re_scaling/steady_re$R/steady_aoa$A
  [ -f $C/constant/polyMesh/boundary ] && exit 0
  N=6000; [ "$A" = 0 ] && N=10000
  BAFFLES=$Z/createBafflesDict nice -n 10 ./make_steady.sh $C $R $A $(ls $Z/meshes/ell_scaled_nnu*_nmu210_re$R.msh) $N > /dev/null 2>&1 || echo "BUILD FAIL $R $A"
  grep -q "Mesh OK" $C/log.checkMesh || echo "CHECKMESH NOT OK $R $A"
  sed -i "s/^writePrecision .*/writePrecision  12;/" $C/system/controlDict'
echo "$(date '+%F %T') BUILD DONE"
