#!/bin/bash
# Re = 500 sweep, 0-6 deg at 0.25 deg, on the elliptic grid E800 (one mesh for every angle).
# 6 deg is the tip-check case (12000 iterations); the rest run 10000 (C_L moved 1e-6 from 10000 to
# 12000 at 6 deg) and are extended if not converged.
cd "$(dirname "$0")/.." || exit 1
python3 -c "print('\n'.join(f'{i*0.25:g}' for i in range(24)))" |
  xargs -P 3 -I{} sh -c 'C=zero_thickness_plate/ell_re500/steady_aoa{}; [ -f $C/constant/polyMesh/boundary ] && exit 0;
     BAFFLES=zero_thickness_plate/createBafflesDict nice -n 19 ./make_steady.sh $C 500 {} zero_thickness_plate/meshes/ell_nnu800_nmu210_re500.msh 10000 > /dev/null 2>&1 || echo "BUILD FAIL {}";
     grep -q "faces" $C/log.checkMesh || echo "CHECKMESH MISSING {}";
     sed -i "s/^writePrecision .*/writePrecision  12;/" $C/system/controlDict'
echo "$(date '+%F %T') BUILD DONE"
