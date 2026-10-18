#!/bin/bash
# Build the zero-thickness plate Re = 500 sweep, 0-6 deg at 0.25 deg, on ONE mesh
# (tip study: T3, tip cell 5.6e-4c).  a3 is a curvature, so every angle shares the mesh.
cd "$(dirname "$0")/.." || exit 1
MSH=zero_thickness_plate/meshes/ztp_T3_tip0.000559017_np240_nsp140_re500.msh
python3 -c "print('\n'.join(f'{i*0.25:g}' for i in range(25)))" |
  xargs -P 3 -I{} sh -c 'C=zero_thickness_plate/steady_re500/steady_aoa{}; [ -f $C/constant/polyMesh/boundary ] && exit 0;
     BAFFLES=zero_thickness_plate/createBafflesDict nice -n 19 ./make_steady.sh $C 500 {} '"$MSH"' 6000 > /dev/null 2>&1 || echo "BUILD FAIL {}";
     grep -q "Mesh OK" $C/log.checkMesh || echo "CHECKMESH FAIL {}";
     sed -i "s/^writePrecision .*/writePrecision  12;/" $C/system/controlDict'
echo "$(date '+%F %T') BUILD DONE"
