#!/bin/bash
# Repo copy: host names, logins and paths replaced by $HOST_B, $REMOTE_B_DIR, $CFD_USER (set them to run).
# make_case.sh <caseDir> <Re> <AoA_deg> <mesh.msh>
#
# Builds a complete static-airfoil case.  AoA is imposed by ROTATING THE
# FREESTREAM, not the mesh, so one mesh serves every angle.  The angle appears
# in exactly three places and all three are written from the same numbers here:
#   0/U          internalField + freestreamValue
#   functions    dragDir / liftDir
# Getting those out of sync is the classic way to produce a plausible-looking
# but wrong C_L, so they are never edited by hand.
set -e
CASE=$1; RE=$2; AOA=$3; MSH=$4
[ -n "$MSH" ] || { echo "usage: make_case.sh <caseDir> <Re> <AoA_deg> <mesh.msh>"; exit 1; }

read CA SA NU <<< $(python3 -c "
import math
a=math.radians($AOA); print(f'{math.cos(a):.9f} {math.sin(a):.9f} {1.0/$RE:.9e}')")

echo "$CASE: Re=$RE  AoA=$AOA  nu=$NU  dir=($CA $SA 0)"
rm -rf $CASE; mkdir -p $CASE/{0,constant,system}

# ---------------------------------------------------------------- 0/U, 0/p ---
cat > $CASE/0/U <<EOF
FoamFile { version 2.0; format ascii; class volVectorField; object U; }

dimensions  [0 1 -1 0 0 0 0];

internalField  uniform ($CA $SA 0);

boundaryField
{
    farfield      { type freestreamVelocity; freestreamValue uniform ($CA $SA 0); }
    airfoil       { type noSlip; }
    frontAndBack  { type empty; }
    #includeEtc "caseDicts/setConstraintTypes"
}
EOF

cat > $CASE/0/p <<'EOF'
FoamFile { version 2.0; format ascii; class volScalarField; object p; }

dimensions  [0 2 -2 0 0 0 0];

internalField  uniform 0;

boundaryField
{
    farfield      { type freestreamPressure; freestreamValue uniform 0; }
    airfoil       { type zeroGradient; }
    frontAndBack  { type empty; }
    #includeEtc "caseDicts/setConstraintTypes"
}
EOF

# ------------------------------------------------------------- constant/ ---
cat > $CASE/constant/physicalProperties <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object physicalProperties; }

viscosityModel  constant;

// Re = $RE at U = 1, c = 1
nu  $NU [m^2/s];
EOF

cat > $CASE/constant/momentumTransport <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object momentumTransport; }

simulationType  laminar;

laminar { model Stokes; }
EOF

# --------------------------------------------------------------- system/ ---
# controlDict holds STAGE-1 settings: writeInterval == endTime means exactly one
# write, at the end, which is the restart point.  run_staged.sh rewrites these
# four keys per stage.  adjustTimeStep stays 'no' everywhere -- the variable-dt
# controller is what corrupted the first attempt.
cat > $CASE/system/controlDict <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }

solver          incompressibleFluid;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.5;
deltaT          4e-05;

writeControl    runTime;
writeInterval   0.5;

purgeWrite      0;
writeFormat     binary;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

adjustTimeStep  no;
maxCo           0.8;
maxDeltaT       0.01;
EOF

cat > $CASE/system/fvSchemes <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes      { default backward; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; div(phi,U) Gauss linearUpwind grad(U); }
laplacianSchemes{ default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }
EOF

cat > $CASE/system/fvSolution <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
    p     { solver GAMG; tolerance 1e-07; relTol 0.01; smoother GaussSeidel; }
    pFinal{ $p; tolerance 1e-08; relTol 0; }
    U     { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-08; relTol 0.01; }
    UFinal{ $U; tolerance 1e-09; relTol 0; }
}
PIMPLE
{
    // PISO (nOuterCorrectors 1) on purpose.  Measured on the 110k-cell mesh:
    // outer2 = 0.416 s/step (159 GAMG it/step), outer1 = 0.176 s/step (53).
    // At Co < 1 with backward ddt the outer loop buys no accuracy here, only
    // 2.4x the cost.
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors      2;

    // Keep this at 1.  The O-grid has max non-orthogonality ~64 deg at the wake
    // block corners; setting it to 0 made the pressure solve converge WORSE and
    // the step 3.7x slower (0.650 s/step), not faster.
    nNonOrthogonalCorrectors 1;

    pRefCell 0;
    pRefValue 0;
}
EOF

cat > $CASE/system/decomposeParDict <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }

// 6, not 12: this box is an i5-12400 with 6 PHYSICAL cores / 12 threads.
// nproc reports hyperthreads and mpirun -np 8 fails with "not enough slots";
// two concurrent 4-rank jobs oversubscribe and cost 1.7x per step.
numberOfSubdomains  6;

method  scotch;
EOF

cat > $CASE/system/functions <<EOF
/*--------------------------------*- C++ -*----------------------------------*\\
    Post-processing for the 2D NACA0012 at Re = $RE, AoA = $AOA deg.

    Reference values follow from c = 1, U = 1, span = 0.1:
        magUInf 1;  lRef 1;  Aref = c*span = 0.1;  rhoInf 1 (kinematic p).

    STAGING CONTRACT -- read before changing a writeControl:
      forceCoeffs is on 'timeStep'.  It appends to a single .dat and NEVER
      creates a time directory, so per-step C_L/C_D through the transient is
      free.  The other three are on 'writeTime', so they fire only when the
      case writes fields.  That is what lets stage 1 run field-free purely by
      setting controlDict's writeInterval == endTime -- this file is IDENTICAL
      in both stages and must stay that way.
\\*---------------------------------------------------------------------------*/

FoamFile { version 2.0; format ascii; class dictionary; object functions; }

// Order matters: Cp and wallShearStress must exist before "surfaces" samples them.

forceCoeffs
{
    type            forceCoeffs;
    libs            ("libforces.so");

    writeControl    timeStep;
    writeInterval   1;
    log             yes;

    patches         (airfoil);

    // Kinematic pressure, so rho is the placeholder name and rhoInf scales it
    rho             rhoInf;
    rhoInf          1;

    magUInf         1;
    lRef            1;      // c
    Aref            0.1;    // c * span (span = 0.1) -- NOT 1.0

    CofR            (0 0 0);
    dragDir         ($CA $SA 0);
    liftDir         (-$SA $CA 0);
    pitchAxis       (0 0 1);
}

Cp
{
    type            pressure;
    libs            ("libfieldFunctionObjects.so");

    // writeTime, NOT timeStep. A field-writing FO creates a time directory for
    // whatever time it fires on, so "timeStep" during a live foamRun mints a new
    // dir every deltaT containing only this field and no U/p.
    executeControl  writeTime;
    writeControl    writeTime;

    field           p;
    result          Cp;

    calcTotal       no;
    calcCoeff       yes;

    // All three are required lookups once calcCoeff is on.
    // p here is kinematic, so it is scaled by rhoInf: Cp = (p - pInf)/(0.5*rhoInf*|UInf|^2)
    pInf            0;
    UInf            ($CA $SA 0);
    rhoInf          1;
}

wallShearStress
{
    type            wallShearStress;
    libs            ("libfieldFunctionObjects.so");

    executeControl  writeTime;
    writeControl    writeTime;

    // Keep this explicit. An explicit list bypasses the isA<wallPolyPatch> filter,
    // so this still works if the mesh is regenerated (gmshToFoam resets every
    // patch to type "patch", losing the "wall" type on airfoil).
    patches         (airfoil);
}

surfaces
{
    type            surfaces;
    libs            ("libsampling.so");

    executeControl  writeTime;
    writeControl    writeTime;

    surfaceFormat   raw;
    interpolationScheme cell;

    fields          (Cp p wallShearStress U);

    surfaces
    {
        airfoil
        {
            type        patch;
            patches     (airfoil);

            // Must be false. wallShearStress has a zero internal field and carries
            // data only on the boundary, so interpolating would smear zeros in.
            interpolate false;
        }
    }
}
EOF

# ------------------------------------------------------------------ mesh ---
# openfoam-dev was removed 2026-10-19; v14 (build 14-a4867fc914c1, same as machineB) replaces it.
source ${FOAM_BASHRC:-/opt/openfoam14/etc/bashrc} >/dev/null 2>&1
gmshToFoam $MSH -case $CASE > $CASE/log.gmshToFoam 2>&1
grep -qi "^End" $CASE/log.gmshToFoam || { echo "gmshToFoam FAILED"; tail -20 $CASE/log.gmshToFoam; exit 1; }

# Opt-in, for meshes whose wall is an INTERNAL surface (zero-thickness plate,
# zero_thickness_plate/ztp_mesher.py): BAFFLES=<createBafflesDict> splits the
# faceZone gmshToFoam wrote into two wall patches in group "airfoil", which every
# "airfoil" entry below (BCs, forces, Cp, wallShearStress, surfaces) matches.
if [ -n "${BAFFLES:-}" ]; then
  cp "$BAFFLES" $CASE/system/createBafflesDict
  createBaffles -case $CASE -overwrite > $CASE/log.createBaffles 2>&1
  grep -qi "^End" $CASE/log.createBaffles || { echo "createBaffles FAILED"; tail -20 $CASE/log.createBaffles; exit 1; }
fi

# gmshToFoam writes every patch as plain "patch".  frontAndBack MUST become
# empty or the run is 3D with one cell in z, and airfoil should be wall so
# wall-distance and wallShearStress behave.
python3 - "$CASE/constant/polyMesh/boundary" <<'PY'
import re, sys
p = sys.argv[1]; s = open(p).read()
def settype(name, t):
    global s
    s2 = re.sub(r'(\b%s\s*\{[^}]*?\btype\s+)\w+;' % name, r'\g<1>%s;' % t, s, count=1, flags=re.S)
    assert s2 != s or ('type            %s;' % t) in s, "patch %s not retyped" % name
    s = s2
settype('frontAndBack', 'empty')
import os
baffled = bool(os.environ.get('BAFFLES'))
if not baffled:      # baffle patches are created as walls by createBaffles
    settype('airfoil', 'wall')
settype('farfield', 'patch')
open(p,'w').write(s)
print("  boundary: frontAndBack->empty, " + ("baffles (group airfoil) are walls" if baffled else "airfoil->wall"))
PY

checkMesh -case $CASE > $CASE/log.checkMesh 2>&1 || true
echo "  cells: $(grep -m1 'cells:' $CASE/log.checkMesh | awk '{print $2}')   $(grep -m1 'Mesh OK\|FAILED' $CASE/log.checkMesh)"
