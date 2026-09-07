#!/bin/bash
# make_steady.sh <caseDir> <Re> <AoA_deg> <mesh.msh> [nIter]
#
# Steady (SIMPLE) counterpart of make_case.sh: the same case in every respect
# except that the time derivative is deleted.
#
# WHY THE SCHEMES ARE NOT CHANGED.  The upstream steady cases used
# `div(phi,U) Gauss upwind`, which is more robust for SIMPLE.  We deliberately
# keep the TRANSIENT case's `Gauss linearUpwind grad(U)` instead.  The whole
# point of this solve is the difference
#
#     C_L,mean - C_L,steady  =  -c2 R^2      (the O(R^2) mean-flow distortion)
#
# and first-order upwind here against second-order upwind there would inject a
# DISCRETISATION difference into exactly that quantity, indistinguishable from
# the physical one being measured.  Same schemes, only ddt removed.
#
# nNonOrthogonalCorrectors is raised 1 -> 2.  That changes only the convergence
# RATE, not the fixed point: at convergence the non-orthogonal correction is
# itself converged however many correctors were used per iteration.
#
# A CONVERGED STEADY SOLVE IS NOT EVIDENCE THE FLOW IS STEADY -- `steadyState`
# deletes the ddt term, so SIMPLE finds a fixed point whether or not it is
# stable.  That is the POINT here: above the onset angle we WANT the unstable
# base flow, because it is the state the Stuart-Landau amplitude is measured
# relative to.  See BIFURCATION_THEORY.md section 1.
set -e
CASE=$1; RE=$2; AOA=$3; MSH=$4; NITER=${5:-6000}
[ -n "$MSH" ] || { echo "usage: make_steady.sh <caseDir> <Re> <AoA> <mesh.msh> [nIter]"; exit 1; }
ROOT=$(cd "$(dirname "$0")" && pwd)

"$ROOT/make_case.sh" "$CASE" "$RE" "$AOA" "$MSH" > "$CASE.build.log" 2>&1 || {
  echo "make_case FAILED for $CASE"; tail -20 "$CASE.build.log"; exit 1; }

cd "$CASE"

sed -i '/ddtSchemes/s/backward/steadyState/' system/fvSchemes
grep -q "steadyState" system/fvSchemes || { echo "ddtSchemes patch FAILED"; exit 1; }

sed -i 's/nNonOrthogonalCorrectors *1/nNonOrthogonalCorrectors 2/' system/fvSolution
cat >> system/fvSolution <<'EOF'

// Under-relaxation is what makes SIMPLE converge; it has no effect on WHICH
// fixed point is found, only on whether the iteration reaches it.
relaxationFactors
{
    equations { U 0.7; }
    fields    { p 0.3; }
}
EOF

sed -i -e 's/^startFrom .*/startFrom       startTime;/' \
       -e 's/^startTime .*/startTime       0;/' \
       -e "s/^endTime .*/endTime         $NITER;/" \
       -e 's/^deltaT .*/deltaT          1;/' \
       -e 's/^writeControl .*/writeControl    timeStep;/' \
       -e "s/^writeInterval .*/writeInterval   $NITER;/" \
       -e 's/^adjustTimeStep .*/adjustTimeStep  no;/' \
       -e 's/^writeFormat .*/writeFormat     ascii;/' \
       system/controlDict

echo "$CASE: steady, AoA=$AOA, $NITER SIMPLE iterations, ascii output"
