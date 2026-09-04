#!/bin/bash
# warm_start.sh <srcCase> <srcTime|latestTime> <dstCase> [nProcs]
#
# Seed <dstCase> with the converged field from <srcCase> so its transient starts
# near the attractor instead of from uniform flow.
#
# Why this is exact here: AoA is imposed by ROTATING THE FREESTREAM, not the
# mesh, so every case at a given Re shares a bit-identical mesh.  There is a 1:1
# cell correspondence and mapFields -consistent is a direct copy, not an
# interpolation.  The mesh check below enforces that assumption rather than
# trusting it.
#
# The velocity vectors are deliberately NOT rotated.  The airfoil occupies the
# same cells in both cases, so the expensive thing to develop -- the boundary
# layer and near wake -- is already right.  Only the far field is off by the
# angle difference, and that corrects itself by advection within a few convective
# times.  Rotating the field would do the opposite: fix the cheap part and
# scramble the expensive part.
#
# The target keeps its OWN 0/U and 0/p, so the freestreamValue entries written
# by make_case.sh for the target angle survive; mapFields overwrites values, not
# boundary-condition definitions.  Verified at the end.
set -e
SRC=$(cd "$1" && pwd); STIME=$2; DST=$(cd "$3" && pwd); NP=${4:-4}
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1

sig() { cat "$1"/constant/polyMesh/{points,faces,owner,neighbour} 2>/dev/null | md5sum | cut -d' ' -f1; }
[ "$(sig $SRC)" = "$(sig $DST)" ] || {
  echo "ABORT: source and target meshes differ -- -consistent would be invalid."
  echo "  Use plain mapFields (interpolating) instead, and expect a longer transient."
  exit 1; }
echo "mesh identity: OK (same polyMesh)"

cd "$DST"
[ -d processor0 ] || { decomposePar > log.decomposePar 2>&1; echo "decomposed onto $NP"; }
sed -i -e 's/^startFrom .*/startFrom       startTime;/' -e 's/^startTime .*/startTime       0;/' system/controlDict

# NOT an MPI application: mapFields runs serially and reads the decomposed
# cases itself via -parallelSource/-parallelTarget.  Launching it under mpirun
# fails with "Invalid option: -parallel".
mapFields "$SRC" -consistent -sourceTime $STIME \
       -parallelSource -parallelTarget > log.mapFields 2>&1
grep -q "^End" log.mapFields || { echo "mapFields FAILED"; tail -30 log.mapFields; exit 1; }
echo "mapped $SRC @ $STIME -> $DST @ 0"

# The whole point of keeping the target's own 0/U: confirm the freestream
# direction is still the TARGET's angle and was not overwritten by the source's.
# Check EVERY processor and report only those that actually own farfield faces --
# interior subdomains have nFaces=0 there and print an empty list, which looks
# like a failure but is not.
echo "--- farfield freestreamValue after mapping (processors owning faces) ---"
_seen=0
for _d in processor*; do
  _nf=$(awk '/^    farfield/{f=1} f&&/nFaces/{print $2; exit}' "$_d/constant/polyMesh/boundary" 2>/dev/null | tr -d ';')
  [ "${_nf:-0}" -gt 0 ] 2>/dev/null || continue
  _seen=$((_seen+1))
  echo "  $_d (nFaces=$_nf): $(foamDictionary $_d/0/U -entry boundaryField/farfield/freestreamValue -value 2>/dev/null | head -1)"
done
[ "$_seen" -gt 0 ] || echo "  WARNING: no processor reported farfield faces -- check the decomposition."
