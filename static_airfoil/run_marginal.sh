#!/bin/bash
# run_marginal.sh <aoa> <srcAoA> <eps> [TMAX]
#
#   ./run_marginal.sh 11.5 12.0 0.1 120.5
#
# Perturbation run for an angle NEAR the bifurcation point, where the usual
# recipe breaks in two ways.
#
# 1. THE SEED IS TOO BIG.  Near onset the saturated amplitude collapses, so
#    "map the base flow from half a degree up" stops being a small perturbation:
#    at 11.5 deg that is dCl ~ 0.0053 against an expected saturated amplitude of
#    0.0018-0.0030.  Starting ABOVE the attractor makes the perturbation decay
#    ONTO it, which reads as sigma < 0 and would wrongly call the angle stable.
#    So the seed is blended: x_s(a) + eps*(x_s(src) - x_s(a)).
#
# 2. THE STARTUP TRANSIENT IS TOO BIG.  mapFields drops phi, the solver rebuilds
#    it, and the first pressure solve threw C_L to -9 at AoA 12 -- fine there
#    because it decayed to ~1e-3, far below that angle's 7.2e-3 attractor, but
#    at 11.5 the attractor may be 1.8e-3 and the startup would swamp the seed.
#    REDUCED, NOT REMOVED, by starting from the STEADY CASE ITSELF so its own
#    converged phi is kept instead of rebuilt.  That is not sufficient on its
#    own: BLENDING U breaks the U-phi consistency (phi is the base flow's flux,
#    U is the blended one), so a divergence error of order eps remains and the
#    first pressure solve still has to absorb it.  Run straight at the
#    measurement dt this diverged in 180 steps with C_L jumping -0.10 on step
#    one.  The low-dt startup below is therefore still required -- keeping phi
#    only makes it gentler, and with eps=0.1 the residual error is ~10x smaller
#    than a full mapFields would leave.
set -e
cd "$(dirname "$0")"
source /opt/openfoam-dev/etc/bashrc >/dev/null 2>&1
A=$1; SRC=$2; EPS=$3; TMAX=${4:-120.5}
REFCASE=${REFCASE:-NACA0012/re500_aoa12}   # known-good TRANSIENT case to take fvSolution from
NP=6; DT=0.0005572582892; CHUNK=20; D=NACA0012/bi_point_search
C=$D/marg_$A

if [ ! -d "$C/processor0" ]; then
  [ -d "$D/steady_$A" ]   || { echo "ABORT: no $D/steady_$A"; exit 1; }
  [ -d "$D/steady_$SRC" ] || { echo "ABORT: no $D/steady_$SRC"; exit 1; }
  rm -rf "$C"; cp -a "$D/steady_$A" "$C"
  rm -rf "$C"/processor* "$C"/postProcessing "$C"/log.steady "$C"/log.decomposePar \
         "$C"/log.reconstructPar "$C"/0
  TA=$(ls -d "$C"/[0-9]* | xargs -n1 basename | sort -g | tail -1)
  TB=$(ls -d "$D/steady_$SRC"/[0-9]* | xargs -n1 basename | sort -g | tail -1)
  echo "  base $D/steady_$A @ $TA   seed source $D/steady_$SRC @ $TB   eps=$EPS"
  mv "$C/$TA" "$C/0"
  for F in U p; do
    python3 analysis/blend_fields.py "$C/0/$F" "$D/steady_$SRC/$TB/$F" "$EPS" "$C/0/$F.new"
    mv "$C/0/$F.new" "$C/0/$F"
  done
  echo "  phi kept from the steady solution (reduces, does NOT remove, the startup jolt)"
  rm -f "$C"/0/Cp "$C"/0/wallShearStress
  # steadyState -> backward, and restore the TRANSIENT fvSolution.
  #
  # DO NOT regex the relaxationFactors block out.  `[^}]*` stops at the FIRST
  # closing brace -- the one ending `equations { U 0.7; }` -- so it deletes
  # `relaxationFactors {` and leaves `fields { p 0.3; }` plus a stray `}` at top
  # level.  A `grep -q relaxationFactors` guard PASSES on that wreckage because
  # the word itself is gone, which is how a malformed dictionary reached the
  # solver and diverged the run.  fvSolution carries no angle dependence, so
  # copy the known-good transient one wholesale.
  sed -i '/ddtSchemes/s/steadyState/backward/' "$C/system/fvSchemes"
  cp "$REFCASE/system/fvSolution" "$C/system/fvSolution"
  rm -rf "$C/0/uniform"      # SIMPLE's time state (value 6000, deltaT 1)
  foamDictionary "$C/system/fvSolution" -entry PIMPLE/nOuterCorrectors >/dev/null 2>&1 \
    || { echo "ABORT: fvSolution does not parse"; exit 1; }
  # `grep -q X && { ...; }` returns 1 under `set -e` when X is ABSENT, which is
  # the case we want -- use an if, not a short-circuit.
  if grep -q relaxationFactors "$C/system/fvSolution"; then
    echo "ABORT: relaxationFactors present"; exit 1; fi
  grep -q backward "$C/system/fvSchemes" || { echo "ABORT: ddtSchemes still steadyState"; exit 1; }
  NNO=$(foamDictionary "$C/system/fvSolution" -entry PIMPLE/nNonOrthogonalCorrectors -value)
  echo "  fvSolution copied from $REFCASE (parses OK), nNonOrthogonalCorrectors=$NNO"
  ( cd "$C" && decomposePar > log.decomposePar 2>&1 )
fi

CLS=$(tail -1 $D/steady_$A/postProcessing/forceCoeffs/0/*.dat | awk '{print $4}')
echo "  Cl_steady($A) = $CLS"
T=$(ls -d $C/processor0/[0-9]* | xargs -n1 basename | sort -g | tail -1)
if awk -v t="$T" 'BEGIN{exit !(t<0.5-1e-9)}'; then
  ( cd "$C" && sed -i -e 's/^startFrom .*/startFrom       latestTime;/' \
      -e 's/^deltaT .*/deltaT          4e-05;/' -e 's/^endTime .*/endTime         0.5;/' \
      -e 's/^writeControl .*/writeControl    runTime;/' -e 's/^writeInterval   .*/writeInterval   0.5;/' \
      -e 's/^adjustTimeStep .*/adjustTimeStep  no;/' system/controlDict
    mpirun -np $NP foamRun -parallel > log.startup 2>&1 || true )
  grep -qa "^End" "$C/log.startup" || { echo "  STARTUP FAILED"; tail -5 "$C/log.startup"; exit 1; }
  CO=$(grep -a "^Courant" "$C/log.startup" | awk '{print $6}' | sort -g | tail -1)
  CLN=$(tail -1 $C/postProcessing/forceCoeffs/*/*.dat | awk '{print $4}')
  echo "  startup to t=0.5 at dt=4e-5: Co_max=$CO  Cl=$CLN  q=$(awk -v a=$CLN -v b=$CLS 'BEGIN{printf "%+.3e", a-b}')"
  T=0.5
fi
while awk -v t="$T" -v m="$TMAX" 'BEGIN{exit !(t<m-1e-9)}'; do
  T=$(awk -v t="$T" -v c=$CHUNK 'BEGIN{printf "%.10g", (int((t-0.5)/c)+1)*c+0.5}')
  ( cd "$C" && sed -i -e 's/^startFrom .*/startFrom       latestTime;/' \
      -e "s/^deltaT .*/deltaT          $DT;/" -e "s/^endTime .*/endTime         $T;/" \
      -e 's/^writeControl .*/writeControl    runTime;/' -e "s/^writeInterval   .*/writeInterval   $T;/" \
      -e 's/^adjustTimeStep .*/adjustTimeStep  no;/' system/controlDict
    mpirun -np $NP foamRun -parallel > log.chunk_$T 2>&1 || true )
  grep -qa "^End" "$C/log.chunk_$T" || { echo "  FAILED at t=$T"; break; }
  echo "  t=$T  $(python3 analysis/eigenvalue.py "$C" --cl-steady "$CLS" 2>/dev/null \
      | awk -F= '/^SIGMA_OSC|^R2_OSC|^FITSPAN_OSC|^OSC_FRAC|^FREQ|^A_OSC_LAST/{printf "%s=%s  ", $1, $2}')"
done
echo "$(date '+%F %T')  MARGINAL RUN DONE: AoA $A"
python3 analysis/eigenvalue.py "$C" --cl-steady "$CLS" --table
