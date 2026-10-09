# AoA 12.0 -- eigenvalue measurement, then promoted to a campaign case

## The measurement (this directory)

Seeded with `x_s(12.5)` mapped onto the 12.0 case, i.e. a perturbation
`dCl = +0.004248` about the 12.0 fixed point.  Startup at dt=4e-5 to t=0.5 to
rebuild `phi` (skipping it diverged the first attempt -- see
`pert_diverged_attempt1.log`), then dt frozen at 5.5726e-4.

**AoA 12.0 is LINEARLY UNSTABLE, and the mode OSCILLATES.**

| quantity | value |
|---|---|
| `sigma_osc` (linear, t=10-16) | **+0.1097** |
| frequency | **0.600011** |
| oscillatory fraction | 0.35-0.40 |

The oscillatory amplitude grows 7x over t=10-36 (8.195e-4 -> 5.703e-3).  The
rate falls as it saturates -- +0.1097 (t=10-16), +0.0962 (16-22), +0.0746
(22-28), +0.0370 (28-34) -- so **the eigenvalue is the EARLY rate**; the tail is
nonlinear.  (`analysis/eigenvalue.py` defaults to fitting the LAST 60% of the
record, which is right for a decaying case and wrong for a growing one.  It
reported 0.053, half the true value, until the head was used instead.)

A complex eigenvalue pair means the LINEAR onset is a Hopf, even though the
SATURATED states at 14-24 deg obey none of the weakly-nonlinear Hopf scalings
(`FINDINGS.md` §17-18).  Those are consistent: a Hopf can be the linear
mechanism while the attractor sits far from the fixed point.

sigma = +0.11 is not marginal (e-folding time 9 t), so **the bifurcation point is
well below 12 deg.**

## Preserved here

- `eigen_record/pert_12.0_postProcessing/` -- the force history the sigma was
  measured from, as it stood at t=41.9
- `eigen_record/pert_12.0_eigenvalue.txt` -- the full readout and local-rate table
- `eigen_record/pert_12.0_log.startup` -- the phi-rebuild startup
- `pert_12.0` -> symlink to `../re500_aoa12`

## Promotion

The case was moved to `re500_aoa12` and handed to `run_staged.sh` with
`DT_FIXED=5.5726e-4`, which finishes phase 1 until the limit cycle converges and
then records **20 orbits at 32 samples per carrier period** -- the campaign
protocol, so the snapshots join the POD set alongside 14-30 deg.
