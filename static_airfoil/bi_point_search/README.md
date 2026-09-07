# Bifurcation-point search, Re=500, 12-14 deg

## Why this is a separate directory and a different method

The amplitude route to the onset angle is dead (`FINDINGS.md` §17-18). Three
independent observables reject the weakly-nonlinear description over 14-24 deg,
and all three follow from one geometric fact: **at every angle measured, the
steady C_L lies BELOW the minimum instantaneous unsteady C_L**. The limit cycle
does not orbit the fixed point, so nothing that extrapolates a perturbation
expansion about that fixed point can work, however well it fits.

This search does not extrapolate anything. It measures the **linear growth rate
sigma directly, at each angle, starting AT the fixed point**, and finds where it
crosses zero. That is the definition of the bifurcation point and it needs no
assumption about the nature of the bifurcation -- super/subcritical, Hopf or
otherwise.

## Method

For each angle alpha in 12.0, 12.5, 13.0, 13.5, 14.0:

1. **Steady solve** (SIMPLE) at alpha -> the fixed point `x_s(alpha)`. Converged
   on the FORCE trace, not residuals (`analysis/steady_summary.py`).
2. **Seed the transient run with `x_s(alpha + 0.5)`.** The initial condition is
   then exactly `x_s(alpha) + [x_s(alpha+0.5) - x_s(alpha)]`: a small, smooth,
   physical perturbation about the alpha fixed point. `mapFields -consistent` is
   a 1:1 copy here because AoA is imposed by rotating the freestream, so every
   case shares a bit-identical mesh.
3. **Run time-accurate**, dt frozen at 5.5726e-4 (the onset-sweep value, so
   sigma is comparable across angles and against `re500_aoa14/16/18`). No field
   snapshots -- `forceCoeffs` only.
4. **Fit `sigma` to the C_L envelope.** `sigma > 0` unstable, `sigma < 0` stable.

## Why seed from a neighbouring STEADY state and not a limit cycle

The AoA 14 attempt (§17) warm-started from AoA 16's *saturated limit cycle* and
the sigma was uninterpretable: the mean-flow displacement was **1.9x the
oscillation amplitude**, so the trajectory was travelling between two attractors
while the oscillation evolved, and no slope of log A was an eigenvalue.

Seeding from a neighbouring *steady* state removes that entirely. The whole
departure from `x_s(alpha)` is `x_s(alpha+0.5) - x_s(alpha)`, which is small
(dC_L ~ 0.003, about 10% of the saturated amplitude at 14 deg) and evolves under
the linearised operator `L(alpha)`. Mean and oscillatory parts are then governed
by the same eigenvalues, which is the regime the measurement assumes.

sigma is independent of the perturbation's size, so the seed magnitude varying a
little between angles does not matter -- only that it stays linear.

## What would falsify the whole approach

If C_L grows straight to saturation with no clean exponential stretch, the
perturbation was too large or the mode is not isolated. Watch the local
`d(ln A)/dt`: a constant value over several e-folds is an eigenvalue; a value
that drifts is not. **Check the SPAN of the fit window before the R^2** -- an R^2
of 0.999 over 6 time units taught nothing (§17).
