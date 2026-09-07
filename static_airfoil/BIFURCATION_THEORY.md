# Bifurcation theory behind the Re=500 measurements

A companion to `ANALYSIS_GLOSSARY.md`. That file explains the **signal-processing**
vocabulary in the verdict lines (carrier, fundamental, subharmonic ratio, the torus
guard). This one explains the **dynamical-systems** theory those diagnostics are
standing in for: where `A^2 = k(alpha - alpha_c)` comes from, what `sigma` actually is,
why period doubling is a different kind of event from shedding onset, and -- the part
that matters most -- which of our numbers the theory licenses and which it does not.

Every formula below is tied to a measured number from this campaign, and the
assumptions are stated as assumptions.

---

## 1. The setting: fixed points and linear stability

Write the incompressible Navier-Stokes system abstractly as

```
    dx/dt = N(x; alpha)
```

`x` is the entire flow field (velocity and pressure on the mesh) and `alpha` is the
angle of attack, our control parameter. Re is held at 500 throughout, so `alpha` is
the only knob.

A **steady solution** (fixed point, base flow) `x_s(alpha)` satisfies `N(x_s) = 0`.
It exists at every angle we ran -- *existence is not stability*. Perturb it,
`x = x_s + eps*u`, and to first order in `eps`

```
    du/dt = L(alpha) u,        L = dN/dx evaluated at x_s
```

`L` is the Jacobian. Its eigenvalues `lambda = sigma + i*omega` govern small
perturbations: each mode behaves like `exp(sigma t) * exp(i omega t)`, an
oscillation at frequency `omega` whose envelope grows or decays at rate `sigma`.

- `sigma < 0` for **every** eigenvalue: the steady flow is linearly stable.
- `sigma > 0` for **any** eigenvalue: it is unstable and will not be observed.

The **leading** eigenvalue is the one with the largest `sigma`. It is the only one
that survives long enough to see, which is why a single scalar `sigma` is a
meaningful thing to measure at all.

> This is exactly the trap recorded in the upstream findings and in `CLAUDE.md`:
> *a converged steady solve is not evidence the flow is steady.* `steadyState`
> deletes the `ddt` term, so SIMPLE converges to `N(x_s) = 0` regardless of the
> sign of `sigma`. It finds the fixed point and tells you nothing about `L`.

## 2. The Hopf bifurcation

As `alpha` increases, the eigenvalues move. A **Hopf bifurcation** is the moment a
complex-conjugate pair crosses the imaginary axis:

```
    sigma(alpha_c) = 0,     omega(alpha_c) = omega_c != 0
```

The `omega_c != 0` is what makes it a Hopf rather than a saddle-node: the mode that
goes unstable is *oscillatory*, so what emerges is a limit cycle -- periodic vortex
shedding -- not a second steady state.

Generically the crossing is transverse, so near `alpha_c`

```
    sigma(alpha)  ~  sigma'(alpha_c) * (alpha - alpha_c)          [LINEAR in alpha]
```

**That linearity is the whole basis of extrapolating to the onset.** It is a
statement about the eigenvalue, and it holds only near `alpha_c`.

## 3. The Stuart-Landau equation

Linear theory says a perturbation grows without bound once `sigma > 0`, which is
false. The nonlinearity saturates it. Near the bifurcation, centre-manifold
reduction collapses the entire infinite-dimensional system onto the single complex
amplitude `A(t)` of the critical mode:

```
    dA/dt = (sigma + i*omega) A  -  (l_r + i*l_i) |A|^2 A
```

This is the **Stuart-Landau equation**. Splitting `A = R e^{i theta}` into a real
amplitude and phase decouples them:

```
    dR/dt      =  sigma R  -  l_r R^3                 (amplitude)
    dtheta/dt  =  omega    -  l_i R^2                 (phase)
```

`l_r` is the **Landau constant**: the leading nonlinear feedback of the mode on its
own growth. Everything we measure follows from these two lines.

### Supercritical vs subcritical -- the sign of `l_r`

| | `l_r > 0` (supercritical) | `l_r < 0` (subcritical) |
|---|---|---|
| saturation | cubic term opposes growth; a stable cycle exists for `sigma > 0` | cubic term *reinforces* growth; the bifurcating cycle is unstable |
| below onset | `R = 0` is the only attractor | a large-amplitude cycle can coexist with the stable fixed point |
| behaviour | amplitude grows continuously from zero | jump to finite amplitude, **hysteresis / bistability** |
| our evidence | `A^2` linear in `alpha` over five angles; AoA 14 decays to zero from a finite seed | not observed |

**This distinction is the single most important caveat in the whole campaign**, and it
is why AoA 14 was worth running. Under a subcritical bifurcation you can be below
`alpha_c` and still shed, if you are handed a large enough perturbation -- and every
case here was warm-started from a shedding state, i.e. handed exactly that.

## 4. The three predictions we actually use

### (a) Amplitude scaling -- how `alpha_c` was measured

Set `dR/dt = 0` in the amplitude equation. The non-trivial fixed point is

```
    R_sat^2 = sigma / l_r
```

and substituting the linear crossing `sigma = sigma'(alpha_c)(alpha - alpha_c)`:

```
    R_sat^2 = [sigma'(alpha_c) / l_r] * (alpha - alpha_c)  =  k (alpha - alpha_c)
```

**`R_sat^2` is linear in `alpha`, and its zero crossing is the onset angle.** This is
the entire content of `analysis/plot_onset.py`. It is also why we plot `A^2` and not
`A`: the square is the quantity theory says is straight. Equivalently
`R_sat ~ (alpha - alpha_c)^{1/2}` -- the square-root law.

**This prediction failed in practice, and the failure is instructive** (`FINDINGS.md`
§17). Fitted over 16-24 deg it gave `alpha_c = 14.83` and so predicted no shedding at
14 deg; AoA 14 sheds, saturating at `A = 0.03413`. Freeing the exponent,
`A = C(alpha - alpha_c)^p`, fits all six angles at `R^2 = 0.9966` with **`p = 0.786`,
not 0.5**, and imposing `p = 0.5` triples the rms error. The square-root law is
therefore disfavoured by the data over this range -- unsurprising, since §4c below
shows the linear regime has already been left. `alpha_c` and `p` also trade off
smoothly, so six points spanning 14-24 deg determine neither. Current position:
**onset below 14 deg, loosely 12.5-13.5**.

The lesson is about *where* the normal form applies, not whether it is right: it is
asymptotic as `alpha -> alpha_c`, and every angle we have may simply be too far
above onset for it.

### (b) Relaxation rate -- why a saturated case cannot show you `sigma`

Perturb about the saturated state, `R = R_sat + r`, and linearise:

```
    dr/dt = (sigma - 3 l_r R_sat^2) r = (sigma - 3 sigma) r = -2 sigma r
```

A converged limit cycle returns to itself at rate **`2 sigma`**, not `sigma`. Two
consequences:

- A saturated case has `sigma_measured ~ 0` **by construction** -- the envelope is
  flat because it is *at* the attractor. Our converged period-1 cases run
  `|sigma| <= 3.6e-6`. That is a null result, not a small growth rate, and it is why
  the amplitude method exists at all: the rate is unmeasurable exactly where the
  amplitude is easiest to measure.
- In principle the warm-start relaxation of AoA 16 and 18 (both seeded above their own
  `R_sat`) carries `2 sigma`. **In practice this failed**: it gave
  `sigma(16) = +0.209 > sigma(18) = +0.159`, which is backwards, at `R^2 ~ 0.87` on
  5-6 points. The relaxation window is short and confounded with the far field
  swinging round to the new angle, an unrelated transient on a similar timescale. The
  numbers were discarded rather than reported.

### (c) Frequency shift -- a check that the normal form is even applicable

The phase equation gives the *observed* frequency of the saturated cycle:

```
    omega_obs = omega - l_i R_sat^2 = omega_c + [omega'(alpha_c) - l_i k] (alpha - alpha_c)
```

so `f_carrier` should also be **linear** in `alpha` near onset. Measured slope
`df/dalpha`:

| interval | slope |
|---|---|
| 16 -> 18 | -0.01373 |
| 18 -> 20 | -0.01832 |
| 20 -> 22 | -0.02485 |
| 22 -> 24 | -0.02668 |

**It steepens monotonically -- it is not linear.** The normal form's linear regime is
therefore narrow, and only 16-18 sit convincingly inside it. This is an honest strike
against using the 22 and 24 points, and it is the reason the `alpha_c` fits are quoted
as a *range* over subsets rather than as one number with a standard error -- and, as
§4a records, that range was itself an artefact: the subsets disagree systematically
because the fitted form is wrong. **The frequency non-linearity here was the warning
sign, recorded and then not acted on.** It said the normal form's linear regime had
been left; the amplitude fit was quoted anyway for another day.

## 5. Critical slowing down

`sigma -> 0` as `alpha -> alpha_c`, so the timescale `1 / |sigma|` **diverges** at the
onset. This is not a numerical annoyance, it is the physics:

- A case just below onset takes forever to visibly decay; just above, forever to grow.
- Measured directly: AoA 14 has e-folding time **79-96 t**, and it sits only ~0.8 deg
  below onset. Reaching the noise floor would take `t ~ 490`, twice the driver's
  `T1_MAX = 240`.
- It is why bracketing the onset by running angles is the *expensive* way to find it,
  and why AoA 16 was the slowest-settling case in the campaign (`D_P2P = 0.00167` at
  the gate, still creeping, versus `1e-5` for AoA 18).

## 6. What our observable actually is

Theory is written in terms of `R`, the amplitude of the critical eigenmode of the
velocity field. We measure `C_L`, a scalar functional of that field. The connection:

```
    C_L(t) = C_L,steady  +  c1 * R cos(omega t + phi)  +  c2 * R^2  +  O(R^3)
```

- The **`O(R)` term** is the projection of the critical mode onto the lift functional.
  So `A(f_carrier) ~ R` to leading order -- this is what makes the `A^2` fit legitimate.
- The **`O(R^2)` term** is the mean-flow distortion plus the second harmonic, produced
  by the quadratic nonlinearity of the convective term. It is why `C_L` mean shifts
  with amplitude, and why harmonics appear at `2f`.

Two consequences for practice:

1. **Use the spectral carrier amplitude, not `C_L` peak-to-peak.** p2p sums the
   fundamental and every harmonic, so it mixes `O(R)` with `O(R^2)` in a ratio that
   drifts across the sweep. `A(f_carrier)` isolates the `O(R)` part.
2. The proportionality constant `c1` may itself drift slowly with `alpha`. That is an
   uncontrolled systematic in the `A^2` fit -- second-order, but real, and another
   reason to trust the near-onset subsets over the full range.

An untested prediction sitting right there: `C_L,mean - C_L,steady` should be
proportional to `R^2`, hence **linear in `alpha`** with the same zero crossing. It
needs the steady base-flow `C_L` at each angle, which `steadyState` can supply cheaply
(the unstable fixed point is exactly what SIMPLE converges to -- the trap in §1 becomes
useful here). That would be a third independent route to `alpha_c`.

## 7. Bifurcations OF the limit cycle: Floquet theory

Once shedding exists, the object whose stability matters is no longer a fixed point but
a **periodic orbit**. Linearise about it over one period `T`; the resulting map has a
monodromy matrix whose eigenvalues `mu_i` are the **Floquet multipliers**. A
perturbation is multiplied by `mu` each period, so:

- `|mu| < 1` for all non-trivial multipliers: the cycle is stable.
- One multiplier leaves the unit circle: the cycle bifurcates.

*Where* it leaves determines what happens, and the three cases are physically distinct:

| multiplier crosses | name | result | our diagnostic |
|---|---|---|---|
| `mu = +1` | fold / transcritical | cycle collides with another, or vanishes | not seen |
| `mu = -1` | **period doubling (flip)** | new cycle of period `2T` | `ORDER = 2` |
| `mu = e^{+-i psi}`, `|mu| = 1` | **Neimark-Sacker** | invariant torus, second incommensurate frequency | `QUASIPERIODIC` |

**Period doubling is not a Hopf bifurcation.** Hopf destabilises a *fixed point* and
creates an oscillation; the flip destabilises an *oscillation* and doubles its period.
Different objects, different linear operator (Jacobian vs monodromy matrix), different
eigenvalue condition. They share only the word "bifurcation". Our two events:

- shedding onset **below 14 deg**, loosely 12.5-13.5 and not yet pinned (§4a) -- Hopf,
  fixed point -> limit cycle;
- period doubling in **(24, 25]** -- flip, limit cycle -> period-2 cycle.

### Why `HARMONIC_DEV` is the right discriminator

A flip produces a **closed orbit**: its spectrum is `f/2` and integer multiples
thereof, because the motion is genuinely periodic with period `2T`. A Neimark-Sacker
produces a **torus**: two frequencies with an irrational ratio, so peaks appear at
`m*f1 + n*f2` and are *not* integer multiples of any single fundamental.

`harmonic_deviation()` measures exactly this -- the largest distance from an integer
over every significant peak, given a candidate fundamental. AoA 30 measured
**0.0000 across eight peaks**: a closed period-2 orbit, not a torus. That is a
Floquet-theoretic statement made with an FFT.

The same normal-form logic applies to the flip: near it, the period-2 component's
amplitude should grow like `(alpha - alpha_flip)^{1/2}`. **We have not tested that**,
and the measured `SUBHARM_RATIO` is non-monotonic (0.949 at 26, 1.458 at 28, 0.632 at
30), so the flip's own amplitude scaling is an open question, not a result.

## 8. What we measure, and which theory licenses it

| measured | tool | theoretical object | status |
|---|---|---|---|
| `SIGMA` from a **decaying** envelope | `growth_rate.py` | `Re(lambda)` of the leading eigenvalue of `L` | valid -- amplitude is small, dynamics linear. AoA 14: `-0.01266`, `R^2 = 0.9988` |
| `SIGMA` from a **saturated** envelope | `growth_rate.py` | nothing; the state is *at* the attractor | null by construction, `<= 3.6e-6`. Used only as a saturation test |
| `SIGMA` from a **growing** envelope | `growth_rate.py` | a nonlinear transient rate | **not** an eigenvalue; finite-amplitude seed cannot separate super- from subcritical |
| `A(f_carrier)` on a converged run | `growth_rate.py` | `c1 * R_sat`, i.e. `R_sat` up to a constant | valid near onset; `c1` drift is an uncontrolled systematic |
| `alpha_c` from `A^2` vs `alpha` | `plot_onset.py` | zero crossing of `sigma(alpha)` | **falsified by AoA 14**; fit is form-biased, see §4a. Onset below 14 deg, unpinned |
| `CARRIER_FREQ` | `limit_cycle.py` | `omega_obs / 2pi` | measured; its *linearity* in `alpha` fails, see §4c |
| `ORDER`, `SUBHARM_RATIO` | `limit_cycle.py` | Floquet multiplier through `-1` | flip confirmed in (24, 25] |
| `HARMONIC_DEV`, `QUASIPERIODIC` | `limit_cycle.py` | flip vs Neimark-Sacker | closed orbit, not a torus (dev 0.0000 at AoA 30) |

## 9. What the theory does **not** license here

1. **Supercriticality is inferred, not proved.** The linearity of `A^2` over five
   angles is strong evidence (a subcritical branch turns back and could not give
   `R^2 = 0.9988`), and AoA 14 decaying to zero *from a finite-amplitude seed* is
   direct evidence at that angle. Neither is a hysteresis loop. The decisive test is
   to approach an angle in (14, 16) from below and from above and show no gap.
2. **One clean `sigma` cannot give a rate-based `alpha_c`.** A line needs two points.
   AoA 12 would supply the second. What one point gives is a bound: the decay at
   AoA 14 is *purely* exponential over a 30% amplitude drop, so `l_r A^2 << |sigma|`
   there, hence `l_r < 3.8`; feeding that through `sigma(16) = l_r R_sat(16)^2` puts
   the zero crossing at `alpha_c` between **14.8 and 15.5 deg** -- consistent with the
   amplitude fit, but a bound and not a measurement.
3. **The normal form is asymptotic as `alpha -> alpha_c`.** AoA 24 sits ~9 deg above
   onset. The frequency data (§4c) shows the linear regime has already been left well
   before there, so the full-range `A^2` fit is being used outside its domain of
   validity and survives only empirically.
4. **`Re sin(alpha)` is not a criterion.** Already recorded in `CLAUDE.md`: it is a
   bluff-body scaling, and for an airfoil the wake width is set by the separation
   point. Nothing above depends on it.

## 10. Symbols

| symbol | meaning | where it appears |
|---|---|---|
| `alpha` | angle of attack, the control parameter | everywhere |
| `alpha_c` | critical angle, the shedding onset | below 14 deg, unpinned |
| `x_s` | steady base flow (fixed point) | §1 |
| `L` | Jacobian about `x_s` | §1 |
| `lambda = sigma + i omega` | eigenvalue of `L` | §1 |
| `sigma` | growth rate; `> 0` unstable | `SIGMA` |
| `omega` | angular frequency; `f = omega / 2pi` | `CARRIER_FREQ` |
| `A`, `R` | complex / real amplitude of the critical mode | `A(f_carrier)` |
| `l_r`, `l_i` | Landau constants, real and imaginary | §3, bounded in §9 |
| `R_sat` | saturated amplitude, `= sqrt(sigma / l_r)` | `A_LAST` on a converged run |
| `k` | slope of `A^2` vs `alpha`, `= sigma' / l_r` | `plot_onset.py` |
| `mu` | Floquet multiplier of the limit cycle | §7 |
| `T` | period of the limit cycle | `PERIOD`, `SAMPLE_PERIOD` |

---

See `ANALYSIS_GLOSSARY.md` for the spectral vocabulary, `FINDINGS.md` §9 and §12-13 for
the period-doubling measurements, and `FINDINGS.md` §16 for the onset.
