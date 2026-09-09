# The shedding bifurcation point at Re=500

**Answer: `alpha_c = 11.46 +- 0.1 deg`**, measured from linear growth rates about
the steady base flow. Bracketed: the flow is **stable at 11.0 deg** and
**unstable at 11.5 deg**.

---

## 1. What "e-folding time" means

`sigma` is a **rate**, in units of 1/time. A small perturbation of amplitude `A`
about the steady flow evolves as

```
    A(t) = A(0) * exp(sigma * t)
```

The **e-folding time** is `1/|sigma|`: the time for the amplitude to change by a
factor of `e = 2.718...`

- `sigma > 0`: e-folding time is how long the disturbance takes to **grow** 2.72x.
- `sigma < 0`: how long it takes to **decay** to 1/2.72 = 37% of its size.

Time here is convective time units (chord / freestream speed), so "e-fold 14 t"
means the disturbance shrinks by 63% in the time the flow travels 14 chord
lengths.

It is just a readable way to say how fast: `sigma = -0.07093` and "e-folding time
14.1 t" are the same statement.

| AoA | `sigma` | e-fold | in words |
|---|---|---|---|
| 11.0 | -0.07093 | 14.1 t | dies out; 63% gone every 14 chords |
| 11.5 | +0.00674 | 148.4 t | grows, but very slowly |
| 12.0 | +0.10970 | 9.1 t | grows briskly |

**The e-folding time diverges at the bifurcation point.** At `alpha_c`,
`sigma = 0` and the e-folding time is infinite: a disturbance neither grows nor
decays. That is *critical slowing down*, and the 16x jump in e-folding time
between 12.0 and 11.5 deg is it being measured directly.

---

## 2. The three measurements

Each: solve for the steady base flow `x_s(alpha)` (SIMPLE), perturb it slightly,
integrate time-accurately, and fit the amplitude of the `f = 0.600` mode.

| AoA | `sigma` | R2 | fit span | evidence | source |
|---|---|---|---|---|---|
| 11.0 | **-0.07093** | 0.99999 | 50 t | decay over **2.2 decades** | collaborator's machine |
| 11.5 | **+0.00674** | 1.00000 | 63 t | growth **+118%**, 21 monotone windows | this machine |
| 12.0 | **+0.10970** | ~0.95 | 6 t (early) | growth 7x, then saturates | this machine |

All three are **the same mode**: `f = 0.600`, decaying at 11.0 and growing at
11.5 and 12.0. A complex-conjugate eigenvalue pair, so the **linear onset is a
Hopf bifurcation**.

`sigma(12.0)` is the weakest of the three: the case saturates, so only the early
window t=10-16 is linear and the fit span is short. It is not used for the
crossing (see below).

---

## 3. How `11.46` is obtained -- and how much it rests on `sigma(11.5)`

Linear interpolation between the two angles that bracket the crossing:

```
    alpha_c = 11.0 + 0.5 * |sigma(11.0)| / (|sigma(11.0)| + sigma(11.5))
            = 11.0 + 0.5 * 0.07093 / (0.07093 + 0.00674)
            = 11.457 deg
```

**The number is not sensitive to `sigma(11.5)`, and that is the point.** Rewrite
the same formula as

```
    alpha_c ~ 11.5 - 0.5 * sigma(11.5) / |sigma(11.0)|
            = 11.5 - 0.0475
```

`sigma(11.5)` enters only through a **small correction to the 11.5 endpoint**,
because it is 10.5x smaller than `|sigma(11.0)|`. Measured sensitivity:

| `sigma(11.5)` assumed | | `alpha_c` | shift |
|---|---|---|---|
| +0.00169 | 4x too small | 11.488 | +0.032 |
| +0.00337 | 2x too small | 11.477 | +0.021 |
| **+0.00674** | **as measured** | **11.457** | -- |
| +0.01348 | 2x too large | 11.420 | -0.036 |
| +0.02696 | 4x too large | 11.362 | -0.094 |
| +0.06740 | **10x** too large | 11.256 | -0.200 |

Even a **tenfold** error in `sigma(11.5)` moves the answer by 0.2 deg. A 20%
error in `sigma(11.0)` -- the number carrying most of the weight -- moves it by
0.008 deg.

So the confidence structure is:

1. **The bracket needs only the SIGNS.** Stable at 11.0, unstable at 11.5, so the
   crossing lies in (11.0, 11.5). This holds no matter what the magnitudes are.
   Both signs are unambiguous: 2.2 decades of monotone decay at 11.0, and +118%
   monotone growth across 21 windows at 11.5.
2. **The position inside that bracket** comes mostly from `sigma(11.0)`, the
   best-measured of the three.
3. **`sigma(11.5)` refines it by 0.05 deg** and is the least load-bearing input,
   despite being the run that took all night.

An independent check: before 11.5 was run at all, interpolating between 11.0 and
12.0 predicted `alpha_c = 11.393`. The 11.5 run then landed the answer at 11.457,
**0.06 deg away** -- an out-of-sample confirmation, not a refit. Fitting a line
through all three gives 11.416 (R2 = 0.9935).

---

## 4. Where the +-0.1 comes from

Not from the fits, which are far tighter than that:

- **dt differs between machines.** AoA 11.0 ran at dt = 7.13e-4, ours at
  5.5726e-4. The scheme's numerical damping is dt-dependent, so the two sigmas
  carry slightly different artificial dissipation. This is the largest known
  systematic and it is not quantified.
- **`sigma(alpha)` is assumed linear** across the half-degree. Standard near a
  transversal crossing, and supported by the three points (R2 = 0.9935), but the
  slope does curve: 0.155/deg over 11.0-11.5 against 0.205/deg over 11.5-12.0.
- **Single mesh, single dt, no independence study.** M2 only.

Not included: statistical error on the fits, which is negligible by comparison.

---

## 5. Why this required a different method

Three earlier estimates of the onset were published in this repo and all three
were wrong: **14.8**, then **12.5-13.5**, then **"below 14"**. Every one came
from extrapolating a property of the *saturated* state -- amplitude, frequency,
or mean-flow shift -- to where it would vanish.

That cannot work here, for a reason that is now measured rather than argued:

- `sigma` falls **16x** between 12.0 and 11.5 deg. The saturated amplitude
  collapses toward zero exactly where the timescale diverges, so the saturated
  state carries almost no information about the crossing.
- Three independent observables reject the weakly-nonlinear (Stuart-Landau)
  description over 14-24 deg: amplitude exponent **0.79** where theory needs 0.5;
  `f_carrier` **not linear** in alpha; mean-flow distortion `~ A^0.88` where
  symmetry **forces** `A^2`.
- All three follow from one geometric fact: **at every angle measured, the steady
  `C_L` lies BELOW the minimum instantaneous unsteady `C_L`**, by 0.44-0.94 of a
  peak-to-peak. The limit cycle does not orbit the fixed point -- it is displaced
  bodily from it -- so perturbation theory about that fixed point never applied.

The eigenvalue is immune to all of this: `sigma` is a property of the
linearisation about the fixed point and does not care where the attractor sits.

See `FINDINGS.md` sections 17-19 and `BIFURCATION_THEORY.md`.

---

## 6. Reproducing it

```bash
./run_bisearch_steady.sh 11.5              # steady base flow (SIMPLE)
./run_marginal.sh 11.5 12.0 0.1 120.5      # blended seed, low-dt startup, 120 t
python3 analysis/eigenvalue.py bi_point_search/marg_11.5 \
        --cl-steady 0.4472496 --freq 0.600011 --table
python3 analysis/plot_eigenvalue.py        # -> analysis/figs/eigenvalue_onset.png
```

Five things that had to be right, each of which failed first:

1. **Start AT the fixed point**, never from a limit cycle. A warm start from a
   saturated state left AoA 14 with a mean-flow displacement **1.9x** the
   oscillation amplitude and no interpretable sigma.
2. **Size the seed against the attractor.** Near onset the attractor collapses:
   the standard `x_s(a+0.5)` seed at 11.5 deg would have been 2-3x too big, so
   the perturbation would decay ONTO the cycle and read `sigma < 0` -- the wrong
   sign. Blended to `eps = 0.1` (`analysis/blend_fields.py`), landing the
   oscillatory seed at 1.4e-5, under 1% of the attractor, so the whole 115 t
   record stayed linear.
3. **Keep `phi`**, by starting from the steady case itself rather than mapping.
   This REDUCES the startup jolt 9x -- it does not remove it, and skipping the
   low-dt startup diverged the run in 180 steps.
4. **Always probe one known frequency** (`--freq`). Tracking the tallest peak
   read `sigma = -0.10` where the mode was at `+0.0067`: a decaying seed
   transient at f=0.2 reached 2.47e-4 against the eigenmode's 1.9e-5.
5. **Check the SPAN before the R2.** An R2 of 0.999 over 6 t and 2-3 independent
   points is worthless; 63 t and 21 windows is a measurement.

## 7. Related files

- `AOA12_NOTE.md` -- the AoA 12 eigenvalue measurement and its promotion to a
  campaign case
- `aoa12_limit_cycle/` -- the saturated field at 12 deg
- `eigen_record/` -- force history and readouts for AoA 12
- `marg_11.5/` -- the decisive run
- `steady_1*/` -- the six steady base flows, all converged to a C_L range of
  0.000e+00 over the last 10% of iterations
