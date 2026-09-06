# Glossary: the terms used in the convergence and bifurcation diagnostics

Written for someone who knows signals but has not been living inside
`analysis/limit_cycle.py`. Every number quoted is measured from this campaign.

---

## 0. What is being analysed

One signal: **`C_L(t)`**, the lift coefficient, written every single time step from
t=0. Nothing here is computed from the flow fields -- all of it comes from that one
column of numbers.

Time is in **convective times**, `t* = t U_inf / c`. With chord c=1 and U_inf=1 they
are the same number. So "t=100" means the freestream has travelled 100 chord lengths,
and a frequency of `f=0.4` means 0.4 oscillations per convective time. That is a
Strouhal number based on chord, `St = f c / U_inf`.

The **spectrum** is the magnitude of the FFT of `C_L(t)` after subtracting its mean.
A sharp peak at frequency `f` means "there is a periodic component at this frequency";
the peak's height is that component's **amplitude**, in the same units as `C_L`.

---

## 1. Fundamental vs carrier -- the distinction the whole analysis turns on

These are two different frequencies that happen to coincide in the simple case, which
is why it took a bug to notice they are not the same thing.

### Fundamental frequency, `f_fund`

**The reciprocal of the true period of the flow.** The time after which the flow
repeats *exactly*, so that every measurable quantity returns to the same value.

A signal that is periodic with period `T` has spectral content **only** at integer
multiples of `1/T` -- nothing in between, and crucially nothing *below*. So the
fundamental is the lowest frequency of which every peak in the spectrum is an integer
multiple. Reported as `FREQ`; `PERIOD = 1/f_fund`.

This is what sets **how long phase 2 runs**: "20 orbits" means 20 x `PERIOD`, twenty
genuine repeats of the flow.

### Carrier frequency, `f_carrier`

**The frequency of the vortex shedding itself** -- the fastest, visually obvious
oscillation you would see in a `C_L` trace or in the wake. Reported as
`CARRIER_FREQ`; `SAMPLE_PERIOD = 1/f_carrier`.

This is what sets **how often snapshots are written**: 32 samples per carrier period,
because the carrier is the fastest thing that has to be resolved.

### How they relate

| state | relationship | example |
|---|---|---|
| period-1 | `f_fund = f_carrier` | AoA 24: both 0.41819 |
| period-2 | `f_fund = f_carrier / 2` | AoA 28: 0.18782 and 0.37564 |

`ORDER` is their ratio: 1, or 2, or in principle higher.

Across the campaign the carrier is smooth and monotonic in angle -- 0.5213, 0.4716,
0.4182, 0.3756, 0.3197 for AoA 20, 22, 24, 28, 30 -- because it is a physical
quantity, the shedding rate, which just slows as the wake widens. The *fundamental*
jumps discontinuously (halves) at the bifurcation. That is why a 0.189 reading at
AoA 28 was obviously wrong: nothing physical drops by half between 24 and 28 degrees.

---

## 2. "Carrier" is a borrowed word, and it does not import all of its meaning

You are right to find this one uncomfortable. The analogy to AM radio is real but
partial, and the differences matter.

**Where it fits.** Past the bifurcation, successive shedding cycles alternate strong
and weak. At AoA 30 the `C_L` maxima go 1.16, 1.58, 1.16, 1.58, ... If you draw an
envelope through the peaks, that envelope itself oscillates -- at exactly half the
shedding frequency. Drawn on a screen it looks precisely like a 0.3196 signal being
amplitude-modulated by a 0.1598 envelope. Hence "carrier" for the fast oscillation
being shaped, and "subharmonic" for the envelope.

**Where it breaks down, and why it matters here:**

1. **There is no modulating source.** In AM, a carrier and an independent message
   signal are multiplied; you could switch the message off and the carrier would
   remain. Here there is one flow with one closed periodic orbit, and the two spectral
   lines are two Fourier components of that single orbit. They are exactly phase-locked
   forever, not because anything locks them, but because they are the same object.

2. **The carrier is not necessarily the largest component.** In AM the carrier normally
   dominates. Here the ordering flips with angle:

   | | AoA 28 | AoA 30 |
   |---|---|---|
   | A at fundamental | **0.2660** | 0.1646 |
   | A at carrier | 0.1825 | **0.2603** |

   Same physical structure, opposite dominance. This is exactly what broke the
   detector -- it assumed "largest peak = carrier" and so read AoA 28 as period-1 at
   0.189. See `FINDINGS.md` section 12.

3. **The frequency ratio is locked at exactly 2, not free.** An AM carrier and its
   modulation have no required relationship. A period-doubled orbit is at *exactly*
   half, to measurement precision, and that exactness is itself the diagnostic
   (section 4).

If you want a term with no radio baggage: `f_carrier` is the **shedding frequency**
and `f_fund` is the **orbit frequency**. Those are what they mean physically.

---

## 3. Harmonics, subharmonics, and the amplitudes

### Harmonics: above, expected, uninformative

Any periodic signal that is not a pure sinusoid has energy at `2f, 3f, 4f, ...`. That
is just Fourier's theorem -- the harmonics encode the *waveform shape*, nothing more.
A perfectly healthy period-1 wake has them: AoA 24 carries 19.5% of its peak amplitude
at `2f`.

So a peak above the shedding frequency tells you nothing about bifurcations.

### Subharmonics: below, and diagnostic

A component *below* the shedding frequency, at `f/2`, cannot be produced by waveform
shape. If the flow really repeated every `1/f`, there would be no content below `1/f`
at all. Its presence means the true period is longer than the shedding period -- the
flow needs two shedding cycles to come back to itself.

This asymmetry is what the detector exploits, with deliberately different thresholds:

- searching **below** the dominant peak: sensitive, 5% threshold, because anything
  down there is meaningful. Period-1 angles read 0.03%, 0.12%, 0.45% -- a real
  subharmonic stands out by orders of magnitude.
- searching **above**: conservative, 35%+ needed, because ordinary harmonics live
  there and must not be mistaken for a doubling.

### The amplitudes

- **`A(f_carrier)`, "carrier amplitude"** -- height of the spectral peak at the
  shedding frequency. How strong the basic vortex shedding is.
- **`A(f_fund)`, "fundamental amplitude"** -- height of the peak at the orbit
  frequency. In a period-1 state this *is* the carrier amplitude, same peak. In a
  period-2 state it is the strength of the cycle-to-cycle alternation.

Both are in `C_L` units, so 0.26 means that component swings `C_L` by +/-0.26.

### `SUBHARM_RATIO`

```
SUBHARM_RATIO = A(f_fund) / A(f_carrier)
```

"How strong is the alternation compared with the shedding it alternates." Dimensionless
so it can be compared across angles.

- **0** -- period-1, no alternation. (Reported as exactly 0 when `ORDER=1`.)
- **0.632** at AoA 30 -- the alternation is 63% as strong as the shedding.
- **1.459** at AoA 28 -- the alternation is *larger* than the shedding oscillation.

It genuinely can exceed 1, and that is not a bug. It was defined as
`A(f_fund)/A(max)` at first, which silently returns 1.000 for **every**
fundamental-dominant case regardless of true strength -- useless as a comparison across
angles. Dividing by the carrier instead is what makes it a bifurcation measure.

`D_SUB` is the same quantity compared between the two convergence windows: it asks
"has the alternation stopped changing?" It matters because the subharmonic saturates
*later* than everything else -- mean, amplitude and frequency can all agree while the
orbit is still splitting.

---

## 4. The torus guard: why a non-integer ratio is disqualifying

This is the sharpest idea in the whole diagnostic, so it is worth the space.

### Two ways to grow a new low frequency

A limit cycle losing stability can produce a new, lower frequency by two quite
different routes, and **they look similar in a spectrum if you only check amplitudes.**

**(a) Period doubling** (flip bifurcation). The trajectory remains a single **closed
loop**, but it now takes two laps of the old cycle to close. The motion is still
exactly periodic, with period `2T`. Therefore all spectral content sits at integer
multiples of `1/(2T)`: you see `f/2, f, 3f/2, 2f, ...`, and the ratio of the shedding
peak to the new one is **exactly 2.000000**.

**(b) Neimark-Sacker** (secondary Hopf, "torus" bifurcation). A second, genuinely
**independent** oscillation appears, at some frequency `f2` that has no reason to be
rationally related to `f1`. The trajectory now moves on the surface of a torus -- it
winds around, and because `f2/f1` is irrational it **never closes and never exactly
repeats**. This is *quasi-periodic* motion. The spectrum fills with combination
frequencies `m*f1 + n*f2`, which are not integer multiples of any single frequency.

### The test

Take every significant peak, divide by the candidate fundamental, and ask whether the
results are integers:

```
AoA 30, converged:
  f        f/f_fund     deviation from integer
  0.15983    1.0000        0.0000
  0.31965    2.0000        0.0000
  0.47948    3.0000        0.0000
  0.63931    4.0000        0.0000
  ... through 8f          max deviation 0.0000
```

Eight peaks, all exact integer multiples. That is one closed period-2 orbit, not a
torus. Reported as `HARMONIC_DEV` (max deviation from an integer, in units of
`f_fund`); `QUASIPERIODIC=yes` fires above 0.05. AoA 28 reads `HARMONIC_DEV=0.0001`.

The same check appears in `spectral_structure()` as a gate on `ORDER`: the ratio
`f_carrier/f_fund` must be within 0.15 of an integer or the state is not called
period-m at all.

### Why this is not academic

**It decides whether the entire snapshot strategy is valid.** Everything downstream
assumes the flow is periodic:

- "32 samples per period" needs there to *be* a period;
- "20 orbits" needs the orbit to *close*;
- POD on a periodic ensemble is a well-posed thing; POD on a truncated slice of
  never-repeating motion is not the same object, and a ROM built on it would be
  extrapolating from an arbitrary window.

A torus would invalidate that premise while still producing a plausible-looking
`C_L` trace, a plausible dominant frequency, and a plausible snapshot set. Checking
commensurability is the cheap test that rules it out.

### It fired for real, once

During AoA 28's transient the early comparison window carried a spurious peak at
0.0992 alongside the fundamental at 0.1848. The ratio `0.371 / 0.0992 = 3.72` is not
close to an integer, so the guard refused to call it period-4 and fell back to
`ORDER=1`, which produced `D_SUB=1.0` and blocked convergence. That was correct: the
0.0992 content was a decaying transient (0.173 -> 0.113 -> 0.079 -> 0.047 of peak
amplitude over successive windows) and it disappeared, after which the run converged
with `D_SUB=0.00063`.

---

## 5. Shorter definitions for the rest

| term | meaning |
|---|---|
| **limit cycle** | A closed, self-sustaining oscillation the flow settles onto and returns to after a disturbance. What phase 1 is waiting for. |
| **attractor** | The state a system settles into. Here: a period-1 or period-2 limit cycle. |
| **orbit** | One complete traversal of the limit cycle -- one *fundamental* period, i.e. two shedding cycles in a period-2 state. |
| **bistability** | Two stable attractors coexisting at the same angle, so the initial condition decides which you find. The reason warm-start direction matters. |
| **p2p** (`CL_P2P`) | Peak-to-peak: 99th minus 1st percentile of `C_L` in the window. Percentiles rather than min/max so one outlier cannot set it. |
| **prominence** | Peak height divided by the *median* of the spectrum. A sharp real peak scores in the thousands; broadband noise with no real peak scores ~1. This is the guard against calling a quiet, nearly-steady plateau "shedding". |
| **noise floor** | The background level of the spectrum where no real signal exists -- about 1e-5 in `C_L` units here. A "subharmonic" at the noise floor is not a subharmonic. |
| **`D_MEAN`, `D_P2P`, `D_FREQ`, `D_SUB`** | Relative disagreement in mean, amplitude, frequency and subharmonic ratio between the first and second halves of the test window. All four must be under `--tol` (0.01) to declare convergence. |
| **`FULL_WINDOW`** | Whether the test window was cut to its full requested length or clamped by the start of the series. Replaced a cycle-recount that failed by rounding. |
| **growth rate `sigma`** | Fitted exponent in `A ~ exp(sigma*t)`. Negative = the perturbation dies (that state is stable); positive = it grows. AoA 24 measured -0.137/t. |
| **e-folding time** | `1/|sigma|`, the time for a perturbation to grow or shrink by a factor of e (~2.718). AoA 24: 7.3 convective times. |
| **critical slowing down** | Near a bifurcation the growth rate passes through zero, so perturbations decay ever more slowly and runs take ever longer to converge. AoA 22 e-folds in 3.5 t, AoA 24 in 7.3 t -- the signature of approaching onset. |
| **spectral leakage** | A finite window smears a sharp peak across neighbouring bins, so a strong peak can put false energy into nearby ones. Why short windows are untrustworthy at the low-frequency end, and why a Hann window is applied before every FFT. |
| **bin width** | Frequency resolution of an FFT, `1/(window length)`. A 32-convective-time window resolves 0.031 -- so two frequencies closer than that cannot be told apart. |
| **sawtooth check** | Unrelated to any of the above: a test for a solver pathology where the solution alternates every other time step. Reports % of steps that reverse direction; clean runs sit near 10-27%, corrupted ones peg at 100%. |

---

## 6. Reading a verdict line

```
t=100.5  shedding=yes  period-2  f_orbit=0.18768707 f_carrier=0.37537413
         Cl=0.903820+-0.808947  d(mean,p2p,f,sub)=(0.00050,0.00123,0.00025,0.00063)
         converged=yes
```

- `shedding=yes` -- there is a real oscillation, not a steady solution or a quiet plateau.
- `period-2` -- the orbit closes after two shedding cycles.
- `f_orbit=0.1877` -- the flow repeats every 5.33 convective times. Sets run length.
- `f_carrier=0.3754` -- vortices shed every 2.66 convective times. Sets snapshot rate.
- `Cl=0.9038+-0.8089` -- mean and peak-to-peak.
- the four `d` values -- window-to-window disagreement, all needing to be under 0.01.
- `converged=yes` -- phase 2 may begin.
