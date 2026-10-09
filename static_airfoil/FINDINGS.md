# FINDINGS

> **Path note (2026-10-09):** case folders have since moved under `NACA0012/` (e.g. `re500_aoa30/` is now
> `NACA0012/re500_aoa30/`, `steady_re500/` is now `NACA0012/steady_re500/`). The text below keeps the paths
> as they were when it was written. See CLAUDE.md, "Folder layout".

Measurement log for the **Re=500 static-airfoil** campaign. Every number here was
measured on this machine unless explicitly marked as an estimate.

**This continues, and does not replace, the August Re=1000 campaign log at**
`~/OpenFOAM/raji-dev/run/Aug/unsteadyRe1000/FINDINGS.md`. That document is
authoritative for the `adjustTimeStep` pathology, shedding onset, mesh independence at
alpha=10, and the physics traps. Read it first.

Last updated 2026-09-04.

---

## 1. Measure the Courant number at t=0.5, not during the impulsive start

**The single largest error made so far.** Sizing dt from a Courant number sampled
early in the impulsive start inflates every wall-clock estimate by ~7x.

Co on the 110k mesh (Re=500, AoA 20), `dt = 4e-5`:

| t | Co_max |
|---|---|
| 4e-5 | 0.120 |
| **8e-4** | **0.515  ← peak** |
| 8e-3 | 0.282 |
| 0.18 | 0.072 |
| 0.5 | **0.0729  ← settled** |

The 49.5k mesh behaves identically in shape: peak 0.384 at t=8e-4, plateau 0.0282
reached by t≈0.24.

Sampling at t=0.008 gave Co=0.282 and therefore a dt **7.3x too small**, which turned
a ~15 h run into a projected 104 h. `run_fixeddt.sh` and `run_staged.sh` both sample
`log.phase0 | tail -1` at t=0.5 and are correct; the error was in an ad-hoc benchmark
built alongside them.

Consequence for the driver: Co is measured while the flow is still largely attached
and **climbs as the separated shear layers develop at AoA 20-40**. `run_staged.sh`
therefore targets Co=0.35 rather than 0.6 for the initial dt, and re-tunes dt between
phase-1 chunks from the Co each chunk actually produced (shrink above 0.80, grow below
0.15, never inside a running `foamRun`).

## 2. Solver settings

Re=500, 110k cells, 4 cores, measured over 50 steps from an identical restart state.

| variant | s/step | GAMG it/step |
|---|---|---|
| A `nOuterCorrectors 2, nCorrectors 2, nNonOrth 1` | 0.416 | 159 |
| B `nOuterCorrectors 2, nCorrectors 2, nNonOrth 0` | 0.650 | 143 |
| **C `nOuterCorrectors 1, nCorrectors 2, nNonOrth 1`** | **0.176** | **53** |
| D `nOuterCorrectors 1, nCorrectors 3, nNonOrth 0` | 0.402 | 84 |
| E = C + GAMG `nCellsInCoarsestLevel 200`, `pFinal 1e-7` | 0.176 | 50 |

- **PISO (C) is 2.4x faster than the 2-outer-corrector setting** with no accuracy cost
  at Co<1 with `backward` differencing. The baseline was doing 8 pressure solves per
  time step; PISO does 4, and the expensive final solve of each set (42-52 GAMG
  iterations against `pFinal tolerance 1e-8, relTol 0`) dominates.
- **Removing the non-orthogonal corrector makes it slower, not faster** (B, D). Max
  mesh non-orthogonality is ~64°, and without the correction the pressure solve
  converges worse. Keep `nNonOrthogonalCorrectors 1`.
- **GAMG tuning bought nothing** (E vs C). Don't spend time there.

The `aoa10/` case's `fvSolution` still has `nOuterCorrectors 2`. Its ~6 h run to t=30
(88k cells, dt=2.5e-4, 120k steps) implies ~0.18 s/step, which normalises to ~0.45 on
the 110k/4-core basis — matching variant A, not PISO. That case has a free 2.4x
available.

## 3. What actually sets the Courant limit: `nte`, not the leading edge

Intuition said the leading-edge clustering (`bump`) drives dt. It does not.

| mesh | cells | nte | TE cell | bump | Co @ 4e-5 (early sample) |
|---|---|---|---|---|---|
| ns200 nte40 nr250 | 110 000 | 40 | 2.06e-4 | 0.08 | 0.2746 |
| ns200 nte40 nr250 | 110 000 | 40 | 2.06e-4 | **0.25** | 0.2705 |
| ns160 nte32 nr180 | 63 360 | 32 | 2.57e-4 | 0.25 | 0.2205 |
| ns150 nte30 nr150 | 49 500 | 30 | 2.74e-4 | 0.30 | 0.2057 |
| ns180 nte20 nr180 | 68 400 | 20 | 4.11e-4 | 0.15 | 0.1666 |
| ns160 nte16 nr160 | 53 760 | 16 | 5.14e-4 | 0.15 | 0.1538 |

Relaxing `bump` 3x changed Co by 1.5%. Changing `nte` tracked Co almost exactly
(0.2746/0.2205 = 1.245 vs. spacing ratio 1.25). **The cells across the truncated
trailing edge are the smallest in the mesh** — the TE gap is 0.0082c, so `nte=40`
gives 2.06e-4 cells.

Caveat: pushing `nte` down raises non-orthogonality, because
`theta_w = pi*nte/(2*ns+nte)` narrows the wake block into a thin wedge. Max
non-orthogonality goes 64° (nte 30-40) → **71°** (nte 16-20). The speed from small
`nte` is partly bought with mesh quality.

(The Co column above is from the *early* sample and is ~7x inflated per §1. The
*relative* comparison is still valid — all rows were sampled identically — but do not
size dt from these numbers.)

## 4. Settled Courant and cost, corrected

Measured at t=0.5:

| mesh | cells | Co @ dt=4e-5 | dt @ Co=0.35 | s/step (4 cores) | BL cells | max non-orth |
|---|---|---|---|---|---|---|
| ns200 nte40 nr250 | 110 000 | **0.0729** | 1.92e-4 | 0.176 | 52 | 63.9° |
| ns150 nte30 nr150 | 49 500 | **0.0282** | 4.96e-4 | 0.059 | 39 | 64.1° |

Two concurrent 4-core jobs on this 12-core box cost **1.7x per step** vs. one job alone
(49.5k mesh: 0.059 → 0.102 s/step). Sequential 8-core runs are preferable, and are
required anyway if the second case is warm-started from the first.

**Projected** total for both runs (AoA 20 cold + AoA 30):

| mesh | sequential, AoA30 warm, 8 cores | concurrent 4+4, both cold |
|---|---|---|
| 110k | ~35 h | ~78 h |
| 49.5k | ~4.5 h | ~10 h |

Still estimated, not measured: **transient length** (assumed ~60-80 convective times
cold, ~20 warm — the dominant term and the reason the convergence gate exists),
**4→8 core scaling** (assumed 1.7x), and **shedding period** (assumed St-based,
T≈1.7 at AoA 20 and ≈2.5 at AoA 30).

## 5. Two-stage run design

Requirement: keep per-step `C_L`/`C_D` through the transient for convergence tracking,
but write **no** fields until the limit cycle is reached, then ≥32 snapshots/cycle for
30+ cycles.

This needs no change to `system/functions` between stages — see CLAUDE.md for the
contract. Phase 1 uses `writeInterval == endTime` to get exactly one write per chunk
(the restart point); this build has no `onEnd` writeControl.

The shedding period is **measured** by FFT of the `C_L` series, not taken from a
Strouhal correlation, and both `writeInterval` and `endTime` derive from it. The
convergence gate (`analysis/limit_cycle.py`) splits the trailing series into two
consecutive windows and requires mean, peak-to-peak and dominant frequency to agree to
1%. Mean alone settles long before amplitude does; amplitude alone misses the slow mean
drift noted in `experiment_plan.md`.

Per-step sampling of `C_p` and `tau_w` as *surface distributions* is not tractable —
`surfaces` writes a directory per sample, so per-step would be ~250k directories and
~1M files per run. They are sampled at field cadence instead. Per-step remains correct
and free for the *integrated* coefficients, which append to a single file.

## 6. Warm start

`mapFields -consistent` between AoA cases at the same Re is an exact 1:1 copy, because
AoA is imposed by rotating the freestream and the meshes are bit-identical. Verified:
after mapping AoA 20 → AoA 30, the target's `internalField` becomes
`nonuniform List<vector>` (the mapped solution) while `freestreamValue` remains
`(0.866025 0.5 0)` — the target's own angle.

Gotchas, all hit in practice:
- `mapFields` is **not** an MPI app; `mpirun ... mapFields ... -parallel` fails with
  `Invalid option: -parallel`. Run it serially with `-parallelSource -parallelTarget`.
- Volume fields only — `phi` is not mapped and is rebuilt from `U`.
- Keep the target's own `0/U`/`0/p`; copying the source time directory carries the
  source's far-field direction with it.
- Do not rotate the velocity field. The airfoil is in the same cells; rotating fixes
  the cheap part (far field) and scrambles the expensive part (boundary layer).

For the full matrix this chains: march AoA upward at fixed Re, warm-starting each run
from the previous, so only the first run at each Re pays a cold transient.

## 7. Open issues

- `aoa10/` has no `constant/polyMesh` and no time directories — the case is stripped,
  so its results are not reproducible from what remains. Its `functions` header claims
  Re=100 but `nu = 1e-3` means it ran at **Re=1000**.
- `aoa10/run_fixeddt.sh` calls `../analysis/sawtooth_check.sh`, which does not exist;
  the script aborts at that gate. The 2dt-sawtooth check it performed has no
  replacement yet in `run_staged.sh`.
- `meshing/test.msh` and `test1.msh` were both built at **Re=100** spacing
  (`h1 = 5e-3`); they under-resolve the boundary layer by 2.24x at Re=500.


## 8. Validating the convergence gate against the documented runs

`analysis/limit_cycle.py` was checked against the four upstream Re=1000 runs, whose
states are independently documented. It reproduces the published numbers:

| case | my f_shed | upstream f_shed | my Cl mean | upstream | verdict |
|---|---|---|---|---|---|
| AOA_10 (t=60) | 0.8549 | 0.857 | 0.4181 | 0.4179 | CONVERGED, SHEDDING |
| AOA_10fine (t=30) | 0.8656 | 0.869 | 0.4197 | 0.4186 | SHEDDING, not yet converged |
| AOA_15 (t=25) | 0.6989 | 0.667 | 0.6735 | 0.6640 | SHEDDING, not yet converged |
| AOA_8 | - | - | 0.3262 | 0.3264 | STEADY |
| AOA_5 | - | - | 0.2482 | 0.2482 | STEADY |

The AOA_15 frequency differs by 4.8% because upstream `shedding_analysis.py` takes a
bare `argmax` over 0.0667-wide bins (0.667 is exactly bin 10); `limit_cycle.py` adds
parabolic interpolation, so the difference is sub-bin, not a disagreement. This matters
because `writeInterval` is derived from the period.

The "not yet converged" verdicts are correct, not failures: AOA_10fine saturates at
t~22 and stops at t=30, AOA_15 saturates at t~17 and stops at t=25, so neither has the
two full saturated windows the gate requires.

Two bugs were found and fixed by this exercise:

- **The window was sized at exactly `2*min_cycles/f0`**, where `f0` is a rough
  first-pass frequency. The half-split then left window B holding marginally fewer than
  `min_cycles` (measured 2.99 for a requested 3), so the `>=` test failed permanently
  and **the gate would never have fired on any case** -- every run would have gone to
  `T1_MAX`. Fixed with a 1.15x margin.
- **The shedding test was `p2p > 1e-4`**, which classified the documented-steady
  alpha=5 and alpha=8 runs as shedding: slow drift in a steady trace produces a
  spurious low-frequency peak with high spectral prominence. Now uses the upstream
  criterion `amp > 5e-3`, which separates the documented cases cleanly (alpha 5/8 at
  8e-4 and 1.5e-3, alpha=10 at 2.5e-2).


## 9. Period doubling between AoA 20 and 30 at Re=500

**The Re=500 wake is not period-1 at all angles, and the analysis chain assumed it
was.** Found while verifying the AoA 30 run, which reported `CONVERGED=no` at the end
of a run that was in fact fully saturated.

### The measurement

Both runs, converged windows only (`analysis/plot_cl.py`):

| | AoA 20 | AoA 30 |
|---|---|---|
| carrier f | 0.5213 | 0.3196 |
| orbit f (fundamental) | 0.5213 | **0.1598** |
| subharmonic amplitude A(f/2) | 1.1e-5 | **0.1646** |
| A(f/2) / A(f) | 0.0001 | **0.632** |
| C_L mean | 0.7867 | 1.0126 |
| C_L peak-to-peak | 0.2813 | 0.8259 |
| per-cycle mean | constant | alternates 0.9639 / 1.0623 |
| per-cycle p2p | constant | alternates 0.41 / 0.83 |

At AoA 30 successive shedding cycles differ; the orbit closes only after **two**
oscillations of C_L, so the true period is 6.257, not 3.129.

### It is a closed period-2 orbit, not a torus

Both bifurcations put a new spectral line below the carrier, so amplitude alone does
not identify which happened. The discriminator is commensurability: a period-doubling
produces one closed orbit, so every peak is an exact integer multiple of the orbit
frequency; a Neimark-Sacker produces a second, incommensurate frequency and the ratios
are irrational.

All eight resolved peaks at AoA 30 sit at f/f0 = 1, 2, 3, 4, 5, 6, 7, 8 with
**maximum deviation 0.0000**. `harmonic_deviation()` in `analysis/limit_cycle.py`
performs this test and reports `HARMONIC_DEV`; `QUASIPERIODIC=yes` fires above 0.05.

### It is saturated, not still growing

The question that decides whether the recorded snapshots are usable. Two independent
measures, one time-domain and one spectral, over the converged window
(`analysis/figs/subharmonic_growth.png`):

- sliding-window A(f/2) = **0.16455 +/- 0.00001** (0.006% scatter);
- fitted exponential growth rate **sigma = +1.6e-6 /t**, i.e. an e-folding time of
  642 000 convective times against a 28 t window -- +0.00% change across it;
- it was within 0.1% of its saturated value from **t = 112.5**, before the gate fired
  at t = 120.5;
- peak-height alternation |A_n - A_{n+1}| flat at 0.40 from t ~ 40 onward.

For contrast the same measures at AoA 20 decay 4.5 decades to a 1e-5 floor -- that
number is the noise floor of the measure, not a weak subharmonic.

The subharmonic saturates **later than the carrier**: at AoA 30 the carrier was within
1% of its final amplitude well before the subharmonic was. A gate that watches only
mean/amplitude/frequency can therefore pass while the orbit is still splitting, which
is why `limit_cycle.py` now also requires `D_SUB <= tol`.

### What it broke

1. **Cycle counting.** The driver took the largest spectral peak as "the" period, which
   is the carrier. AoA 30's phase 2 ran `NCYCLES=20` carrier periods = **9.97 actual
   orbits** -- half the requested data. Sampling was unaffected and in fact better than
   asked (32/carrier period = 64 per orbit).
2. **The convergence gate.** Alternate cycles genuinely differ, so a window holding a
   non-integer number of fundamental periods -- or split into halves that do -- compares
   different phases of the orbit and reports the alternation as a difference between
   windows. That is what produced the spurious `CONVERGED=no`; the linear trend across
   the whole window is -0.15%.

### The fix

`analysis/limit_cycle.py` now separates two periods and reports both:

- `PERIOD` -- the fundamental, `ORDER` x carrier. Sets the phase-2 **run length** and
  the convergence window, which is cut to a whole number of fundamental periods and
  split at the exact midpoint so each half holds the same whole number.
- `SAMPLE_PERIOD` -- the carrier. Sets the field **write interval**, so `PER_CYCLE`
  samples per carrier period are kept whatever the order.

New keys: `ORDER`, `SUBHARM_RATIO`, `HARMONIC_DEV`, `QUASIPERIODIC`, `CARRIER_FREQ`,
`SAMPLE_PERIOD`, `D_SUB`. `spectrum()` also now resamples onto a uniform grid before
the FFT -- the series is only piecewise uniform because `run_staged.sh` re-tunes dt
between chunks, and a window straddling a boundary smeared the peak, which matters most
exactly when the peak being looked for is a subharmonic.

Regression: verdicts are **unchanged on every period-1 case** (upstream AOA_5, AOA_8,
AOA_10, AOA_10fine, AOA_15 and re500_aoa20 -- same SHEDDING and CONVERGED, frequencies
within 0.4%, all sub-bin). Only AoA 30 changes: `ORDER=2`, and `CONVERGED` flips
no -> yes with `D_MEAN=0.00003, D_P2P=0.00001, D_FREQ=0.00002, D_SUB=0.00020`.

### The cycle-count gate, regressed and fixed properly

The window-sizing fix above dropped the 1.15x margin recorded in section 8, on the
reasoning that cutting whole fundamental periods made it unnecessary. **It did not**,
and the same bug class came back on AoA 22.

The window is cut to `2*n_per` fundamental periods using the first-pass frequency
`f_fund0`; `cycles` then recounts that same span with the window-B frequency `f_fund`.
So `cycles = n_per * f_fund/f_fund0`, and it lands just under `n_per` whenever the
refined estimate comes in slightly low. AoA 22 at t=60.5:

```
CYCLES_IN_WINDOW=5.99   (required 6.0)
D_MEAN=0.00025  D_P2P=0.00094  D_FREQ=0.00067  D_SUB=0.00000   (tolerance 0.01)
CONVERGED=no
```

Every agreement test passed by two orders of magnitude; the run was blocked by one
hundredth of a cycle. Which side of 6.0 a case lands on is a coin flip -- AoA 20 and
AoA 30 got 6.00 and 6.01. Unlike the section 8 version this is intermittent rather
than permanent (each chunk re-rolls it), so the cost is wasted chunks, ~35 min each,
not an infinite hang.

Fixed by asking the question the check is actually for -- *was the window cut to full
length, or clamped by the start of the series?* -- which removes the disagreement
between the two frequency estimates instead of papering over it with a margin. New key
`FULL_WINDOW`. Regression across all eight available cases (upstream AOA_5, AOA_8,
AOA_10, AOA_10fine, AOA_15 and re500_aoa20/22/30): **every verdict unchanged**. AOA_15
sits at 5.99 cycles but is rejected on `D_MEAN=0.0153` either way, and AOA_10fine on
`D_P2P=0.2246`, so nothing marginal was let through.

### What is NOT established

- **The onset angle.** Bracketed only to (20, 30]. The 1-degree sweep localises it.
- **Mesh independence of the bifurcation.** The only mesh study in the campaign was at
  alpha=10 on a period-1 state. Period doubling is a delicate bifurcation and nobody
  has shown this one survives refinement. The subharmonic carrying 63% of the carrier
  is an argument that it is not marginal; it is not evidence.
- **Whether AoA 30 is bistable.** It was warm-started from AoA 20 and took 120 t to
  settle. If period-1 and period-2 attractors coexist there, the initial condition
  selects which one is found. A cold start at AoA 30 would settle this.
- dt changed twice during the AoA 30 transient (5.09e-4 -> 1.86e-4 -> 4.34e-4), though
  the whole converged window is at a single dt. Halving dt from t=120.5 would rule out
  a temporal-discretisation artefact for ~1-2 h.
- This is a 2D result. In 3D a spanwise instability may pre-empt or modify it.

### Extension to 20 orbits, and a counting trap

AoA 30 was extended in place from t=182.8768 to 245.6497 (`extend_phase2.sh`,
2 h 12 m on 6 ranks, 0.0548 s/step): **9.97 -> 20.00 orbits**, 1288 snapshots on a
uniform grid, 16 GB. Sawtooth 10.4%. Final gate: `ORDER=2, SUBHARM_RATIO=0.632326,
HARMONIC_DEV=0.0000, CONVERGED=yes` with D_MEAN 6e-5 and D_P2P, D_FREQ, D_SUB all
0.00000 -- ten more orbits changed C_L mean by 1.7e-4 (1.012583 -> 1.012750) and p2p
by 1.2e-5, so the state was as saturated as the growth-rate fit said it was.

`dt` and `writeInterval` were deliberately **not** recomputed. The fresh period
estimate would have moved the interval from 224 to 225 x dt, putting a seam in the
sampling rate partway through what is meant to be one uniformly weighted POD
ensemble. Keeping it byte-identical continues the same grid exactly; `writeInterval /
dt = 224.000000`, so writes still land on time steps.

**The trap:** OpenFOAM names a time directory from its *value*, so a write landing on
a whole number is `213`, not `213.0`. The glob `[0-9]*.[0-9]*` -- used in the driver's
own snapshot count -- silently skips those. Here it undercounted by one and made a
perfectly uniform series appear to have a 2x gap at t=213. Scripts now match on the
value (`count_snaps`); **any downstream loader must do the same** or it will drop
snapshots and mis-space its time axis.

Also note `t = t_conv` itself (120.5) is a phase-1 chunk restart point, not a phase-2
sample: it sits 0.0898 before the first grid point instead of 0.0972. It is a valid
field, but **drop it from the POD set** -- the uniform grid starts at 120.589797.

### Onset bracket after the odd-angle runs

The collaborator's odd-angle runs (8 ranks, same mesh and scripts) report **AoA 23
period-1** and **AoA 25 period-2**. With our AoA 20 (period-1) and AoA 30 (period-2)
that brackets the onset to **(23, 25]** -- from a 10-degree window down to 2 degrees.

**AoA 24 is therefore the decisive angle**, and it is the one most exposed to the
initial-condition question raised above: if period-1 and period-2 attractors coexist
near onset, a warm start from a period-1 state biases toward finding period-1. The
precedent says the bias is not fatal -- AoA 30 was warm-started from period-1 AoA 20
and still grew and saturated its subharmonic -- but the argument only runs one way.
A period-2 result at 24 from a period-1 seed is trustworthy; a period-1 result is not,
and would need a confirmation run seeded from AoA 30.

## 10. Interruption recovery, and the phase-2 resume gap

A power cut took the machine down at 22:36 on 2026-09-05, mid-chunk in AoA 22's phase
1; power returned at 23:32. Nothing was lost beyond the partial chunk. The staged
design is what made that true: **phase 1 writes one field snapshot per chunk purely as
a restart point**, so the worst case is re-running the chunk in progress. The t=20.5
restart point was verified intact on all six processors (all fields present, sizes
proportional to per-rank cell counts, `uniform/time` present) and all 2315
reconstructed snapshots in the two finished cases checked clean.

`last -x reboot` shows the *boot*, not the loss, so 23:32 was power returning; the mtime
of the interrupted chunk log gives the actual moment the machine died.

**The gap this exposed:** phase-1 chunk boundaries and phase-2 snapshots are both just
time directories -- nothing on disk distinguishes them. The resume logic assumed any
restart point was a chunk boundary, so an interruption during **phase 2** would have
resumed with `T_CONV` set to the interruption time and silently recorded the wrong
number of orbits. `run_staged.sh` now stamps a `.phase2` marker (`T_CONV`, `dt`,
`writeInterval`, `endTime`, `ORDER`) as phase 2 begins, and the resume branch reads it
to finish phase 2 on the *original* dt and write interval -- so the snapshots either
side of an interruption stay on one uniform grid.

`nohup` alone did not survive the event; use `setsid nohup ... &` so the job is in its
own session.

## 11. Even-angle runs, and the onset localised to one degree

AoA 22 and 24 were run on this machine (6 ranks, warm-started up the ladder
20 -> 22 -> 24), alongside the collaborator's odd angles. Both are **period-1**.

| AoA | order | f (orbit) | C_L | C_D | t_conv | snapshots | source |
|---|---|---|---|---|---|---|---|
| 20 | 1 | 0.5213 | 0.7867 +- 0.2813 | 0.4157 | 60.5 | 1025 | here |
| 22 | 1 | 0.4716 | 0.8399 +- 0.3282 | 0.4742 | 80.5 | 644 | here |
| 23 | 1 | - | - | - | - | - | collaborator |
| 24 | 1 | 0.4182 | 0.8790 +- 0.3967 | 0.5320 | 60.5 | 641 | here |
| 25 | 2 | - | - | - | - | - | collaborator |
| 30 | 2 | 0.1598 | 1.0128 +- 0.8259 | 0.7354 | 120.5 | 1288 | here |

**Onset is now bracketed to (24, 25].** Frequency falls and both coefficients rise
monotonically with angle across the whole set.

### Why the period-1 result at AoA 24 is trustworthy despite a period-1 seed

A period-1 verdict from a period-1 warm start is the one case that cannot be taken at
face value, so the subharmonic was measured directly rather than read off the detector
flag. It **decays**:

| t | A(f0/2) at AoA 24 | ratio to carrier |
|---|---|---|
| 31.9 | 0.002895 | 0.0157 |
| 41.5 | 0.000586 | 0.0031 |
| 60.5 (phase 2 starts) | 0.000262 | 0.0014 |
| >= 85 | 0.000071 +- 0.000005 | 0.0004 |

Fitted sigma = **-0.137/t** on the transient (e-fold 7.3 t), settling to 7.1e-5, which
is the same noise floor AoA 22 reaches (6.6e-5). The period-2 perturbation the warm
start injected died. That is a measurement of **linear stability of the period-1 orbit
at 24 degrees**, which a period-1 seed cannot manufacture -- a decaying perturbation is
evidence wherever it started.

It does **not** exclude a coexisting period-2 attractor at 24 (bistability); that still
needs a run seeded from AoA 30. But the question is now about coexistence, not about
whether the transient was integrated long enough.

### Critical slowing down, and why the onset is not extrapolated

The decay rate roughly halves from AoA 22 (sigma = -0.282, e-fold 3.5 t) to AoA 24
(sigma = -0.137, e-fold 7.3 t) -- the expected signature of approaching a bifurcation.

**Do not extrapolate it.** A linear fit of sigma against AoA over 20/22/24/30 returns
onset at 18.3 degrees, which contradicts every run: at AoA 20 the subharmonic is
already on the 1e-5 noise floor so the fit is fitting noise, and at AoA 30 it is
saturating so the measured rate is nonlinearly suppressed, far below the true linear
growth rate. Two clean points cannot locate the crossing. The bracket comes from the
runs.

### A blind spot in the convergence gate

AoA 24 converged at t=60.5, the earliest `T1_MIN` allows, with A(f0/2) still 4x above
the floor and falling. `D_SUB` is only computed when `ORDER > 1`, so a subharmonic that
is present but **decaying** is invisible to the gate -- the run can be declared
converged while that mode is still settling.

The cost here was small and was measured, not assumed: contamination at phase-2 start
is 0.14% of carrier amplitude, decaying to 0.04% and reaching AoA 22's floor by t~85.
An order of magnitude below any discretisation error, so the AoA 24 snapshot set was
left as it is.

It will get worse closer to onset, where critical slowing down lengthens the decay.
**Proposed:** report `SUBHARM_RATIO` at `ORDER=1` as well, and require it either below
a floor or stationary between the two windows before declaring convergence. Not yet
implemented.

## 12. The carrier is not always the tallest peak

Found at AoA 28, which reported `period-1  f_orbit=0.189` at its first chunk verdict --
against a carrier trend of 0.5213, 0.4716, 0.4182, ..., 0.3202 that leaves no room for
0.189 at 28 degrees.

The spectrum is unambiguous and is the same structure as AoA 30, with the **dominance
flipped**:

| | AoA 28 | AoA 30 |
|---|---|---|
| fundamental | 0.18931 **A=3006** | 0.15932 A=8580 |
| carrier (2x) | 0.37813 A=2018 | 0.32015 **A=12806** |
| ratio | carrier/sub = 0.671 | sub/carrier = 0.670 |
| higher harmonics | 3x, 4x, 5x | 3x ... 7x |

`subharmonic_order()` searched only BELOW the largest peak, so at AoA 28 it found
nothing and called the state period-1 at the orbit frequency. Consequences, had it
reached phase 2:

- `SAMPLE_PERIOD` = 5.29 instead of 2.65, so **32 samples per orbit = 16 per shedding
  cycle**, half the intended density, silently;
- `ORDER` reported 1, mislabelling the physics in the results table;
- `D_SUB` not computed, so the gate lost its subharmonic-stationarity check exactly
  where it is needed.
- `PERIOD` came out right, but by luck.

**`A(2f)/A(f)` cannot be used to detect this**: it reads 0.270 at AoA 30 (period-2) and
0.195 at AoA 24 (period-1) -- no threshold separates them.

### The fix: find the two frequencies independently

`spectral_structure()` replaces `subharmonic_order()`:

- **`f_carrier`** = the HIGHEST peak still carrying `carrier_frac`=0.5 of the maximum
  amplitude. "Resolve everything with at least half the peak's amplitude" is exactly
  the condition the sampling rate has to meet. Verified to pick 0.3767 at AoA 28,
  0.3197 at AoA 30, and the single peak at every period-1 angle.
- **`f_fund`** = the LOWEST peak above `min_ratio`=0.05, subject to the existing
  local-maximum and prominence guards. The looser threshold matters near onset, where
  the subharmonic emerges small and `carrier_frac` would miss it.
- **`ORDER`** = round(`f_carrier`/`f_fund`), rejected unless within 0.15 of an integer
  (an incommensurate ratio is a torus, not a subharmonic).
- `f_fund` is then re-derived as `f_carrier/ORDER`, since the carrier is the sharper
  peak and its refined location is the more precise of the two.

`SUBHARM_RATIO` is now `A(f_fund)/A(f_carrier)`, not `A(f_fund)/A(max)`. The old form
reads 1.0 for **every** fundamental-dominant case regardless of subharmonic strength;
the new one is comparable across angles, which is what makes it a bifurcation diagram.
AoA 28 reads 1.50, AoA 30 reads 0.63.

Regression: all nine pre-existing cases (upstream AOA_5/8/10/10fine/15 and
re500_aoa20/22/24/30) produce **byte-identical output**. Only AoA 28 changes, from
`ORDER=1, f=0.1885` to `ORDER=2, f_orbit=0.1884, f_carrier=0.3767`.

This was caught while AoA 28 was still in phase 1, so no snapshots were written at the
wrong rate.

## 13. The dominance crossover at AoA 26, and a metric that compared two different things

AoA 26 sits almost exactly at `SUBHARM_RATIO = 1` -- the angle where the subharmonic
and the carrier are the same height. That is the worst possible place for any code that
identifies a frequency by "which peak is tallest", because the answer flips on noise.

`window_stats()` still used `dominant_freq` (the largest peak) for its `freq` field,
even after `spectral_structure()` fixed the main path. At t=80.5:

```
window A (t=20.1-50.5)  ratio=1.0724  ->  tallest = FUNDAMENTAL  0.19922
window B (t=50.5-80.9)  ratio=0.9478  ->  tallest = CARRIER      0.39532

D_FREQ  = |0.19922 - 0.39532| / 0.39532 = 0.49604
```

So `D_FREQ` was comparing a fundamental against a carrier -- two different physical
quantities -- and reporting the factor of two as disagreement. On the carriers, which
is what the test means to compare, the two windows agreed to **0.00709**, well inside
the 0.01 tolerance.

It was pinned at ~0.5 across every chunk (0.50564 at t=60.5, 0.49636 at t=80.5) while
`D_MEAN`, `D_P2P` and `D_SUB` all improved normally. **The run could never have
converged**; it would have burned to `T1_MAX=240` and then proceeded with a warning on
whatever state it happened to be in.

Fixed by having `window_stats()` report the structure-aware carrier, and while there,
`order` and `ratio` too -- so `D_SUB` now reuses window A's already-computed ratio
instead of recomputing it. Regression across all ten other cases: **byte-identical**.
Only AoA 26 changes, `D_FREQ` 0.49564 -> 0.00685.

**The pattern worth remembering:** this is the third bug in the same family. Every one
came from identifying a frequency by peak height rather than by structure --
`subharmonic_order` searching only below the tallest peak (section 12),
`carrier_freq()` in `plot_cl.py` doing the same, and now `window_stats`. Amplitude
ordering is not a reliable label. Anywhere the code needs to know *which* frequency it
is looking at, it must ask `spectral_structure()`.

## 14. The even-angle set, complete

All six even angles 20-30 at Re=500 are run and converged on the M2 mesh (63,360
cells, 6 ranks), warm-started along the ladder 20 -> 22 -> 24 for the period-1 branch
and 30 -> 28 -> 26 for the period-2 branch.

| AoA | order | f_carrier | f_orbit | C_L mean | C_L p2p | C_D | SUBHARM_RATIO | snapshots |
|---|---|---|---|---|---|---|---|---|
| 20 | 1 | 0.52126 | 0.52126 | 0.7867 | 0.2813 | 0.4157 | 0 | 1025 |
| 22 | 1 | 0.47155 | 0.47155 | 0.8399 | 0.3282 | 0.4742 | 0 | 644 |
| 23 | 1 | - | - | - | - | - | - | collaborator |
| 24 | 1 | 0.41819 | 0.41819 | 0.8790 | 0.3967 | 0.5320 | 0 | 641 |
| 25 | 2 | - | - | - | - | - | - | collaborator |
| 26 | 2 | 0.39527 | 0.19764 | 0.9090 | 0.6694 | 0.5896 | 0.949 | 1289 |
| 28 | 2 | 0.37541 | 0.18770 | 0.9039 | 0.8090 | 0.6299 | 1.458 | 1282 |
| 30 | 2 | 0.31965 | 0.15983 | 1.0128 | 0.8259 | 0.7355 | 0.632 | 1288 |

Every case: `CONVERGED=yes`, `QUASIPERIODIC=no`, `HARMONIC_DEV <= 0.0002`, all four
window-agreement metrics at or below 1e-4, snapshot grids uniform with **zero** off-grid
intervals and cumulative drift below 6e-7, ~32 samples per carrier period, 20 orbits.
6179 snapshots, 77 GB.

### What the numbers say

**The carrier is monotonic and smooth** through the bifurcation -- 0.5213, 0.4716,
0.4182, 0.3953, 0.3754, 0.3197. It is a physical shedding rate and knows nothing about
the orbit splitting. The *orbit* frequency halves discontinuously at onset. That
contrast is what exposed the AoA 28 detection bug: no physical quantity drops by half
between 24 and 28 degrees.

**The subharmonic branch is single-humped**, not monotonic: 0 below onset, 0.949 at 26,
peaking near 1.458 at 28, falling to 0.632 by 30. No structural change is needed
between 28 and 30 to explain it -- an earlier guess that a second bifurcation might sit
there is not supported.

**C_L peak-to-peak jumps 69% across onset** (0.397 at 24 -> 0.669 at 26) against 17%
and 21% for the two steps before it. The orbit does not merely split; it grows sharply.

**C_L mean is non-monotonic** -- 0.8790, 0.9090, 0.9039, 1.0128 for 24, 26, 28, 30 --
with a shallow dip at 28 that coincides with the subharmonic maximum. It is the only
break in an otherwise clean trend. C_D is monotonic throughout.

### Still open

- **Onset is bracketed to (24, 25]** by AoA 24 (period-1, subharmonic decaying at
  sigma = -0.137/t) and the collaborator's AoA 25 (period-2). Not narrowed further.
- **Supercritical vs subcritical** is not settled. The curve between 24 and 26 is where
  a subcritical branch would show, and one point at 26 cannot resolve it. Nothing seen
  so far looks like hysteresis.
- **Bistability, mesh independence and dt independence** are all still untested, by
  decision, deferred to a later pass. Note that the period-2 angles were warm-started
  from the period-2 branch, so finding period-2 there is the seed confirming itself --
  the onset bracket rests on AoA 24 and 25, not on 26/28/30.

## 15. The dt controller diverged AoA 18, and the Courant number could not have caught it

**Date:** 2026-09-07. **Cost:** ~35 min of compute, no data lost.

The first attempt at AoA 18 (warm-started from the saturated AoA 20 state, to
measure how a shedding oscillation evolves below the shedding onset) blew up in
its second phase-1 chunk. `foamRun` died with SIGFPE inside `GAMGSolver::scale`,
and `analysis/limit_cycle.py` reported `CL_MEAN=2.5e206`.

### What happened

| | |
|---|---|
| phase 0 | `Co_max=0.0251` at t=0.5, sawtooth 22.8% (clean) |
| dt chosen | `5.5726e-4` (target Co 0.35) |
| chunk 1, t=0.5→20.5 | ran clean. `Co_max=0.417`, sawtooth 7.8%, `C_L=0.8014` at the write |
| re-tune | `0.417 < CO_LO=0.45` → grow. `dt × 0.75/0.417 = 1.80×` → **clamped at `DT_MAX=1.0e-3`** |
| chunk 2, t=20.5→… | diverged after **191 steps**, at t=20.69 |

### Why the Courant number was blind to it

At the moment of divergence the Courant number was **0.60** — below `RETUNE_CO`
(0.75) and below the `CO_HI` abort threshold (0.80). The controller's own
diagnostic said the run was comfortably inside its band while the solution was
already exploding. It was not a marginal miss; Co was never going to see this.

The reason is in the ratio: **mean Co was 0.0036 against a max of 0.60, a factor
of 166.** The maximum Courant number in this mesh is set by a handful of
trailing-edge cells (`nte` is the Courant limit — the table in §3 already said
so) and reports on nothing else. The failure is somewhere else entirely: the
mesh has 64° maximum non-orthogonality and runs `nNonOrthogonalCorrectors 1`, so
the non-orthogonal part of the pressure Laplacian is carried as a single lagged
deferred correction. That correction's weight grows with dt while Co, dominated
by cells that are small but well-shaped, barely moves.

**A Courant target is not a stability criterion on this mesh.** It is a
resolution criterion that happens to correlate with stability over the range
where it was calibrated. The ceiling has to be an absolute dt.

### The fix

1. **`DT_MAX` 1.0e-3 → 6.0e-4.** Chosen from evidence, not theory: `6.331e-4`
   (AoA 26) is the largest dt that has ever run to convergence on this mesh, and
   `1.0e-3` is the only one that has ever diverged. Every case in the completed
   campaign ran in 3.6e-4 – 6.3e-4.
2. **`DT_GROW_MAX=1.5`**, a cap on one re-tune's growth ratio. Secondary — AoA 26
   grew by 1.76× and survived, so the ratio alone does not discriminate — but it
   makes the controller walk up to a bad dt with a recoverable overshoot instead
   of leaping past every value it has evidence for.
3. **Bounded retry.** A diverged chunk loses nothing: the restart point is the
   previous chunk's write and is untouched. `try_run` now halves dt and goes
   again, up to twice, instead of ending the run and waiting for a human.
4. **`DT_FIXED`** env override — hold dt and disable re-tuning entirely (below).

### The trap inside the trap: a diverged chunk poisons the force history

`forceCoeffs` opens `postProcessing/forceCoeffs/<restartTime>/` on restart, and
`analysis/limit_cycle.py` concatenates **every** directory it finds (§ Restart
caveat). The diverged chunk had already written 191 steps of garbage rising to
`4.75e204` into `forceCoeffs/20.5/`. A retry restarts from the same t=20.5 and
therefore **appends to that same directory** — so the blow-up would have
survived the retry, silently, and contaminated C_L for the rest of the case.

`try_run` moves the directory to `.diverged/` *before* retrying. Quarantine, not
delete: it is evidence.

### Why dt is now frozen for the onset sweep

`run_onset.sh` exports `DT_FIXED=5.5726e-4` — AoA 18's own validated value — for
every angle in the sweep, which is a different trade-off from the POD runs:

- `analysis/growth_rate.py` fits an exponential to the C_L envelope across the
  whole phase-1 record. The scheme's numerical damping is dt-dependent, so a
  re-tune between chunks puts a **step into the quantity being fitted**, and that
  artefact is indistinguishable from the physical decay it exists to measure.
- Comparing σ *between* angles only means something if the numerical damping is
  identical at both. One shared dt makes σ(18) − σ(16) physics rather than
  discretisation.

Lower angles separate less and run at lower Co, so reusing AoA 18's dt is
conservative at 16.

For a POD run the re-tune stays on and stays worth it: dt only has to keep Co in
band there, and a change between chunks costs one Euler step out of ~100k.


## 16. The shedding onset, located from saturated amplitudes rather than a bracket

**Both AoA 18 and AoA 16 shed.** Neither brackets the onset from below -- but they
were not needed for that. The onset is `alpha_c ~ 14.8 deg`, and the method that
locates it was **tested out of sample before AoA 16 ran**.

| AoA | order | f_carrier | C_L mean | C_L p2p | C_D | A(f_carrier) | snaps |
|---|---|---|---|---|---|---|---|
| 16 | 1 | 0.58537 | 0.630309 | 0.139020 | 0.304484 | 0.069449 | 645 |
| 18 | 1 | 0.55790 | 0.711513 | 0.217655 | 0.356784 | 0.108713 | 643 |

Both converged at t=60.5, `FULL_WINDOW=yes`, every window metric at 1e-5 or below,
19.96 orbits each, ~32.2 samples per carrier period, zero off-grid intervals, 7.9 GB
each.

### Saturation was measured, not assumed

§11 records the convergence gate's blind spot for a slowly-settling amplitude, so
`analysis/growth_rate.py` was run on both. On the converged record every period-1 case
now sits in a narrow band:

| AoA | sigma | verdict |
|---|---|---|
| 16 | +3.58e-6 | saturated |
| 18 | -7.31e-7 | saturated |
| 20 | -1.96e-7 | saturated |
| 22 | +2.12e-6 | saturated |
| 24 | +3.47e-6 | saturated |

A live rate for comparison is `+2.1e-2` (AoA 30's subharmonic, §9) -- four orders
larger. **Anything under ~1e-5 is saturated**, and that threshold is now calibrated
rather than asserted.

AoA 16 is the case where this mattered. At the t=60.5 verdict its `D_P2P` was 0.00167,
two orders looser than AoA 18's 1e-5 and the only metric anywhere near tolerance; on
the full phase-2 record it fell to 0.00003. It really was still creeping at the gate,
and it really had stopped by the end. Close to onset everything slows down, so this is
the region where trusting the gate alone would eventually cost a wrong answer.

### The method

Near a supercritical Hopf the saturated amplitude obeys the normal form
`A^2 = k (alpha - alpha_c)`. **Every period-1 case in the campaign is a saturated
amplitude**, so the onset is already latent in runs paid for long ago; it only has to
be read out.

Use the spectral carrier amplitude, not `C_L` peak-to-peak: p2p folds in every harmonic
and the harmonic content is not constant across the sweep. Exclude the period-2 cases
-- the scaling describes a single Hopf mode and means nothing once the orbit doubles.

| fit | `alpha_c` | R2 |
|---|---|---|
| 18, 20, 22 -- **the fit that made the prediction** | 14.63 | 0.99657 |
| 16, 18, 20 (nearest onset) | 14.74 | 0.99884 |
| 16, 18, 20, 22 | 14.65 | 0.99863 |
| all five (16-24) | 14.83 | 0.99631 |

Spread 14.63-14.83, mean **14.71 deg**. `analysis/plot_onset.py` draws it;
`analysis/figs/onset_amplitude.png`.

### The out-of-sample test, which is the part that matters

Fitted on 18/20/22 alone, the scaling predicted `A(16) = 0.0701` for an angle it had
never seen. Measured: **0.0694, an error of -0.9%.** C_L peak-to-peak as an independent
basis predicted 0.1330 against a measured 0.1390.

A fit that merely looks tidy cannot do that. R2=0.997 on four points is weak evidence
on its own -- almost anything is straight over a narrow range -- but the amplitudes span
a factor of 2.7 (0.069 to 0.189) and the extrapolation held at an angle outside the
fitted set. That is what makes ~14.8 deg a measurement rather than a curve-fitting
artefact.

Two independent corroborations, neither used in the fit: per-chunk `Co_max` at AoA 16
fell 0.4129 -> 0.3320 -> 0.3151 as the oscillation relaxed to its weaker saturated
level; and `f_carrier` rises smoothly and monotonically toward the onset across the
whole period-1 range (0.4182 at 24 -> 0.5854 at 16) with no kink.

### What is still open

The scaling **assumes supercriticality**, and its own linearity in `A^2` is evidence
for that -- a subcritical branch turns back on itself and could not produce R2=0.9988
across five angles -- but it is not proof. The direct test is a genuinely steady case
below onset, which needs an angle at or under ~14 deg. Two related caveats:

- Every case in this sweep was warm-started from a shedding state, so a shedding
  result never proves the steady state is unstable, only that a shedding attractor
  exists there. The asymmetry in §16 of `growth_rate.py`'s docstring applies: a
  DECAYING case measures a linear eigenvalue, a saturating one does not.
- This is the same supercritical-vs-subcritical question left open for the
  period-doubling near 25 deg (§11). One steady case at 14 deg would settle the
  shedding half of it.

Cheapest decisive test: AoA 14, warm-started from AoA 16 exactly as here. If the
oscillation decays to a steady state, `growth_rate.py` returns a genuine linear
`sigma < 0`, which both confirms supercriticality and gives an independent estimate of
`alpha_c` from the RATE rather than the amplitude -- two different observables agreeing
on one number.


## 17. AoA 14 sheds -- the onset extrapolation was wrong, and how it fooled me

**AoA 14 is a saturated period-1 limit cycle**, not a decaying transient:
`C_L = 0.5513 +- 0.06835`, `A(f_carrier) = 0.034126`, `f = 0.608322`,
`C_D = 0.26190`. The envelope spans **0.034122-0.034126 over t=48-67, a spread of
0.010%**, slope `+6.5e-6`, below the 1e-5 saturation floor.

That falsifies §16. The `A^2` fit over 16-24 deg predicted `alpha_c = 14.83` and
therefore **no shedding at 14 deg**. There is shedding at 14 deg.

### The error, and the two chances I had to catch it

I read the record at t~41 -- mid-relaxation, where the local decay rate genuinely was
varying 6x -- and concluded "not converged". When the gate then passed the case at
t=60.5 I overrode it and killed the run 12 snapshots into phase 2. The gate was right.

What I was watching was the warm-start relaxation from AoA 16's amplitude (0.0694)
down to AoA 14's own saturated amplitude (0.0341). It finished around t~45. The
local rate tells the story unambiguously once the record is long enough:

| t | A | d(ln A)/dt |
|---|---|---|
| 34.65 | 0.035238 | -0.04105 |
| 38.40 | 0.034292 | -0.01179 |
| 42.15 | 0.034075 | -0.00194 |
| 47.14 | 0.034020 | -0.00030 |
| 52.14 | 0.034019 | -0.00000 |
| 55.89 | 0.034022 | +0.00002 |

**A decay rate that itself decays to zero is an approach to an asymptote, not a decay
to zero amplitude.** The first `sigma` I quoted, `-0.01266 at R2=0.99876`, was fitted
over **t=11.67-18.06 -- 6.4 time units**, about 3.8 carrier periods, containing 2-3
independent envelope points. R2 on 3 points is not evidence. *Check the span before
the R2.*

### Why the extrapolation failed, diagnosed

Not a stray point -- a wrong functional form. Fitted `alpha_c` rises monotonically
with the highest angle included:

| subset | alpha_c | R2 |
|---|---|---|
| 14, 16 | 13.36 | -- |
| 14, 16, 18 | 13.77 | 0.96834 |
| 14-20 | 14.00 | 0.97635 |
| 14-24 (all six) | 14.26 | 0.98717 |
| 16, 18, 20 | 14.74 | 0.99884 |
| 16-24 | 14.83 | 0.99631 |

A correct model would give the same answer from every subset. This is lever-arm bias.

Freeing the exponent, `A = C (alpha - alpha_c)^p`, fits all six at **R2 = 0.9966**
with **p = 0.786** and `alpha_c = 12.76` -- against the normal form's `p = 0.5`.
Imposing `p = 0.5` gives a best rms of 0.0092 (at `alpha_c = 13.5`) versus 0.0027
free, so the square-root law is **disfavoured**, not merely unconstrained.

But `alpha_c` and `p` trade off smoothly, so six points over 14-24 deg fix neither:

| alpha_c | p | rms in A |
|---|---|---|
| 10.0 | 1.359 | 0.01033 |
| 12.0 | 0.962 | 0.00579 |
| 13.0 | 0.722 | **0.00270** |
| 13.5 | 0.566 | 0.00606 |

### Current position

**Onset is below 14 deg, loosely 12.5-13.5, and the amplitude method cannot pin it
from these angles.** Consistent with upstream's 8-10 deg at Re=1000 (lower Re is more
stable), though that bound is too weak to discriminate.

The §16 out-of-sample success at AoA 16 (-0.9%) was weaker evidence than it looked:
16 deg sits inside the well-behaved range, so it tested interpolation dressed as
extrapolation. And §4c of `BIFURCATION_THEORY.md` had already recorded that
`f_carrier` is *not* linear in alpha -- i.e. the normal form's linear regime had been
left -- which I wrote down as a caveat and then kept quoting the amplitude fit anyway.

### What settles it

1. **AoA 12**, warm-started from 14. Sheds -> onset below 12 and the exponent gains
   its most constraining point. Decays to steady -> onset bracketed in (12, 14)
   directly, plus the first genuine eigenvalue, since that far below onset `sigma` is
   well separated from the readjustment timescale.
2. **Base-flow perturbation runs.** `steadyState` converges to the unstable fixed
   point (the §1 trap in `BIFURCATION_THEORY.md`, used deliberately); seed the
   unsteady solver with `x_s + eps(x_shed - x_s)`. Starting AT the fixed point makes
   the mean-flow displacement zero by construction -- and that displacement, measured
   at **1.9x the oscillation amplitude** during AoA 14's transient, is exactly what
   makes a warm-started `sigma` uninterpretable.

### Side effect

AoA 14 is a valid converged period-1 case, so its snapshots join the POD set. The
usable sweep is now **seven angles, 14-30**.

## 18. Steady base flows, 0-24 deg -- and why the weakly-nonlinear framework does not apply

13 SIMPLE solves at 2-degree spacing, 6000 iterations each, 6 MPI ranks, 48 min
total. Same case as the transient runs in every respect **except that the time
derivative is deleted** -- `gradSchemes`, `divSchemes`, `laplacianSchemes` all
unchanged. Upstream's steady cases used `div(phi,U) Gauss upwind`; copying that
would have put a first-vs-second-order discretisation difference straight into
`C_L,mean - C_L,steady`, which is the quantity being measured.

Validation: **`C_L = 0.000000` at AoA 0**, exactly, as a symmetric section must
give -- confirms mesh symmetry and the `dragDir`/`liftDir` wiring. And at AoA 10
this gives `C_L = 0.426072` against upstream's `0.458842` at Re=1000, 7.1% lower,
the right direction and size for thicker boundary layers at lower Re.

### SIMPLE converges up to 20 deg and limit-cycles above it

| AoA | `C_L,steady` | `C_D,steady` | C_L range, last 10% | verdict |
|---|---|---|---|---|
| 0 | 0.000000 | 0.174132 | 2.2e-10 | converged |
| 2 | 0.112979 | 0.175460 | 0 | converged |
| 4 | 0.218670 | 0.179485 | 0 | converged |
| 6 | 0.309705 | 0.186303 | 0 | converged |
| 8 | 0.379699 | 0.195962 | 0 | converged |
| 10 | 0.426072 | 0.208307 | 0 | converged |
| 12 | 0.452294 | 0.222876 | 0 | converged |
| 14 | 0.465639 | 0.239060 | 0 | converged |
| 16 | 0.472710 | 0.256362 | 0 | converged |
| 18 | 0.477540 | 0.274533 | 0 | converged |
| 20 | 0.482076 | 0.293524 | 0 | converged |
| 22 | 0.494523 | 0.316797 | **6.046e-03** | LIMIT-CYCLE |
| 24 | 0.505358 | 0.340897 | **1.258e-02** | LIMIT-CYCLE |

**Convergence must be judged on the force trace, not the residuals.** At AoA 24
the C_L range over the last 50%, 30%, 20% and 10% of iterations is 1.258e-02 in
*every* window. An identical range across shrinking windows is a sustained
oscillation in iteration space -- the physical instability leaking into the
SIMPLE iteration -- not slow convergence. A converging trace's range must fall as
the window shrinks. `analysis/steady_summary.py` tests exactly that.

This is a numerical property of the iteration, **not** a stability criterion. It
must not be read as "the base flow becomes unstable above 20 deg"; that would be
the §1 circularity of `BIFURCATION_THEORY.md` run backwards.

### The result that matters: the limit cycle does not orbit the fixed point

| AoA | `C_L,steady` | limit-cycle range | gap below the cycle | in p2p |
|---|---|---|---|---|
| 14 | 0.465639 | [0.51749, 0.58571] | +0.0519 | 0.76 |
| 16 | 0.472710 | [0.56080, 0.69982] | +0.0881 | 0.63 |
| 18 | 0.477540 | [0.60269, 0.82034] | +0.1251 | 0.57 |
| 20 | 0.482076 | [0.64605, 0.92735] | +0.1640 | 0.58 |
| 22 | 0.494523 | [0.67577, 1.00403] | +0.1812 | 0.55 |
| 24 | 0.505358 | [0.68065, 1.07735] | +0.1753 | 0.44 |

**At every angle the steady lift lies BELOW the minimum instantaneous unsteady
lift**, by 0.44-0.76 of a full peak-to-peak. The saturated flow is not an
oscillation *about* the base flow; it is displaced bodily away from it.

The mean-flow distortion measures this directly. Stuart-Landau requires it to be
`O(R^2)` -- forced by symmetry, since the linear term averages to zero over a
cycle. Measured:

| AoA | shift | `shift/A^2` | `shift/A` |
|---|---|---|---|
| 14 | 0.085962 | 74.14 | 2.525 |
| 16 | 0.157599 | 32.68 | 2.269 |
| 18 | 0.233973 | 19.80 | 2.152 |
| 20 | 0.304624 | 15.47 | 2.171 |
| 22 | 0.345377 | 13.23 | 2.137 |
| 24 | 0.373642 | 10.48 | 1.979 |

`shift/A^2` varies by 7x. `shift/A` is flat at ~2.2. Global fit
**`shift ~ A^0.879`, R2 = 0.99823** -- linear in the amplitude, off by a whole
power from the theory.

### Consequence

Three independent observables now reject the weakly-nonlinear description over
14-24 deg:

1. amplitude exponent `p = 0.79`, not 0.5 (§17);
2. `f_carrier` not linear in alpha (`BIFURCATION_THEORY.md` §4c);
3. mean-flow distortion `~ A^0.88`, not `A^2` (here).

All three follow from the one geometric fact above. **Extrapolating `alpha_c`
from any amplitude scaling over this range is dead** -- not imprecise; the
functional forms the extrapolation rests on do not hold, so no refit and no extra
angle rescues it.

The eigenvalue route is unaffected: `sigma` is a property of the linearisation
about the fixed point and does not care that the attractor is far away.

### Physical result, independent of any bifurcation question

The steady base flow is **stalled from about 12 deg on** -- `C_L,steady` runs
0.4523 (12) to 0.4821 (20), +6.6% over 8 degrees, against a near-zero-incidence
slope of 0.056/deg which is itself only 51% of the inviscid `2*pi` (heavy viscous
decambering at Re=500). Over the same 12-20 deg the time-accurate mean climbs
0.55 to 0.79.

**Essentially all the lift rise with incidence above ~12 deg comes from the
vortex shedding**, and the shed-vortex state out-lifts the steady state at every
instant of every cycle. That inverts the upstream Re=1000 picture at alpha=10,
where the steady solve *over*-predicted lift by 9.8%.

Figure: `analysis/figs/steady_vs_mean.png` (`analysis/plot_steady.py`).

## 19. The bifurcation point, measured: alpha_c = 11.46 deg

Three linear growth rates, measured directly about the steady base flow -- not
extrapolated from any property of the saturated state, which is what failed over
14-24 deg (sections 17-18).

| AoA | `sigma` | R2 | e-fold | verdict | how |
|---|---|---|---|---|---|
| 11.0 | **-0.07093** | 0.99999 | 14.1 t | STABLE | decay over 2.2 decades |
| 11.5 | **+0.00674** | 1.00000 | 148.4 t | UNSTABLE | growth, 63 t span, 118% |
| 12.0 | **+0.10970** | -- | 9.1 t | UNSTABLE | growth, early linear window |

**Bracketed: stable at 11.0, unstable at 11.5.** Zero crossing at
**alpha_c = 11.46 deg** (11.42 from a linear fit to all three; 11.39 from the
11.0/12.0 bracket alone, which was the prediction made BEFORE the 11.5 run).
`dsigma/dalpha = 0.181` per degree.

All three sit on the same mode, **f = 0.600** -- the same frequency decaying at
11.0 and growing at 11.5 and 12.0. A complex pair, so the LINEAR onset is a Hopf
even though the saturated states at 14-24 deg obey none of the weakly-nonlinear
Hopf scalings. Those are consistent: a Hopf can be the linear mechanism while the
attractor sits far from the fixed point (section 18).

AoA 11.0 came from a collaborator's machine. Verified before use: directions and
`Aref` correct for 11 deg, mesh identical to M2 (non-orthogonality 64.0766,
skewness 0.821208), and its converged `C_L = 0.441572` lands **+0.12%** off our
own steady lift curve interpolated to 11 deg. Its dt (7.13e-4) differs from ours
(5.5726e-4), so the two sigmas carry slightly different numerical damping -- the
main reason to quote **11.46 +- 0.1** rather than three digits.

### Critical slowing down, measured

`sigma` falls by **16x between 12.0 and 11.5 deg** (e-fold 9.1 t -> 148.4 t)
across half a degree. That is why no amplitude-based method could ever have
found this point: the saturated amplitude collapses toward zero exactly where
the timescale diverges.

### What made the 11.5 measurement work

Three choices, each fixing a failure recorded earlier:

1. **Start AT the fixed point.** Removes the mean-flow displacement that made
   AoA 14's warm-started sigma uninterpretable (section 17): there the
   displacement was 1.9x the oscillation amplitude.
2. **Blend the seed, eps=0.1.** Near onset the attractor is tiny; the standard
   `x_s(a+0.5)` seed would have been dCl~0.0053 against an attractor of
   0.0018-0.0030 -- 2-3x TOO BIG, so the perturbation would have decayed ONTO the
   cycle and read as sigma<0. Blending put the oscillatory seed at 1.4e-5, under
   1% of the attractor, and the whole 115 t record stayed linear. No saturating
   tail to discard, unlike AoA 12 where only t=10-16 was usable.
3. **Track ONE known frequency.** See below.

### Two self-inflicted failures worth recording

**A malformed dictionary diverged the run twice, and the guard passed.** Stripping
`relaxationFactors` with `re.sub(r"relaxationFactors\s*\{[^}]*\}")` deletes the
opening brace and the FIRST closing one -- which belongs to the nested
`equations { U 0.7; }` -- leaving `fields { p 0.3; }` and an orphan `}` at top
level. `grep -q relaxationFactors` then passes, because the word is gone. **A
guard that checks for the absence of the symptom you just removed proves nothing
about what you left behind.** Fixed by copying the known-good transient
`fvSolution` wholesale (it has no angle dependence) and validating with
`foamDictionary`, which fails on any malformation, not just the anticipated one.

**`Y.max()` hijacked the fit -- the peak-height rule again.** `eigenvalue.py`
tracked the tallest spectral peak. A slowly-decaying non-oscillatory seed mode
made one large excursion around t=25-35; its content at f=0.2 was **2.47e-4**
against the eigenmode's **1.9e-5**, so between consecutive windows the tracker
jumped from the transient to the eigenmode and the fit read **sigma = -0.10**
where the mode was at **+0.0067** -- the wrong SIGN. Adding an `f_min` floor did
not fix it; the excursion reached exactly 0.2. Only `--freq 0.600011`, probing
one known frequency, is correct.

`CLAUDE.md` already said **never identify a frequency by peak height**, in bold,
after three earlier bugs. It was written into a brand-new tool anyway. The rule
is not "prefer the right peak" -- it is *know which frequency you want and probe
there*.

Figure: `analysis/figs/eigenvalue_onset.png` (`analysis/plot_eigenvalue.py`).
