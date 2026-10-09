# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

2D laminar CFD campaign on a static NACA0012 (OpenFOAM-dev), producing field/force
data for a POD-based reduced-order-model thesis. Scope and matrix live in
`experiment_plan.md`; the wider campaign (pitching, step-change, actuated) is in
`../plan/simulation_plan.md`.

**Read the prior campaign's findings before changing solver settings.** The August
Re=1000 work lives at `~/OpenFOAM/raji-dev/run/Aug/unsteadyRe1000/` and its
`FINDINGS.md` is the authoritative record: shedding onset, a silent solver pathology,
and a mesh-independence study. `NACA0012/aoa10/` in this repo is a hand-edited *derivative* of
`unsteadyRe1000/AOA_10fine/` and disagrees with it (`nOuterCorrectors 2` vs 1,
`numberOfSubdomains 4` vs 6, `writeInterval 0.025` vs 5) -- treat the upstream case as
the source of truth, not the copy here.

## Environment

- **OpenFOAM**: `source /opt/openfoam-dev/etc/bashrc` (OpenFOAM.org dev, build `dev-84333c37c7aa`).
  Solver is `foamRun` with `solver incompressibleFluid;` in `controlDict`.
- **gmsh Python API lives only in the `research` conda env.** System `python3` has no
  `gmsh` module. Use `/home/raji/miniconda3/envs/research/bin/python` for anything
  under `meshing/`; plain `python3` is fine for `analysis/`.
- **This box has 6 PHYSICAL cores, not 12.** i5-12400, 6 cores / 12 threads; `nproc`
  reports hyperthreads. `mpirun -np 8` fails with "not enough slots". Two concurrent
  4-rank jobs oversubscribe the 6 physical cores and cost ~1.7x per step (measured:
  0.059 -> 0.102 s/step). Budget **6 MPI ranks total** across all running jobs.
- 23 GB RAM.

## Commands

```bash
# mesh  (h1 = 0.05/sqrt(Re) -> a mesh is Re-SPECIFIC)
cd meshing && /home/raji/miniconda3/envs/research/bin/python mesher.py OUT.msh \
    --re 500 -r 30 --nr 250 --ns 200 --nte 40 --bump 0.08

./make_case.sh <caseDir> <Re> <AoA_deg> <mesh.msh>   # full case + gmshToFoam + checkMesh
./run_staged.sh <caseDir> [nCores]                   # two-stage run, default 4
./warm_start.sh <src> <srcTime|latestTime> <dst> [np]  # seed dst from a converged src
./bench_mesh.sh <mesh.msh> <label>                   # throwaway case -> Co + s/step + ETA

./run_onset.sh 18 16      # shedding-onset sweep: every angle seeded from AoA 20,
                          # dt frozen, so sigma is comparable across angles

python3 analysis/limit_cycle.py <caseDir> [--tol 0.01]   # KEY=VALUE, eval-able from sh
python3 analysis/growth_rate.py <caseDir> [--tmin T]     # C_L envelope sigma; SIGMA<0 => below onset
./make_steady.sh <caseDir> <Re> <AoA> <mesh.msh> [nIter]  # SIMPLE base flow
./run_steady.sh 24 22 20 ... 0    # sequential, 6 MPI ranks each (NOT 6 serial
                                  # jobs -- that was 4.2x slower per case)
                                  # cases live in NACA0012/steady_re500/ (default) or
                                  # STEADY_DIR=NACA0012/steady_re1000 ./run_steady.sh ...
python3 analysis/steady_summary.py [--root NACA0012/steady_re1000]  # Cl_steady + convergence,
                                         # judged on the FORCE trace, never on residuals
/home/raji/miniconda3/envs/research/bin/python analysis/lasso_re1000.py
                                         # Re500 vs Re1000 Lasso -> STEADY_LASSO_Re1000.md
python3 analysis/plot_steady.py          # lift curve + mean-flow distortion
                                         # -> analysis/figs/steady_vs_mean.png

./run_marginal.sh 11.5 12.0 0.1 120.5    # eigenvalue run near onset (blended seed)
python3 analysis/blend_fields.py <A> <B> <eps> <out>     # out = A + eps*(B-A)
python3 analysis/eigenvalue.py <case> --cl-steady <Cl> --freq <f>  # ALWAYS pass --freq
python3 analysis/plot_eigenvalue.py                      # sigma(alpha) -> onset
                                                         # -> analysis/figs/eigenvalue_onset.png

python3 analysis/campaign_data.py         # one table: Cl/Cd mean, min, max, order,
                                          # f_carrier -- every figure reads from here
python3 analysis/plot_psd.py              # C_L/C_D amplitude spectra, one figure per
                                          # angle -> analysis/figs/PSD/psd_aoa<NN>.png
python3 analysis/plot_coeffs_vs_aoa.py    # Cl, Cd, Cl/Cd vs AoA with the limit-cycle
                                          # min/max and both bifurcations marked
                                          # -> analysis/figs/coeffs_vs_aoa.png (+ panels)

python3 analysis/plot_onset.py                           # A^2 vs AoA (FALSIFIED, kept as a record)
                                                         # -> analysis/figs/onset_amplitude.png
python3 analysis/plot_cl.py                              # C_L histories + period-2 saturation
                                                         # -> analysis/figs/*.png

# post-processing: plain `postProcess` CANNOT construct the momentum transport model
foamPostProcess -solver incompressibleFluid -time 90:100
```

## Architecture

### Folder layout: cases are grouped by section (since 2026-10-09)

Scripts, meshes, `analysis/`, `lsb/` and the write-ups stay at the root. **Case folders live under the
section they belong to**, and every driver and loader resolves its paths there:

- `NACA0012/`: `steady_re{100,500,1000,2000}/`, `steady_re_scaling/steady_re<RE>/`, the unsteady
  Re=500 campaign `re500_aoa{12..30}/`, `aoa10/`, `unsteady_re2000/`, `bi_point_search/` and
  `mesh_convergence_M1/`.
- `NACA0006/`: `steady_re500/` (formerly `steady_naca0006_re500/`), `steady_re500_fine/`, and
  `steady_re_scaling/steady_re<RE>/`.
- `collect(re)` (`analysis/lasso_re1000.py`) reads NACA0012. `collect(re, "NACA0006/steady_re500")`
  reads any other folder, relative to this directory.
- `run_steady.sh` defaults to `STEADY_DIR=NACA0012/steady_re500`.
- New drivers must create cases under the section folder too, e.g. `NACA0012/re500_aoa$A`.
- **Older text uses the old paths.** `FINDINGS.md`, `INCIDENT_*.md`, run logs and `claude_md/`
  predate the move and keep the old root-level paths: prefix `NACA0012/` when reading them.
- The plate lives in `../plate_case/re<RE>/`.

### `meshing/` — layered so the arithmetic is testable without gmsh

`airfoil.py` (pure numpy) → `mesh_math.py` (pure arithmetic) → `mesher.py` (gmsh).
Neither of the first two imports gmsh, so geometry and spacing can be diffed and
plotted standalone. Contracts `mesher.py` depends on, documented in `airfoil.py`:
`upper`/`lower` both run LE→TE and **share an identical first row** — the LE point
must be created once and its tag reused, or the splines never close the curve loops.

`mesher.py.__main__` is the production path: a 3-block transfinite O-ring
(upper/lower/wake) extruded 0.1 in z, tagged `internal`/`airfoil`/`farfield`/`frontAndBack`.
`wake_half_angle` picks `theta_w = pi*nte/(2*ns+nte)` so the far-field circle is
discretised uniformly; any other value leaves a visible kink at the cut lines.

**`outer_domain.py` is not used by the production path.** It is a separate,
self-contained alternative (rectangular outer domain + wake refinement boxes wrapped
around the ring) with its own `__main__`. Don't assume it participates.

### AoA is imposed by rotating the freestream, not the mesh

This is the single most important thing to know. One mesh serves every angle at a
given Re. The angle appears in **three** places and `make_case.sh` writes all of them
from one computation: `0/U` `internalField` + `freestreamValue`, `functions`
`dragDir`/`liftDir`, and `Cp`'s `UInf`. Editing one by hand yields a plausible but
wrong `C_L`. It is also what makes `warm_start.sh` exact.

Reference values follow from c=1, U=1, span=0.1: `lRef 1`, **`Aref 0.1`** (c*span, not
1.0), `magUInf 1`, `rhoInf 1` (pressure is kinematic).

`gmshToFoam` resets every patch to type `patch`. `make_case.sh` retypes
`frontAndBack`→`empty` (mandatory — otherwise it is a 3D run one cell thick) and
`airfoil`→`wall`.

### The staging contract — spans `run_staged.sh` and `system/functions`

`system/functions` is **byte-identical in both stages**; field output is controlled
entirely by `controlDict`. Do not "fix" this by maintaining two functions files.

- `forceCoeffs` is on `writeControl timeStep`. It appends to a single
  `postProcessing/forceCoeffs/<startTime>/*.dat` and **never creates a time
  directory**, so per-step `C_L`/`C_D` through the transient costs nothing.
- `Cp`, `wallShearStress`, `surfaces` are on `executeControl writeTime` — they fire
  only when the case writes fields. A field-writing FO on `timeStep` mints a new time
  directory every `deltaT` containing only that field and no `U`/`p`.

Stages:
1. **Phase 0** — impulsive start, `dt=4e-5` to `t=0.5`. **Measure Co at t=0.5, never
   earlier.** Co peaks ~14x higher around `t=1e-3`; sampling inside that spike gives a
   dt roughly 7x too small and an ETA ~7x too long.
2. **Phase 1** — transient, field-free, in chunks. `writeInterval == endTime` gives
   exactly one write per chunk, which is the restart point. (There is no `onEnd`
   writeControl in this build; valid values are `timeStep, writeTime, outputTime,
   adjustableRunTime, runTime, runTimes, clockTime, cpuTime, none`.) Loop exits when
   `analysis/limit_cycle.py` reports `CONVERGED=yes`.
3. **Phase 2** — fields at ≥32 samples per **measured** carrier period, for
   `NCYCLES` **orbits**. Those are two different periods once the wake has
   period-doubled; see below.

**The Courant number is not a stability criterion on this mesh — `DT_MAX` is.**
Mean Co runs ~166x below max Co, because the max is set by a few trailing-edge
cells and reports on nothing else. AoA 18 diverged at **Co = 0.60**, below both
`RETUNE_CO` (0.75) and the `CO_HI` abort (0.80), after the re-tune grew dt onto
the old `DT_MAX=1.0e-3` clamp; the real limit is the lagged non-orthogonal
pressure correction on 64-degree cells, whose weight grows with dt while Co
barely moves. `DT_MAX` is now **6.0e-4**, below the largest dt ever run to
convergence here (6.331e-4, AoA 26). A diverged chunk is retried at dt/2 (twice,
then it gives up) and its `postProcessing/forceCoeffs/<restart>/` is quarantined
to `.diverged/` **before** the retry -- a retry restarts at the same time and so
appends to that same directory, which would otherwise leave 1e204-magnitude
samples in the series forever. See `FINDINGS.md` section 15.

**`DT_FIXED=<dt> ./run_staged.sh ...`** holds dt and disables re-tuning. Use it
for anything that measures a *rate* rather than a converged state: the scheme's
numerical damping is dt-dependent, so a re-tune mid-record puts a step into the
envelope that `analysis/growth_rate.py` is fitting, and cross-angle sigma
comparisons are only meaningful at one shared dt. `run_onset.sh` sets it.

`adjustTimeStep` is `no` everywhere — the variable-dt controller corrupted an earlier
attempt. dt changes only *between* chunks, so it is rigorously fixed within any single
`foamRun`. `writeInterval` is snapped **down** to an integer multiple of dt, which
lands writes exactly on time steps without needing `adjustableRunTime` (that option
perturbs dt and would reintroduce the problem). Each dt change costs one Euler step,
since `backward` needs two old time levels and a restart has one.

### Two periods, not one — the wake period-doubles

Past some angle in (20°, 30°] at Re=500 the wake is no longer period-1: successive
shedding cycles alternate and the orbit closes only after **two** oscillations of
`C_L`. Measured at AoA 30: carrier f=0.3196, saturated subharmonic at f/2=0.1598
carrying **63%** of the carrier's amplitude, all eight resolved peaks exact integer
multiples of f/2 (max deviation 0.0000 — a closed period-2 orbit, not a torus).
AoA 20 is period-1, with f/2 power at the 1e-5 noise floor.

**The carrier is not always the tallest peak.** At AoA 30 the carrier dominates
(subharmonic 0.67 of it); at AoA 28 the subharmonic dominates (carrier 0.67 of it) --
same structure, flipped. `A(2f)/A(f)` cannot tell them apart (0.270 at period-2 AoA 30
vs 0.195 at period-1 AoA 24), so `spectral_structure()` finds the two independently:
`f_carrier` is the highest peak still carrying 50% of the maximum amplitude, `f_fund`
the lowest peak above 5%, and `ORDER` their ratio, rejected unless within 0.15 of an
integer. `SUBHARM_RATIO` is `A(f_fund)/A(f_carrier)` so it stays comparable across
angles. See `FINDINGS.md` sections 12 and 13.

**Never identify a frequency by peak height.** Three separate bugs came from it, the
worst at AoA 26 where `SUBHARM_RATIO ~ 1` makes the tallest peak flip between windows
on noise -- `D_FREQ` then compared a fundamental against a carrier, sat pinned at 0.496,
and the run could never have converged. Anywhere the code needs to know *which*
frequency it has, it must call `spectral_structure()`.

Taking the largest spectral peak as "the" period gets **both** downstream decisions
wrong, so `analysis/limit_cycle.py` reports two:

| key | meaning | drives |
|---|---|---|
| `SAMPLE_PERIOD` | the **carrier** — fastest oscillation present | `writeInterval` |
| `PERIOD` | the **fundamental** = `ORDER` × carrier — when the orbit closes | phase-2 `endTime`, convergence window |

Using the carrier for run length records `NCYCLES/ORDER` orbits — AoA 30's "20 cycles"
was 9.97. Using a window holding a non-integer number of fundamental periods makes the
cycle-to-cycle alternation look like drift, which is what made a fully saturated AoA 30
report `CONVERGED=no`. The window is now cut to a whole number of fundamental periods
and split at the exact midpoint. `ORDER=1` recovers the old behaviour exactly, and
every period-1 case was regression-checked as unchanged.

The convergence gate tests `FULL_WINDOW` (was the window cut to length, or clamped by
the start of the series?), **not** a recount of cycles inside it. The window is
`2*n_per` fundamental periods by construction, so recounting it with a second, slightly
different frequency estimate lands a hair under the requirement about half the time and
blocks an otherwise-converged case. Measured on AoA 22: 5.99 cycles against 6.0 with
every agreement test two orders of magnitude inside tolerance.

Also reported: `SUBHARM_RATIO`, `HARMONIC_DEV` + `QUASIPERIODIC` (period-m vs
Neimark–Sacker torus — a torus has incommensurate peaks), and `D_SUB`, the
window-to-window agreement of the subharmonic. `D_SUB` matters because the subharmonic
saturates **later** than the carrier: mean, amplitude and frequency can all agree while
the orbit is still splitting.

**Unresolved:** the onset angle is bracketed only to (20°, 30°]; the bifurcation has
not been shown mesh- or dt-independent; and AoA 30 was warm-started, so bistability is
not excluded. `FINDINGS.md` §9 has the numbers and the cheapest tests.

`reconstructPar` runs only for `t >= t_conv`: phase-1 chunk writes are restart points,
not snapshots, and must not enter the POD set. `t_conv` itself is one of those writes
and sits **off** the phase-2 grid (0.0898 before the first sample, not 0.0972) -- drop
it too.

**Counting time directories:** OpenFOAM names them from the value, so a write on a
whole number is `213`, not `213.0`, and `[0-9]*.[0-9]*` silently skips it. Use the
`count_snaps` helper in `run_staged.sh`; a loader that globs on the decimal point will
drop snapshots and mis-space its time axis.

`extend_phase2.sh <case> <t_conv> <nOrbits> [nCores]` lengthens a finished phase 2 in
place, holding `dt` and `writeInterval` fixed so the new snapshots continue the same
uniform grid.

### Restart caveat

On restart `forceCoeffs` opens a **new** `postProcessing/forceCoeffs/<restartTime>/`
directory, so the series is split across stage directories with duplicated times at the
seams. `analysis/limit_cycle.py` concatenates, sorts by time, and keeps the last sample
at each time. Any new force analysis must do the same.

### Warm start

`mapFields -consistent` is a direct 1:1 copy here, not an interpolation, because the
meshes are bit-identical across AoA. `warm_start.sh` enforces that with an md5 of
`points/faces/owner/neighbour` rather than trusting it — across different Re the
meshes differ (`h1` is Re-dependent) and `-consistent` would be silently invalid.

- `mapFields` is **not** an MPI app. Under `mpirun` it fails with
  `Invalid option: -parallel`. It runs serially and reads decomposed cases itself via
  `-parallelSource`/`-parallelTarget`.
- It maps **volume fields only**, so `phi` is not transferred; the solver rebuilds it
  from the mapped `U` and the first pressure solve absorbs the divergence error.
- Keep the target's own `0/U`/`0/p`. `mapFields` overwrites values but preserves BC
  definitions, so `freestreamValue` stays at the target angle. Copying the source time
  directory instead would drag the source's far-field direction along with it.
- **Do not rotate the velocity vectors.** The airfoil occupies the same cells in both
  cases; the boundary layer and near wake are already correct, and the far field
  self-corrects by advection within a few convective times.

## Settings that were measured, not guessed

Re=500, 110k-cell mesh, 4 cores. Re-measure before changing any of these.

| change | result |
|---|---|
| `nOuterCorrectors 1` (PISO) vs `2` | **0.176 vs 0.416 s/step** — 2.4x, no accuracy cost at Co<1 with `backward` |
| `nNonOrthogonalCorrectors 0` vs `1` | 0.650 vs 0.176 s/step — **slower**; max non-orthogonality is ~64°, keep it at 1 |
| GAMG `nCellsInCoarsestLevel`/`pFinal 1e-7` | no measurable gain |
| `bump` (LE clustering) 0.08→0.25 | no effect on Co; the LE is *not* the Courant limit |
| `nte` (TE cells) | **this** sets the Courant limit — TE cells are the smallest in the mesh |

Settled Courant at `dt=4e-5`, t=0.5: **0.0729** (110k, ns200/nte40/nr250) and **0.0282**
(49.5k, ns150/nte30/nr150).

See `FINDINGS.md` for the full measurement log and the reasoning behind each choice,
`ANALYSIS_GLOSSARY.md` defines the diagnostic vocabulary (carrier vs fundamental,
subharmonic ratio, the torus guard) for anyone reading the verdict lines,
`BIFURCATION_THEORY.md` is the companion on the dynamical-systems side -- Hopf,
Stuart-Landau, Floquet multipliers, where `A^2 = k(alpha - alpha_c)` comes from, and a
table of which measurements the theory licenses and which it does not,
and `INCIDENT_2026-09-05_power_outage.md` covers interruption recovery --
including why phase 1 writes one field per chunk and why that write is not in the POD
set.

### The shedding onset is alpha_c = 11.46 deg, measured from eigenvalues

**Do not quote 14.8, 13.0, or "below 14".** Those all came from extrapolating
properties of the SATURATED state and were falsified in turn. The onset is now
measured directly, by perturbing the steady base flow and fitting the growth rate
of the f=0.600 mode:

| AoA | `sigma` | R2 | verdict |
|---|---|---|---|
| 11.0 | **-0.07093** | 0.99999 | STABLE (decay over 2.2 decades) |
| 11.5 | **+0.00674** | 1.00000 | UNSTABLE (growth, 63 t span) |
| 12.0 | **+0.10970** | -- | UNSTABLE |

Bracketed and crossing at **alpha_c = 11.46 +- 0.1 deg**, `dsigma/dalpha = 0.181`
per degree. All three are the same mode, f=0.600 -- a complex pair, so the LINEAR
onset is a Hopf even though the saturated states obey none of the weakly-nonlinear
Hopf scalings. See `FINDINGS.md` section 19; figure
`analysis/figs/eigenvalue_onset.png`.

**Amplitude-based extrapolation cannot find this point, in principle.** `sigma`
falls 16x between 12.0 and 11.5 deg (e-fold 9.1 t -> 148 t): the saturated
amplitude collapses exactly where the timescale diverges. Three observables
reject the weakly-nonlinear description over 14-24 deg -- amplitude exponent 0.79
not 0.5, `f_carrier` not linear in alpha, and mean-flow distortion `~A^0.88` not
`A^2` -- and all three follow from one geometric fact: **at every angle measured,
the steady C_L lies BELOW the minimum instantaneous unsteady C_L**, by 0.44-0.94
of a peak-to-peak. The limit cycle does not orbit the fixed point.

### Measuring an eigenvalue: what the campaign learned the hard way

`./run_marginal.sh <aoa> <srcAoA> <eps> [TMAX]` implements all of this.

1. **Start AT the fixed point** (steady solve), never from a limit cycle. A
   warm start from a saturated state left AoA 14 with a mean-flow displacement
   **1.9x** the oscillation amplitude and no interpretable sigma.
2. **Size the seed against the ATTRACTOR, not against habit.** Near onset the
   attractor collapses: at 11.5 deg the standard `x_s(a+0.5)` seed is 2-3x too
   big, so the perturbation decays ONTO the cycle and reads sigma<0. Blend with
   `analysis/blend_fields.py` (eps~0.1) to land an order of magnitude below it.
3. **Keep `phi`** by starting from the steady case itself -- but that only
   REDUCES the startup jolt (9x, measured), it does not remove it. The low-dt
   startup at 4e-5 is still required; skipping it diverged the run in 180 steps.
4. **Always pass `--freq`** to `analysis/eigenvalue.py`. See below.
5. **Check the SPAN before the R2.** An R2 of 0.999 over 6 t and 2-3 independent
   points taught nothing at AoA 14; 63 t and 21 windows is a measurement.

**Never identify a frequency by peak height -- this now has five instances.**
The newest: `eigenvalue.py` tracked `Y.max()`, a decaying seed transient reached
2.47e-4 at f=0.2 against the eigenmode's 1.9e-5, the tracker jumped between them,
and the fit reported **sigma = -0.10 where the mode was at +0.0067** -- wrong
sign. An `f_min` floor did NOT fix it (the transient reached exactly 0.2). Probe
the frequency you want: `--freq 0.600011`.

**A guard that checks for the absence of what you just deleted proves nothing.**
A regex strip of `relaxationFactors` left `fields { p 0.3; }` and an orphan `}`
at top level; `grep -q relaxationFactors` passed because the word was gone, and
the malformed dictionary diverged the solver twice. Copy known-good files whole
and validate with `foamDictionary`, which fails on any malformation rather than
the one you anticipated.

## Physics traps documented upstream

- **A converged steady solve is not evidence the flow is steady.** `steadyState`
  deletes the ddt term, so SIMPLE converges to a fixed point of the steady equations
  whether or not it is stable. At alpha=10, Re=1000 it overpredicts Cl by 9.8%.
- **`Re*sin(alpha)` does not predict shedding onset.** It is a bluff-body criterion
  applied to a streamlined body; for an airfoil the wake width is set by where the
  boundary layer separates. Do not use it to estimate shedding frequency either.
  Measured at Re=1000: onset between 8 and 10 deg; f_shed 0.867 (alpha=10), 0.667
  (alpha=15).
- **`wallShearStress` returns `-devTau & n`** (`wallShearStress.C:82`), the negative of
  the traction on the wall. Physical skin friction is the negation of the sampled
  value; getting this backwards makes an attached boundary layer look 100% separated.
- **The transient can go quiet and then restart.** At alpha=10 the impulsive-start
  disturbance decays to a near-constant Cl~0.353 over t=8-12 before the real
  instability re-emerges and saturates near t~22. A run stopped at t=10 reports
  "steady, no shedding" and is wrong by 15%. `analysis/limit_cycle.py` guards against
  this with a spectral-prominence test, and the driver's `T1_MIN` is the backstop.

## Other agent configs

A Gemini CLI config exists at `~/.gemini/settings.json`. To import user-level items
(MCP servers, slash commands, subagents, skills, instructions), reply `/import` to scan
and list what's importable, then `/import --yes=<digest>` to apply. If `/import` isn't
available on this surface, run `claude import` from a terminal.
