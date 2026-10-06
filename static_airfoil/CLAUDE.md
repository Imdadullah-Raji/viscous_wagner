# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Snapshot of `static_airfoil/` as of 2026-10-06. The authoritative, longer guide is the root
`CLAUDE.md`; this file condenses it and adds the steady-lift fitting layer built on 2026-10-06,
which the root file only partly describes.

## What this repository is

A 2D laminar CFD campaign on a static NACA0012 (OpenFOAM-dev), producing force and field data for
a POD reduced-order-model thesis. It has two strands:

1. **Unsteady shedding at Re = 500**: staged transient runs, limit-cycle detection, period
   doubling, and the eigenvalue-measured onset α_c = 11.46°. Records live in `FINDINGS.md`,
   `BIFURCATION_THEORY.md` and `ANALYSIS_GLOSSARY.md`.
2. **Steady-branch lift fitting at Re = 100, 500, 1000, 2000**: SIMPLE sweeps, odd-polynomial,
   Lasso and Legendre fits, and an Re-scaling of the cubic term. Records live in
   `STEADY_LASSO_Re1000.md` and `STEADY_LASSO_ReSweep.md`, and the reader-facing PDF is
   `reports/steady_lift_fits/`.

Task requests sometimes arrive as files in `instructions/` (for example `oct6_reSweep.md`); read
that folder when asked.

## Environment

- **OpenFOAM:** `source /opt/openfoam-dev/etc/bashrc`. The solver is `foamRun` with
  `solver incompressibleFluid`.
- **Post-processing:** use `foamPostProcess -solver incompressibleFluid`. Plain `postProcess`
  cannot build the transport model.
- **Python:** `/home/raji/miniconda3/envs/research/bin/python`. Only this env has gmsh, sklearn,
  pandas and scipy, and system `python3` lacks pandas. Plain `python3` works only for the
  numpy-only scripts (`steady_summary.py`, `limit_cycle.py`, `growth_rate.py`).
- **6 physical cores.** Budget 6 MPI ranks in total across all jobs, and run cases sequentially:
  one 6-rank job is faster than 3+3 or 2+2+2 (measured several times).
- **LaTeX:** pdflatex/latexmk (TeX Live 2022) is installed; `siunitx` is not.

## Commands

```bash
# mesh: the M2 recipe; only --re changes (first cell h1 = 0.05/sqrt(Re), so a mesh is Re-specific)
cd meshing && /home/raji/miniconda3/envs/research/bin/python mesher.py \
    ../meshes/M2_ns160_nr180_b025_re<RE>.msh --re <RE> -r 30 --nr 180 --ns 160 --nte 32 --bump 0.25
#   (Re=500's file is meshes/M2_ns160_nr180_b025.msh; the recipe reproduces it byte for byte)

# steady (SIMPLE) sweep: cases live in steady_re<RE>/steady_aoa<A>
./make_steady.sh steady_re<RE>/steady_aoa<A> <RE> <A> <mesh.msh> [nIter=6000]
STEADY_DIR=steady_re<RE> ./run_steady.sh 0 0.5 1 ...    # sequential, 6 ranks each; default dir steady_re500
steady_re<RE>/sweep.sh                                   # the per-Re driver that built + ran a sweep
python3 analysis/steady_summary.py --root steady_re<RE>  # verdict per angle, judged on the C_L trace

# unsteady
./make_case.sh <caseDir> <Re> <AoA> <mesh.msh>
./run_staged.sh <caseDir> 6                  # phase 0 -> phase 1 (forces) -> phase 2 (fields)
NO_FIELDS=1 ./run_staged.sh <caseDir> 6      # stop after phase 1: "does it shed?", forces only
DT_FIXED=<dt> ./run_staged.sh <caseDir> 6    # freeze dt for growth-rate (sigma) measurements
python3 analysis/limit_cycle.py <caseDir>    # KEY=VALUE: SHEDDING, ORDER, FREQ, CONVERGED, ...
python3 analysis/growth_rate.py <caseDir>    # envelope sigma; < 0 means decaying

# steady-lift fits (research env)
PY=/home/raji/miniconda3/envs/research/bin/python
$PY analysis/lasso_re1000.py [--re 100 500 1000 2000] [--amax 8] [--odd]
(cd analysis && $PY legendre_odd_fit.py --amax 8)   # odd Legendre fit + rewrite into powers of theta
(cd analysis && $PY cubic_vs_re.py)                 # OLS b1, b3 vs Re, power law
(cd analysis && $PY explain_figs.py)                # explainer figures for the PDF report
cd reports/steady_lift_fits && latexmk -pdf steady_lift_fits.tex && latexmk -c
```

There is no test suite. Correctness is checked by regression: rerun a script with its default
arguments and `diff` the `lasso_summary.txt` it writes against the previous copy. Every change to
`lasso_re1000.py` so far has been checked this way.

## Architecture: the big picture

### Angle of attack is a property of the inflow, not the mesh

One mesh serves every angle at a given Re. `make_case.sh` writes the angle into three places from
one computation: `0/U`, the `forceCoeffs` lift and drag directions, and `Cp`'s `UInf`. Never edit
one of them by hand. This is also why `warm_start.sh` can map fields 1:1 between angles; it checks
mesh md5 because meshes differ across Re.

Reference values: `Aref 0.1` (chord × span), `lRef 1`, `magUInf 1`.

`gmshToFoam` resets patch types, so `make_case.sh` retypes `frontAndBack` → `empty` and
`airfoil` → `wall`.

`make_steady.sh` wraps `make_case.sh` and changes only three things: ddt becomes `steadyState`,
2 non-orthogonal correctors, and relaxation U 0.7 / p 0.3. The schemes are kept identical on
purpose, so that steady and unsteady differ only in the time derivative.

### Staged unsteady runs (`run_staged.sh` + `system/functions`)

- **Phase 0:** impulsive start at dt = 4e-5 to t = 0.5. Courant is measured *at* t = 0.5.
- **Phase 1:** field-free chunks of 20 time units, with dt re-tuned only between chunks.
  `adjustTimeStep` is `no` everywhere. The loop exits when `limit_cycle.py` reports
  `CONVERGED=yes`.
- **Phase 2:** fields at ≥ 32 samples per carrier period for 20 orbits.
- **`DT_MAX = 6e-4` is the real stability limit, not Courant.** A diverged chunk is retried at
  dt/2, and its `forceCoeffs/<restart>/` is quarantined to `.diverged/` first.
- **Stopping a run early:** kill the `run_staged.sh` driver *before* `mpirun`. Killing the solver
  alone looks like a divergence and triggers the quarantine-and-retry.
- **`system/functions` is byte-identical across stages.** `forceCoeffs` runs every step into one
  `.dat`; field function objects fire only at write times.
- **Restarts split the force series** across `postProcessing/forceCoeffs/<restartTime>/` with
  duplicated times at the seams. Concatenate, sort, and keep the last sample per time.
- **Count time directories by value** (`213`, not `213.0`).

### Two periods once the wake period-doubles

`limit_cycle.py` reports both the carrier (`SAMPLE_PERIOD`, which drives `writeInterval`) and the
fundamental (`PERIOD`, which drives run length and the convergence window). Never identify a
frequency by peak height. Always go through `spectral_structure()`, and pass `--freq` to
`eigenvalue.py`.

### Steady-branch data and its traps

- **SIMPLE convergence is judged on the force trace** (`steady_summary.verdict`):
  - `CONVERGED`: the range over the last 10 % is below tolerance.
  - `CREEPING`: the range is still shrinking.
  - `LIMIT-CYCLE`: the range is flat across windows, so the base flow is strongly unstable and
    the last iterate is *not* a fixed point.

  Only `CONVERGED` cases enter any fit.
- **A converged steady solve is not evidence the flow is steady.** Above the onset the shedding
  mean sits far above the steady value; even the minimum instantaneous C_L is above it.
- **Where SIMPLE first oscillates:** Re = 100 never (to 20°), Re = 500 at 22°, Re = 1000 at 14°,
  Re = 2000 at 9.5°. The physical onset lies at or below these angles. It is measured only at
  Re = 500 (11.46°); Re = 1000 is bracketed at 8–10°; Re = 2000 is known only to shed at 9.5°.

### Fitting layer (`analysis/`)

- **`lasso_re1000.py`** is the shared engine. Despite the name it takes any list of Re: `collect()`
  loads one Re's steady cases with their verdicts, `analyse()` runs the LARS path, Lasso at
  λ = 1e-4, LOO-CV, a 500-resample bootstrap, and OLS cubic/quintic fits.
  - **Output folders:** the defaults (Re 500 vs 1000, 0–12°) write to `figs/re1000_steady/`.
    Other choices go to `re1000_steady_0to<amax>deg[_odd]/` or `resweep_steady_...`, so earlier
    figures are never overwritten.
  - **Other scripts import from it** (`collect`, `ROOT`).
- **Symmetry:** C_L(−α) = −C_L(α), so only odd powers are physical. A Lasso that selects an even
  power, such as θ⁴, is showing collinearity, not physics.
- **Raw powers of θ are nearly collinear on 0–8°** (correlations 0.91–0.99, condition number
  ~1e12 for θ¹..⁹). As a result, Lasso selection at fixed λ depends on the fit range, and
  monomial coefficients above θ³ are unstable.
- **The robust quantities** are the OLS b1 and b3, and the odd-Legendre coefficients on
  x = θ/θ_max mirrored to [−1, 1] (condition number 2.7).
- **Rewriting in θ:** a Legendre fit truncated at P_k equals the OLS fit in θ…θ^k. Quote θ-form
  coefficients only at the truncation actually used.
- **Data floor:** `forceCoeffs` writes ~7 significant digits, roughly 1e-7 in C_L. Coefficients
  or residuals below that are rounding, not physics.

### Plot conventions

Each Re keeps one colour in every figure: Re 500 `#2a78d6`, 1000 `#eb6834`, 100 `#1baf7a`,
2000 `#eda100`. Multi-series lift/drag figures use one shared legend below the panels, because a
legend inside the axes covers data.

## Directory map (non-obvious parts only)

| path | contents |
|---|---|
| `steady_re{100,500,1000,2000}/` | SIMPLE cases `steady_aoa<A>/`, build logs, `run_steady.log`, `sweep.sh`. All steady Re = 500 cases were moved here on 2026-10-06. |
| `unsteady_re2000/re2000_aoa9.5/` | forces-only shedding check, stopped at t ≈ 30 once shedding was confirmed |
| `re500_aoa*/` | Re = 500 staged unsteady cases (POD snapshots) |
| `bi_point_search/`, `mesh_convergence_M1/` | self-contained studies with their own steady cases and scripts |
| `lift_coeffs.csv` | the original Re = 500 lift table (23 angles, 0–12°). Several fit scripts read it from the repo root, so do not move it. |
| `bizarre_experiments/` | one-off fits (A·x/√(1+Bx²) saturating model, Re = 500) |
| `reports/steady_lift_fits/` | LaTeX report for a non-specialist reader; figures via `\graphicspath{../../analysis/figs/}` |
| `aoa10/` | hand-edited derivative of an upstream August case. The upstream case is the source of truth. |
| `~/OpenFOAM/raji-dev/run/Aug/unsteadyRe1000/FINDINGS.md` | prior Re = 1000 campaign. Read it before changing solver settings. |
