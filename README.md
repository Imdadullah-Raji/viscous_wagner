# The Viscous Wagner Problem

How does an airfoil at ultra-low Reynolds number respond to a step change in angle of attack?
In potential flow the answer is the Wagner function, which relaxes to the steady lift 2πα. This
programme asks what replaces it when viscosity dominates. It was prompted by S. L. Brunton and
S. T. M. Dawson (2021), *Improved Approximation to the Wagner Function by Sparse Identification of
Nonlinear Dynamics*.

2D laminar CFD with OpenFOAM (Foundation), on a NACA0012, a NACA0006 and a thick flat plate,
at Re = 100–2000.

## Where it stands: the steady baseline

A step response has to be measured against the steady state it settles to, so the first task is
the steady lift curve.

- **The lift curve is an odd cubic from the very onset:** C_L = a1·θ + a3·θ³. Higher odd terms are
  small near the origin. This is the first clear departure from the linear 2πα of potential flow.
  The slope a1 is about half of 2π, and a3 is negative and measurable down to 1°.
- **a3 depends on Re, and on the body:**
  - NACA0012: |a3| ~ Re^0.68.
  - NACA0006: |a3| ~ Re^0.42.
  - Both are local power laws whose exponents drift slowly.
  - The plate's a3 is not monotonic: it peaks near Re 400–450. This is most likely an effect of its
    sharp leading edge.

  See `static_airfoil/analysis/A3_RE_SCALING.md` and `plate_case/analysis/A3_RE_SCALING.md`.

## How the analysis got there

The commit history follows the research decisions in order.

1. **September: shedding.** The NACA0012 campaign at Re = 500 started with no aim beyond watching
   it shed.
   - A Claude agent found that the wake period-doubles past 20°. That led to Kurtulus (2016),
     Rossi, Durante et al., and Gupta et al., *Two- and three-dimensional wake transitions of a
     NACA0012 airfoil* (JFM, doi:10.1017/jfm.2022.958).
   - These runs save fields over 20 shedding cycles, and twice as long for period-2 flows.
   - The shedding onset was measured from eigenvalues: α_c = 11.46°.
   - Notes: `static_airfoil/FINDINGS.md`, `BIFURCATION_THEORY.md`, `ANALYSIS_GLOSSARY.md`.
2. **Late September: large-angle pitching.** This validation work against the literature is not
   in this repository. See the predefense report repository (`Thesis/reports/predefense`).
3. **October: the steady baseline**, avoiding large-angle steps.
   - LASSO over all powers of α, on the pre-shedding Re = 500 and 1000 curves. This was probably a
     misstep, but it is kept as a record (`STEADY_LASSO_Re1000.md`).
   - Odd powers only, because C_L is odd in α for a symmetric section.
   - Powers of α are nearly collinear on a small window (condition number about 1e12 for θ…θ⁹), so
     the coefficients were unstable. Odd Legendre polynomials on the rescaled angle are well
     conditioned (about 2.7) and give stable results (`STEADY_LASSO_ReSweep.md`,
     `reports/steady_lift_fits/`).
   - The cubic coefficient showed power-law scaling with Re on NACA0012. The whole process was then
     repeated for NACA0006 and the plate.
   - A laminar-separation-bubble explanation of the curvature was tested and rejected
     (`static_airfoil/lsb/HYPOTHESIS_TEST.md`).

Where each idea and hypothesis came from (the user, ChatGPT, Claude) is in
[`HYPOTHESIS_PROVENANCE.md`](HYPOTHESIS_PROVENANCE.md).

## Layout

```
plan/                 campaign plan (September)
static_airfoil/       NACA sections: drivers (*.sh), meshing/, analysis/, lsb/, notes (*.md),
                      instructions/ (the research brief for each step), reports/ (LaTeX),
                      NACA0012/, NACA0006/ (per-sweep driver scripts and notes), sample_cases/
plate_case/           thick flat plate: plate_mesher.py, run_sweep.py, run_batch.sh, analysis/,
                      sample_cases/
```

`static_airfoil/CLAUDE.md` is the working reference: commands, the staging contract, why the angle
is set by rotating the freestream, and which settings were measured rather than guessed.

**Not included:** fields, force histories, meshes (they regenerate from `mesher.py` and
`plate_mesher.py`), logs and PDFs. These stay in the working tree at `/home/raji/Research/Thesis`.

## Machines and OpenFOAM versions

| machine | cores used | OpenFOAM |
|---|---|---|
| Intel i5-12400 (main box; 6 physical cores) | 6 ranks | dev (build `dev-e4f852e7fd6a` in the October logs) |
| AMD Ryzen 7 7700 (8 cores) | 6 ranks | v14 |
| MacBook Air M2, Ubuntu in a 6-core Multipass VM | 6 ranks | 13 |
| Google Cloud c3-standard-8 (4 physical cores) | 4 ranks | v14 (build `14-7b05503f98a8` in its case logs) |

All are Foundation releases. On the plate's 0.25° lift curve, dev, 13 and v14 agree to about 4e-7
in C_L.

## About this history

The repository was assembled on 2026-10-10 from a working tree that had not been under version
control.

- **Commit dates are the dates of the research step,** recovered from file creation and
  modification times, run-log headers and agent session logs.
- **Where an earlier version of a file survived, it was committed first.** Sources: `.bak` and
  `.pre_*` copies, a September copy of the mesher, and an October 6 copy of `CLAUDE.md`.
- **Otherwise a file enters with its final content** at the date it was created, so a few contain
  later edits.
- **The 2026-10-09 folder move was replayed.** Files it edited are committed first with their
  original paths; the replacements were reversed exactly and checked by re-applying them. The move
  is then one commit.

**Authorship.** Every commit message ends with an "Authorship" list.
- *Claude-generated:* written by Claude Code agents running Claude Opus 5.5. This covers most drivers, analysis scripts
  and notes.
- *Codex-generated:* written by an OpenAI Codex session on 2026-10-05 (odd-cubic fits, literature
  note, M1 mesh study, the saturating-fit experiment), plus the AoA 30 flow render.
- *Written by the user:* the plans, the instructions and the first fitting scripts
  (`preshedding_slope.py`, `fit_quad.py`, `lasso_reg.py`, `extract_lift_coeffs.sh`).

Attribution comes from the Claude and Codex session logs. September has no surviving Claude logs;
those scripts are attributed to Claude because that campaign was run through Claude Code (see
`CLAUDE.md`).
