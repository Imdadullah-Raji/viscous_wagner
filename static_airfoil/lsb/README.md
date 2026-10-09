# Upper-surface separation: plate, NACA0006, NACA0012 at Re = 500

These are steady (SIMPLE) solutions at each case's latest iteration. The script is
`streamlines_lsb.py`; run it with the research env and one or more `--aoa` values. Figures go to
`figs/aoa<A>/`:

- `streamlines_full.png`: the whole section.
- `streamlines_aft.png`: the rear, x = 0.35–1.2.
- `streamlines_le.png`: the leading edge, x = −0.03–0.4.
- `cf_upper.png`: upper-surface skin friction.
- `separation.txt`: the separation and reattachment points.

How to read the figures and data:

- Colour is the chord-wise velocity u_x, on a colour scale diverging about 0, so blue means
  reversed flow. The black line is u_x = 0.
- The geometry sits at zero incidence, because the freestream is rotated rather than the body.
  The mesh x-axis is therefore the chord.
- Skin friction is C_f = 2τ_x. It is the *negated* `wallShearStress` output, because that function
  object returns −devTau·n.
- ▼ marks separation and ▲ marks reattachment, both found from zero crossings of C_f for x > 0.02.
- The small u_x < 0 patch just below each leading edge is the stagnation region, not separation.
  The skin friction there is positive.

Prerequisites for each case:

- `foamPostProcess -func writeCellCentres -latestTime` must have been run.
- The case must be reconstructed. The NACA0012 cases at 7.5°, 8.5°, 9° and 9.5° were run
  decomposed; `reconstructPar -latestTime` was added on 2026-10-08.

## Separation on the upper surface (x/c)

| α | plate | NACA0006 | NACA0012 |
|---:|---|---|---|
| 7 | attached (min C_f +0.0047) | attached (+0.0025) | 0.715 → 0.990 |
| 7.5 | attached (+0.0030) | *no case* | 0.655 → 0.991 |
| 8 | **0.137 → 0.424, a closed bubble** (min C_f −0.0008) | 0.832 → 0.972, a closed bubble | 0.599 → 0.991 |
| 8.5 | 0.079 → 0.979 | *no case* | 0.548 → 0.991 |
| 9 | 0.061 → 0.988 | 0.542 → 0.985 | 0.501 → 0.991 |
| 9.5 | 0.051 → 0.991 | *no case* | 0.459 → 0.991 |
| 10 | 0.045 → 0.992 | 0.356 → 0.988 | 0.420 → 0.991 |

A reversed-flow end at about 0.99 is the trailing-edge closure, not a reattachment. Only the
plate at 8° and the NACA0006 at 8° reattach upstream of the trailing edge. The NACA0006 sweep has
integer angles only.

## Plate LSB formation and growth, fine mesh (2026-10-09, `instructions/oct9_workflow.md`)

Steady SIMPLE, Re = 500, on the fine plate mesh
(`plate_case/mesh_conv_re500/meshes/plate_ns160_nr260_n195_h025.msh`, h1 = 0.025/√Re). Every
case is a cold start with 6000 iterations; the cases are in `plate_fine_re500/steady_aoa<A>`.
They cover 7.5–8.5° in steps of 0.1°, and 7.91–7.99° in steps of 0.01°. The 0.01° cases also
write at iterations 3000, 4000 and 5000 so iteration convergence can be checked.

Machines: this machine ran 7.5–8.0, 8.4, 8.5, 7.91–7.96, 7.98 and 7.99. banana-cfd
(OpenFOAM 13, aarch64 QEMU VM) ran 8.1–8.3 and 7.97, using `run_remote.sh`. That VM's speed
varied by up to 7×, so cases were moved back here as needed. Its logs are copied as
`run_remote.log`.

**Formation: α_f = 7.9087°.**
- The smallest upper-surface C_f (0.02 < x < 0.9) falls exactly linearly with angle:
  −0.00752 per degree, with a maximum residual of 1.2e-6 over 7.90–8.00.
- Its first zero is at x/c = 0.225.
- Iteration error is ±2e-6 in C_f,min, which is ±0.0003° in α_f.
- The production mesh gives 7.892°, so the formation angle is mesh-insensitive.

**Growth near formation.**
- L ∝ (α − α_f)^n, with local n between 0.53 and 0.61 over 7.91–8.00
  (`figs/plate_fine/formation.png`).
- This is the √ law of a parabolic C_f minimum pushed through zero.
- n then rises to 0.65, 0.78 and 0.87 over 8.0–8.3: growth turns nearly linear, because the
  C_f minimum is flat-bottomed.
- The bubble grows mostly downstream:

  | α | separation x/c | reattachment x/c |
  |---|---|---|
  | 7.91 | 0.216 | 0.241 |
  | 8.0 | 0.143 | 0.401 |
  | 8.3 | 0.094 | 0.841 |

**Reattachment reaches the trailing edge.**

| α | 8.3 | 8.4 | 8.5 |
|---|---|---|---|
| reattachment x/c | 0.841 | 0.967 | 0.979 |

- The reattachment point is pinned just ahead of the rounded trailing edge, which starts at
  x ≈ 0.988.
- The bubble is still closed at 8.5°: its ψ = 0 line returns to the wall, and there is no
  reversed flow in the near wake.
- **"Burst" is therefore not yet defined or observed on the plate.** The "0.079 → 0.979" entry
  for 8.5° in the table above is still a closed bubble.

**Lift does not register formation.** C_L increments per 0.1° fall smoothly from 0.0047 to
0.0033 across 7.5–8.5°, with no kink at α_f.

**Resolution caveat.** At 8.0° the bubble is 0.0019c tall, and its reversed flow occupies only
the first cell row (cell centre 0.00056c; magenta ticks in `bubble_zoom_*.png`). The formation
angle depends only on the sign of the first-cell velocity, so it is robust. Bubble lengths near
formation are not: the production mesh gives a length 10 % different at 8°.

**Scripts and outputs.**
- `plate_sweep.py` → `figs/plate_sweep/`. Its 0.95 "burst" cut-off is wrong for the plate; see
  above.
- `formation_fine.py` → `figs/plate_fine/formation.{txt,csv,png}`.
- `plate_bubble_zoom.py` → `figs/plate_fine/bubble_zoom_*.png`: the stream function with the
  wall-normal axis stretched.
- `streamlines_lsb.py --only plate --plate-root plate_fine_re500 --tag plate_fine` →
  `figs/plate_fine/aoa<A>/`.
