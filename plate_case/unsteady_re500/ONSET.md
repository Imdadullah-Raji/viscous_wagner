# Plate shedding onset, Re = 500 (oct8_workflow.md)

2026-10-08/09. Unsteady runs with `static_airfoil/run_staged.sh`, `NO_FIELDS=1` (forces only,
no snapshots). Mesh `../re500/meshes/plate_ns160_nr180.msh` (63,360 cells), started
impulsively from uniform flow. Each run was stopped by hand once its state was clear: the
driver was killed first, then `mpirun`. Each log ends with a `MANUAL STOP` note.

## Result

**Onset at α_c ≈ 11.3–11.4°.** It is bracketed by steady flow at 10.75° and shedding at 12°,
with 11.375° marginal.

| α | outcome | evidence | stopped at |
|---:|---|---|---:|
| 8.5 | **steady** | C_L → 0.60240; the shedding-band residual decays to ~1e-7, the data floor; `shedding=no` at t = 60.5 | t ≈ 61 |
| 9.5 | **steady** | C_L → 0.62280; residual ~1e-7; `shedding=no` at t = 60.5 | t ≈ 61 |
| 10.75 | **steady** | C_L → 0.62780; residual 5e-4 at t ≈ 7 → ~1e-6 to 1e-7 for t ≥ 38; no f ≈ 0.63 mode | t = 44.1 |
| 11.375 | **marginal** | f ≈ 0.63 oscillation (period ~1.5) holds at a peak-to-peak of 6.4e-5 to 7.2e-5 from t = 14 to 33, σ ≈ +0.009 | t = 35.5 |
| 12 | **sheds** | f ≈ 0.63 mode grows from t ≈ 10 at σ ≈ 0.15 and saturates at a peak-to-peak of ~0.023; mean C_L 0.653, still rising | t = 51.6 |
| 15 | **sheds** | saturated from t ≈ 12: C_L ≈ 0.80 ± 0.10, f ≈ 0.62 (period 1.60–1.65); mean still creeping | t = 52.6 |

- Linear interpolation of σ between 11.375° (+0.009) and 12° (+0.15) gives σ = 0 at about
  **11.33°**. The NACA0012's eigenvalue onset at the same Re is 11.46°
  (`static_airfoil/FINDINGS.md` §19). At Re = 500 the two sections lose stability at nearly the
  same angle, despite very different steady lift curves and separation patterns.
- The steady LSB at 8°, which bursts by 8.5° (`static_airfoil/lsb/`), lies **well below** the
  shedding onset. The bubble forms and bursts on a stable steady branch.
- The time-averaged lift rises above the steady value once the flow sheds. At 12° it is
  0.653 against the steady 0.621; at 15° it is about 0.80. This matches the NACA0012.

## Caveats

1. **σ here is not a clean eigenvalue measurement.** It comes from an impulsive start with no
   deliberate seed. The time step was changed at t = 20.5 and t = 40.5, and the mean lift was
   still drifting. The 11.375° value (+0.009) is at the resolution limit of this method. Whether
   that angle is just above or just below onset is **not resolved**. A precise α_c needs
   `static_airfoil/run_marginal.sh`: a fixed dt, a start from the steady solution, a small blended
   seed, and σ fitted at `--freq ~0.63`.
2. **No run reached a converged limit cycle**, because all were stopped early by request. The
   12° and 15° means and amplitudes are good to a few percent at best.
3. **One mesh only** (the base plate mesh). The 8° bubble moved about 10 % on the finer mesh
   (`../mesh_conv_re500/RESULTS.md`). The onset has not been checked for mesh independence.
4. The 15° case was started only because 9.5° was steady (user request). 12°, 10.75° and
   11.375° then bisect the bracket.

## Files

- **Cases:** `aoa{8.5,9.5,10.75,11.375,12,15}/`.
- **Driver logs:** `run_aoa<A>.log`.
- **Build logs:** `aoa<A>.build.log`.
- **Chain scripts:** `run_aoa9.5_after_meshconv.sh`, `run_aoa15_if_95_steady.sh` and
  `run_aoa12_after_15.sh`. These wait on log markers, not on `pgrep` of command lines, which
  matched their own launching shells.
