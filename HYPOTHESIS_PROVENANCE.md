# Hypothesis provenance

Where each idea and hypothesis in this programme came from: the user, ChatGPT, or Claude. Who wrote
the *code* for a step is a separate question, recorded in each commit's "Authorship" list.

Claude here means Claude Opus 5.5, used either in a chat conversation or through Claude Code
during a working session.

## The steady lift-curve fits (October 2026)

As stated by the user in `static_airfoil/instructions/oct10_gitVersionControl.md`:

| idea | source | where it enters | files |
|---|---|---|---|
| LASSO regression on the lift curve | the user | 2026-10-05 | `analysis/lasso_reg.py` |
| Inspecting the LASSO path | suggested by ChatGPT | 2026-10-05 | `analysis/lasso_path.py` (code by Codex) |
| Q-Q and residual plots (concepts new to the user) | suggested by ChatGPT | 2026-10-05 | `analysis/preshedding_slope.py` (code by the user) |
| A cubic is enough to describe the nonlinearity at the origin; fit odd powers only, by symmetry | the user | 2026-10-05 | `analysis/fit_odd_cubic*.py`, `check_odd_cubic.py` |
| Use an orthogonal basis instead of monomials | suggested by Claude Opus 5.5 in a conversation | 2026-10-06 | — |
| Legendre polynomials as that basis | the user | 2026-10-06 | `analysis/legendre_odd_fit.py` |
| Try a saturating lift law | suggested by Claude Opus 5.5 in a conversation | 2026-10-06 | — |
| The specific saturating form A·x/√(1+B·x²) | the user | 2026-10-06 | `bizarre_experiments/` |
| a3 follows a scaling law in Re (NACA0012) | suggested by Claude Code (Claude Opus 5.5) mid-session. **It computed the power-law fit on its own, without being asked.** | 2026-10-06 | `analysis/cubic_vs_re.py`; tested further in `legendre_vs_re.py`, `STEADY_LASSO_ReSweep.md` |

Paths are under `static_airfoil/`.

## Other hypotheses, as recorded in the project files

| idea | source | record |
|---|---|---|
| The NACA0012 wake period-doubles past 20° at Re 500 | reported by a Claude agent; this led the user to Kurtulus (2016), Rossi, Durante et al., and Gupta et al. (JFM, doi:10.1017/jfm.2022.958) | `oct10_gitVersionControl.md`; `FINDINGS.md` §9 |
| A laminar separation bubble explains why cubic dominance wanes from NACA0012 to NACA0006 to the plate | the user (tested and rejected) | `instructions/lsb_hypothesis.md`; `lsb/HYPOTHESIS_TEST.md` |
| dC_L/dα is quadratic at small angles for all three bodies, with a3 following a local scaling law in Re | the user's claim | `instructions/Quadratic_hypothesis.md` |
| The plate's non-monotonic a3(Re) comes from its sharp leading edge | an inference written by Claude Code (Claude Opus 5.5) in the a3 analysis; drawn from three bodies, not tested by varying the leading edge alone | `plate_case/analysis/A3_RE_SCALING.md`, `static_airfoil/analysis/A3_RE_SCALING.md` |

## The zero-thickness plate and the stress test (October 2026)

| idea | source | where it enters | record |
|---|---|---|---|
| A zero-thickness plate, to test the effect of nose radius on a3 and to set up step responses against the classical Wagner function | the user | 2026-10-18 | `zero_thickness_plate/zeroThickPlate.md`, `zeroThickPlate.geo` |
| The thick plate's leading edge is blunter than NACA0006's, so "sharp leading edge" is the wrong explanation for its a3(Re) | found by Claude Code (Claude Opus 5.5) | 2026-10-17 | `SHARP_LE_CORRECTIONS.md` |
| An orthogonal elliptic-coordinate grid for the plate; the H-grid's θ^3.4 is a mesh artefact | suggested by Claude Code (Claude Opus 5.5) | 2026-10-18 | `zero_thickness_plate/ell_mesher.py`, `A3_RE500.md` |
| Stress-test "cubic at the origin": do higher coefficients vanish, how a3 varies with a5, how far truncations extrapolate, and whether a singularity such as a bifurcation explains where they fail | the user | 2026-10-19 | `analysis/CUBIC_STRESS_TEST.md` |
| The cubic range shrinks with Re because a singularity moves toward α = 0, as the thick plate's bubble formation angle does | an inference written by Claude Code (Claude Opus 5.5); withdrawn for Re 500–900 on re-converged data | 2026-10-19 | `CUBIC_STRESS_TEST.md` §5, §7 |
| Pursue stable a5 and a7 for the zero-thickness plate at Re 500, 700, 900 | the user | 2026-10-19 | `SESSION_HANDOFF_2026-10-20.md` |
| The production solver tolerances stall SIMPLE short of the fixed point; that, not noise, limits a5/a7 | found by Claude Code (Claude Opus 5.5) | 2026-10-19 | `zero_thickness_plate/tighten_case.sh`, `CUBIC_STRESS_TEST.md` §7 |
| One-signed coefficients a3…a11 point to a real-axis singularity near 9–11° (e.g. a fold) | an inference written by Claude Code (Claude Opus 5.5); a hypothesis, not tested | 2026-10-20 | `CUBIC_STRESS_TEST.md` §7 |

Paths are under `static_airfoil/`.
