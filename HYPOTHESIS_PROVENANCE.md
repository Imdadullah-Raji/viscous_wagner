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
