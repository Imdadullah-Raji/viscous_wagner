# Session hand-off, 2026-10-20

The session ran 2026-10-09 to 2026-10-20. Most of the context is in the transcript; this file
points to where things are.

## Where to look

| what | where |
|---|---|
| **Full transcript** (everything, verbatim) | `~/.claude/projects/-home-raji-Research-Thesis-static-airfoil/8fe1c66e-f6eb-468e-918b-0b0d34827803.jsonl` |
| Persistent memory (machines, rules, ZTP) | `~/.claude/projects/-home-raji-Research-Thesis-static-airfoil/memory/` (`MEMORY.md` is the index) |
| Current plan (a5/a7 tight convergence) | `~/.claude/plans/flickering-prancing-iverson.md` |
| Project instructions | `CLAUDE.md` (still says OpenFOAM dev; see "Environment" below) |
| **Cubic stress test** (a3 robustness, a5, resolution, R) | `analysis/CUBIC_STRESS_TEST.md` / `.pdf`, script `analysis/cubic_stress_test.py` (run with the `research` env python), outputs `analysis/cubic_stress_test/` |
| Summary report (7 pp) | `reports/cubic_baseline/cubic_baseline.tex` / `.pdf` |
| Zero-thickness plate (ZTP), Re 500 + mesh history | `zero_thickness_plate/A3_RE500.md`, `zero_thickness_plate/zeroThickPlate.md` |
| ZTP mesher (elliptic grid, used) | `zero_thickness_plate/ell_mesher.py`; meshes in `zero_thickness_plate/meshes/` |
| ZTP Re sweep (Re 200-900, 0-6 deg at 0.5 deg; production tolerances) | `zero_thickness_plate/ell_re500/` (0.25-deg grid), `zero_thickness_plate/ell_re_scaling/steady_re<RE>/` |
| **Tight runs** (in progress) | `zero_thickness_plate/tight/re{500,700,900}/`, Step 0 test in `zero_thickness_plate/tight_test/` |
| Docs-to-fix list for the sharp-LE corrections (user edits) | `SHARP_LE_CORRECTIONS.md` |
| Git repo for the viscous Wagner programme | `~/Research/Repositories/viscous-wagner` (notes: `GIT_REPO.md`, `HYPOTHESIS_PROVENANCE.md`); nothing from the ZTP work committed yet |

## Status, 2026-10-20 06:45

**Goal:** stable a5 and, if possible, a7 for the ZTP at Re 500, 700 and 900.

**Finding:**
- Production SIMPLE runs **stall** at the linear-solver tolerance: p 1e-7, U 1e-8, with 0 iterations
  per solve near the end.
- The stall leaves C_L 2.5e-7 to 3e-6 from the fixed point, acting as noise between angles. a5 and a7
  need about 1e-9.
- `fvSolution.tight` (tolerance 1e-12) plus about 3000 iterations fixes this.
- Write precision (12) is not a limit.
- Extending the angle range does not help, because the series radius is only about 6–8° at Re 900.

**Runs:**

| machine | queue | state |
|---|---|---|
| this machine (6 ranks) | Re 900 quarter angles (9 new), `tight/run_quarters.sh`: production 6000, then tight 3000 | production phase 8/9 done |
| machine B (4 + 4 ranks) | Re 500 tight continuations (19); Re 900 tight continuations (8) | Re 500: 14 done; 4–5° (5 cases) re-pushed and relaunched (`run_re500_b.log`). Re 900: 7/8 |
| cfd-test (4 ranks) | Re 700 tight continuations (9), then Re 700 quarter angles (9 new, `run_re700q.log`) | continuations done; quarters being pushed and launched |

**Check:** run `zero_thickness_plate/tight/status.sh` (read-only; all three machines).

**Remote folders:**
- machineB: `<machine-B work dir>/ztp_tight/`
- cfd-test: `/home/raji/cases/ztp_tight/`

**Pull back:**
- machineB: `rsync -a --exclude 'processor*' machineB:<dir>/ <local>/`
- cfd-test: the same, plus `--rsync-path="sudo rsync" -e "ssh -l <cfd-test login>"`

## Next steps

1. When the queues finish, pull `re500`/`re900` from machine B and `re700` from cfd-test into
   `zero_thickness_plate/tight/`.
2. Check that the tight phase converged: the C_L change over the last 1000 iterations is ≲ 1e-10
   (`analysis/steady_summary.py --root zero_thickness_plate/tight/re<RE>`; rng(10%)). Extend
   with `zero_thickness_plate/tighten_case.sh <iters> <case>` if not.
3. Analysis (plan Step 2):
   - Odd Legendre fits P1–P9 and P1–P11 on windows A = 2.5…≈R/2 (0.25° steps), using `fit_window()`
     in `analysis/cubic_stress_test.py`.
   - Call a5 or a7 resolved only on a plateau with spread/se ≲ 3.
   - Confirm that a3 does not move.
   - Compare sqrt|a5/a7| with sqrt|a3/a5|.
4. Update `analysis/CUBIC_STRESS_TEST.md`/`.pdf` and `reports/cubic_baseline`. Flag that the a5 mesh
   sensitivity (E560 vs E800) is unmeasured.

## Rules and pitfalls

- **machine B is a friend's machine.** Run there only when the user explicitly asks; this job was
  authorised. Work only in `<machine-B work dir>/`.
- **cfd-test:** log in as `-l <cfd-test login>` (sudo), run as `sudo -u raji`. Never touch
  the user's `eldredgeRamp*` cases.
- **Environment:**
  - OpenFOAM dev was **removed** from this machine on 2026-10-19. v14 is at `/opt/openfoam14`, the same
    build as machine B (`14-a4867fc914c1`); continuations here and on machine B agree to 6e-12.
  - cfd-test runs a different v14 build (`14-7b05503f98a8`), so keep each Re's tight runs on one build.
  - `make_case.sh` now sources v14 (`FOAM_BASHRC=` overrides). Other root scripts still source the
    missing dev path; the user has not yet decided.
- **Tight from a cold or warm seed is slow** (about 3.5 s/iteration). Always run production tolerances
  first, then tight (`run_quarters.sh`).
- **Never `pgrep`/`pkill` a pattern that appears in your own command line.** It kills the tool shell.
  Put kill logic in a script file and filter by process cwd.
- Never rewrite or force-push git history; no commits until the user asks.
