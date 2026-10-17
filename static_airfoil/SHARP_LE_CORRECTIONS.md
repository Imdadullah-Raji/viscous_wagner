# Documents that call the plate's leading edge "sharp"

2026-10-17. This is a pointer list only: the wording is left for you to change.

## The issue

The plate's leading edge is **not sharper** than the airfoils'. By radius it is blunter than
NACA0006's.

| body | t/c | leading-edge radius |
|---|---:|---:|
| plate (stadium, `plate_mesher.py --th 0.023`) | 0.023 | 0.0115c (semicircle, r = t/2) |
| NACA0006 | 0.06 | 0.0040c (1.1019·t²·c) |
| NACA0012 | 0.12 | 0.0159c |

What sets the plate apart is the **shape**, not the radius. Its semicircular nose joins the flat
sides with a jump in curvature at the shoulder, and a flat upper surface follows. On a NACA section
the curvature changes smoothly from the nose aft.

## Where it appears

### Workspace (`/home/raji/Research/Thesis`), also copied into the repo

| file | line | text |
|---|---:|---|
| `static_airfoil/analysis/A3_RE_SCALING.md` | 87 | table header `plate (sharp LE)` |
| same | 102 | "The plate has no scaling law for a3. The most likely cause is its sharp leading edge." |
| same | 104 | "Neither rounded-nose airfoil shows that bubble before 6°": implies the plate's nose is not rounded |
| same | 105 | "most likely a leading-edge effect, not a thickness effect" |
| same | 106 | "not tested by varying the LE shape alone" (consistent with the issue; check that it still reads right) |
| `static_airfoil/instructions/oct10_gitVersionControl.md` | 75 | your provenance line: "the plate's a3 peak comes from its sharp leading edge" |

### Repo only (`/home/raji/Research/Repositories/viscous-wagner`)

| file | line | text |
|---|---:|---|
| `README.md` | 29–30 | "This is most likely an effect of its sharp leading edge." |
| `HYPOTHESIS_PROVENANCE.md` | 34 | row "The plate's non-monotonic a3(Re) comes from its sharp leading edge" |

### Commit messages (already on GitHub)

Rewriting history was ruled out, so these stay. A later commit could note the correction.

| commit | date | text |
|---|---|---|
| `adc5b68` | 2026-10-06 | subject "Add a thick flat plate as the sharp-leading-edge counterpart" |
| `3cd3de8` | 2026-10-10 | body "non-monotonic a3 is most likely a sharp-leading-edge effect" |

## Checked and fine

- `plate_case/analysis/A3_RE_SCALING.md`: talks about the approach to leading-edge separation, not
  the edge being sharp.
- `lsb/HYPOTHESIS_TEST.md` and `reports/lsb_hypothesis/lsb_hypothesis.tex`: no "sharp".

Search used:
`grep -rn -i -E "sharp|rounded-nose|LE shape" --include=*.md --include=*.tex static_airfoil plate_case`.
