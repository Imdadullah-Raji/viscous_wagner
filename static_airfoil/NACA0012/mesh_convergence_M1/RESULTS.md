# M1 steady mesh-refinement comparison

Re = 500. Baseline: 63,360 cells; M1: 110,000 cells. Same steady equations and spatial schemes.
This is a two-mesh sensitivity study, not a three-grid Richardson/GCI assessment.

Convergence requires at least 6,000 iterations, and both Cl and Cd must satisfy:
last-1,000-iteration range and mean shift from the preceding 1,000 iterations
each ≤ 1e-5 × max(|tail mean|, 0.01).

| Angle | M1 status | Baseline Cl | M1 Cl | ΔCl (%) | Baseline Cd | M1 Cd | ΔCd (%) |
|---|---|---:|---:|---:|---:|---:|---:|
| 4 | CONVERGED | 0.2186696 | 0.21859337 | -0.03486008 | 0.1794853 | 0.17946875 | -0.0092181461 |
| 6 | CONVERGED | 0.3097054 | 0.30961701 | -0.028539817 | 0.1863028 | 0.18628115 | -0.011623312 |
| 8 | CONVERGED | 0.3796988 | 0.3796168 | -0.021596768 | 0.1959624 | 0.19593176 | -0.015637753 |
| 10 | CONVERGED | 0.4260719 | 0.42599456 | -0.018152405 | 0.2083067 | 0.20826254 | -0.02119849 |
| 12 | CONVERGED | 0.4522936 | 0.45218288 | -0.024479444 | 0.2228761 | 0.2228119 | -0.028806068 |

Differences from unconverged cases are provisional. Inspect comparison.csv for force ranges and drift.
Physical stability and steady-solver convergence are distinct; no shedding-onset conclusion follows from this comparison.
