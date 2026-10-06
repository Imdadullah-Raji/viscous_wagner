# Saturating lift fit

Unweighted nonlinear least squares on the 17 original CSV observations from 0–8°. Angle x is in radians; no intercept.

`Cl = 3.31190786059 * x / sqrt(1 + 24.0441764718 * x^2)`

- A = 3.31190786059
- B = 24.0441764718
- Training R² = 0.999949075885
- Positive-angle asymptote A/√B = 0.675419029632
- Four initial guesses converged to the same objective value.

| Set (converged cases only) | n | Model RMSE | Cubic RMSE |
|---|---:|---:|---:|
| training | 17 | 0.000845534128 | 0.000102738583 |
| interpolation | 14 | 0.000776683524 | 8.26862838e-05 |
| extrapolation_8_12 | 12 | 0.0172643587 | 0.0117707461 |
| extrapolation_12_24 | 4 | 0.0809788501 | 0.330788182 |

The comparison cubic includes an intercept and is independently fitted to the same 17 training points.
The 22° and 24° saved SIMPLE traces are unconverged; their final-iteration values are provisional and excluded from RMSE summaries.
Saturation is a model assumption, not a demonstrated large-angle physical limit.

| Angle | Observed | Predicted | Observed − predicted | Status |
|---|---:|---:|---:|---|
| 8 | 0.37969880 | 0.38156728 | -0.00186848 | CONVERGED |
| 8.5 | 0.39347780 | 0.39732455 | -0.00384675 | CONVERGED |
| 9 | 0.40576240 | 0.41214871 | -0.00638631 | CONVERGED |
| 9.25 | 0.41136040 | 0.41922372 | -0.00786332 | CONVERGED |
| 9.5 | 0.41660560 | 0.42608116 | -0.00947556 | CONVERGED |
| 9.75 | 0.42150520 | 0.43272653 | -0.01122133 | CONVERGED |
| 10 | 0.42607190 | 0.43916536 | -0.01309346 | CONVERGED |
| 10.25 | 0.43031980 | 0.44540325 | -0.01508345 | CONVERGED |
| 10.5 | 0.43426190 | 0.45144578 | -0.01718388 | CONVERGED |
| 10.75 | 0.43790950 | 0.45729855 | -0.01938905 | CONVERGED |
| 11 | 0.44127850 | 0.46296711 | -0.02168861 | CONVERGED |
| 11.5 | 0.44724960 | 0.47377355 | -0.02652395 | CONVERGED |
| 12 | 0.45229360 | 0.48390835 | -0.03161475 | CONVERGED |
| 14 | 0.46563930 | 0.51854258 | -0.05290328 | CONVERGED |
| 16 | 0.47270970 | 0.54545080 | -0.07274110 | CONVERGED |
| 18 | 0.47754050 | 0.56652067 | -0.08898017 | CONVERGED |
| 20 | 0.48207580 | 0.58318367 | -0.10110787 | CONVERGED |
| 22 | 0.49452320 | 0.59650458 | -0.10198138 | LIMIT-CYCLE |
| 24 | 0.50535760 | 0.60727061 | -0.10191301 | LIMIT-CYCLE |

Run: `conda run -n research python bizarre_experiments/fit_saturating.py`

[Diagnostics through 12°](saturating_diagnostics_12deg.png) · [Extrapolation through 24°](saturating_extrapolation_24deg.png) · [CSV](saturating_predictions.csv)
