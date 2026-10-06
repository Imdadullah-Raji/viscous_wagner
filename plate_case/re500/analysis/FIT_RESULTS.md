# Odd-cubic lift fit, plate vs NACA0012, Re = 500

C_L = c0 + a*alpha + b*alpha^3, alpha in radians, OLS. Primary window 0-8 deg.

| | NACA0012 | plate | plate / NACA |
|---|---:|---:|---:|
| a (per rad), 0–8° | 3.2674 | 4.9519 | 1.516 |
| b (per rad³), 0–8° | -28.2359 | -37.0137 | 1.311 |
| intercept c0 | 8.65959e-05 | -0.00234494 |  |
| C_L/α at 0.5° (per rad) | 3.26892 | 4.82314 | 1.475 |
| b/a | -8.64171 | -7.47464 | 0.865 |
| cubic maximum (deg) | 11.2528 | 12.0995 | 1.075 |
| RMSE 0–8° (fit) | 0.000102739 | 0.00178797 |  |
| max rel. error 0.5–8° | 0.00191105 | 0.0296011 |  |
| RMSE 8.5–12° (extrapolated) | 0.0127989 | 0.0440449 |  |

Window sensitivity (intercept fit, upper limit varied):

| window | NACA a | NACA b | plate a | plate b |
|---|---:|---:|---:|---:|
| 0–4° | 3.2715 | -28.57 | 4.8298 | -22.12 |
| 0–6° | 3.2720 | -28.67 | 4.8574 | -27.45 |
| 0–8° | 3.2674 | -28.24 | 4.9519 | -37.01 |
| 0–10° | 3.2438 | -26.80 | 5.1155 | -47.39 |
