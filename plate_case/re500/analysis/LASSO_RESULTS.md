# Degree-9 Lasso, plate vs NACA0012 (Re = 500, 0–12°, θ in rad)

PolynomialFeatures(9) → StandardScaler → Lasso(alpha=0.0001); same pipeline as static_airfoil/analysis/lasso_reg.py.

| term | NACA0012 standardized | NACA0012 raw | plate standardized | plate raw |
|---|---:|---:|---:|---:|
| intercept | | +0.00207083 | | +0.000423183 |
| θ^1 | +0.192151 | +3.21236 | +0.287588 | +4.80788 |
| θ^2 | 0 | 0 | 0 | 0 |
| θ^3 | -0.0641861 | -26.0362 | 0 | 0 |
| θ^4 | 0 | 0 | -0.117712 | -240.891 |
| θ^5 | 0 | 0 | 0 | 0 |
| θ^6 | 0 | 0 | 0 | 0 |
| θ^7 | 0 | 0 | 0 | 0 |
| θ^8 | 0 | 0 | 0 | 0 |
| θ^9 | +0.00338461 | +19915.8 | +0.0162792 | +95790.8 |

| | NACA0012 | plate |
|---|---:|---:|
| nonzero terms at λ=1e-4 | 3 | 3 |
| training RMSE | 0.00105 | 0.00169 |

Path note: θ^1 enters first for both; after that the high powers SWAP rather than
accumulate (θ^9 is replaced by θ^8, then θ^7, ...) because the powers of θ over
0–0.21 rad are nearly collinear.  The selected set at a given λ is therefore not a
unique physical model; see lasso_path_*.csv for the knots.
