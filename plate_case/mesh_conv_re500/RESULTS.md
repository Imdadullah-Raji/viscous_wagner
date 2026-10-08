# Mesh convergence, Re = 500: plate and NACA0006 (oct8_workflow.md)

2026-10-08. Each section is run on a finer mesh at 0°, 2°, 6° and 8° with steady SIMPLE
(`make_steady.sh`, 6000 iterations) and compared with the production mesh.

## Meshes

| | base (production) | fine |
|---|---|---|
| first cell h1 | 0.05/√Re = 2.24e-3 | **0.025/√Re = 1.12e-3** |
| radial cells | 180 (radial growth 1.0347) | **260 (radial growth 1.0254)** |
| cells inside the Blasius layer, 5/√Re = 0.224 | 44 | **72** |
| tangential cells | ns 160, nte 32 | unchanged |
| total cells | 63,360 | 91,520 |

- **Plate:** `meshes/plate_ns160_nr260_n195_h025.msh`, built with
  `plate_mesher.py --re 500 --nr 260 --n1 95 --h1-coef 0.025`. The inner ring has 95 cells so its
  offset stays at 0.43c, against 0.434c on the base mesh.
- **NACA0006:** `static_airfoil/meshes/M2f_ns160_nr260_b025_h025_naca0006_re500.msh`, built with
  `mesher.py --naca 0006 --re 500 -r 30 --nr 260 --ns 160 --nte 32 --bump 0.25 --h1-coef 0.025`.
- Both meshers gained `--h1-coef` (default 0.05). With the defaults, both reproduce the
  production meshes byte for byte.
- checkMesh passes on all 8 cases. Plate fine mesh: non-orthogonality 54°, skewness 1.52.
  NACA0006 fine mesh: non-orthogonality 68.3°, skewness 0.73.
- Only the wall-normal resolution was refined, as the instructions specify. Tangential
  resolution is not tested.

## Results

Every fine case is CONVERGED on its force trace (`steady_summary.py`). All base cases were
already converged.

| section | α | C_L base | C_L fine | ΔC_L | C_D base | C_D fine | ΔC_D |
|---|---:|---:|---:|---:|---:|---:|---:|
| plate | 0 | 0 | 0 | – | 0.153513 | 0.153581 | +0.044 % |
| plate | 2 | 0.167585 | 0.167577 | −0.005 % | 0.155056 | 0.155130 | +0.048 % |
| plate | 6 | 0.476049 | 0.476031 | −0.004 % | 0.168616 | 0.168721 | +0.062 % |
| plate | 8 | 0.584495 | 0.584821 | +0.056 % | 0.182149 | 0.182254 | +0.058 % |
| NACA0006 | 0 | 0 | 0 | – | 0.156942 | 0.156946 | +0.003 % |
| NACA0006 | 2 | 0.146938 | 0.146956 | +0.012 % | 0.158346 | 0.158349 | +0.002 % |
| NACA0006 | 6 | 0.411298 | 0.411449 | +0.037 % | 0.170259 | 0.170257 | −0.001 % |
| NACA0006 | 8 | 0.506829 | 0.507261 | +0.085 % | 0.181514 | 0.181500 | −0.008 % |

- **Forces are converged in wall-normal resolution to < 0.1 %** for both sections at every angle.
- **The change grows with angle.** For C_L it rises from ≤ 0.012 % at 2° to 0.06–0.09 % at 8°,
  which is where the separated regions are. These differences are far below the fit residuals
  that matter. The odd-cubic domain results (`static_airfoil/STEADY_NACA0006_DOMAIN.md`) rest on
  differences of 1e-3 and larger, while the largest change here is 4.3e-4 (NACA0006 at 8°).

### Separated regions at 8° (upper-surface C_f zero crossings, `static_airfoil/lsb/`)

| section | mesh | separation x/c | reattachment x/c | length | min C_f |
|---|---|---:|---:|---:|---:|
| plate | base | 0.137 | 0.424 | 0.288 | −0.00083 |
| plate | fine | 0.143 | 0.401 | 0.258 | −0.00068 |
| NACA0006 | base | 0.832 | 0.972 | 0.139 | −0.00088 |
| NACA0006 | fine | 0.836 | 0.971 | 0.135 | −0.00084 |

- **Both bubbles survive refinement.**
- **The NACA0006's trailing-edge bubble is converged:** its length changes by 3 % and its
  minimum C_f by 5 %.
- **The plate's leading-edge bubble is not converged:** it is 10 % shorter and 18 % weaker on
  the fine mesh, while C_L changes by only 0.06 %. A weak bubble near its formation angle is
  sensitive to the near-wall grid. **The LSB formation and burst study should use the fine plate
  mesh**, or at least check against it, because the formation angle may shift with resolution.

## Files

- Plate cases are in `steady_aoa{0,2,6,8}/`. NACA0006 cases are in
  `static_airfoil/steady_naca0006_re500_fine/`.
- The driver is `run_after_unsteady.sh` and the log is `run_steady.log`.
- Base values are from `../re500/RESULTS.md` and `static_airfoil/steady_naca0006_re500/run_steady.log`.
