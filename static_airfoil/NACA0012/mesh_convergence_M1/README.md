# Re = 500 steady mesh-refinement study

## Status

Five M1 cases are prepared for **4°, 6°, 8°, 10°, and 12°**. All pass `checkMesh`.
No parallel solver has started: this Codex session's sandbox denies MPI/PMIx socket creation
(`socket() failed with errno=1`; PMIx listener cannot start).

## Meshes and preserved settings

- Baseline cases: `../steady_aoa{4,6,8,10,12}`, **63,360 cells**.
- Fine mesh: `../meshes/M1_ns200_nr250_b025.msh`, **110,000 cells**, not 1.2 million.
- M1 maximum non-orthogonality: 64.232°; maximum skewness: 1.05303; `Mesh OK`.
- Cell-count ratio: 1.7361; nominal 2-D spacing ratio `sqrt(110000/63360) ≈ 1.3176`.
  This is not a uniformly scaled mesh sequence; the nominal ratio alone does not establish an order of accuracy.
- Same Re, geometry, domain, boundary conditions, force normalization, steady equations,
  spatial schemes, solver tolerances and under-relaxation as the existing steady cases.
- Uniform-flow initialization, following the original cases. No cross-mesh field mapping.
- Baseline post-processing, dictionaries and mesh-check logs are copied into `baseline/`.
  The original cases remain untouched.
- Production writes restart fields every 2,000 iterations, keeps two decomposed checkpoints,
  and uses 12-digit field output. These output settings do not change the discrete equations.

This is a **two-mesh sensitivity comparison**. It cannot independently determine observed
order or a three-grid GCI. It concerns steady-branch lift and drag, not shedding means or
mesh independence of the Hopf onset.

## Run

From a normal local terminal in the workspace, with no other OpenFOAM solver running:

```bash
bash mesh_convergence_M1/run_study.sh > mesh_convergence_M1/run.log 2>&1
```

For a run that survives closing the terminal:

```bash
nohup bash mesh_convergence_M1/run_study.sh > mesh_convergence_M1/run.log 2>&1 < /dev/null &
```

The runner sources OpenFOAM-dev and uses the `research` Python environment. It checks MPI
before starting the benchmark and refuses to launch while another `foamRun` is visible.

## Requested core scheduling

1. **3+3 first:** cold-start copies of the 4° and 6° cases each run 1,000 iterations
   concurrently on three physical cores. They use disjoint CPU-affinity sets.
2. **Then 6:** fresh copies of those same two cases run sequentially on all six physical
   cores for the same 1,000 iterations. This compares the same work and initial states.
3. Select the schedule with the smaller total wall time for the pair. Save measurements
   and the choice in `benchmark.json`. Decomposition is excluded from these timings;
   solver startup and final output are included.
4. Run all five production cases using that schedule, with at most six MPI ranks active.
   With 3-rank jobs the final unpaired case uses three cores.

On this PC, physical-core representatives are Linux CPUs `0,2,4,6,8,10`. The runner
discovers them from topology and affinity rather than assuming 12 logical CPUs means
12 physical cores. Single-threaded BLAS/OpenMP settings prevent extra worker threads.

Benchmark cases live under `benchmarks/`, separate from production cases.
Benchmark results are reused if restarting an interrupted study. A directory lock prevents
two copies of this driver from running concurrently.

## Convergence and comparison

Each production case runs at least 6,000 SIMPLE iterations. Both Cl and Cd must have:

- A range over the last 1,000 iterations below `1e-5 * max(abs(tail mean), 0.01)`.
- A change in mean from the preceding 1,000-iteration window below the same threshold.

If needed, the runner extends by 2,000 iterations up to 18,000. Failure to satisfy these
checks is reported, not silently classified as convergence. The final fields are reconstructed.
The force loader sorts restart segments numerically and keeps the later record at duplicate
iteration numbers.

Outputs:

- `comparison.csv`: baseline and M1 forces, convergence status, ranges, drift and percentage changes.
- `RESULTS.md`: readable comparison, updated after completed cases/pairs.
- `mesh_comparison.png`: lift/drag curves and percentage differences, once results exist.
- `force_convergence.png`: M1 force histories after iteration 1,000.
- `steady_aoa*/log.steady_to_*`: solver logs.
- `steady_aoa*/convergence.json`: final convergence diagnostics.

Refresh the report without launching solvers:

```bash
bash mesh_convergence_M1/run_study.sh --report-only
```

At 12° a converged steady solution can be physically unstable. The steady force comparison
must not be interpreted as a test for whether the actual flow sheds vortices.
