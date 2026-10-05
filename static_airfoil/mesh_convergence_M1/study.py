"""Benchmark 3+3 then 6 ranks; solve five M1 cases; compare steady forces.

Only writes under this study directory. Original cases are not modified.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
ANGLES = (4, 6, 8, 10, 12)
TOL = 1e-5


def command(args, log):
    with Path(log).open("w") as stream:
        subprocess.run([str(a) for a in args], stdout=stream,
                       stderr=subprocess.STDOUT, check=True)


def set_entry(case, file, entry, value):
    result = subprocess.run(["foamDictionary", str(case / "system" / file),
                             "-entry", entry, "-set", str(value)],
                            capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)


def physical_cpus():
    seen, cpus = set(), []
    for cpu in sorted(os.sched_getaffinity(0)):
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        key = tuple((topology / name).read_text().strip()
                    for name in ("physical_package_id", "core_id"))
        if key not in seen:
            seen.add(key)
            cpus.append(cpu)
    if len(cpus) != 6:
        raise RuntimeError(f"Expected 6 available physical cores, found {cpus}")
    return cpus


def mpi(case, cpus, log):
    # Disjoint Linux CPU affinity sets prevent the two MPI jobs from sharing
    # physical cores; no hyperthread sibling is included. MPI inherits affinity.
    start = time.monotonic()
    command(["taskset", "-c", ",".join(map(str, cpus)), "mpirun",
             "--bind-to", "none", "-np", len(cpus), "foamRun",
             "-case", case, "-parallel"], log)
    if not re.search(r"^End\s*$", Path(log).read_text(), re.M):
        raise RuntimeError(f"Solver did not finish normally: {log}")
    return time.monotonic() - start


def forces(case):
    records = {}
    files = sorted(case.glob("postProcessing/forceCoeffs/*/*.dat"),
                   key=lambda p: (float(p.parent.name), p.name))
    for path in files:
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = [float(v) for v in line.split()]
            records[values[0]] = (values[2], values[3])  # later restart wins
    return np.array([(t, *records[t]) for t in sorted(records)]).reshape(-1, 3)


def stats(case):
    a = forces(case)
    if not len(a):
        return {"status": "NOT_RUN"}
    if not np.isfinite(a).all():
        raise RuntimeError(f"Non-finite forces in {case}")
    tail = a[a[:, 0] >= a[-1, 0] - 1000, 1:]
    prev = a[(a[:, 0] >= a[-1, 0] - 2000) & (a[:, 0] < a[-1, 0] - 1000), 1:]
    spread = np.ptp(tail, axis=0)
    mean = tail.mean(axis=0)
    drift = np.abs(mean - prev.mean(axis=0)) if len(prev) else np.full(2, np.inf)
    threshold = TOL * np.maximum(np.abs(mean), 0.01)
    converged = (a[-1, 0] >= 6000 and len(tail) >= 1000 and len(prev) >= 1000
                 and np.all(spread <= threshold) and np.all(drift <= threshold))
    return dict(status="CONVERGED" if converged else "NOT_CONVERGED",
                iteration=float(a[-1, 0]), cd=float(a[-1, 1]), cl=float(a[-1, 2]),
                cd_range=float(spread[0]), cl_range=float(spread[1]),
                cd_drift=float(drift[0]), cl_drift=float(drift[1]))


def cells(case):
    return int(re.search(r"cells:\s+(\d+)", (case / "log.checkMesh").read_text())[1])


def report():
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for angle in ANGLES:
        fine = ROOT / f"steady_aoa{angle}"
        coarse = ROOT / "baseline" / fine.name
        c, f = stats(coarse), stats(fine)
        row = {"angle_deg": angle, "coarse_cells": cells(coarse), "fine_cells": cells(fine)}
        row.update({f"coarse_{key}": value for key, value in c.items()})
        row.update({f"fine_{key}": value for key, value in f.items()})
        for coefficient in ("cl", "cd"):
            if coefficient in f:
                row[f"delta_{coefficient}_pct"] = 100 * (f[coefficient]-c[coefficient]) / abs(c[coefficient])
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "comparison.csv", index=False)
    lines = ["# M1 steady mesh-refinement comparison", "",
             "Re = 500. Baseline: 63,360 cells; M1: 110,000 cells. Same steady equations and spatial schemes.",
             "This is a two-mesh sensitivity study, not a three-grid Richardson/GCI assessment.", "",
             "Convergence requires at least 6,000 iterations, and both Cl and Cd must satisfy:",
             "last-1,000-iteration range and mean shift from the preceding 1,000 iterations",
             "each ≤ 1e-5 × max(|tail mean|, 0.01).", "",
             "| Angle | M1 status | Baseline Cl | M1 Cl | ΔCl (%) | Baseline Cd | M1 Cd | ΔCd (%) |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        def value(key):
            return f"{row[key]:.8g}" if key in row else "—"
        lines.append(f"| {row['angle_deg']} | {row['fine_status']} | " + " | ".join(
            value(k) for k in ("coarse_cl", "fine_cl", "delta_cl_pct", "coarse_cd", "fine_cd", "delta_cd_pct")) + " |")
    lines += ["", "Differences from unconverged cases are provisional. Inspect comparison.csv for force ranges and drift.",
              "Physical stability and steady-solver convergence are distinct; no shedding-onset conclusion follows from this comparison."]
    (ROOT / "RESULTS.md").write_text("\n".join(lines) + "\n")
    if not any("fine_cl" in row for row in rows):
        return
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for j, coefficient in enumerate(("cl", "cd")):
        for prefix, marker in (("coarse", "o"), ("fine", "s")):
            axes[0, j].plot(df.angle_deg, df[f"{prefix}_{coefficient}"], marker + "-", label=prefix)
        axes[0, j].set_ylabel(coefficient.upper())
        axes[0, j].legend()
        axes[1, j].plot(df.angle_deg, df[f"delta_{coefficient}_pct"], "o-")
        axes[1, j].axhline(0, color="black", lw=0.7)
        axes[1, j].set(xlabel="Angle (degrees)", ylabel="(M1 − baseline) / |baseline| (%)")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "mesh_comparison.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for angle in ANGLES:
        a = forces(ROOT / f"steady_aoa{angle}")
        if len(a):
            for ax, column in zip(axes, (2, 1)):
                mask = a[:, 0] >= 1000
                ax.plot(a[mask, 0], a[mask, column], label=f"{angle}°")
    axes[0].set_ylabel("Cl")
    axes[1].set(ylabel="Cd", xlabel="SIMPLE iteration (startup below 1,000 omitted)")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend()
    fig.tight_layout()
    fig.savefig(ROOT / "force_convergence.png", dpi=180)
    plt.close(fig)


def benchmark(cpus):
    destination = ROOT / "benchmark.json"
    if destination.exists():
        return json.loads(destination.read_text())
    results = {}
    for ranks in (3, 6):  # requested ordering: 3+3 first, then 6
        jobs = []
        for angle in (4, 6):
            case = ROOT / "benchmarks" / f"aoa{angle}_np{ranks}"
            case.mkdir(parents=True, exist_ok=True)
            for name in ("0", "constant", "system"):
                shutil.copytree(ROOT / f"steady_aoa{angle}" / name, case / name, dirs_exist_ok=True)
            set_entry(case, "decomposeParDict", "numberOfSubdomains", ranks)
            set_entry(case, "controlDict", "startFrom", "startTime")
            set_entry(case, "controlDict", "startTime", 0)
            set_entry(case, "controlDict", "endTime", 1000)
            set_entry(case, "controlDict", "writeInterval", 1000)
            command(["decomposePar", "-case", case, "-force"], case / "log.decomposePar")
            jobs.append(case)
        print(f"Benchmark starting: {'3+3 concurrent' if ranks == 3 else '6 sequential'}", flush=True)
        start = time.monotonic()
        if ranks == 3:
            with ThreadPoolExecutor(max_workers=2) as pool:
                futures = [pool.submit(mpi, case, cpus[3*i:3*i+3], case / "log.benchmark")
                           for i, case in enumerate(jobs)]
                durations = [future.result() for future in futures]
        else:
            durations = [mpi(case, cpus, case / "log.benchmark") for case in jobs]
        elapsed = time.monotonic() - start
        results[str(ranks)] = {"pair_seconds": elapsed, "individual_seconds": durations,
                               "aggregate_iterations_per_second": 2000 / elapsed}
        print(f"Benchmark {ranks}: {results[str(ranks)]}", flush=True)
    results["selected_ranks"] = min((3, 6), key=lambda n: results[str(n)]["pair_seconds"])
    results["physical_cpu_ids"] = cpus
    results["note"] = "Both schedules run the same cold-start 4° and 6° cases for 1,000 iterations each; wall times include solver startup and final write. Decomposition excluded."
    destination.write_text(json.dumps(results, indent=2) + "\n")
    return results


def solve(angle, cpus):
    case = ROOT / f"steady_aoa{angle}"
    current = stats(case)
    if current["status"] == "CONVERGED":
        return current
    existing = list(case.glob("processor[0-9]*"))
    if existing and len(existing) != len(cpus):
        raise RuntimeError(f"Rank count changed for {case}; preserve/reconstruct restart before changing decomposition")
    set_entry(case, "decomposeParDict", "numberOfSubdomains", len(cpus))
    if not existing:
        command(["decomposePar", "-case", case], case / "log.decomposePar")
    checkpoints = [int(p.name) for p in (case / "processor0").iterdir()
                   if p.is_dir() and p.name.isdigit()]
    latest = max(checkpoints, default=0)
    set_entry(case, "controlDict", "startFrom", "latestTime")
    set_entry(case, "controlDict", "writePrecision", 12)
    set_entry(case, "controlDict", "writeInterval", 2000)
    set_entry(case, "controlDict", "purgeWrite", 2)
    for end in range(max(6000, (latest // 2000 + 1) * 2000), 18001, 2000):
        set_entry(case, "controlDict", "endTime", end)
        print(f"AoA {angle}: {len(cpus)} ranks, CPUs {cpus}, target {end}", flush=True)
        duration = mpi(case, cpus, case / f"log.steady_to_{end}")
        current = stats(case)
        print(f"AoA {angle}: {current}; stage wall time {duration:.2f}s", flush=True)
        if current["status"] == "CONVERGED":
            break
    command(["reconstructPar", "-case", case, "-latestTime"], case / "log.reconstructPar")
    (case / "convergence.json").write_text(json.dumps(current, indent=2) + "\n")
    return current


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    if args.report_only:
        report()
        return
    with (ROOT / ".run.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if subprocess.run(["pgrep", "-x", "foamRun"], capture_output=True).returncode == 0:
            raise RuntimeError("Another foamRun is active; the study requires all six physical cores")
        cpus = physical_cpus()
        command(["mpirun", "--bind-to", "none", "-np", 2, "/bin/hostname"], ROOT / "log.mpi_preflight")
        result = benchmark(cpus)
        ranks = result["selected_ranks"]
        print(f"Selected {ranks} ranks per case from measured pair throughput", flush=True)
        if ranks == 3:
            with ThreadPoolExecutor(max_workers=2) as pool:
                for pair in ((4, 6), (8, 10), (12,)):
                    futures = [pool.submit(solve, angle, cpus[3*i:3*i+3]) for i, angle in enumerate(pair)]
                    for future in futures:
                        future.result()
                    report()
        else:
            # The prepared AoA 4 decomposition is disposable only while at zero.
            case = ROOT / "steady_aoa4"
            if len(list(case.glob("processor[0-9]*"))) == 3 and not len(forces(case)):
                for part in case.glob("processor[0-9]*"):
                    shutil.rmtree(part)
            for angle in ANGLES:
                solve(angle, cpus)
                report()
        report()
        if any(stats(ROOT / f"steady_aoa{a}")["status"] != "CONVERGED" for a in ANGLES):
            raise RuntimeError("Some cases did not converge within 18,000 iterations; see RESULTS.md")
        print("Study complete: all five M1 cases meet force convergence criteria", flush=True)


if __name__ == "__main__":
    main()
