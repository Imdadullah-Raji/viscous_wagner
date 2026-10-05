"""Steady (SIMPLE) C_L sweep for the thick flat plate.  Usage: run_sweep.py --re 1000

Cases are built by static_airfoil/make_steady.sh, so equations, schemes, BCs,
force normalisation and the forceCoeffs setup are identical to the airfoil
steady cases; only the mesh differs.  AoA set matches static_airfoil/lift_coeffs.csv.

Steps:
  1. build every steady_aoa<A> case that does not exist yet
  2. benchmark 2+2+2 vs 3+3 vs 6 ranks (cold start, 1000 iterations, pinned to
     the six physical cores, decomposition excluded) -> benchmark.json
  3. run all angles with the fastest layout, one job per CPU slot pulled from a
     queue; convergence gate as in mesh_convergence_M1/study.py
  4. lift_coeffs.csv + RESULTS.md

Re-running resumes: converged cases are skipped, benchmark.json is reused.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
AIRFOIL = ROOT.parent / "static_airfoil"
# Set by configure(): everything for one Reynolds number lives in re<RE>/,
# including its own mesh (h1 = 0.05/sqrt(Re), so meshes are Re-specific).
RE = WORK = MESH = None


def configure(re):
    global RE, WORK, MESH
    RE = int(re) if float(re).is_integer() else re
    WORK = ROOT / f"re{RE}"
    MESH = WORK / "meshes" / "plate_ns160_nr180.msh"
ANGLES = [a / 2 for a in range(0, 21)] + [11.0, 12.0]
TOL = 1e-5
BENCH_ITERS = 1000


def name(angle):
    return f"{angle:g}"


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
        key = tuple((topology / n).read_text().strip()
                    for n in ("physical_package_id", "core_id"))
        if key not in seen:
            seen.add(key)
            cpus.append(cpu)
    if len(cpus) != 6:
        raise RuntimeError(f"Expected 6 available physical cores, found {cpus}")
    return cpus


def mpi(case, cpus, log):
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
        return {"status": "DIVERGED", "iteration": float(a[-1, 0])}
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


def build(angle):
    case = WORK / f"steady_aoa{name(angle)}"
    if (case / "constant" / "polyMesh" / "boundary").exists():
        return case
    command([AIRFOIL / "make_steady.sh", case, RE, name(angle), MESH], WORK / f"{case.name}.make.log")
    if "Mesh OK" not in (case / "log.checkMesh").read_text():
        raise RuntimeError(f"checkMesh failed for {case}")
    return case


def benchmark(cpus):
    destination = WORK / "benchmark.json"
    if destination.exists():
        return json.loads(destination.read_text())
    results = {}
    for ranks in (2, 3, 6):
        n_jobs = 6 // ranks
        jobs = []
        for i in range(n_jobs):
            src = WORK / "steady_aoa4"
            case = WORK / "benchmarks" / f"np{ranks}_job{i}"
            if case.exists():
                shutil.rmtree(case)
            case.mkdir(parents=True)
            for d in ("0", "constant", "system"):
                shutil.copytree(src / d, case / d)
            set_entry(case, "decomposeParDict", "numberOfSubdomains", ranks)
            set_entry(case, "controlDict", "endTime", BENCH_ITERS)
            set_entry(case, "controlDict", "writeInterval", BENCH_ITERS)
            command(["decomposePar", "-case", case, "-force"], case / "log.decomposePar")
            jobs.append(case)
        label = "+".join([str(ranks)] * n_jobs)
        print(f"Benchmark {label}: {n_jobs} concurrent job(s) x {BENCH_ITERS} iterations", flush=True)
        start = time.monotonic()
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            futures = [pool.submit(mpi, c, cpus[ranks*i:ranks*(i+1)], c / "log.benchmark")
                       for i, c in enumerate(jobs)]
            durations = [f.result() for f in futures]
        elapsed = time.monotonic() - start
        results[label] = {"ranks": ranks, "wall_seconds": elapsed, "job_seconds": durations,
                          "aggregate_iterations_per_second": n_jobs * BENCH_ITERS / elapsed}
        print(f"Benchmark {label}: {results[label]}", flush=True)
    best = max((k for k in results), key=lambda k: results[k]["aggregate_iterations_per_second"])
    results["selected"] = best
    results["selected_ranks"] = results[best]["ranks"]
    results["physical_cpu_ids"] = cpus
    results["note"] = (f"Cold-start AoA 4 copies, {BENCH_ITERS} SIMPLE iterations per job; "
                       "metric = total iterations / wall time with all jobs concurrent. "
                       "Decomposition excluded, solver startup and final write included.")
    destination.write_text(json.dumps(results, indent=2) + "\n")
    return results


def solve(angle, cpus):
    case = WORK / f"steady_aoa{name(angle)}"
    current = stats(case)
    if current["status"] in ("CONVERGED", "DIVERGED"):
        return current
    existing = list(case.glob("processor[0-9]*"))
    if existing and len(existing) != len(cpus):
        raise RuntimeError(f"Rank count changed for {case}; reconstruct before changing decomposition")
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
        print(f"{time.strftime('%H:%M:%S')} AoA {name(angle)}: {len(cpus)} ranks, CPUs {cpus}, target {end}", flush=True)
        try:
            duration = mpi(case, cpus, case / f"log.steady_to_{end}")
        except (subprocess.CalledProcessError, RuntimeError) as err:
            current = stats(case) | {"status": "SOLVER_FAILED", "error": str(err)}
            print(f"AoA {name(angle)}: {current}", flush=True)
            break
        current = stats(case)
        print(f"{time.strftime('%H:%M:%S')} AoA {name(angle)}: {current}; stage wall time {duration:.1f}s", flush=True)
        if current["status"] in ("CONVERGED", "DIVERGED"):
            break
    command(["reconstructPar", "-case", case, "-latestTime"], case / "log.reconstructPar")
    (case / "convergence.json").write_text(json.dumps(current, indent=2) + "\n")
    return current


def report():
    rows = [(a, stats(WORK / f"steady_aoa{name(a)}")) for a in ANGLES]
    with (WORK / "lift_coeffs.csv").open("w") as f:
        f.write("alpha cl\n")
        for a, s in rows:
            if "cl" in s:
                f.write(f"{a} {s['cl']:.6e}\n")
    lines = [f"# Plate steady sweep, Re = {RE}", "",
             f"Mesh `{MESH.relative_to(WORK)}` (63,360 hexes). Cases built by `static_airfoil/make_steady.sh`.",
             "Convergence: >= 6000 iterations; Cl and Cd range over the last 1000 iterations and mean shift",
             "from the preceding 1000 each <= 1e-5 * max(|mean|, 0.01); extended to at most 18000.", "",
             "| AoA | status | iter | Cl | Cd | Cl range | Cl drift |", "|---:|---|---:|---:|---:|---:|---:|"]
    for a, s in rows:
        g = lambda k: f"{s[k]:.6g}" if k in s else "—"
        lines.append(f"| {name(a)} | {s['status']} | {g('iteration')} | {g('cl')} | {g('cd')} | "
                     f"{g('cl_range')} | {g('cl_drift')} |")
    (WORK / "RESULTS.md").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--re", type=float, required=True)
    parser.add_argument("--ranks", type=int, choices=(2, 3, 6),
                        help="skip the benchmark and use this many ranks per job")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--build-only", action="store_true",
                        help="build every case (mesh conversion, checkMesh) and stop; for "
                             "handing cases to other machines")
    parser.add_argument("--step", type=float,
                        help="angle grid 0..--max in this step instead of the default "
                             "0..10 by 0.5 plus 11, 12 (lsb_hypothesis.md: 0.25)")
    parser.add_argument("--max", type=float, default=12.0)
    parser.add_argument("--skip", type=float, nargs="+", default=[],
                        help="angles solved elsewhere (e.g. on the remote machine): kept in the "
                             "report, never queued here")
    args = parser.parse_args()
    configure(args.re)
    if args.step:
        global ANGLES
        n = int(round(args.max / args.step))
        ANGLES = [round(i * args.step, 6) for i in range(n + 1)]
    if args.report_only:
        report()
        return
    with (WORK / ".run.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if subprocess.run(["pgrep", "-x", "foamRun"], capture_output=True).returncode == 0:
            raise RuntimeError("Another foamRun is active; the sweep requires all six physical cores")
        cpus = physical_cpus()
        for a in ANGLES:
            build(a)
        print(f"Built {len(ANGLES)} cases", flush=True)
        if args.build_only:
            return
        ranks = args.ranks or benchmark(cpus)["selected_ranks"]
        print(f"Selected {ranks} ranks per job ({6 // ranks} concurrent)", flush=True)

        # Hardest (highest AoA) first so the long tails do not run alone at the end.
        todo = queue.Queue()
        skip = {round(a, 6) for a in args.skip}
        for a in sorted(ANGLES, reverse=True):
            if round(a, 6) not in skip:
                todo.put(a)

        def worker(slot):
            while True:
                try:
                    a = todo.get_nowait()
                except queue.Empty:
                    return
                solve(a, slot)
                report()

        slots = [cpus[ranks*i:ranks*(i+1)] for i in range(6 // ranks)]
        with ThreadPoolExecutor(max_workers=len(slots)) as pool:
            for f in [pool.submit(worker, s) for s in slots]:
                f.result()
        report()
        bad = [name(a) for a in ANGLES if stats(WORK / f"steady_aoa{name(a)}")["status"] != "CONVERGED"]
        print("Sweep complete" + (f"; NOT converged: {bad}" if bad else ": all cases converged"), flush=True)


if __name__ == "__main__":
    main()
