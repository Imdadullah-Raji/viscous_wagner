"""Plate LSB formation and burst, Re = 500, 7.5-8.5 deg in 0.1 deg (oct9_workflow.md).

Reads the steady cases in lsb/plate_fine_re500/steady_aoa<A> (fine wall-normal mesh,
h1 = 0.025/sqrt(Re)) and, for comparison, the production-mesh cases in
plate_case/re500/ at whatever angles in the range exist there.

For each angle: upper-surface C_f zero crossings (C_f = 2 tau_x, tau the NEGATED
wallShearStress), the bubble if one closes upstream of the trailing edge, min C_f, and
the final C_L, C_D.  Classification:
    attached : C_f > 0 everywhere aft of x = 0.02
    closed   : separation followed by reattachment at x < XTE_CLOSE
    burst    : reversed flow reaching the trailing edge (no reattachment before XTE_CLOSE)

Writes figs/plate_sweep/{sweep.csv, sweep.txt, sweep.png, cf_upper.png}.
Run with the research env:  python plate_sweep.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from streamlines_lsb import ROOT, THESIS, latest, read_patch, upper_cf, crossings  # noqa: E402

FINE = ROOT / "lsb" / "plate_fine_re500"
BASE = THESIS / "plate_case" / "re500"
OUT = ROOT / "lsb" / "figs" / "plate_sweep"
XTE_CLOSE = 0.95          # a reversed-flow end aft of this is the TE closure, not reattachment


def angles(root):
    out = []
    for p in root.glob("steady_aoa*"):
        m = re.fullmatch(r"steady_aoa([\d.]+)", p.name)
        if m and p.is_dir() and 7.45 <= float(m.group(1)) <= 8.55:
            out.append((float(m.group(1)), p))
    return sorted(out)


def forces(case):
    """Last C_L, C_D and the relative C_L change over the last 1000 iterations."""
    rows = []
    for f in sorted((case / "postProcessing" / "forceCoeffs").glob("*/*.dat")):
        a = np.loadtxt(f, comments="#", usecols=(0, 2, 3))
        rows.append(np.atleast_2d(a))
    a = np.vstack(rows)
    a = a[np.argsort(a[:, 0], kind="stable")]
    _, keep = np.unique(a[::-1, 0], return_index=True)
    a = a[::-1][keep]
    cl, cd = a[-1, 2], a[-1, 1]
    j = np.searchsorted(a[:, 0], a[-1, 0] - 1000)
    return cl, cd, abs(cl - a[j, 2]) / max(abs(cl), 1e-12), int(a[-1, 0])


def analyse(case):
    d = latest(case)
    Cw = read_patch(d / "C")[:, :2]
    tau = -read_patch(d / "wallShearStress")[:, :2]
    x, cf = upper_cf(Cw, tau)
    cr = crossings(x, cf)
    aft = (x > 0.02) & (x < 0.9)          # the LE bubble region; excludes the TE corner
    seps = [a for a, k in cr if k == "separation"]
    reats = [a for a, k in cr if k == "reattachment" and a < XTE_CLOSE]
    if not seps:
        kind, xs, xr = "attached", np.nan, np.nan
    elif reats and reats[0] > seps[0]:
        kind, xs, xr = "closed", seps[0], reats[0]
    else:
        kind, xs, xr = "burst", seps[0], np.nan
    return dict(kind=kind, x_sep=xs, x_reat=xr, cf_min=float(cf[aft].min()),
                x_cfmin=float(x[aft][np.argmin(cf[aft])]), crossings=cr, x=x, cf=cf,
                time=d.name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    for mesh, root in (("fine", FINE), ("base", BASE)):
        for a, case in angles(root):
            try:
                r = analyse(case)
                r["cl"], r["cd"], r["dcl"], r["iter"] = forces(case)
            except (FileNotFoundError, ValueError, AttributeError) as e:
                print(f"skip {mesh} {a}: {type(e).__name__} {e}")
                continue
            res[(mesh, a)] = r

    lines = [f"{'mesh':5} {'aoa':>5} {'kind':9} {'x_sep':>7} {'x_reat':>7} {'length':>7} "
             f"{'Cf_min':>9} {'x(Cfmin)':>8} {'C_L':>9} {'C_D':>9} {'dCL/1000it':>10} iter"]
    csv = ["mesh,aoa,kind,x_sep,x_reat,length,cf_min,x_cfmin,cl,cd,dcl_rel_1000,iter"]
    for (mesh, a), r in sorted(res.items()):
        L = r["x_reat"] - r["x_sep"]
        lines.append(f"{mesh:5} {a:5.1f} {r['kind']:9} {r['x_sep']:7.3f} {r['x_reat']:7.3f} "
                     f"{L:7.3f} {r['cf_min']:9.5f} {r['x_cfmin']:8.3f} {r['cl']:9.6f} "
                     f"{r['cd']:9.6f} {r['dcl']:10.1e} {r['iter']}")
        csv.append(f"{mesh},{a},{r['kind']},{r['x_sep']},{r['x_reat']},{L},{r['cf_min']},"
                   f"{r['x_cfmin']},{r['cl']},{r['cd']},{r['dcl']},{r['iter']}")

    # Brackets, per mesh: formation = last attached -> first closed; burst = last closed ->
    # first burst.  Formation is also interpolated on min C_f crossing zero.
    for mesh in ("fine", "base"):
        rs = sorted((a, r) for (m, a), r in res.items() if m == mesh)
        if not rs:
            continue
        lines.append("")
        for k0, k1, name in (("attached", "closed", "formation"), ("closed", "burst", "burst")):
            br = [(a0, a1) for (a0, r0), (a1, r1) in zip(rs, rs[1:])
                  if r0["kind"] == k0 and r1["kind"] == k1]
            lines.append(f"{mesh}: {name} bracket " +
                         (", ".join(f"({a0:.1f}, {a1:.1f}]" for a0, a1 in br) or "not in range"))
        for (a0, r0), (a1, r1) in zip(rs, rs[1:]):
            if r0["cf_min"] > 0 >= r1["cf_min"]:
                ac = a0 + r0["cf_min"] * (a1 - a0) / (r0["cf_min"] - r1["cf_min"])
                lines.append(f"{mesh}: min C_f = 0 (linear interp.) at {ac:.3f} deg")

    (OUT / "sweep.txt").write_text("\n".join(lines) + "\n")
    (OUT / "sweep.csv").write_text("\n".join(csv) + "\n")
    print("\n".join(lines))

    # --- figures
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
    for mesh, mk, ls in (("fine", "o", "-"), ("base", "s", "--")):
        rs = sorted((a, r) for (m, a), r in res.items() if m == mesh)
        if not rs:
            continue
        A = np.array([a for a, _ in rs])
        xs = np.array([r["x_sep"] for _, r in rs])
        xr = np.array([r["x_reat"] if r["kind"] == "closed" else
                       (1.0 if r["kind"] == "burst" else np.nan) for _, r in rs])
        ax[0].plot(A, xs, mk + ls, color="#c2410c", mfc="none" if mesh == "base" else None,
                   label=f"separation ({mesh})")
        ax[0].plot(A, xr, mk + ls, color="#2a78d6", mfc="none" if mesh == "base" else None,
                   label=f"reattachment ({mesh}; 1 = reaches TE)")
        ax[1].plot(A, [r["cf_min"] for _, r in rs], mk + ls, color="k",
                   mfc="none" if mesh == "base" else None, label=mesh)
    ax[0].set_ylabel("x/c"); ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=7)
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set_ylabel(r"min upper $C_f$, 0.02<x<0.9"); ax[1].set_xlabel(r"$\alpha$ (deg)")
    ax[1].legend(fontsize=7)
    fig.suptitle("Plate, Re = 500, steady: upper-surface separation")
    fig.tight_layout(); fig.savefig(OUT / "sweep.png", dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    rs = sorted((a, r) for (m, a), r in res.items() if m == "fine")
    cm = plt.get_cmap("viridis")
    for i, (a, r) in enumerate(rs):
        ax.plot(r["x"], r["cf"], color=cm(i / max(len(rs) - 1, 1)), lw=1, label=f"{a:.1f}")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlim(0, 1); ax.set_ylim(-0.01, 0.03)
    ax.set_xlabel("x/c"); ax.set_ylabel(r"upper $C_f$"); ax.legend(fontsize=6, ncol=2, title="α")
    ax.set_title("Plate, Re = 500, fine mesh")
    fig.tight_layout(); fig.savefig(OUT / "cf_upper.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
