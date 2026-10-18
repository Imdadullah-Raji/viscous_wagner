"""Zero-thickness plate, Re = 500: small-angle a1, a3, against the other three bodies.

Same method as analysis/a3_re_scaling_airfoil.py (and the plate's a3_re_scaling.py): odd Legendre
fits on windows 0-A (x = alpha/A, mirrored);
  a1     P1+P3 on 0-1 deg;
  a3     mean of the P1-P7 a3 over windows A in [2.25, 4] deg, std as the error bar;
  a3_01  P1+P3 on 0-1 deg, an independent check.
Data: steady_re500/ (0-6 deg at 0.25 deg, one mesh: tip-study T3, tip cell 5.6e-4c).
Only cases with 6000 iterations and a CONVERGED force trace are used.

  python3 zero_thickness_plate/a3_zero_thickness.py   (from static_airfoil/)
  -> zero_thickness_plate/a3_re500/
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SA = HERE.parent
sys.path.insert(0, str(SA / "analysis"))
from lasso_re1000 import collect                      # noqa: E402
from a3_re_scaling_airfoil import coeffs, PLATEAU     # noqa: E402

OUT = HERE / "a3_re500"          # set in main() from --mesh
MESHES = {"ell": ("zero_thickness_plate/ell_re500", "a3_re500_ell", "elliptic grid E800"),
          "h": ("zero_thickness_plate/steady_re500", "a3_re500", "H-type grid T3 (superseded: tip-cell distortion)")}
OTHERS = {  # Re = 500 rows of the existing a3(Re) tables
    "plate t/c=0.023": SA.parent / "plate_case/analysis/a3_re_scaling/a3_vs_re.csv",
    "NACA0006": SA / "analysis/a3_re_scaling/NACA0006/a3_vs_re.csv",
    "NACA0012": SA / "analysis/a3_re_scaling/NACA0012/a3_vs_re.csv",
}


def main():
    global OUT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", choices=sorted(MESHES), default="ell")
    mesh = ap.parse_args().mesh
    base, outdir, label = MESHES[mesh]
    OUT = HERE / outdir
    OUT.mkdir(exist_ok=True)
    df = collect(500, base)
    bad = df[(df.verdict != "CONVERGED") | (df.iters < 6000)]
    ok = df.drop(bad.index)
    r, T = coeffs(ok[["aoa", "cl"]].reset_index(drop=True))
    df.to_csv(OUT / "lift_curve.csv", index=False)
    T.to_csv(OUT / "windows.csv", index=False)

    rows = [dict(body=f"zero-thickness plate ({mesh})", t_c=0.0, a1=r["a1"], a3=r["a3"], a3_sd=r["a3_sd"], a3_01=r["a3_01"])]
    for name, f in OTHERS.items():
        s = pd.read_csv(f).set_index("Re").loc[500]
        rows.append(dict(body=name, t_c={"plate t/c=0.023": 0.023, "NACA0006": 0.06, "NACA0012": 0.12}[name],
                         a1=s.a1, a3=s.a3, a3_sd=s.a3_sd, a3_01=s.a3_01))
    C = pd.DataFrame(rows)
    C.to_csv(OUT / "compare_bodies.csv", index=False)

    f = lambda v: f"{v:+.5f}"
    lines = [f"Zero-thickness plate, Re = 500, {label}: small-angle coefficients (odd Legendre, x = alpha/A, mirrored)",
             f"a3 = mean P1-P7 over windows {PLATEAU[0]}-{PLATEAU[1]} deg (+- std); a1, a3_01 from P1+P3 on 0-1 deg",
             f"cases: {len(df)}, used {len(ok)}; excluded (not converged or < 6000 iterations): "
             f"{' '.join(f'{a:g}' for a in bad.aoa) or 'none'}", "",
             ok[["aoa", "cl", "cd", "iters", "verdict"]].to_string(index=False, float_format=lambda v: f"{v:.6f}"), "",
             "P1-P7 a3 per window:", T[["A", "n", "a1_P13", "a3_P13", "a3_P1357"]].to_string(index=False, float_format=f), "",
             "Re = 500, all bodies (thin-airfoil theory: a1 = 2*pi = 6.283):",
             C.to_string(index=False, float_format=lambda v: f"{v:.4f}")]
    txt = "\n".join(lines)
    (OUT / "a3_re500.txt").write_text(txt + "\n")
    print(txt)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    th = np.radians(ok.aoa)
    ax[0].plot(ok.aoa, ok.cl, "o", ms=3, color="#c2410c", label="zero-thickness plate")
    ax[0].plot(ok.aoa, r["a1"] * th + r["a3"] * th**3, "-", lw=1, color="k",
               label=rf"$a_1\theta + a_3\theta^3$ ($a_1$={r['a1']:.3f}, $a_3$={r['a3']:.2f})")
    ax[0].plot(ok.aoa, 2 * np.pi * th, ":", color="0.5", label=r"$2\pi\theta$")
    ax[0].set(xlabel="α (deg)", ylabel="$C_L$", title="lift curve, Re = 500")
    ax[0].legend(fontsize=8)
    ax[1].plot(T.A, T.a3_P1357, "o-", ms=3, color="#c2410c")
    ax[1].axvspan(*PLATEAU, color="#fde68a", alpha=0.35, zorder=0)
    ax[1].set(xlabel="window 0–A (deg)", ylabel=r"$a_3$ from P1–P7 (rad$^{-3}$)", title="P1–P7 a3 per window (yellow: plateau)")
    ax[2].errorbar(C.t_c, -C.a3, yerr=C.a3_sd, fmt="o", color="#c2410c", label=r"$|a_3|$")
    for _, row in C.iterrows():
        ax[2].annotate(row.body, (row.t_c, -row.a3), fontsize=7, xytext=(4, 4), textcoords="offset points")
    a2 = ax[2].twinx(); a2.plot(C.t_c, C.a1, "s", color="#2a78d6"); a2.set_ylabel(r"$a_1$ (rad$^{-1}$)", color="#2a78d6")
    ax[2].set(xlabel="thickness t/c", ylabel=r"$|a_3|$ (rad$^{-3}$)", title="all bodies, Re = 500")
    for a_ in ax:
        a_.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(OUT / "a3_re500.png", dpi=160); plt.close(fig)
    deficit_figure()


BODIES = {  # Re = 500 lift curves, folders relative to static_airfoil/
    "NACA0012": ("NACA0012/steady_re500", "#2a78d6"),
    "NACA0006": ("NACA0006/steady_re500", "#16a34a"),
    "plate t/c=0.023": ("../plate_case/re500", "#9333ea"),
    "zero-thickness plate, H grid (T3)": ("zero_thickness_plate/steady_re500", "#f59e0b"),
    "zero-thickness plate, elliptic grid (E800)": ("zero_thickness_plate/ell_re500", "#c2410c"),
}


def deficit_figure():
    """Is the nonlinearity cubic at the origin?  For C_L = a1 th + a3 th^3 + ..., the deficit
    d = a1 th - C_L goes as th^3 and a3_eff = d / th^3 tends to a constant as th -> 0.
    a1 = Richardson limit of C_L/th from the 0.25 and 0.5 deg points.  0.25-deg grid, 0-6 deg."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    lines = ["Deficit d = a1*th - C_L, Re = 500 (a1 from C_L/th at 0.25, 0.5 deg, theta^2-extrapolated)",
             f"{'body':22s} {'a1':>8s} {'exp 1-4deg':>10s}  a3_eff = d/th^3 at 0.5, 1, 2, 4, 6 deg"]
    for name, (base, col) in BODIES.items():
        df = collect(500, base)
        df = df[(df.verdict == "CONVERGED") & (df.aoa <= 6 + 1e-9) & np.isclose((df.aoa * 4) % 1, 0) & (df.aoa > 0)]
        th, cl = np.radians(df.aoa.to_numpy()), df.cl.to_numpy()
        s = cl / th
        a1 = (4 * s[0] - s[1]) / 3
        d = a1 * th - cl
        eff = d / th**3
        k = (df.aoa >= 1) & (df.aoa <= 4)
        p = np.polyfit(np.log(th[k]), np.log(d[k]), 1)[0]
        deg = df.aoa.to_numpy()
        m = deg >= 0.5            # 0.25 deg is used to fix a1
        ax[0].plot(deg[m], eff[m], "o-", ms=3, color=col, label=f"{name}")
        ax[1].loglog(deg[m], d[m], "o", ms=3, color=col, label=f"{name} (slope {p:.2f} on 1–4°)")
        loc = np.gradient(np.log(d[m]), np.log(th[m]))
        ax[2].plot(deg[m][1:-1], loc[1:-1], "o-", ms=3, color=col, label=name)
        at = lambda x: eff[np.argmin(abs(deg - x))]
        lines.append(f"{name:22s} {a1:8.4f} {p:10.3f}  " + "  ".join(f"{at(x):6.2f}" for x in (0.5, 1, 2, 4, 6)))
    g = np.array([0.5, 6])
    ax[1].loglog(g, 3e-6 * (g / 0.5) ** 3, "k--", lw=0.8, label=r"$\propto\theta^3$")
    ax[2].axhline(3, color="k", ls="--", lw=0.8)
    ax[0].set(xlabel="α (deg)", ylabel=r"$(a_1\theta - C_L)/\theta^3$ (rad$^{-3}$)",
              title=r"effective $|a_3|$: flat for a cubic, drifting otherwise")
    ax[1].set(xlabel="α (deg)", ylabel=r"$a_1\theta - C_L$", title="lift deficit")
    ax[2].set(xlabel="α (deg)", ylabel=r"$d\,\ln(a_1\theta - C_L)/d\,\ln\theta$", title="local exponent (3 = cubic)")
    for a_ in ax:
        a_.grid(alpha=0.25, which="both"); a_.legend(fontsize=7)
    fig.suptitle("Re = 500: is the lift nonlinearity cubic at the origin?")
    fig.tight_layout(); fig.savefig(OUT / "deficit_four_bodies.png", dpi=160); plt.close(fig)
    txt = "\n".join(lines)
    (OUT / "deficit_four_bodies.txt").write_text(txt + "\n")
    print("\n" + txt)


if __name__ == "__main__":
    main()
