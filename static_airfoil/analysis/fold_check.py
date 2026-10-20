"""Fold test, zero-thickness plate, Re 500: the steady chains beyond 6 deg against the models fitted on
the re-converged 0-4.25 deg data.

Models (all fitted on tight 0.25-4.25 deg C_L): odd polynomials to theta^9 and theta^11; the rational
fit [3/2] (pole near 10.8 deg, CUBIC_STRESS_TEST.md section 8); the square-root branch fit (branch point
near 8.9 deg, section 10).  Chain data: zero_thickness_plate/fold/re500_{asc,desc}, a case counted once
it ended and its geometrically extrapolated remainder is below 1e-5 (production tolerances, so ~1e-5 is
the precision; differences below that mean nothing).

  python analysis/fold_check.py     (research env)  -> analysis/fold_check/
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lasso_re1000 import collect, load                              # noqa: E402
from cubic_stress_test import fit_window                            # noqa: E402
from drag_even_fit import cases                                      # noqa: E402
from residue_test import sqrt_model                                  # noqa: E402

warnings.filterwarnings("ignore")
OUT = HERE / "fold_check"
FOLD = HERE.parent / "zero_thickness_plate/fold"
COL = {"poly θ⁹": "#2a78d6", "poly θ¹¹": "#eb6834", "rational [3/2]": "#1baf7a", "√ branch": "#eda100"}
LS = {"poly θ⁹": "-", "poly θ¹¹": "--", "rational [3/2]": "-.", "√ branch": ":"}


def chain(name):
    """Converged chain cases: alpha -> (C_L, C_D, remainder estimate)."""
    out = {}
    for c in sorted((FOLD / name).glob("steady_aoa*")):
        lg = c / "log.steady"
        if not c.is_dir() or not lg.exists() or b"\nEnd" not in lg.read_bytes():
            continue
        try:
            it, cd, cl = load(str(c))
        except Exception:
            continue
        it = it.astype(int); i = {k: j for j, k in enumerate(it)}; tN = it[-1]
        if tN < 3000:
            continue
        v = np.array([cl[i[k]] for k in range(tN - 1500, tN + 1, 500)]); d = np.diff(v)
        q = d[-1] / d[-2] if d[-2] else 0
        rem = abs(d[-1] * q / (1 - q)) if 0 < q < 0.95 else abs(d[-1])
        out[float(c.name[10:])] = (cl[-1], cd[-1], rem)
    return out


def models():
    df = cases(500); d = df[df.ok_l & (df.aoa > 0)]
    th, cl = np.radians(d.aoa.to_numpy()), d.cl.to_numpy()
    m = {}
    for degs, name in (((1, 3, 5, 7, 9), "poly θ⁹"), ((1, 3, 5, 7, 9, 11), "poly θ¹¹")):
        a, *_ = fit_window(d, d.aoa.max(), degs)
        m[name] = (lambda a: lambda t: sum(a[k] * t ** k for k in a))(a)
    s = th.max(); z = (th / s) ** 2; y = cl / th; q = np.zeros(2); L, M = 3, 2
    for _ in range(40):
        Qp = 1 + sum(q[k] * z ** (k + 1) for k in range(M))
        A = np.column_stack([z ** k for k in range(L + 1)] + [-y * z ** (k + 1) for k in range(M)])
        W = (th / np.abs(Qp))[:, None]; sol, *_ = np.linalg.lstsq(A * W, y * W[:, 0], rcond=None); q = sol[L + 1:]
    p = sol[:L + 1]
    m["rational [3/2]"] = lambda t: t * np.polyval(p[::-1], (t / s) ** 2) / (1 + sum(q[k] * ((t / s) ** 2) ** (k + 1) for k in range(M)))
    dd = df[df.ok_l & df.ok_d & (df.aoa > 0)]
    sm = sqrt_model(np.radians(dd.aoa.to_numpy()), dd.cl.to_numpy(), dd.cd.to_numpy(), 8e-9, 6e-10)
    def sq(t):
        t = np.atleast_1d(t); z = t * t; out = np.full_like(t, np.nan)
        ok = z < sm["zs"]
        out[ok] = t[ok] * (sm["cL"][0] + sm["cL"][1] * z[ok] + sm["cL"][2] * z[ok] ** 2
                           + sm["cL"][3] * np.sqrt(1 - z[ok] / sm["zs"]))
        return out if len(out) > 1 else out[0]
    m["√ branch"] = sq
    return m, d, sm["tdeg"]


def main():
    OUT.mkdir(exist_ok=True)
    asc, desc = chain("re500_asc"), chain("re500_desc")
    m, tight, ts = models()
    prod = collect(500, "zero_thickness_plate/ell_re500"); prod = prod[(prod.aoa > 4.3) & (prod.verdict == "CONVERGED")]
    lines = [f"Fold test, Re 500: models fitted on tight 0.25-{tight.aoa.max():g} deg; square-root branch point {ts:.2f} deg",
             "Chain cases counted once ended with extrapolated remainder < 1e-5 (production tolerances).", ""]
    for nm, ch in (("ascending (this machine)", asc), ("descending (cfd-test)", desc)):
        lines.append(f"== {nm}: " + (", ".join(f"{a:g}" for a in sorted(ch)) if ch else "none finished yet"))
        prev = None
        for a in sorted(ch):
            cl, cd, rem = ch[a]
            slope = "" if prev is None else f"  dC_L/dalpha {(cl - prev[1]) / (a - prev[0]):.4f} /deg"
            lines.append(f"  {a:5g}  C_L {cl:.7f}  C_D {cd:.7f}  remainder {rem:.0e}{slope}" +
                         "".join(f"  {k}: {float(f(np.radians(a))) - cl:+.1e}" for k, f in m.items()))
            prev = (a, cl)
    both = sorted(set(asc) & set(desc))
    lines.append("")
    lines.append("== hysteresis (ascending - descending C_L at common angles): " +
                 (", ".join(f"{a:g}: {asc[a][0] - desc[a][0]:+.1e}" for a in both) if both else "no common angle finished yet"))
    txt = "\n".join(lines); print(txt); (OUT / "fold_check_re500.txt").write_text(txt + "\n")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    tt = np.radians(np.linspace(0.25, 12, 400))
    for k, f in m.items():
        yy = np.array([float(f(t)) for t in tt]); yy[(yy < 0) | (yy > 0.8)] = np.nan   # no spike through the pole
        ax[0].plot(np.degrees(tt), yy, LS[k], color=COL[k], lw=2, label=k)
    ax[0].plot(tight.aoa, tight.cl, "o", color="#333333", ms=4, label="tight data (fit)")
    ax[0].plot(prod.aoa, prod.cl, "s", mfc="white", mec="#333333", ms=5, label="production 4.5–6°")
    for ch, mk, lab in ((asc, "^", "chain ↑"), (desc, "v", "chain ↓")):
        if ch:
            ax[0].plot(list(ch), [v[0] for v in ch.values()], mk, color="#111111", ms=8, label=lab)
    ax[0].axvline(ts, color="#999999", lw=1); ax[0].annotate(f"√ branch {ts:.1f}°", (ts, 0.76), xytext=(4, 0),
                                                            textcoords="offset points", fontsize=8, color="#555555")
    ax[0].set_xlabel("α (deg)"); ax[0].set_ylabel("$C_L$"); ax[0].set_ylim(0, 0.8)
    ax[0].set_title("Models fitted on 0–4.25°, against the chains", fontsize=11)
    pts = sorted(set(asc) | set(desc) | set(prod.aoa))
    obs = {**{a: v for a, v in zip(prod.aoa, prod.cl)}, **{a: v[0] for a, v in desc.items()}, **{a: v[0] for a, v in asc.items()}}
    for k, f in m.items():
        e = [abs(float(f(np.radians(a))) - obs[a]) for a in pts]
        ax[1].semilogy(pts, e, LS[k], marker="o", color=COL[k], lw=2, ms=6, mec="white", mew=1.2, label=k)
        last = [(a, v) for a, v in zip(pts, e) if np.isfinite(v)][-1]
        dy = {"poly θ⁹": -7, "poly θ¹¹": 6}.get(k, 0)
        ax[1].annotate(k, last, xytext=(5, dy), textcoords="offset points", fontsize=8, color="#333333", va="center")
    ax[1].axhspan(1e-7, 1e-5, color="#e5e5e2", lw=0); ax[1].annotate("data precision", (pts[0], 3e-6), fontsize=8, color="#555555")
    ax[1].set_xlabel("α (deg)"); ax[1].set_ylabel("|model − data|")
    ax[1].set_title("Extrapolation error (lower is better)", fontsize=11)
    for a in ax:
        a.grid(color="#e5e5e2", lw=0.8); a.spines[["top", "right"]].set_visible(False)
    ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fold_extrapolation.png", dpi=150)


if __name__ == "__main__":
    main()
