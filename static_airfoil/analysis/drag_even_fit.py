"""The drag side of the small-angle expansion, zero-thickness plate, re-converged ("tight") data.

For a symmetric body C_D is even in alpha: C_D = d0 + d2 th^2 + d4 th^4 + ...  It is a functional of
the same steady solution as C_L, so it must have the same singularities: an independent check of the
real-axis singularity near 10.3-10.8 deg (Re 500) that analysis/pade_singularity.py found in C_L.

Per Re:
  1. convergence of C_D on the tight tail (same criterion as C_L), and tight - production for both;
  2. even Legendre fits P0..P8 / P0..P10 on windows 0-A (alpha = 0 included: the production 0-deg case,
     where the flow is symmetric), converted to powers of theta;
  3. body axes: C_N = C_L cos th + C_D sin th (odd), C_A = C_D cos th - C_L sin th (even).  Exactly to
     second order: n1 = a1 + d0 and d2 = n1 + A2 - A0/2.  The NET suction fraction eta = -A2/n1: 1 is the
     potential flat plate (C_A = -2 pi th^2 against C_N = 2 pi th, zero drag rise); in viscous flow A2
     also holds the change of skin friction with alpha, so eta is net of it;
  4. singularity: rational least squares C_D = P_L(z)/Q_M(z), z = th^2, and AAA, with the doublet tests;
  5. drag polar: (C_D - d0)/C_L^2 against C_L^2.

  python analysis/drag_even_fit.py 500 700 900 [--prod]   (research env)  -> analysis/drag_even_fit/
  --prod: the same analysis on the biased production values at the same angles
  --bodies: the four bodies at Re 500 on production data.  On the zero-thickness plate production C_D
            gives d0 and d2 to 3e-5 and d4 to 0.3-0.6 %, d6 not at all, so only d0, d2, d4, the body-axis
            quantities and the polar ratio are reported, and no singularity analysis.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import AAA

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lasso_re1000 import collect, load                              # noqa: E402
from cubic_stress_test import fit_window                            # noqa: E402
from pade_singularity import nearest, fmt, rational_fit             # noqa: E402

warnings.filterwarnings("ignore")
OUT = HERE / "drag_even_fit"
ROOT = HERE.parent / "zero_thickness_plate"
TIGHT = ROOT / "tight"
ZERO = {500: "ell_re500", 700: "ell_re_scaling/steady_re700", 900: "ell_re_scaling/steady_re900"}
COLOR = {500: "#2a78d6", 700: "#eb6834", 900: "#1baf7a"}            # categorical slots 1-3, validated
MARK = {500: "o", 700: "s", 900: "^"}
OK = 3e-9                                                           # usable: tail change and remainder


def tail(y, it, i, tN):
    """Last-500 change and the Aitken-extrapolated remainder of one force series."""
    v = np.array([y[i[k]] for k in range(tN - 2500, tN + 1, 500)])
    d = np.diff(v); q = d[-1] / d[-2] if d[-2] else np.nan
    rem = d[-1] * q / (1 - q) if 0 < q < 0.9 else np.nan
    return d[-1], rem


def usable(d500, rem):
    return abs(d500) <= OK and (abs(rem) <= OK if np.isfinite(rem) else abs(d500) <= 1e-9)


def cases(re):
    rows = []
    for c in sorted((TIGHT / f"re{re}").glob("steady_aoa*")):
        if not c.is_dir() or not list(c.glob("log.steady_to*")):
            continue
        lg = c / "log.steady"
        if not lg.exists() or b"\nEnd" not in lg.read_bytes():
            continue
        it, cd, cl = load(str(c)); it = it.astype(int); i = {k: j for j, k in enumerate(it)}
        t0 = sorted(int(d.name) for d in (c / "postProcessing/forceCoeffs").iterdir() if d.name != "0")[0]
        tN = it[-1]
        dl, rl = tail(cl, it, i, tN); dd, rd = tail(cd, it, i, tN)
        rows.append(dict(aoa=float(c.name[10:]), cl=cl[-1], cd=cd[-1], cl_prod=cl[i[t0]], cd_prod=cd[i[t0]],
                         dl500=dl, rl=rl, dd500=dd, rd=rd, ok_l=usable(dl, rl), ok_d=usable(dd, rd)))
    df = pd.DataFrame(rows).sort_values("aoa").reset_index(drop=True)
    # alpha = 0: the production case (no tight run exists), held to the same tail test.  The Re 900 one
    # stopped at 6000 iterations still moving 2e-7 per 500, so it fails and the fit sets d0 itself.
    it, cd, cl = load(str(ROOT / ZERO[re] / "steady_aoa0")); it = it.astype(int)
    dd, rd = tail(cd, it, {k: j for j, k in enumerate(it)}, it[-1])
    zero = dict(aoa=0.0, cl=0.0, cd=cd[-1], cl_prod=0.0, cd_prod=cd[-1], dl500=0, rl=0, dd500=dd, rd=rd,
                ok_l=True, ok_d=usable(dd, rd))
    return pd.concat([pd.DataFrame([zero]), df], ignore_index=True)


def axes(df):
    th = np.radians(df.aoa)
    return df.assign(cn=df.cl * np.cos(th) + df.cd * np.sin(th), ca=df.cd * np.cos(th) - df.cl * np.sin(th))


def rational_even(th, y, L, M, iters=30):
    """y = P_L(z)/Q_M(z), Q(0) = 1, z = (th/s)^2; Sanathanan-Koerner.  Poles/zeros in unscaled z."""
    s = th.max(); z = (th / s) ** 2; q = np.zeros(M)
    for _ in range(iters):
        Qp = 1 + sum(q[k] * z ** (k + 1) for k in range(M))
        A = np.column_stack([z ** k for k in range(L + 1)] + [-y * z ** (k + 1) for k in range(M)])
        W = (1 / np.abs(Qp))[:, None]
        sol, *_ = np.linalg.lstsq(A * W, y * W[:, 0], rcond=None)
        q = sol[L + 1:]
    p = sol[:L + 1]
    Q = 1 + sum(q[k] * z ** (k + 1) for k in range(M))
    rms = np.sqrt(np.mean((np.polyval(p[::-1], z) / Q - y) ** 2))
    return np.roots(np.r_[q[::-1], 1.0]) * s ** 2, (np.roots(p[::-1]) * s ** 2 if L else np.array([])), rms


def analyse(re, prod=False):
    df = cases(re)
    if prod:
        df = df.assign(cl=df.cl_prod, cd=df.cd_prod)
    out = [f"######## Re {re}, {'PRODUCTION' if prod else 'tight'} data (alpha = 0 from the production 0-deg case, if converged)"]
    out.append(" aoa   C_D(final)        C_D tight-prod  dC_D(last500)  remain   ok_D | C_L tight-prod  ok_L")
    for _, r in df.iterrows():
        out.append(f" {r.aoa:4g}  {r.cd:.12f}  {r.cd - r.cd_prod:+.2e}      {r.dd500:+.1e}     {r.rd:+.1e}  "
                   f"{'Y' if r.ok_d else '-'}    | {r.cl - r.cl_prod:+.2e}       {'Y' if r.ok_l else '-'}")
    D = df[df.ok_d].reset_index(drop=True)
    B = axes(df[df.ok_d & df.ok_l]).reset_index(drop=True)
    amax = D.aoa.max()
    res = {}
    for name, data, col, fam in (("C_D", D, "cd", "even"), ("C_A", B, "ca", "even"), ("C_N", B, "cn", "odd")):
        out.append(f"\n== {name} ({fam}), {len(data)} angles 0-{data.aoa.max():g}")
        orders = ((0, 2, 4, 6, 8), (0, 2, 4, 6, 8, 10)) if fam == "even" else ((1, 3, 5, 7, 9), (1, 3, 5, 7, 9, 11))
        for degs in orders:
            for A in sorted(set(list(np.arange(2.0, data.aoa.max() + 1e-9, 0.5)) + [data.aoa.max()])):
                n = ((data.aoa <= A + 1e-9) & (data.aoa > (-1 if fam == "even" else 0))).sum()
                if n < len(degs) + 2:
                    continue
                a, se, rms, n = fit_window(data, A, degs, col=col, include_zero=(fam == "even"))
                k = degs[1]; k2 = degs[2]
                ratio = np.degrees(np.sqrt(abs(a[k] / a[k2]))) if a[k2] else np.nan
                out.append(f"  P{degs[0]}-P{max(degs):<2} 0-{A:<4g} n={n:2d} rms={rms:.1e} " +
                           " ".join(f"c{j}={a[j]:+.7g}±{se[j]:.1g}" for j in degs) +
                           f" | sqrt|c{k}/c{k2}|={ratio:.1f}°")
                res[(name, max(degs), A)] = a
    # body-axis bookkeeping at the widest window, higher order
    top_e = max(k for k in res if k[0] == "C_D")[1]
    aD, aA = res[("C_D", top_e, amax)], res[("C_A", top_e, B.aoa.max())]
    aN = res[("C_N", top_e + 1, B.aoa.max())]
    eta = -aA[2] / aN[1]
    out.append(f"\n  body axes (P0-P{top_e} / P1-P{top_e + 1}, widest window): d0 = {aD[0]:.7f}, d2 = {aD[2]:.5f}, "
               f"n1 = {aN[1]:.5f}, A0 = {aA[0]:.7f}, A2 = {aA[2]:.5f}")
    aL, *_ = fit_window(B[B.aoa > 0], B.aoa.max(), tuple(range(1, top_e + 2, 2)), col="cl")
    out.append(f"  checks: n1 - (a1 + d0) = {aN[1] - (aL[1] + aD[0]):+.1e};  d2 - (n1 + A2 - A0/2) = "
               f"{aD[2] - (aN[1] + aA[2] - aA[0] / 2):+.1e}")
    out.append(f"  net suction fraction eta = -A2/n1 = {eta:.4f};  drag-polar leading ratio d2/a1^2 = {aD[2] / aL[1] ** 2:.5f}")
    out.append(f"  flat-plate Blasius scale 2 x 1.328/sqrt(Re) = {2 * 1.328 / np.sqrt(re):.4f}  (d0 = {aD[0]:.4f})")
    # singularity, C_D and C_A
    for name, data, col in (("C_D", D, "cd"), ("C_A", B, "ca"), ("C_N", B, "cn"), ("C_L", B, "cl")):
        odd = name in ("C_N", "C_L")
        data = data[data.aoa > 0] if odd else data
        th, y = np.radians(data.aoa.to_numpy()), data[col].to_numpy()
        out.append(f"\n== singularity from {name}: rational least squares {'th P/Q' if odd else 'P/Q'} in z, "
                   "nearest GENUINE pole (doublets skipped)")
        for L, M in ((1, 1), (2, 1), (1, 2), (2, 2), (3, 2), (2, 3), (3, 3)):
            if len(th) < L + M + 3:
                continue
            p, zz, rms = rational_fit(th, y, L, M) if odd else rational_even(th, y, L, M)
            g = nearest(p, zz, genuine=True)
            zero = ""
            if g is not None and len(zz):
                zr = zz[np.argmin(np.abs(zz - g[0]))]
                zero = f"   nearest zero {np.degrees(abs(np.sqrt(complex(zr)))):.2f}°"
            out.append(f"  [{L}/{M}] rms={rms:.1e}  {fmt(g[0]) if g else 'no genuine pole'}{zero}")
        out.append(f"   AAA on {name}{'/theta' if odd else ''}(z)  (genuine: residue > 1e-4 x largest near AND no zero within 1 % in z):")
        for rtol in (1e-5, 1e-6, 1e-7):
            r = AAA(th ** 2, y / th if odd else y, rtol=rtol)
            pol, rs, zer = r.poles(), np.abs(r.residues()), r.roots()
            near = [(v, q) for v, q in zip(pol, rs) if abs(np.sqrt(complex(v))) < np.radians(40)]
            big = max((q for _, q in near), default=0)
            far = lambda v: (np.min(np.abs(zer - v)) / abs(v) if len(zer) else np.inf) >= 0.01
            gen = sorted({float(round(np.degrees(abs(np.sqrt(complex(v)))), 2)) for v, q in near
                          if q > 1e-4 * big and far(v)})
            zs = sorted({float(round(np.degrees(abs(np.sqrt(complex(v)))), 2)) for v in zer
                         if abs(np.sqrt(complex(v))) < np.radians(40)})
            out.append(f"    rtol {rtol:.0e}: m={len(r.support_points):2d}  genuine poles {gen}   zeros {zs[:4]}")
    return df, D, aD, out


def figure(results):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for re, (df, D, aD, _) in results.items():
        P = D[D.aoa > 0]
        c2 = P.cl.to_numpy() ** 2
        ax[0].plot(c2, (P.cd - aD[0]) / c2, MARK[re] + "-", color=COLOR[re], lw=2, ms=8, mec="white", mew=1.5,
                   label=f"Re {re}")
        ax[0].annotate(f"Re {re}", (c2[-1], ((P.cd - aD[0]) / c2).iloc[-1]), xytext=(6, 0),
                       textcoords="offset points", va="center", color="#333333", fontsize=9)
        th = np.radians(P.aoa.to_numpy())
        ax[1].plot(np.degrees(th), (P.cd - aD[0]) / th ** 2, MARK[re] + "-", color=COLOR[re], lw=2, ms=8,
                   mec="white", mew=1.5, label=f"Re {re}")
    ax[0].set_xlabel("$C_L^2$"); ax[0].set_ylabel("$(C_D - d_0)/C_L^2$")
    ax[0].set_title("Drag polar: departure from a parabola", fontsize=11)
    ax[1].set_xlabel("α (deg)"); ax[1].set_ylabel("$(C_D - d_0)/\\theta^2$")
    ax[1].set_title("Drag rise per $\\theta^2$ (tends to $d_2$ as α → 0)", fontsize=11)
    for a in ax:
        a.grid(color="#e5e5e2", lw=0.8); a.spines[["top", "right"]].set_visible(False)
        a.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "drag.png", dpi=150)


def bodies():
    from cubic_stress_test import BODIES
    out = ["Four bodies, Re 500, PRODUCTION data: d0, d2, d4 (P0-P8 even fit), a1 and n1 (P1-P9 odd fit),",
           "net suction fraction eta = -A2/n1, polar ratio d2/a1^2.  Windows 0-4 and 0-6 deg (spread = resolution).", "",
           f"{'body':22s} {'A':>3s}  {'d0':>10s}  {'d2':>9s}  {'d4':>9s}  {'a1':>8s}  {'n1':>8s}  {'A2':>8s}  {'eta':>6s}  {'d2/a1^2':>8s}  rms(C_D)"]
    for name, (base, _) in BODIES.items():
        df = collect(500, base)
        df = axes(df[df.verdict == "CONVERGED"].sort_values("aoa").reset_index(drop=True))
        for A in (4, 6):
            if df.aoa.max() < A - 1e-9:
                continue
            aD, _, rms, _ = fit_window(df, A, (0, 2, 4, 6, 8), col="cd", include_zero=True)
            aA, *_ = fit_window(df, A, (0, 2, 4, 6, 8), col="ca", include_zero=True)
            aN, *_ = fit_window(df, A, (1, 3, 5, 7, 9), col="cn")
            aL, *_ = fit_window(df, A, (1, 3, 5, 7, 9), col="cl")
            out.append(f"{name:22s} {A:3g}  {aD[0]:10.7f}  {aD[2]:9.5f}  {aD[4]:9.3f}  {aL[1]:8.5f}  {aN[1]:8.5f}  "
                       f"{aA[2]:8.4f}  {-aA[2] / aN[1]:6.4f}  {aD[2] / aL[1] ** 2:8.5f}  {rms:.1e}")
    txt = "\n".join(out)
    print(txt)
    (OUT / "drag_four_bodies_re500.txt").write_text(txt + "\n")


def main():
    if "--bodies" in sys.argv:
        OUT.mkdir(exist_ok=True); bodies(); return
    prod = "--prod" in sys.argv
    res_list = [int(a) for a in sys.argv[1:] if not a.startswith("--")] or [500, 700, 900]
    OUT.mkdir(exist_ok=True)
    results = {}
    for re in res_list:
        df, D, aD, out = analyse(re, prod)
        results[re] = (df, D, aD, out)
        txt = "\n".join(out)
        print(txt + "\n")
        (OUT / f"drag_re{re}{'_prod' if prod else ''}.txt").write_text(txt + "\n")
    if not prod:
        figure(results)


if __name__ == "__main__":
    main()
