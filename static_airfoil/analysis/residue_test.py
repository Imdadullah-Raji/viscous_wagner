"""Residue-cancellation test: does the singular part of the force point along the lift?

Section 9 of CUBIC_STRESS_TEST.md: C_L (Re 500) and C_A carry a real-axis singularity, C_D shows none.
Hypothesis: near theta_s the singular parts cancel in C_D, i.e. the singular force is along the lift.

Separate Padé fits put the poles in different places, so their residues cannot be compared.  Here both
forces share ONE denominator:  C_L = theta P_L(z)/Q(z),  C_D = P_D(z)/Q(z),  z = theta^2,  Q(0) = 1.
Partial fractions then give each force's residue at the same pole z_k:
    r_L = theta_k P_L(z_k)/Q'(z_k)   (C_L as a function of z),     r_D = P_D(z_k)/Q'(z_k),
and the direction of the singular force from the lift axis, phi = atan2(r_D, r_L).
Hypothesis: |r_D/r_L| << 1.  A singular force along the body normal would give tan(theta_s) (~0.19 at
10.8 deg); a random direction gives order 1.

  1. joint fit, Sanathanan-Koerner, each block weighted by 1/sigma of its own noise;
  2. doublet test on the shared poles: a zero of BOTH numerators within 1 % -> spurious;
  3. orders [L/M], windows, 200 noise refits;
  4. model comparison: with Q fixed, refit P_D with and without the factor (z - z_s);
  5. synthetic calibration first: smooth parts plus a shared sqrt branch with known c_D/c_L.
     First version (B_L tuned to a3 only) FAILED: its singularity was ~4x weaker than the data's (a5 -99
     vs -367), the pole came out at 14 deg, and nonzero ratios were recovered ~30 % low.  Now the truth is
     a fit of the REAL data:  C_L/theta = poly2(z) + B_L h(z),  h = sqrt(1 - z/z_s), z_s scanned;
     C_D = poly3(z) + B_D h(z) with the same z_s.  That fit is also a third, model-based route to the ratio:
     B_D/(theta_s B_L).
  Verdict rule.  As fixed in advance it pooled every order; the low orders do not fit the data (C_L rms
  1e2-1e4 x the noise) and their ratios drove it.  AMENDED after seeing the data, and reported as such:
  only fits reaching the noise floor (rms_L <= 1.5 sigma_L, rms_D <= 2 sigma_D) count.  Both are printed.

  python analysis/residue_test.py [RE]    (research env)  -> analysis/residue_test/
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
from cubic_stress_test import fit_window                            # noqa: E402
from drag_even_fit import cases                                      # noqa: E402

warnings.filterwarnings("ignore")
OUT = HERE / "residue_test"
RNG = np.random.default_rng(1)
ORDERS = [(L, M) for M in (1, 2, 3) for L in (1, 2, 3)]


def joint_fit(th, cl, cd, L, M, sL, sD, iters=40):
    """Shared-denominator fit.  Returns (pL, pD, q) in scaled z (z/s2) and the scale s2."""
    s2 = th.max() ** 2; z = th ** 2 / s2; yL = cl / th; n = len(th)
    q = np.zeros(M)
    for _ in range(iters):
        Qp = np.abs(1 + sum(q[k] * z ** (k + 1) for k in range(M)))
        Z = np.column_stack([z ** k for k in range(L + 1)])
        zero = np.zeros_like(Z)
        AL = np.hstack([Z, zero, np.column_stack([-yL * z ** (k + 1) for k in range(M)])])
        AD = np.hstack([zero, Z, np.column_stack([-cd * z ** (k + 1) for k in range(M)])])
        wL = (th / sL / Qp)[:, None]; wD = (1 / sD / Qp)[:, None]       # C_L residual in C_L units
        A = np.vstack([AL * wL, AD * wD]); b = np.r_[yL * wL[:, 0], cd * wD[:, 0]]
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        q_new = sol[2 * (L + 1):]
        if np.allclose(q_new, q, rtol=1e-13, atol=1e-16):
            q = q_new; break
        q = q_new
    return sol[:L + 1], sol[L + 1:2 * (L + 1)], q, s2


def poles_and_residues(pL, pD, q, s2):
    """Shared poles (unscaled z), residues of C_L(z) and C_D(z), and the doublet flags."""
    if len(q) == 0:
        return []
    Q = np.poly1d(np.r_[q[::-1], 1.0]); PL = np.poly1d(pL[::-1]); PD = np.poly1d(pD[::-1])
    zl, zd = PL.r, PD.r
    out = []
    for zk in Q.r:                                                   # scaled z
        dQ = Q.deriv()(zk)
        tk = np.sqrt(complex(zk * s2))
        rL = tk * PL(zk) / dQ * s2                                   # residue in unscaled z
        rD = PD(zk) / dQ * s2
        near = lambda r: (np.min(np.abs(r - zk)) / abs(zk)) if len(r) else np.inf
        out.append(dict(z=zk * s2, rL=rL, rD=rD, dbl_L=near(zl) < 0.01, dbl_D=near(zd) < 0.01))
    return out


def lead(res):
    """Nearest pole that is not a doublet of BOTH numerators (a zero in P_D alone is the hypothesis)."""
    for r in sorted(res, key=lambda r: abs(r["z"])):
        if not (r["dbl_L"] and r["dbl_D"]):
            return r
    return None


def theta_deg(z):
    t = np.sqrt(complex(z)); t = -t if t.real < 0 else t
    return np.degrees(t.real), np.degrees(t.imag)


def ratio(r):
    return (r["rD"] / r["rL"]).real if abs(np.imag(r["z"])) < 1e-12 else np.nan


def run(th, cl, cd, L, M, sL, sD, nboot=200):
    pL, pD, q, s2 = joint_fit(th, cl, cd, L, M, sL, sD)
    r = lead(poles_and_residues(pL, pD, q, s2))
    if r is None:
        return None
    zr = r["z"]; rats = []
    for _ in range(nboot):
        pL2, pD2, q2, _ = joint_fit(th, cl + RNG.normal(0, sL, len(th)), cd + RNG.normal(0, sD, len(th)), L, M, sL, sD)
        rr = lead(poles_and_residues(pL2, pD2, q2, s2))
        if rr is not None and abs(rr["z"] - zr) < 0.3 * abs(zr):     # same pole (within 30 % in z)
            rats.append(ratio(rr))
    rats = np.array([x for x in rats if np.isfinite(x)])
    Q = 1 + sum(q[k] * (th ** 2 / s2) ** (k + 1) for k in range(M))
    rmsL = np.sqrt(np.mean((th * np.polyval(pL[::-1], th ** 2 / s2) / Q - cl) ** 2))
    rmsD = np.sqrt(np.mean((np.polyval(pD[::-1], th ** 2 / s2) / Q - cd) ** 2))
    return dict(r=r, rat=ratio(r), lo=np.percentile(rats, 16) if len(rats) else np.nan,
                hi=np.percentile(rats, 84) if len(rats) else np.nan, found=len(rats) / nboot,
                rmsL=rmsL, rmsD=rmsD, pD=pD, q=q, s2=s2)


def constrained(th, cd, q, s2, L, zs):
    """Q fixed: C_D rms with P_D free (degree L) and with P_D = (z - z_s) R(z), deg R = L - 1."""
    z = th ** 2 / s2; zsc = zs / s2
    Q = 1 + sum(q[k] * z ** (k + 1) for k in range(len(q)))
    Z = np.column_stack([z ** k for k in range(L + 1)]) / Q[:, None]
    free = np.linalg.lstsq(Z, cd, rcond=None)[0]
    Zc = np.column_stack([(z - zsc) * z ** k for k in range(L)]) / Q[:, None]
    cons = np.linalg.lstsq(Zc, cd, rcond=None)[0]
    return np.sqrt(np.mean((Z @ free - cd) ** 2)), np.sqrt(np.mean((Zc @ cons - cd) ** 2))


def sigmas(df, A):
    _, _, sL, _ = fit_window(df[df.aoa > 0], A, (1, 3, 5, 7, 9, 11), col="cl")
    _, _, sD, _ = fit_window(df, A, (0, 2, 4, 6, 8, 10), col="cd", include_zero=True)
    return sL, sD


def sqrt_model(th, cl, cd, sL, sD):
    """Fit C_L/theta = poly2(z) + B_L h, h = sqrt(1 - z/z_s), scanning z_s; then C_D = poly3(z) + B_D h."""
    z = th ** 2; best = None
    for tdeg in np.arange(5.0, 25.0, 0.02):
        zs = np.radians(tdeg) ** 2
        if zs <= z.max():
            continue
        h = np.sqrt(1 - z / zs)
        X = np.column_stack([np.ones_like(z), z, z ** 2, h]) * th[:, None]
        c, *_ = np.linalg.lstsq(X, cl, rcond=None)
        rms = np.sqrt(np.mean((X @ c - cl) ** 2))
        if best is None or rms < best[0]:
            best = (rms, tdeg, c)
    rmsL, tdeg, cL = best
    zs = np.radians(tdeg) ** 2; h = np.sqrt(1 - z / zs)
    XD = np.column_stack([np.ones_like(z), z, z ** 2, z ** 3, h])
    cD, *_ = np.linalg.lstsq(XD, cd, rcond=None)
    rmsD = np.sqrt(np.mean((XD @ cD - cd) ** 2))
    return dict(tdeg=tdeg, zs=zs, cL=cL, cD=cD, rmsL=rmsL, rmsD=rmsD,
                ratio=cD[-1] / (np.radians(tdeg) * cL[-1]))


def adequate(res, sL, sD):
    return res["rmsL"] <= 1.5 * sL and res["rmsD"] <= 2 * sD


def synthetic(th, cl_real, cd_real, sL, sD, lines):
    """Truth = the sqrt-branch fit of the real data, with B_D replaced to give known ratios."""
    sm = sqrt_model(th, cl_real, cd_real, sL, sD)
    z = th ** 2; h = np.sqrt(1 - z / sm["zs"]); ts = np.radians(sm["tdeg"]); BL = sm["cL"][-1]
    lines.append(f"Square-root model of the real data (route 3): C_L/theta = poly2(z) + B_L h, h = sqrt(1 - z/z_s), z_s scanned:")
    lines.append(f"  theta_s = {sm['tdeg']:.2f} deg, B_L = {BL:+.4f}, rms_L = {sm['rmsL']:.1e};  C_D = poly3(z) + B_D h: "
                 f"B_D = {sm['cD'][-1]:+.3e}, rms_D = {sm['rmsD']:.1e};  model ratio B_D/(theta_s B_L) = {sm['ratio']:+.4f}")
    lines.append("Synthetic calibration: that fit as the truth, B_D set for known ratios, noise sigma_L, sigma_D;")
    lines.append("  recovered by every order that reaches the noise floor (orders that do not are not scored)")
    ok_orig, ok_disc, rec = True, True, {}
    for rho in (0.0, 0.01, 0.05, 0.19, 1.0):
        cl = th * (np.column_stack([np.ones_like(z), z, z ** 2]) @ sm["cL"][:3] + BL * h) + RNG.normal(0, sL, len(th))
        BD = rho * ts * BL
        cd = np.column_stack([np.ones_like(z), z, z ** 2, z ** 3]) @ sm["cD"][:4] + BD * h + RNG.normal(0, sD, len(th))
        got, vals = [], []
        for L, M in ORDERS:
            res = run(th, cl, cd, L, M, sL, sD, nboot=60)
            if res is None or not adequate(res, sL, sD) or not np.isfinite(res["rat"]):
                continue
            a, _ = theta_deg(res["r"]["z"])
            got.append(f"[{L}/{M}] {a:5.2f}deg {res['rat']:+.4f}"); vals.append(res["rat"])
            ok_orig &= abs(res["rat"] - rho) < max(0.05, 0.2 * rho)
        rec[rho] = vals
        lines.append(f"  true {rho:4.2f}: " + ("; ".join(got) if got else "no adequate fit"))
    ok_disc = (all(abs(v) < 0.01 for v in rec[0.0]) and bool(rec[0.0]) and
               all(v > 0.1 for v in rec[0.19]) and bool(rec[0.19]))
    lines.append(f"  original criterion (every adequate order within max(0.05, 20 %)): {'PASSED' if ok_orig else 'FAILED'}")
    lines.append(f"  discrimination (0 -> |r| < 0.01, 0.19 -> r > 0.1; amended, after seeing data): {'PASSED' if ok_disc else 'FAILED'}")
    return ok_orig, ok_disc


def main():
    re = int(next((a for a in sys.argv[1:] if not a.startswith("--")), 500))
    OUT.mkdir(exist_ok=True)
    df = cases(re)
    df = df[df.ok_l & df.ok_d & (df.aoa > 0)].reset_index(drop=True)
    th_all, cl_all, cd_all = np.radians(df.aoa.to_numpy()), df.cl.to_numpy(), df.cd.to_numpy()
    sL, sD = sigmas(cases(re).query("ok_l and ok_d"), df.aoa.max())
    lines = [f"Residue test, zero-thickness plate, Re {re}, tight data: {len(df)} angles "
             f"{df.aoa.min():g}-{df.aoa.max():g} deg; block noise sigma_L = {sL:.1e}, sigma_D = {sD:.1e}", ""]
    calib_orig, calib_disc = synthetic(th_all, cl_all, cd_all, sL, sD, lines)
    lines.append("")
    lines.append("Real data: shared-denominator fits; nearest pole that is not a doublet of both numerators")
    lines.append("  r_D/r_L with its 16-84 % range over 200 noise refits; 'P_D zero' = P_D alone vanishes at the pole")
    rows = []
    windows = sorted({w for w in (3.5, 4.0, df.aoa.max()) if w <= df.aoa.max() + 1e-9})
    for A in windows:
        m = df.aoa.to_numpy() <= A + 1e-9
        th, cl, cd = th_all[m], cl_all[m], cd_all[m]
        for L, M in ORDERS:
            if 2 * (L + 1) + M >= 2 * len(th) - 2:
                continue
            res = run(th, cl, cd, L, M, sL, sD)
            if res is None:
                lines.append(f"  0-{A:<4g} [{L}/{M}]  no genuine shared pole"); continue
            a, b = theta_deg(res["r"]["z"])
            fD, fC = constrained(th, cd, res["q"], res["s2"], L, res["r"]["z"]) if abs(b) < 1e-9 else (np.nan, np.nan)
            phi = np.degrees(np.arctan(res["rat"])) if np.isfinite(res["rat"]) else np.nan
            tag = "  P_D zero" if res["r"]["dbl_D"] else ""
            lines.append(f"  0-{A:<4g} [{L}/{M}] rmsL {res['rmsL']:.1e} rmsD {res['rmsD']:.1e}  pole {a:6.2f}{b:+.2f}i deg  "
                         f"r_D/r_L {res['rat']:+.4f} ({res['lo']:+.4f}..{res['hi']:+.4f}, found {res['found']:.0%})  "
                         f"phi {phi:+6.2f} deg  C_D rms free/no-pole {fD:.2e}/{fC:.2e}{tag}")
            rows.append(dict(A=A, L=L, M=M, th=a, im=b, rat=res["rat"], lo=res["lo"], hi=res["hi"], fD=fD, fC=fC,
                             ok=adequate(res, sL, sD)))
    def verdict_of(sel):
        rat = np.array([r["rat"] for r in sel]); lo = np.nanmin([r["lo"] for r in sel]); hi = np.nanmax([r["hi"] for r in sel])
        ts = np.tan(np.radians(np.median([r["th"] for r in sel])))
        span = max(abs(lo), abs(hi))
        v = ("CONFIRMED" if span < 0.05 else
             "REJECTED" if (lo <= ts <= hi) or (lo <= -ts <= hi) or (lo <= 1 <= hi) or (lo <= -1 <= hi) else "INCONCLUSIVE")
        return v, np.median(rat), lo, hi, ts
    allreal = [r for r in rows if abs(r["im"]) < 1e-9 and np.isfinite(r["rat"])]
    real = [r for r in allreal if r["ok"]]
    if allreal:
        v0, m0, lo0, hi0, _ = verdict_of(allreal)
        lines += ["", f"Original rule (all {len(allreal)} real-pole fits): median {m0:+.4f}, intervals {lo0:+.4f}..{hi0:+.4f} -> {v0}"]
    if real:
        rat = np.array([r["rat"] for r in real]); lo = np.nanmin([r["lo"] for r in real]); hi = np.nanmax([r["hi"] for r in real])
        ts = np.tan(np.radians(np.median([r["th"] for r in real])))
        span = max(abs(lo), abs(hi))
        verdict = ("CONFIRMED" if span < 0.05 else
                   "REJECTED" if (lo <= ts <= hi) or (lo <= -ts <= hi) or (lo <= 1 <= hi) or (lo <= -1 <= hi) else "INCONCLUSIVE")
        cost = np.nanmedian([r["fC"] / r["fD"] for r in real])
        lines += [f"Amended rule (the {len(real)} fits at the noise floor): median {np.median(rat):+.4f}, intervals "
                  f"{lo:+.4f}..{hi:+.4f};  tan(theta_s) = {ts:.3f} -> {verdict}",
                  f"  C_D rms cost of forbidding the pole outright: median ratio no-pole/free = {abs(cost):.1f}",
                  f"  calibration: original criterion {'passed' if calib_orig else 'FAILED'}, discrimination "
                  f"{'passed' if calib_disc else 'FAILED'}"]
        fig, ax = plt.subplots(figsize=(9, 3.8))
        x = np.arange(len(real))
        ax.errorbar(x, rat, yerr=[rat - np.array([r["lo"] for r in real]), np.array([r["hi"] for r in real]) - rat],
                    fmt="o", color="#2a78d6", ms=7, mec="white", mew=1.5, capsize=3, label="r_D / r_L (16–84 % noise range)")
        for yv, lab in ((0, "hypothesis: 0"), (ts, f"along the body normal: tan θ_s = {ts:.2f}"), (-ts, None)):
            ax.axhline(yv, color="#888888", lw=1, ls="-" if yv == 0 else "--")
            if lab:
                ax.annotate(lab, (len(real) - 0.5, yv), xytext=(0, 4), textcoords="offset points", ha="right",
                            fontsize=8, color="#555555")
        ax.set_xticks(x, [f"[{r['L']}/{r['M']}]\n0–{r['A']:g}°" for r in real], fontsize=7)
        ax.set_ylabel("$r_D / r_L$"); ax.set_title(f"Direction of the singular force, Re {re}", fontsize=11)
        ax.grid(color="#e5e5e2", lw=0.8); ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(OUT / f"residue_re{re}.png", dpi=150)
    txt = "\n".join(lines); print(txt)
    (OUT / f"residue_re{re}.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
