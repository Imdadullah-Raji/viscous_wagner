"""Tests of instructions/lsb_hypothesis.md on the plate, Re = 500 (steady SIMPLE).

Claims:
  1. Lift is mostly linear when the LSB has not formed.
  2. When the LSB forms, cubic terms explain the curvature well.
  3. When trailing-edge separation happens, the cubic term loses dominance.

Markers (fine mesh, lsb/plate_fine_re500):
  alpha_f  : min upper-surface C_f (0.02 < x < 0.9) crosses zero  (formation_fine.py)
  alpha_TE : "TRAILING-EDGE SEPARATION" AS DEFINED FOR THIS PLATE (2026-10-09, by
             decision): the angle at which the leading-edge bubble's REATTACHMENT POINT
             REACHES THE TRAILING EDGE, i.e. the steady linear growth of x_reat(alpha)
             (8.0-8.3 deg, ~1.47 c/deg) extrapolated to x_flat, the end of the flat upper
             surface where the rounded TE begins.  Value used: ALPHA_TE = 8.40 deg,
             bracket (8.3, 8.4].  It is NOT a separation that starts at the trailing
             edge and moves forward (the plate's TE never separates on its own), and it is
             NOT the bubble opening into the wake: up to 9.0 deg the bubble stays closed,
             its rear stagnation point pinned in the TE corner (x_reat 0.967 -> 0.989).

Lift curve (base mesh, plate_case/re500, every CONVERGED case on the 0.25-deg grid 0-12):
  (A) Legendre fit on alpha >= alpha_TE, interval mapped to [-1, 1], all degrees;
      reported: linear slope, and slope/curvature from the lowest adequate degree.
  (B) the same on alpha_f <= alpha <= alpha_TE, base and fine points separately.
  (C) odd Legendre (P1..P9), mirrored, on 0-A for A = 6, 6.25, ..., 8:
      c1 c3 c5 c7 c9, ratios, P1+P3 collected into a1 theta + a3 theta^3.
  Claim-1 check: on 0-A for A below alpha_f, how much of C_L a straight line leaves
      unexplained, against P1+P3.

Writes figs/hypothesis/: markers.txt, fits.txt, *.csv, lift_drag_markers.png,
odd_window_scan.png, local_fits.png.   Run with the research env.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.polynomial import legendre as Lg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THESIS = ROOT.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "analysis"))
from streamlines_lsb import read_patch, upper_cf, crossings, latest  # noqa: E402
from legendre_odd_fit import ODD, fit, mirror, to_theta_monomial  # noqa: E402

FINE = HERE / "plate_fine_re500"
OUT = HERE / "figs" / "hypothesis"
ALPHA_TE = 8.40        # see the docstring: reattachment reaches the TE (decision, 2026-10-09)
ONSET = 11.33          # plate shedding onset, plate_case/unsteady_re500/ONSET.md (approx.)


# ----------------------------------------------------------------------------- data
def plate_data(re_):
    """Base-mesh sweep: last force sample of each case, verdict from RESULTS.md
    (run_sweep.py's range+drift gate).  Cases not run yet are skipped."""
    root = THESIS / "plate_case" / f"re{re_}"
    status = {}
    for line in (root / "RESULTS.md").read_text().splitlines():
        m = re.match(r"\|\s*([\d.]+)\s*\|\s*(\w+)\s*\|", line)
        if m:
            status[float(m.group(1))] = m.group(2)
    rows = []
    for case in root.glob("steady_aoa*"):
        if not case.is_dir() or not (case / "postProcessing" / "forceCoeffs").is_dir():
            continue
        a = float(case.name.removeprefix("steady_aoa"))
        cl, cd = forces_last(case)
        rows.append((a, cl, cd, status.get(a, "MISSING")))
    return pd.DataFrame(rows, columns=["aoa", "cl", "cd", "verdict"]).sort_values("aoa")


def fine_cases():
    out = []
    for p in FINE.glob("steady_aoa*"):
        m = re.fullmatch(r"steady_aoa([\d.]+)", p.name)
        if m and p.is_dir():
            out.append((float(m.group(1)), p))
    return sorted(out)


def forces_last(case):
    rows = []
    for f in sorted((case / "postProcessing" / "forceCoeffs").glob("*/forceCoeffs.dat")):
        rows.append(np.loadtxt(f, comments="#", usecols=(0, 2, 3), ndmin=2))
    a = np.vstack(rows)
    a = a[np.argsort(a[:, 0], kind="stable")]
    return a[-1, 2], a[-1, 1]


def surface(case):
    d = latest(case)
    Cw = read_patch(d / "C")[:, :2]
    tau = -read_patch(d / "wallShearStress")[:, :2]
    c0 = Cw.mean(axis=0)
    o = np.argsort(np.arctan2(Cw[:, 1] - c0[1], Cw[:, 0] - c0[0]))
    return Cw[o], tau[o]


def x_flat_end(Cw):
    """Chord-wise end of the flat upper surface: last face still within 1e-4 c of the
    flat-surface height (the rounded trailing edge starts there)."""
    up = Cw[Cw[:, 1] > 0]
    yt = np.median(up[(up[:, 0] > 0.3) & (up[:, 0] < 0.7), 1])
    flat = up[np.abs(up[:, 1] - yt) < 1e-4]
    return float(flat[:, 0].max())


def separation_table(cases):
    rows, xflat = [], None
    for a, case in cases:
        try:
            Cw, tau = surface(case)
        except (FileNotFoundError, ValueError, AttributeError):
            continue
        xflat = xflat or x_flat_end(Cw)
        x, cf = upper_cf(Cw, tau)
        m = (x > 0.02) & (x < 0.9)
        cr = crossings(x, cf)
        s = [v for v, k in cr if k == "separation"]
        r = [v for v, k in cr if k == "reattachment" and s and v > s[0]]
        cl, cd = forces_last(case)
        rows.append(dict(aoa=a, cf_min=float(cf[m].min()), x_sep=s[0] if s else np.nan,
                         x_reat=r[0] if r else np.nan, cl=cl, cd=cd))
    return pd.DataFrame(rows).sort_values("aoa").reset_index(drop=True), xflat


def markers(sep, xflat):
    """alpha_f from the C_f,min zero; alpha_TE where reattachment reaches xflat."""
    d = sep[(sep.aoa >= 7.85) & (sep.aoa <= 8.05)]
    p = np.polyfit(d.aoa, d.cf_min, 1)
    af = -p[1] / p[0]
    closed = sep[np.isfinite(sep.x_sep) & (sep.x_reat < xflat)]
    opened = sep[np.isfinite(sep.x_sep) & ~(sep.x_reat < xflat)]
    ate, bracket = np.nan, None
    if len(closed) and len(opened):
        a0 = closed.aoa.max()
        later = opened[opened.aoa > a0]
        if len(later):
            a1 = later.aoa.min()
            r0 = closed.loc[closed.aoa == a0, "x_reat"].item()
            # reattachment is pinned close to xflat, so interpolating in x is only a rough
            # estimate; the bracket is the hard result
            r_prev = closed[closed.aoa < a0].tail(1)
            if len(r_prev):
                ap, rp = r_prev.aoa.item(), r_prev.x_reat.item()
                slope = (r0 - rp) / (a0 - ap)
                ate = min(a0 + (xflat - r0) / slope, a1) if slope > 0 else np.nan
            bracket = (a0, a1)
    return af, ate, bracket


# ----------------------------------------------------------------------------- fits
def local_legendre(a, cl, lo, hi, deg):
    """Legendre fit on [lo, hi] mapped to [-1, 1].  Returns coefficients, slope and
    curvature at mid-interval (per degree), RMSE, max |res|."""
    x = (2 * a - (lo + hi)) / (hi - lo)
    c = Lg.legfit(x, cl, deg)
    r = cl - Lg.legval(x, c)
    s = 2 / (hi - lo)
    d1 = Lg.legval(0.0, Lg.legder(c)) * s
    d2 = Lg.legval(0.0, Lg.legder(c, 2)) * s ** 2 if deg >= 2 else 0.0
    lin = np.polyfit(a, cl, 1)[0]
    return dict(lo=lo, hi=hi, n=len(a), deg=deg, coeffs=c, slope_mid=d1, curv_mid=d2,
                lin_slope=lin, rmse=np.sqrt(np.mean(r ** 2)), max_res=np.abs(r).max())


def fmt_local(name, f):
    c = " ".join(f"c{k}={v:+.6f}" for k, v in enumerate(f["coeffs"]))
    return (f"{name}: [{f['lo']:.3f}, {f['hi']:.3f}] deg, n={f['n']}, degree {f['deg']}: {c}\n"
            f"    linear-LSQ slope {f['lin_slope']:.5f}/deg = {np.rad2deg(f['lin_slope']):.4f}/rad;"
            f"  mid-interval slope {f['slope_mid']:.5f}/deg = {np.rad2deg(f['slope_mid']):.4f}/rad;"
            f"  curvature {f['curv_mid']:+.5f}/deg^2;  RMSE {f['rmse']:.1e}, max|res| {f['max_res']:.1e}")


def odd_window(df, amax):
    """Odd Legendre on 0..amax (amax must be a computed angle), mirrored.  The number of
    odd terms is capped at the number of nonzero angles, so a short window is never
    under-determined; missing higher coefficients are NaN."""
    d = df[(df.aoa <= amax + 1e-9)]
    x0, y0 = d.aoa.to_numpy() / amax, d.cl.to_numpy()
    n = int((x0 > 0).sum())
    degs = ODD[:min(len(ODD), n)]
    x, y = mirror(x0, y0)
    c, rmse, mx = fit(x, y, degs)
    c = np.r_[c, np.full(len(ODD) - len(c), np.nan)]
    p1, r1, m1 = fit(x, y, (1,))
    p13, r13, m13 = fit(x, y, (1, 3)) if n >= 2 else (p1, r1, m1)
    th = to_theta_monomial(p13, amax, (1, 3))
    return dict(amax=amax, n=n, terms=len(degs), **{f"c{k}": v for k, v in zip(ODD, c)},
                c3_c1=c[1] / c[0], c5_c3=c[2] / c[1], c7_c3=c[3] / c[1], c9_c3=c[4] / c[1],
                rmse_fit=rmse, maxres_P1=m1, maxres_P13=m13, rmse_P1=r1, rmse_P13=r13,
                cl_end=y0[-1], a1=th[1], a3=th[3])


# ----------------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = plate_data(500)
    base = base[base.verdict == "CONVERGED"].reset_index(drop=True)
    sep, xflat = separation_table(fine_cases())
    af, _, _ = markers(sep, xflat)
    lin = sep[(sep.aoa >= 7.995) & (sep.aoa <= 8.305)]
    q = np.polyfit(lin.aoa, lin.x_reat, 1)
    ate_lin = (xflat - q[1]) / q[0]
    ate, bracket = ALPHA_TE, (8.3, 8.4)
    sep.to_csv(OUT / "fine_separation.csv", index=False)

    L = [f"x_flat (end of flat upper surface, start of rounded TE) = {xflat:.4f}",
         f"alpha_f  = {af:.4f} deg   (fine mesh, C_f,min = 0)"]
    L += [f"alpha_TE = {ALPHA_TE:.2f} deg  [DEFINITION: reattachment of the LE bubble reaches the "
          "trailing edge -- NOT a TE-initiated separation, NOT an open/burst bubble]",
          f"    x_reat linear on 8.0-8.3 deg: slope {q[0]:.3f} c/deg; reaches x_flat at {ate_lin:.3f} deg; "
          f"bracket ({bracket[0]:g}, {bracket[1]:g}]",
          "    x_reat above: " + ", ".join(f"{r.aoa:g}: {r.x_reat:.4f}" for r in sep[sep.aoa >= 8.35].itertuples())]
    L.append(f"plate shedding onset ~{ONSET} deg (unsteady runs); steady branch used above it")
    (OUT / "markers.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))

    ate_use = ate if np.isfinite(ate) else (bracket[1] if bracket else np.nan)

    # ---- fits
    F = []
    if np.isfinite(ate_use):
        post = base[base.aoa >= ate_use - 1e-9]
        if len(post) >= 3:
            for deg in (1, 2, 3):
                if len(post) > deg + 1:
                    F.append(fmt_local(f"(A) post-TE, base mesh", local_legendre(
                        post.aoa.to_numpy(), post.cl.to_numpy(), post.aoa.min(), post.aoa.max(), deg)))
            post_lo = post[post.aoa <= ONSET]
            if len(post_lo) >= 3:
                F.append(fmt_local(f"(A') post-TE below shedding onset, base mesh", local_legendre(
                    post_lo.aoa.to_numpy(), post_lo.cl.to_numpy(), post_lo.aoa.min(),
                    post_lo.aoa.max(), min(2, len(post_lo) - 2))))
        for label, dset in (("base", base), ("fine", sep)):
            mid = dset[(dset.aoa >= af - 1e-9) & (dset.aoa <= ate_use + 1e-9)]
            for deg in (1, 2, 3):
                if len(mid) > deg + 1:
                    F.append(fmt_local(f"(B) alpha_f..alpha_TE, {label} mesh", local_legendre(
                        mid.aoa.to_numpy(), mid.cl.to_numpy(), af, ate_use, deg)))
        # what a global odd cubic (fit 0..alpha_TE) implies for the slope/curvature mid-(B)
        g = odd_window(base, np.floor(ate_use * 4) / 4)
        am = np.deg2rad(0.5 * (af + ate_use))
        F.append(f"global P1+P3 on 0-{g['amax']:g}: a1={g['a1']:.4f}/rad, a3={g['a3']:.3f}/rad^3 -> at "
                 f"mid-(B) slope {np.deg2rad(g['a1'] + 3 * g['a3'] * am**2):.5f}/deg, "
                 f"curvature {np.deg2rad(1)**2 * 6 * g['a3'] * am:+.5f}/deg^2")
    (OUT / "fits.txt").write_text("\n".join(F) + "\n")
    print("\n".join(F))

    have = set(np.round(base.aoa, 6))
    W = pd.DataFrame([odd_window(base, A) for A in np.arange(6.0, 8.0 + 1e-9, 0.25)
                      if round(A, 6) in have])
    W.to_csv(OUT / "odd_windows_6to8.csv", index=False)
    C1 = pd.DataFrame([odd_window(base, A) for A in np.arange(1.0, 7.76, 0.25)
                       if round(A, 6) in have])
    C1.to_csv(OUT / "claim1_windows.csv", index=False)
    with pd.option_context("display.width", 200, "display.float_format", "{:.4g}".format):
        print(W[["amax", "n", "c1", "c3", "c5", "c7", "c9", "c3_c1", "c5_c3", "a1", "a3",
                 "maxres_P1", "maxres_P13"]].to_string(index=False))
        print(C1[["amax", "n", "cl_end", "maxres_P1", "maxres_P13", "a1", "a3"]]
              .to_string(index=False))

    # ---- figures
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for k, (col, lab) in enumerate((("cl", r"$C_L$"), ("cd", r"$C_D$"))):
        ax[k].plot(base.aoa, base[col], "o-", ms=3.5, color="#2a78d6", lw=1, label="base mesh (fits)")
        ax[k].plot(sep.aoa, sep[col], "s", ms=3, mfc="none", color="#c2410c", label="fine mesh (LSB cases)")
        ax[k].axvline(af, color="k", ls="--", lw=1, label=rf"$\alpha_f$ = {af:.3f}°")
        if np.isfinite(ate_use):
            ax[k].axvline(ate_use, color="m", ls="--", lw=1, label=rf"$\alpha_{{TE}}$ = {ate_use:.2f}° (reattachment reaches TE)")
        ax[k].axvspan(ONSET, base.aoa.max() + 0.2, color="0.9", zorder=0, label="sheds (steady branch)")
        ax[k].set(xlabel=r"$\alpha$ (deg)", ylabel=lab)
        ax[k].grid(alpha=0.25)
    ax[0].legend(fontsize=8, loc="lower right")
    ins = ax[0].inset_axes([0.08, 0.5, 0.36, 0.42])
    z = base[(base.aoa >= 7) & (base.aoa <= 10)]
    ins.plot(z.aoa, z.cl, "o-", ms=3, color="#2a78d6", lw=1)
    zz = sep[(sep.aoa >= 7) & (sep.aoa <= 10)]
    ins.plot(zz.aoa, zz.cl, "s", ms=2.5, mfc="none", color="#c2410c")
    ins.axvline(af, color="k", ls="--", lw=0.8)
    if np.isfinite(ate_use):
        ins.axvline(ate_use, color="m", ls="--", lw=0.8)
    ins.tick_params(labelsize=7)
    fig.suptitle("Plate, Re = 500, steady: LSB formation and trailing-edge separation")
    fig.tight_layout(); fig.savefig(OUT / "lift_drag_markers.png", dpi=170); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for k in ODD:
        ax[0].plot(W.amax, W[f"c{k}"], "o-", ms=3, label=f"c{k}")
    ax[0].set(xlabel="window 0–A (deg)", ylabel="Legendre coefficient", yscale="symlog")
    ax[0].set_yscale("symlog", linthresh=1e-5); ax[0].legend(fontsize=7)
    for k in ("c3_c1", "c5_c3", "c7_c3", "c9_c3"):
        ax[1].plot(W.amax, W[k].abs(), "o-", ms=3, label=k.replace("_", "/"))
    ax[1].set(xlabel="window 0–A (deg)", ylabel="|ratio|", yscale="log"); ax[1].legend(fontsize=7)
    ax[2].semilogy(C1.amax, C1.maxres_P1, "o-", ms=3, label="P1 only (straight line)")
    ax[2].semilogy(C1.amax, C1.maxres_P13, "s-", ms=3, label="P1+P3")
    ax[2].axvline(af, color="k", ls="--", lw=0.8)
    ax[2].set(xlabel="window 0–A (deg)", ylabel=r"max |residual| in $C_L$"); ax[2].legend(fontsize=7)
    for a_ in ax:
        a_.grid(alpha=0.25)
    fig.suptitle("Plate, Re = 500: odd-Legendre fits")
    fig.tight_layout(); fig.savefig(OUT / "odd_window_scan.png", dpi=170); plt.close(fig)

    # ---- claim 3: the odd expansion on 0-A beyond the markers (A up to 12)
    X = pd.DataFrame([odd_window(base, A) for A in np.arange(4.0, 12.0 + 1e-9, 0.25)
                      if round(A, 6) in have])
    X.to_csv(OUT / "odd_windows_4to12.csv", index=False)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for k, col in zip(("c5_c3", "c7_c3", "c9_c3"), ("#eb6834", "#1baf7a", "#7b3fbf")):
        ax[0].plot(X.amax, X[k], "o-", ms=3, color=col, label=k.replace("_", "/"))
    ax[0].axhline(0, color="0.6", lw=0.8)
    ax[0].set(xlabel="window 0–A (deg)", ylabel="ratio to c3", title="higher odd terms relative to the cubic")
    ax[1].plot(X.amax, X.maxres_P13, "s-", ms=3, color="k", label="max |res|, P1+P3")
    ax[1].set(xlabel="window 0–A (deg)", ylabel=r"$C_L$", title="two-term (θ, θ³) fit error", yscale="log")
    for a_ in ax:
        a_.axvline(af, color="k", ls="--", lw=0.9)
        a_.axvline(ALPHA_TE, color="m", ls="--", lw=0.9)
        a_.axvspan(ONSET, 12.2, color="0.9", zorder=0)
        a_.grid(alpha=0.25)
    ax[0].legend(fontsize=8)
    fig.suptitle(r"Plate, Re = 500: odd-Legendre fit on 0–A  (dashed: $\alpha_f$ black, $\alpha_{TE}$ magenta; grey: sheds)")
    fig.tight_layout(); fig.savefig(OUT / "cubic_dominance_0to12.png", dpi=170); plt.close(fig)

    # ---- does C_L register the markers?  slope and curvature on the dense fine-mesh grid
    g = sep[np.isclose(np.round(sep.aoa * 10, 6) % 1, 0) & (sep.aoa >= 7.5 - 1e-9)]
    A_, C_ = g.aoa.to_numpy(), g.cl.to_numpy()
    sl = np.diff(C_) / np.diff(A_); am = 0.5 * (A_[1:] + A_[:-1])
    cv = np.diff(sl) / np.diff(am); ac = 0.5 * (am[1:] + am[:-1])
    pd.DataFrame(dict(aoa=am, slope_per_deg=sl)).to_csv(OUT / "fine_slope.csv", index=False)
    pd.DataFrame(dict(aoa=ac, curvature_per_deg2=cv)).to_csv(OUT / "fine_curvature.csv", index=False)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(am, np.rad2deg(sl), "o-", ms=3, color="#c2410c")
    ax[0].set(xlabel=r"$\alpha$ (deg)", ylabel=r"$dC_L/d\alpha$ (per rad)", title="lift slope (fine mesh, 0.1° differences)")
    ax[1].plot(ac, cv, "o-", ms=3, color="#c2410c")
    ax[1].set(xlabel=r"$\alpha$ (deg)", ylabel=r"$d^2C_L/d\alpha^2$ (per deg²)",
              title="curvature (noise ≈ ±4e-4 from 1e-6 force noise)")
    for a_ in ax:
        a_.axvline(af, color="k", ls="--", lw=0.9, label=rf"$\alpha_f$ = {af:.3f}°")
        a_.axvline(ALPHA_TE, color="m", ls="--", lw=0.9, label=rf"$\alpha_{{TE}}$ = {ALPHA_TE:.2f}°")
        a_.grid(alpha=0.25)
    ax[0].legend(fontsize=8)
    fig.suptitle("Plate, Re = 500: no kink in C_L at either marker")
    fig.tight_layout(); fig.savefig(OUT / "slope_curvature_markers.png", dpi=170); plt.close(fig)


if __name__ == "__main__":
    main()
