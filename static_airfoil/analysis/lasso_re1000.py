"""Steady-lift Lasso analysis at Re=1000, side by side with Re=500.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/lasso_re1000.py

Repeats the analysis of lasso_reg.py / lasso_path.py (degree-9 polynomial in
theta [rad], standardized, Lasso at lambda = 1e-4, LARS path) on the Re=1000
steady (SIMPLE) sweep in NACA0012/steady_re1000/, and on the Re=500 sweep in
NACA0012/steady_re500/ restricted to the SAME angle grid, so the two are comparable
point for point.  Then three checks that the Lasso selection alone cannot give:

  * leave-one-out cross-validated lambda (is 1e-4 a sensible penalty?),
  * bootstrap selection frequency at lambda = 1e-4 (is {theta, theta^3}
    stable, or an arbitrary pick among nearly collinear powers?),
  * OLS odd cubic / odd quintic (does theta^5 buy anything once theta^3 is in?).

Only cases whose FORCE trace is CONVERGED (steady_summary.verdict) enter a fit.
Outputs go to analysis/figs/re1000_steady/.
"""
from pathlib import Path
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import Lasso, LassoCV, LinearRegression, lars_path
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from steady_summary import load, verdict  # noqa: E402

OUT = ROOT / "analysis" / "figs" / "re1000_steady"
LAMBDA = 1e-4
ACTIVE_AT = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 1e-6)
DEGREE = 9
POWERS = list(range(1, DEGREE + 1))              # --odd -> [1, 3, 5, 7, 9]
LABEL = "degree-9"
AMAX = 12.0                                     # upper fit angle (deg); --amax
FIT_GRID = np.arange(0, AMAX + 1e-4, 0.5)       # degrees, both Re
# colour follows the Reynolds number, never its rank (slots 1-4 of the palette)
COLOR = {500: "#2a78d6", 1000: "#eb6834", 100: "#1baf7a", 2000: "#eda100"}
RES = (500, 1000)                                # --re
TAG = "re500_vs_re1000"                          # figure-name suffix
MUTED, INK = "#9a9893", "#52514e"


def collect(re, base=None):
    """Steady NACA0012 cases at one Re.  `base` overrides the folder, relative to static_airfoil/
    (e.g. "NACA0006/steady_re500")."""
    rows = []
    base = ROOT / base if base else ROOT / "NACA0012" / f"steady_re{re}"
    if not base.is_dir():                    # Re-scaling sweep (200..900) lives one level down
        base = ROOT / "NACA0012" / "steady_re_scaling" / f"steady_re{re}"
    for case in base.glob("steady_aoa*"):
        try:
            a = float(case.name.removeprefix("steady_aoa"))
        except ValueError:
            continue
        d = load(str(case))
        if d is None:
            continue
        it, cd, cl = d
        v, _ = verdict(it, cl, 1e-4)
        rows.append((a, cl[-1], cd[-1], v, int(it[-1])))
    return pd.DataFrame(rows, columns=["aoa", "cl", "cd", "verdict", "iters"]
                        ).sort_values("aoa").reset_index(drop=True)


def fit_set(df):
    m = np.isin(np.round(df.aoa, 6), np.round(FIT_GRID, 6)) & (df.verdict == "CONVERGED")
    return df[m]


def design(theta):
    # identical to PolynomialFeatures(DEGREE, include_bias=False) when POWERS
    # is 1..DEGREE; --odd keeps only the powers allowed by C_L(-a) = -C_L(a)
    scaler = StandardScaler()
    X = scaler.fit_transform(np.column_stack([theta ** p for p in POWERS]))
    return X, scaler


def lasso_at(X, y, lam):
    return Lasso(alpha=lam, tol=1e-12, max_iter=1_000_000).fit(X, y)


def analyse(re, df):
    d = fit_set(df)
    theta, y = np.deg2rad(d.aoa.to_numpy()), d.cl.to_numpy()
    X, scaler = design(theta)
    lams, _, path = lars_path(X, y - y.mean(), method="lasso")
    model = lasso_at(X, y, LAMBDA)
    std = model.coef_
    raw = std / scaler.scale_
    intercept = model.intercept_ - scaler.mean_ @ raw
    rmse = np.sqrt(np.mean((model.predict(X) - y) ** 2))
    # The path is not monotone (theta^4 is active at 1e-3, then drops out as
    # theta^3 enters), so report the active set at fixed penalties rather than
    # a single "entry" lambda.
    active = {lam: [POWERS[j] for j in range(len(POWERS))
                    if abs(np.interp(lam, lams[::-1], path[j][::-1])) > 1e-10]
              for lam in ACTIVE_AT}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        cv = LassoCV(alphas=np.logspace(-8, -1, 120), cv=LeaveOneOut(),
                     max_iter=200_000, tol=1e-10).fit(X, y)

    rng = np.random.default_rng(0)
    nboot, sel = 500, np.zeros(len(POWERS))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        for _ in range(nboot):
            i = rng.integers(0, len(y), len(y))
            Xb, _ = design(theta[i])
            sel += np.abs(lasso_at(Xb, y[i], LAMBDA).coef_) > 1e-10
    sel /= nboot

    ols = {}
    for name, powers in (("odd cubic", (1, 3)), ("odd quintic", (1, 3, 5))):
        A = np.column_stack([theta ** p for p in powers])
        lr = LinearRegression().fit(A, y)
        res = y - lr.predict(A)
        ols[name] = dict(powers=powers, b0=lr.intercept_, b=lr.coef_,
                         rmse=np.sqrt(np.mean(res ** 2)),
                         maxres=np.abs(res).max(), r2=lr.score(A, y),
                         resid=res, model=lr)
    return dict(re=re, d=d, theta=theta, y=y, lams=lams, path=path, std=std,
                raw=raw, intercept=intercept, rmse=rmse, active=active,
                cv=cv, sel=sel, ols=ols, scaler=scaler, model=model)


def style(ax):
    ax.grid(alpha=0.25, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def fig_lift(dfs):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for re, df in dfs.items():
        ok = df.verdict == "CONVERGED"
        for ax, col in ((a1, "cl"), (a2, "cd")):
            ax.plot(df.aoa[ok], df[col][ok], "o-", ms=4, lw=2, color=COLOR[re],
                    label=f"Re = {re}")
            if (~ok).any():
                ax.plot(df.aoa[~ok], df[col][~ok], "x", ms=8, mew=2, color=COLOR[re],
                        label=f"Re = {re}, not converged")
    for ax, lab in ((a1, r"$C_L$ (steady)"), (a2, r"$C_D$ (steady)")):
        ax.axvspan(0, AMAX, color="0.5", alpha=0.07, lw=0)
        ax.set(xlabel=r"$\alpha$ (deg)", ylabel=lab)
        style(ax)
    a1.text(0.3, a1.get_ylim()[1] * 0.95, f"fit range 0-{AMAX:g}°", color=INK, va="top")
    a1.set_title("Steady-branch lift", loc="left")
    a2.set_title("Steady-branch drag", loc="left")
    # one shared legend below both panels: inside an axis it covers the data
    h, l = a1.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=min(len(l), 4), frameon=False)
    fig.tight_layout(rect=(0, 0.06 * ((len(l) + 3) // 4) + 0.02, 1, 1))
    fig.savefig(OUT / f"steady_lift_drag_{TAG}.png", dpi=200)
    plt.close(fig)


def fig_path(r):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1]})
    hi = {POWERS.index(1): ("#2a78d6", r"$\theta$"),
          POWERS.index(3): ("#eb6834", r"$\theta^3$")}
    others = [p for p in POWERS if p not in (1, 3)]
    first_other = POWERS.index(others[0])
    for j, row in enumerate(r["path"]):
        if j in hi:
            a1.plot(r["lams"], row, color=hi[j][0], lw=2.2, label=hi[j][1], zorder=3)
        else:
            a1.plot(r["lams"], row, color=MUTED, lw=1.0, alpha=0.8,
                    label=rf"other $\theta^k$ (k = {', '.join(map(str, others))})"
                    if j == first_other else None)
    a1.axhline(0, color="black", lw=0.6)
    a1.set_ylabel("coefficient of standardized feature")
    a1.legend(frameon=False)
    a1.set_title(f"Re = {r['re']}: {LABEL} Lasso coefficient path (0-{AMAX:g}°)", loc="left")
    pred = r["y"].mean() + design(r["theta"])[0] @ r["path"]
    a2.plot(r["lams"], np.sqrt(np.mean((pred - r["y"][:, None]) ** 2, axis=0)),
            color="black", lw=1.5)
    a2.set(ylabel="training RMSE", xlabel=r"Lasso penalty $\lambda$; stronger $\rightarrow$")
    for ax in (a1, a2):
        ax.set_xscale("symlog", linthresh=1e-8)
        ax.axvline(LAMBDA, color=MUTED, ls="--")
        ax.axvline(r["cv"].alpha_, color=MUTED, ls=":")
        style(ax)
    a1.text(LAMBDA, 0.98, r"  $\lambda$=1e-4", transform=a1.get_xaxis_transform(),
            va="top", color=INK)
    a1.text(r["cv"].alpha_, 0.88, r"  $\lambda_{LOO}$", transform=a1.get_xaxis_transform(),
            va="top", color=INK)
    fig.tight_layout()
    fig.savefig(OUT / f"lasso_path_re{r['re']}.png", dpi=200)
    plt.close(fig)


def fig_compare(res):
    k = np.array(POWERS)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
    w = 0.76 / len(res)
    offs = (np.arange(len(res)) - (len(res) - 1) / 2) * w
    for off, r in zip(offs, res):
        a1.bar(k + off, r["std"], w * 0.92, color=COLOR[r["re"]], label=f"Re = {r['re']}")
        a2.bar(k + off, r["sel"], w * 0.92, color=COLOR[r["re"]], label=f"Re = {r['re']}")
    a1.axhline(0, color="black", lw=0.6)
    a1.set(xlabel=r"power $k$ of $\theta^k$", ylabel="standardized coefficient",
           xticks=k)
    a1.set_title(r"Lasso coefficients at $\lambda$ = 1e-4", loc="left")
    a2.set(xlabel=r"power $k$ of $\theta^k$", ylabel="fraction of 500 bootstraps",
           xticks=k, ylim=(0, 1.05))
    a2.set_title(r"Bootstrap selection frequency at $\lambda$ = 1e-4", loc="left")
    for ax in (a1, a2):
        style(ax)
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / f"lasso_coefficients_{TAG}.png", dpi=200)
    plt.close(fig)


def fig_cubic(res):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1]})
    g = np.linspace(0, AMAX, 400)
    tg = np.deg2rad(g)
    for r in res:
        o = r["ols"]["odd cubic"]
        c = COLOR[r["re"]]
        a1.plot(np.rad2deg(r["theta"]), r["y"], "o", ms=5, color=c,
                label=f"Re = {r['re']} steady")
        a1.plot(g, o["model"].predict(np.column_stack((tg, tg ** 3))), color=c, lw=2,
                label=rf"Re = {r['re']}: ${o['b0']:+.4f}{o['b'][0]:+.3f}\theta{o['b'][1]:+.2f}\theta^3$")
        a2.plot(np.rad2deg(r["theta"]), 1e3 * o["resid"], "o-", ms=4, lw=1.2, color=c,
                label=f"Re = {r['re']}")
        o5 = r["ols"]["odd quintic"]
        a2.plot(np.rad2deg(r["theta"]), 1e3 * o5["resid"], "--", lw=1.2, color=c,
                alpha=0.7, label=f"Re = {r['re']}, with $\\theta^5$")
    a2.axhline(0, color="black", lw=0.6)
    a1.set(ylabel=r"$C_L$")
    a1.set_title(r"OLS odd cubic $b_0 + b_1\theta + b_3\theta^3$ ($\theta$ in rad)", loc="left")
    a2.set(xlabel=r"$\alpha$ (deg)", ylabel=r"residual $\times 10^3$")
    for ax in (a1, a2):
        style(ax)
        ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / f"odd_cubic_fit_{TAG}.png", dpi=200)
    plt.close(fig)


def report(r):
    L = [f"=== Re = {r['re']}  (n = {len(r['y'])}, "
         f"{r['d'].aoa.min():g}-{r['d'].aoa.max():g} deg) ==="]
    L.append(f"Lasso lambda = {LAMBDA:g}: intercept {r['intercept']:+.9g}, "
             f"train RMSE {r['rmse']:.3e}")
    for j, pw in enumerate(POWERS):
        L.append(f"  theta^{pw}: std {r['std'][j]:+.6e}  raw {r['raw'][j]:+.6e}  "
                 f"bootstrap {r['sel'][j]:.3f}")
    for lam, act in r["active"].items():
        L.append(f"  active at lambda {lam:.0e}: " + ", ".join(f"theta^{k}" for k in act))
    cvm = r["cv"].mse_path_.mean(axis=1)
    L.append(f"LOO-CV lambda = {r['cv'].alpha_:.3e} (LOO RMSE {np.sqrt(cvm.min()):.3e}); "
             f"active: {[f'theta^{POWERS[j]}' for j in np.nonzero(np.abs(r['cv'].coef_) > 1e-10)[0]]}")
    for name, o in r["ols"].items():
        terms = " ".join(f"b{p}={b:+.6g}" for p, b in zip(o["powers"], o["b"]))
        L.append(f"OLS {name}: b0={o['b0']:+.3e} {terms}  R2={o['r2']:.8f}  "
                 f"RMSE={o['rmse']:.3e}  max|res|={o['maxres']:.3e}")
    c = r["ols"]["odd cubic"]
    L.append(f"  b1/(2 pi) = {c['b'][0] / (2 * np.pi):.4f};  b3/b1 = {c['b'][1] / c['b'][0]:+.4f}")
    return "\n".join(L)


def main():
    global AMAX, FIT_GRID, OUT, POWERS, LABEL, RES, TAG
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--amax", type=float, default=12.0,
                    help="upper end of the fit range in degrees (default 12)")
    ap.add_argument("--odd", action="store_true",
                    help="odd powers only (theta, theta^3, ..., theta^9)")
    ap.add_argument("--re", type=int, nargs="+", default=list(RES),
                    help="Reynolds numbers to compare (default 500 1000)")
    args = ap.parse_args()
    AMAX = args.amax
    RES = tuple(args.re)
    TAG = "_vs_".join(f"re{r}" for r in RES)
    if args.odd:
        POWERS, LABEL = [1, 3, 5, 7, 9], "odd-only degree-9"
    FIT_GRID = np.arange(0, AMAX + 1e-4, 0.5)
    if RES != (500, 1000):
        OUT = OUT.parent / f"resweep_steady_0to{AMAX:g}deg"
    elif AMAX != 12.0:  # keep the 0-12 outputs the report links to
        OUT = OUT.parent / f"re1000_steady_0to{AMAX:g}deg"
    if args.odd:
        OUT = OUT.parent / f"{OUT.name}_odd"
    OUT.mkdir(parents=True, exist_ok=True)
    dfs = {re: collect(re) for re in RES}
    for re, df in dfs.items():
        df.to_csv(OUT / f"steady_coeffs_re{re}.csv", index=False)
    res = [analyse(re, dfs[re]) for re in RES]
    for r in res:
        pd.DataFrame({"lambda": r["lams"], **{f"std_theta^{pw}": r["path"][j]
                                               for j, pw in enumerate(POWERS)}}
                     ).to_csv(OUT / f"lasso_path_re{r['re']}.csv", index=False)
    fig_lift(dfs)
    for r in res:
        fig_path(r)
    fig_compare(res)
    fig_cubic(res)
    text = "\n\n".join(report(r) for r in res)
    (OUT / "lasso_summary.txt").write_text(text + "\n")
    print(text)
    for re, df in dfs.items():
        bad = df[df.verdict != "CONVERGED"]
        if len(bad):
            print(f"\nRe={re} NOT converged (excluded from fits):\n{bad.to_string(index=False)}")


if __name__ == "__main__":
    main()
