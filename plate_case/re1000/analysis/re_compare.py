"""Plate, Re = 500 vs Re = 1000: does the C_L(alpha) shape hold?

Uses CONVERGED cases only (each case's convergence.json, written by
run_sweep.py).  At Re 1000, 8.5-12 deg did not converge (SIMPLE cycles), so
every fit is on 0-8 deg -- the same window as the airfoil cubic -- for BOTH Re.
All forms are odd in alpha with no intercept (symmetric section).
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import Lasso, lars_path
from sklearn.preprocessing import StandardScaler

PLATE = Path(__file__).resolve().parents[2]
RES = {500: "C3", 1000: "C1"}
HI = 8.0
FORMS = {"a·α + b·α³": (1, 3), "a·α + b·α⁵": (1, 5), "a·α + b·α³ + c·α⁵": (1, 3, 5)}
POWERS = np.array([1, 3, 5, 7, 9])


def load(re):
    rows, failed = [], []
    for c in (PLATE / f"re{re}").glob("steady_aoa*"):
        if not c.is_dir():
            continue
        s = json.loads((c / "convergence.json").read_text())
        a = float(c.name.removeprefix("steady_aoa"))
        (rows if s["status"] == "CONVERGED" else failed).append((a, s.get("cl", np.nan)))
    return np.array(sorted(rows)), sorted(failed)


def lsq(deg, cl, powers):
    t = np.deg2rad(deg)
    X = np.column_stack([t**p for p in powers])
    c, *_ = np.linalg.lstsq(X, cl, rcond=None)
    r = cl - X @ c
    nz = deg > 0
    return c, np.sqrt(np.mean(r**2)), np.max(np.abs(r[nz] / cl[nz]))


D = {re: load(re) for re in RES}
L = ["# Plate: does the lift-curve shape hold from Re 500 to Re 1000?", "",
     "Converged cases only; all fits on 0–8°, odd in α, no intercept, α in rad.", ""]
for re, (d, failed) in D.items():
    L.append(f"- Re {re}: {len(d)} converged angles"
             + (f"; NOT converged (excluded): {', '.join(f'{a:g}°' for a, _ in failed)}" if failed else ""))

L += ["", "| quantity | Re 500 | Re 1000 | ratio 1000/500 |", "|---|---:|---:|---:|"]
fits = {}
for re, (d, _) in D.items():
    deg, cl = d[:, 0], d[:, 1]
    m = deg <= HI + 1e-9
    fits[re] = {k: lsq(deg[m], cl[m], p) for k, p in FORMS.items()}
    fits[re]["secant"] = cl[deg > 0][0] / np.deg2rad(deg[deg > 0][0])
    i = np.argmax(cl)
    fits[re]["peak"] = (deg[i], cl[i])


def row(label, f):
    v5, v10 = f(fits[500]), f(fits[1000])
    L.append(f"| {label} | {v5:.5g} | {v10:.5g} | {v10 / v5:.3f} |")


row("slope at origin, C_L/α at 0.5° (per rad)", lambda F: F["secant"])
for k, p in FORMS.items():
    for j, pw in enumerate(p):
        row(f"{k}: coef of α^{pw}", lambda F, k=k, j=j: F[k][0][j])
    row(f"{k}: RMSE 0–8°", lambda F, k=k: F[k][1])
    row(f"{k}: max rel. error", lambda F, k=k: F[k][2])
L += ["", f"Highest converged C_L: Re 500 {fits[500]['peak'][1]:.4f} at {fits[500]['peak'][0]:g}°; "
      f"Re 1000 {fits[1000]['peak'][1]:.4f} at {fits[1000]['peak'][0]:g}° "
      "(Re 1000 values above 8° are unconverged and excluded)."]

# window sensitivity of the cubic and quintic b at each Re
L += ["", "Window sensitivity (two-term fits):", "",
      "| window | Re 500 α³ b | Re 1000 α³ b | Re 500 α⁵ b | Re 1000 α⁵ b |", "|---|---:|---:|---:|---:|"]
for hi in (4, 6, 8):
    vals = []
    for pw in (3, 5):
        for re in RES:
            d = D[re][0]; m = d[:, 0] <= hi + 1e-9
            vals.append(lsq(d[m, 0], d[m, 1], (1, pw))[0][1])
    L.append(f"| 0–{hi}° | {vals[0]:.2f} | {vals[1]:.2f} | {vals[2]:.0f} | {vals[3]:.0f} |")

# odd-powers Lasso path (same pipeline as re500/analysis/lasso_odd.py), on 0-8 deg
L += ["", "Odd-powers Lasso on 0–8° (θ…θ⁹, scaled not centred, no intercept): active sets along the path", ""]
paths = {}
for re in RES:
    d = D[re][0]; m = d[:, 0] <= HI + 1e-9
    t, y = np.deg2rad(d[m, 0]), d[m, 1]
    sc = StandardScaler(with_mean=False)
    X = sc.fit_transform(t[:, None] ** POWERS)
    pen, _, coef = lars_path(X, y, method="lasso")
    paths[re] = (pen, coef)
    seq = []
    for i in range(len(pen)):
        s = "{" + ", ".join(f"θ^{p}" for p, c in zip(POWERS, coef[:, i]) if c != 0) + "}"
        if not seq or seq[-1] != s:
            seq.append(s)
    L.append(f"- Re {re}: " + " → ".join(seq))
text = "\n".join(L)
print(text)
(Path(__file__).parent / "RE_COMPARE.md").write_text(text + "\n")

# ---- figure ---------------------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(18, 10))
g = np.linspace(0, 12, 300); gt = np.deg2rad(g)
for re, c in RES.items():
    d, failed = D[re]
    ax[0, 0].plot(d[:, 0], d[:, 1], "o", ms=4, color=c, label=f"Re {re} converged")
    if failed:
        f = np.array(failed)
        ax[0, 0].plot(f[:, 0], f[:, 1], "x", ms=6, color=c, label=f"Re {re} NOT converged (last iterate)")
    co = fits[re]["a·α + b·α³ + c·α⁵"][0]
    ax[0, 0].plot(g, co[0] * gt + co[1] * gt**3 + co[2] * gt**5, "-", lw=1, color=c,
                  label=f"Re {re}: {co[0]:.3f}α {co[1]:+.1f}α³ {co[2]:+.0f}α⁵")
    nz = d[:, 0] > 0
    ax[0, 1].plot(d[nz, 0], d[nz, 1] / np.deg2rad(d[nz, 0]), "o", ms=4, color=c, label=f"Re {re}")
    a0 = fits[re]["secant"]
    ax[0, 2].plot(d[nz, 0], d[nz, 1] / (a0 * np.deg2rad(d[nz, 0])), "o", ms=4, color=c, label=f"Re {re}")
    for j, (k, p) in enumerate(FORMS.items()):
        co = fits[re][k][0]
        m = d[:, 0] <= HI
        pred = sum(ci * np.deg2rad(d[m, 0]) ** pw for ci, pw in zip(co, p))
        ax[1, j].plot(d[m, 0], d[m, 1] - pred, "o-", ms=3, lw=0.8, color=c, label=f"Re {re}")
        ax[1, j].set(title=f"residual: {k} (0–8°)", xlabel="α (deg)")
ax[0, 0].set(xlabel="α (deg)", ylabel="$C_L$", title="Plate C_L: 3-term odd fit on 0–8°")
ax[0, 1].set(xlabel="α (deg)", ylabel="$C_L/α$ (per rad)", title="Secant slope")
ax[0, 2].set(xlabel="α (deg)", ylabel="$C_L/(a_0 α)$", title="Normalised by each Re's slope at origin")
for x in ax.flat:
    x.grid(alpha=0.25); x.legend(fontsize=7)
for x in ax[1]:
    x.axhline(0, color="k", lw=0.6)
fig.suptitle("Thick flat plate: Re 500 vs Re 1000 (steady SIMPLE)")
fig.tight_layout()
fig.savefig(Path(__file__).parent / "re_compare.png", dpi=150)
