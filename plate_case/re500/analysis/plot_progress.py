"""Interim plate C_L vs the NACA0012 steady data and its 0-8 deg odd cubic.
Reads only converged plate cases (lift_coeffs.csv is written from converged
stats by run_sweep.py)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plate = np.loadtxt(ROOT / "lift_coeffs.csv", skiprows=1, ndmin=2)
foil = np.loadtxt(ROOT.parents[1] / "static_airfoil" / "lift_coeffs.csv", skiprows=1, ndmin=2)
A_F, B_F = 3.267395529, -28.235887546          # airfoil odd cubic, fit 0-8 deg, theta in rad

th = np.deg2rad(plate[:, 0]); cl = plate[:, 1]
secant = cl / th                                   # = a + b*theta^2 for an odd cubic
print("secant Cl/theta per rad:", dict(zip(plate[:, 0], np.round(secant, 4))))
lo = plate[:, 0] <= 9.0
X = np.column_stack((th[lo], th[lo] ** 3))
(a, b), *_ = np.linalg.lstsq(X, cl[lo], rcond=None)
print(f"provisional odd cubic on {plate[lo,0].min():g}-{plate[lo,0].max():g} deg: a={a:.4f} b={b:.4f}")

g = np.linspace(0, 12, 300); gt = np.deg2rad(g)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(foil[:, 0], foil[:, 1], "o", ms=4, color="C0", label="NACA0012 steady (computed)")
ax[0].plot(g, A_F * gt + B_F * gt**3, "-", color="C0", lw=1,
           label=f"NACA0012 fit 0–8°: {A_F:.3f}α {B_F:+.2f}α³")
ax[0].plot(plate[:, 0], cl, "s", ms=5, color="C3", label="plate steady (converged so far)")
ax[0].plot(g, a * gt + b * gt**3, "--", color="C3", lw=1,
           label=f"plate PROVISIONAL fit ≤9°: {a:.2f}α {b:+.1f}α³ (extrapolated <7°)")
ax[0].plot(g, 2 * np.pi * gt, ":", color="gray", lw=1, label="thin-airfoil 2πα")
ax[0].axvspan(0, plate[:, 0].min(), color="C3", alpha=0.06)
ax[0].set(xlabel="α (deg)", ylabel="$C_L$", ylim=(0, 0.9), title="Re = 500, steady SIMPLE")
ax[0].legend(fontsize=8, loc="upper left")
nz = foil[:, 0] > 0
ax[1].plot(foil[nz, 0], foil[nz, 1] / np.deg2rad(foil[nz, 0]), "o", ms=4, color="C0",
           label="NACA0012")
ax[1].plot(plate[:, 0], secant, "s", ms=5, color="C3", label="plate")
ax[1].axhline(A_F, color="C0", lw=1, ls="-", label=f"NACA0012 a = {A_F:.3f}")
ax[1].set(xlabel="α (deg)", ylabel="$C_L/α$ (per rad)", xlim=(0, 12.5),
          title="Secant slope: → a as α → 0 (odd cubic: a + bα²)")
ax[1].legend(fontsize=8)
for x in ax: x.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(ROOT / "analysis" / "cl_progress.png", dpi=150)
