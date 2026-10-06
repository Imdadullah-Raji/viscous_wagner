"""C_L(t) and amplitude spectrum of the unsteady Re=2000, AoA 9.5 run.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/plot_unsteady_re2000.py

Force series is split across postProcessing/forceCoeffs/<restart>/ dirs with
duplicated times at the seams: concatenate, sort, keep the LAST sample at each
time (same rule as analysis/limit_cycle.py).
-> analysis/figs/unsteady_re2000/cl_aoa9.5.png
"""
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "unsteady_re2000" / "re2000_aoa9.5"
OUT = ROOT / "analysis" / "figs" / "unsteady_re2000"
CL_STEADY_9 = 0.246865      # last CONVERGED steady angle, steady_re2000/steady_aoa9
COL, INK, MUTED = "#eda100", "#52514e", "#9a9893"   # Re=2000 slot


def load():
    t, cd, cl = [], [], []
    for f in sorted(glob.glob(str(CASE / "postProcessing/forceCoeffs/*/forceCoeffs.dat")),
                    key=lambda p: float(Path(p).parent.name)):
        a = np.loadtxt(f, comments="#", usecols=(0, 2, 3))
        t.append(a[:, 0]); cd.append(a[:, 1]); cl.append(a[:, 2])
    t, cd, cl = map(np.concatenate, (t, cd, cl))
    # keep last sample per time: reverse, unique keeps first occurrence
    _, idx = np.unique(t[::-1], return_index=True)
    keep = len(t) - 1 - idx
    o = np.argsort(t[keep])
    return t[keep][o], cd[keep][o], cl[keep][o]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t, cd, cl = load()
    w = t >= 20.0                                   # saturated window
    tw, yw = t[w], cl[w] - cl[w].mean()
    dt = np.median(np.diff(tw))
    n = len(yw)
    win = np.hanning(n)
    # The window is only ~10 t long (bin width ~0.1), so zero-pad 16x to
    # interpolate the peak; this refines its location, not the resolution.
    npad = 16 * n
    A = 2 * np.abs(np.fft.rfft(yw * win, npad)) / win.sum()
    f = np.fft.rfftfreq(npad, dt)
    fpk = f[1:][np.argmax(A[1:])]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.4),
                                 gridspec_kw={"width_ratios": [2.2, 1]})
    a1.plot(t, cl, color=COL, lw=0.8, label="unsteady $C_L(t)$")
    a1.axhline(cl[w].mean(), color=INK, lw=1, ls="--",
               label=f"mean, t ≥ 20: {cl[w].mean():.4f}")
    a1.axhline(CL_STEADY_9, color=MUTED, lw=1.2, ls=":",
               label=f"steady (SIMPLE) at 9°: {CL_STEADY_9:.4f}")
    a1.axvspan(20, t[-1], color="0.5", alpha=0.07, lw=0)
    # The impulsive start from uniform flow puts |C_L| up to ~6e3 in the first
    # ~17 steps (t < 1e-3); scale the axis to t > 0.5 so the physics is visible.
    late = cl[t > 0.5]
    pad = 0.08 * np.ptp(late)
    a1.set(xlabel="t (convective units)", ylabel="$C_L$", xlim=(0, t[-1]),
           ylim=(late.min() - pad, late.max() + pad))
    a1.text(0.2, late.max() + 0.5 * pad, f"impulsive-start spike (t < 1e-3, |$C_L$| ≤ "
            f"{np.abs(cl).max():.0f}) off scale", color=INK, fontsize=8, va="center")
    a1.set_title("Re = 2000, α = 9.5°: lift history (stopped once shedding was confirmed)",
                 loc="left")
    a2.plot(f, A, color=COL, lw=1.5)
    a2.set(xlim=(0, 4 * fpk), xlabel="f", ylabel="$C_L$ amplitude")
    a2.set_title(f"Spectrum, t ≥ 20: peak f = {fpk:.3f}", loc="left")
    for ax in (a1, a2):
        ax.grid(alpha=0.25, lw=0.6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    a1.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "cl_aoa9.5.png", dpi=200)
    plt.close(fig)
    print(f"t_end={t[-1]:.3f} n={len(t)} dt={dt:.4g}")
    print(f"t>=20: Cl mean {cl[w].mean():.5f} min {cl[w].min():.5f} max {cl[w].max():.5f} "
          f"p2p {np.ptp(cl[w]):.5f}  Cd mean {cd[w].mean():.5f}  f_peak {fpk:.4f}")


if __name__ == "__main__":
    main()
