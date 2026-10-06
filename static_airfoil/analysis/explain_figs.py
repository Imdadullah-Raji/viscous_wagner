"""Two explanatory figures for the steady-lift fitting report.

Run with: /home/raji/miniconda3/envs/research/bin/python analysis/explain_figs.py

  symmetry.png    C_L(alpha) is odd for a symmetric section: computed points for
                  alpha >= 0 and their mirror images for alpha < 0 (Re = 1000).
  lookalike.png   (a) theta, theta^2, theta^3 over 0-8 deg, each scaled to 1 at
                  8 deg: all rise monotonically, a straight line explains most of
                  each, and their pairwise correlations are 0.91-0.99.
                  (b) P1, P3, P5 on the mirrored, rescaled axis: visibly
                  different shapes.
-> analysis/figs/partner_report/
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import legendre as L

from lasso_re1000 import collect, ROOT

OUT = ROOT / "analysis" / "figs" / "partner_report"
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3
INK, MUTED = "#52514e", "#9a9893"


def style(ax):
    ax.grid(alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def symmetry():
    df = collect(1000)
    d = df[(df.aoa <= 8) & (df.verdict == "CONVERGED")]
    a, cl = d.aoa.to_numpy(), d.cl.to_numpy()
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.6)
    ax.plot(a, cl, "o", color=S1, ms=6, label="computed, α ≥ 0")
    ax.plot(-a[a > 0], -cl[a > 0], "o", mfc="none", mec=S1, mew=1.5, ms=6,
            label="mirror image: $C_L(-α) = -C_L(α)$")
    ax.annotate("same section, flipped:\nlift reverses sign", xy=(-6, -cl[a == 6][0]),
                xytext=(-7.6, 0.12), color=INK, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=MUTED))
    ax.set(xlabel="angle of attack α (deg)", ylabel="$C_L$ (steady)",
           xlim=(-8.8, 8.8))
    ax.set_title("NACA0012 at Re = 1000: lift is an odd function of angle", loc="left")
    style(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "symmetry.png", dpi=200)
    plt.close(fig)


def lookalike():
    deg = np.linspace(0, 8, 200)
    th = np.deg2rad(deg)
    grid = np.deg2rad(np.arange(0, 8.0001, 0.5))
    corr = lambda p, q: np.corrcoef(grid ** p, grid ** q)[0, 1]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, c in zip((1, 2, 3), (S1, S2, S3)):
        y = (th / th[-1]) ** k
        lab = r"$\theta$" if k == 1 else rf"$\theta^{k}$"
        a1.plot(deg, y, color=c, lw=2.2, label=f"{lab}, scaled to 1 at 8°")
    a1.text(0.3, 0.93, f"correlation over the 0–8° data grid:\n"
            rf"$\theta$ vs $\theta^2$: {corr(1, 2):.2f}    "
            rf"$\theta$ vs $\theta^3$: {corr(1, 3):.2f}" "\n"
            rf"$\theta^2$ vs $\theta^3$: {corr(2, 3):.2f}    "
            rf"$\theta^3$ vs $\theta^5$: {corr(3, 5):.2f}",
            transform=a1.transAxes, va="top", fontsize=9, color=INK)
    a1.set(xlabel="angle α (deg)", ylabel="scaled value", ylim=(-0.05, 1.12))
    a1.set_title("Raw powers of θ: all rise together", loc="left")
    a1.legend(frameon=False, loc="lower right", fontsize=9)

    x = np.linspace(-1, 1, 400)
    for k, c in zip((1, 3, 5), (S1, S2, S3)):
        a2.plot(x, L.legval(x, np.eye(k + 1)[k]), color=c, lw=2.2, label=rf"$P_{k}(x)$")
    a2.axhline(0, color="black", lw=0.6)
    a2.set(xlabel=r"$x = \theta/\theta_{max}$, mirrored to [−1, 1]", ylabel="value",
           ylim=(-1.15, 1.15))
    a2.set_title("Odd Legendre polynomials: distinct shapes", loc="left")
    a2.legend(frameon=False, loc="lower right", fontsize=9)
    for ax in (a1, a2):
        style(ax)
    fig.tight_layout()
    fig.savefig(OUT / "lookalike.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    symmetry()
    lookalike()
