#!/usr/bin/env python3
"""C_L history plots and the period-2 saturation test.

Three figures:
  cl_history.png        -- C_L(t) for both angles, whole run, phase boundaries
  cl_zoom.png           -- last ~6 shedding periods, where the period doubling
                           at AoA 30 is visible as alternating peak heights
  aoa30_period2.png     -- is the period-2 component saturated, growing or
                           decaying?  Two independent measures on a log axis.

The period-2 measures are deliberately independent of each other:
  (a) peak-height alternation |A_n - A_{n+1}| between consecutive C_L maxima --
      a purely time-domain quantity, no spectral assumption;
  (b) the amplitude of the f0/2 Fourier component in a sliding window of 6
      fundamental periods -- a purely spectral quantity, insensitive to which
      individual peak is picked.
AoA 20 is carried through both as the period-1 control: whatever it reads is
the noise floor of the measure, not a signal.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from limit_cycle import load_series, dominant_freq
from limit_cycle import spectral_structure as structure

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGS = os.path.join(HERE, "figs")

# t_conv is read from the case's .phase2 marker where one exists; the two runs
# that finished before that marker was introduced carry it here.
CASES = [("re500_aoa14", r"AoA $14^\circ$", 60.5, "#bcbd22"),
         ("re500_aoa16", r"AoA $16^\circ$", 60.5, "#7f7f7f"),
         ("re500_aoa18", r"AoA $18^\circ$", 60.5, "#17becf"),
         ("re500_aoa20", r"AoA $20^\circ$", 60.5, "#1f77b4"),
         ("re500_aoa22", r"AoA $22^\circ$", 80.5, "#2ca02c"),
         ("re500_aoa24", r"AoA $24^\circ$", 60.5, "#ff7f0e"),
         ("re500_aoa26", r"AoA $26^\circ$", 100.5, "#8c564b"),
         ("re500_aoa28", r"AoA $28^\circ$", 100.5, "#9467bd"),
         ("re500_aoa30", r"AoA $30^\circ$", 120.5, "#d62728")]


def t_conv_of(case, fallback):
    p = os.path.join(ROOT, "NACA0012", case, ".phase2")
    if os.path.exists(p):
        for line in open(p):
            if line.startswith("T_CONV="):
                return float(line.split("=", 1)[1])
    return fallback


def uniform(t, y, dt=None):
    """Resample onto a uniform grid.  dt changes between chunks, so the raw
    series is only piecewise-uniform and cannot be handed to an FFT."""
    dt = dt or float(np.median(np.diff(t)))
    tu = np.arange(t[0], t[-1], dt)
    return tu, np.interp(tu, t, y)


def carrier_freq(t, cl, t_conv):
    """SHEDDING frequency measured on the converged part only.

    Must not be `dominant_freq` (the largest peak): at AoA 28 the largest peak is
    the orbit fundamental at 0.1877, not the carrier at 0.3754, so every window
    and sampling decision derived from it would be off by a factor of two.
    See FINDINGS.md section 12.
    """
    sel = t >= t_conv
    return structure(t[sel], cl[sel])[1]


def sliding_subharmonic(t, y, f0, n_per=6.0, step_frac=0.25):
    """Amplitude of the f0/2 and f0 components in a sliding window.

    Window = n_per fundamental (2/f0) periods so the subharmonic bin is well
    resolved; the window is an integer number of subharmonic periods, so a true
    period-2 signal lands on a bin centre and does not leak.
    """
    tu, yu = uniform(t, y)
    dt = tu[1] - tu[0]
    win = int(round(n_per * 2.0 / f0 / dt))
    step = max(1, int(win * step_frac))
    tc, a_sub, a_car = [], [], []
    w = np.hanning(win)
    cg = w.sum() / win                      # coherent gain of the window
    for i in range(0, len(tu) - win, step):
        seg = yu[i:i + win]
        seg = seg - seg.mean()
        Y = np.abs(np.fft.rfft(seg * w)) * 2.0 / (win * cg)
        f = np.fft.rfftfreq(win, dt)
        def amp_at(target):
            k = int(np.argmin(np.abs(f - target)))
            return float(Y[max(1, k - 1):k + 2].max())   # +/-1 bin, ignore DC
        tc.append(tu[i + win // 2])
        a_sub.append(amp_at(f0 / 2.0))
        a_car.append(amp_at(f0))
    return np.array(tc), np.array(a_sub), np.array(a_car)


def peak_alternation(t, cl, f0):
    """|A_n - A_{n+1}| for consecutive C_L maxima, timestamped at the midpoint."""
    tu, yu = uniform(t, cl)
    dt = tu[1] - tu[0]
    idx, _ = find_peaks(yu, distance=int(0.6 / f0 / dt))
    # Drop maxima within half a carrier period of either end: the series is cut
    # mid-cycle, so the last "peak" is often just the highest point of a
    # truncated rise and reads as a spurious alternation.
    edge = int(0.5 / f0 / dt)
    idx = idx[(idx > edge) & (idx < len(tu) - edge)]
    if len(idx) < 3:
        return np.array([]), np.array([]), np.array([]), np.array([])
    tp, ap = tu[idx], yu[idx]
    return tp, ap, 0.5 * (tp[:-1] + tp[1:]), np.abs(np.diff(ap))


def main():
    os.makedirs(FIGS, exist_ok=True)
    data = {}
    for case, label, t_conv, colour in CASES:
        t, cl, cd = load_series(os.path.join(ROOT, "NACA0012", case))
        t_conv = t_conv_of(case, t_conv)
        f0 = carrier_freq(t, cl, t_conv)
        data[case] = dict(t=t, cl=cl, cd=cd, f0=f0, label=label,
                          t_conv=t_conv, colour=colour)
        print(f"{case}: carrier f0={f0:.5f} (T={1/f0:.4f}), "
              f"t_conv={t_conv}, t_end={t[-1]:.2f}")

    # ---- Figure 1: full C_L history -------------------------------------
    # Two columns sharing y: the impulsive start needs a log time axis (the
    # added-mass spike lives at t < 1e-3 and reaches O(1e4)), the limit cycle
    # needs a linear one.  An inset would sit on top of the trace it explains.
    n = len(CASES)
    fig, axes = plt.subplots(n, 2, figsize=(12, 3.4 * n), sharex="col",
                             gridspec_kw=dict(width_ratios=[1, 5], wspace=0.16))
    for (axs, axm), (case, *_) in zip(axes, CASES):
        d = data[case]
        body = d["cl"][d["t"] >= 1.0]
        lo, hi = body.min(), body.max()
        pad = 0.10 * (hi - lo)

        m = (d["t"] > 0) & (d["t"] < 1.0)
        axs.semilogx(d["t"][m], d["cl"][m], lw=0.7, color=d["colour"])
        axs.set_ylabel(r"$C_L$")
        axs.grid(alpha=0.3)
        axs.set_title(f"impulsive start (peak {d['cl'].max():.0f})", fontsize=8)

        axm.plot(d["t"], d["cl"], lw=0.4, color=d["colour"])
        axm.set_ylim(lo - pad, hi + pad)
        axm.axvline(d["t_conv"], color="k", ls="--", lw=1)
        axm.axvspan(d["t_conv"], d["t"][-1], color="k", alpha=0.05)
        axm.text(d["t_conv"], hi + pad, "  limit cycle / fields written",
                 va="top", ha="left", fontsize=8)
        axm.set_xlim(0, max(dd["t"][-1] for dd in data.values()) * 1.01)
        # separate y scale from the startup panel -- keep its own tick labels
        axm.grid(alpha=0.3)
        axm.set_title(f"{d['label']}   Re=500   "
                      rf"$\overline{{C_L}}$="
                      rf"{d['cl'][d['t']>=d['t_conv']].mean():.4f},"
                      rf"  $f_0$={d['f0']:.4f}  ($T$={1/d['f0']:.3f})",
                      fontsize=10)
    axes[-1][0].set_xlabel(r"$t$")
    axes[-1][1].set_xlabel(r"$t$  (convective times, $c/U_\infty$)")
    fig.suptitle("Lift coefficient history, static NACA0012, Re=500", y=0.98)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "cl_history.png"), dpi=150)

    # ---- Figure 2: zoom on the limit cycle -------------------------------
    fig, axes = plt.subplots(len(CASES), 1, figsize=(11, 2.6 * len(CASES)))
    for ax, (case, *_) in zip(axes, CASES):
        d = data[case]
        span = 6.0 / d["f0"]
        sel = d["t"] >= d["t"][-1] - span
        ax.plot(d["t"][sel], d["cl"][sel], lw=1.0, color=d["colour"])
        tp, ap, _, _ = peak_alternation(d["t"][d["t"] >= d["t_conv"]],
                                        d["cl"][d["t"] >= d["t_conv"]], d["f0"])
        m = tp >= d["t"][-1] - span
        ax.plot(tp[m], ap[m], "kv", ms=5, mfc="none")
        ax.set_ylabel(r"$C_L$")
        ax.set_xlabel(r"$t$")
        ax.grid(alpha=0.3)
        ax.set_title(f"{d['label']}: last 6 carrier periods "
                     f"(markers = $C_L$ maxima)", fontsize=10)
    fig.suptitle("Period-1 at 20-24$^\\circ$, period-2 at 30$^\\circ$: "
                 "alternating peak heights", y=0.995)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "cl_zoom.png"), dpi=150)

    # ---- Figure 3: is the period-2 state saturated? ----------------------
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for case, *_ in CASES:
        d = data[case]
        f0 = d["f0"]
        tc, a_sub, a_car = sliding_subharmonic(d["t"], d["cl"], f0)
        axes[0].semilogy(tc, a_sub, color=d["colour"], lw=1.4,
                         label=f"{d['label']}  $A(f_0/2)$")
        axes[0].semilogy(tc, a_car, color=d["colour"], lw=0.9, ls=":",
                         alpha=0.7, label=f"{d['label']}  $A(f_0)$")
        _, _, tm, dA = peak_alternation(d["t"], d["cl"], f0)
        m = dA > 0
        axes[1].semilogy(tm[m], dA[m], ".", ms=2.5, color=d["colour"],
                         label=d["label"])
        axes[0].axvline(d["t_conv"], color=d["colour"], ls="--", lw=0.8)
        axes[1].axvline(d["t_conv"], color=d["colour"], ls="--", lw=0.8)
        d["sub"] = (tc, a_sub, a_car)

    axes[0].set_ylabel("component amplitude")
    axes[0].set_title("Sliding-window Fourier amplitude at the subharmonic "
                      r"$f_0/2$ and the carrier $f_0$ (6-period window)",
                      fontsize=10)
    axes[0].legend(fontsize=8, ncol=2); axes[0].grid(alpha=0.3, which="both")
    axes[1].set_ylabel(r"$|A_n - A_{n+1}|$")
    axes[1].set_xlabel(r"$t$")
    axes[1].set_title(r"Peak-height alternation between consecutive $C_L$ maxima",
                      fontsize=10)
    axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3, which="both")
    fig.suptitle("Growth and saturation of the period-2 component", y=0.99)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "subharmonic_growth.png"), dpi=150)

    # ---- Figure 4: every angle on one orbit-normalised axis ---------------
    # x is time in ORBITS (fundamental periods), not convective times, because
    # the orbit period runs from 1.92 at AoA 20 to 5.33 at AoA 28 -- a factor of
    # 2.8.  On a shared real-time axis the low angles would be a blur.  Two
    # orbits of each therefore occupy the same width, and the period-2 cases
    # visibly need TWO oscillations to repeat while the period-1 cases need one.
    # Each trace is phase-aligned on its tallest C_L maximum so the alternation
    # always starts on the strong half-cycle.
    NORB = 2.2
    fig, ax = plt.subplots(figsize=(12, 7.2))
    for case, *_ in CASES:
        d = data[case]
        t, cl = d["t"], d["cl"]
        sel = t >= d["t_conv"]
        t, cl = t[sel], cl[sel]
        f_fund, f_car, order, ratio = structure(t, cl)
        T_orb = 1.0 / f_fund
        tu, yu = uniform(t, cl)
        dt = tu[1] - tu[0]
        idx, _ = find_peaks(yu, distance=int(0.6 / f_car / dt))
        # tallest peak that still leaves NORB orbits of trace after it
        ok = idx[tu[idx] <= tu[-1] - NORB * T_orb]
        k = ok[int(np.argmax(yu[ok]))] if len(ok) else idx[0]
        t0 = tu[k]
        m = (tu >= t0) & (tu <= t0 + NORB * T_orb)
        lab = (f"{d['label']}  period-{order},  "
               rf"$f_{{shed}}${f_car:.4f},  $T_{{orb}}${T_orb:.3f}")
        ax.plot((tu[m] - t0) / T_orb, yu[m], lw=1.6, color=d["colour"], label=lab)
    for x in (1.0, 2.0):
        ax.axvline(x, color="k", ls=":", lw=0.9, alpha=0.6)
    ax.annotate("orbit closes", xy=(1.0, ax.get_ylim()[0]),
                xytext=(1.02, ax.get_ylim()[0] + 0.03 * np.ptp(ax.get_ylim())),
                fontsize=8, color="0.3")
    ax.set_xlabel("time in orbits  ($t / T_{orbit}$, aligned on the tallest $C_L$ peak)")
    ax.set_ylabel(r"$C_L$")
    ax.set_title("Re=500 static NACA0012: $C_L$ over two orbits at every angle\n"
                 "period-1 repeats every oscillation; period-2 alternates and "
                 "repeats only after two", fontsize=11)
    ax.grid(alpha=0.3)
    # Legend below the axes: at upper-right it sat squarely on AoA 30's peak.
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.11),
              ncol=3, framealpha=0.95)
    ax.set_xlim(0, NORB)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "cl_all_angles.png"), dpi=160)

    # ---- quantitative verdict -------------------------------------------
    print("\n--- period-2 saturation test (converged window only) ---")
    for case, *_ in CASES:
        d = data[case]
        tc, a_sub, a_car = d["sub"]
        m = tc >= d["t_conv"] + 2.0 / d["f0"]     # drop the window straddling t_conv
        if m.sum() < 3:
            print(f"{case}: too few windows past t_conv"); continue
        x, y = tc[m], a_sub[m]
        # exponential fit: log A = log A0 + sigma t  -> sigma is the growth rate
        sigma, logA0 = np.polyfit(x, np.log(y), 1)
        span = x[-1] - x[0]
        ratio = (a_sub[m] / a_car[m]).mean()
        print(f"{case}: A(f0/2) = {y.mean():.6f} +/- {y.std():.6f}   "
              f"ratio A(f0/2)/A(f0) = {ratio:.5f}")
        # "Saturation" only means something for a state that HAS a subharmonic.
        # Applied to a period-1 case it fits an exponential to the noise floor and
        # prints "never settled", which reads like a problem and is not one.
        if ratio < 0.01:
            print(f"    period-1: subharmonic at the noise floor, "
                  f"{1/ratio:.0f}x below the carrier -- saturation test N/A")
            print(f"    decay rate over this window: sigma = {sigma:+.3e} /t"
                  + (f"  (e-fold {abs(1/sigma):.0f} t, still falling)"
                     if sigma < -1e-3 else "  (flat -- at the floor)"))
            continue
        print(f"    growth rate sigma = {sigma:+.3e} /t  -> "
              f"{100*(np.exp(sigma*span)-1):+.2f}% change over the {span:.1f}t window "
              f"(e-fold time {abs(1/sigma):.0f}t)" if sigma != 0 else "")
        # when did it get within 1% / 0.1% of its saturated value, and stay?
        sat = y.mean()
        for frac in (0.05, 0.01, 0.001):
            near = np.abs(a_sub - sat) <= frac * sat
            k = len(near) - 1
            while k > 0 and near[k - 1]:
                k -= 1
            print(f"    within {100*frac:5.1f}% of saturation from t = {tc[k]:.1f}"
                  if near[-1] else f"    never settled to {100*frac:.1f}%")
    print(f"\nfigures written to {FIGS}/")


if __name__ == "__main__":
    main()
