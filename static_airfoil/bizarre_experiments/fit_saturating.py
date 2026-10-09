"""Fit A*x/sqrt(1+B*x²) on the original 0–8° CSV data; x in radians."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
sys.path.insert(0, str(ROOT / 'analysis'))
from steady_summary import load, verdict


def model(x, a, b):
    return a*x / np.sqrt(1+b*x*x)


def main():
    data = pd.read_csv(ROOT/'lift_coeffs.csv', sep=r'\s+')
    train = data[(data.iloc[:, 0] >= 0) & (data.iloc[:, 0] <= 8)]
    degrees = train.iloc[:, 0].to_numpy()
    x = np.deg2rad(degrees)
    y = train.iloc[:, 1].to_numpy()
    # Require real predictions throughout the 0–24° validation interval.
    # Allow negative B; the data, rather than a positivity constraint, select saturation.
    lower_b = -1 / np.deg2rad(24)**2 + 1e-8
    solutions = [least_squares(lambda p: model(x, *p)-y, [3.3, b],
                               bounds=([-np.inf, lower_b], [np.inf, np.inf]),
                               xtol=1e-13, ftol=1e-13, gtol=1e-13)
                 for b in (0, 10, 30, 100)]
    if not all(s.success for s in solutions):
        raise RuntimeError('An optimization failed')
    best = min(solutions, key=lambda s: s.cost)
    np.testing.assert_allclose([s.cost for s in solutions], best.cost, atol=1e-12)
    a, b = best.x
    cubic = np.linalg.lstsq(np.column_stack((np.ones(len(x)), x, x**3)), y, rcond=None)[0]
    rows = [(d, cl, 'training', 'CONVERGED') for d, cl in zip(degrees, y)]
    case_root = ROOT/'NACA0012'/'steady_re500' if (ROOT/'NACA0012'/'steady_re500').is_dir() else ROOT
    for case in case_root.glob('steady_aoa*'):
        if not case.is_dir():
            continue
        d = float(case.name.removeprefix('steady_aoa'))
        if not 0 <= d <= 24 or np.any(np.isclose(d, degrees)):
            continue
        trace = load(str(case))
        if trace is None:
            continue
        it, _, cl = trace
        state, _ = verdict(it, cl, 1e-4)
        rows.append((d, cl[-1], 'interpolation' if d < 8 else
                     'extrapolation_8_12' if d <= 12 else 'extrapolation_12_24', state))
    df = pd.DataFrame(rows, columns=['angle_deg', 'observed', 'group', 'status']).sort_values('angle_deg')
    if df.angle_deg.max() < 24:
        raise RuntimeError('Expected Re=500 validation cases through 24 degrees')
    theta = np.deg2rad(df.angle_deg.to_numpy())
    df['predicted'] = model(theta, a, b)
    df['residual'] = df.observed-df.predicted
    df['cubic_predicted'] = cubic[0]+cubic[1]*theta+cubic[2]*theta**3
    df['cubic_residual'] = df.observed-df.cubic_predicted
    df.to_csv(OUT/'saturating_predictions.csv', index=False)
    lines = ['# Saturating lift fit', '',
             'Unweighted nonlinear least squares on the 17 original CSV observations from 0–8°. Angle x is in radians; no intercept.', '',
             f'`Cl = {a:.12g} * x / sqrt(1 + {b:.12g} * x^2)`', '',
             f'- A = {a:.12g}', f'- B = {b:.12g}',
             f'- Training R² = {1-np.sum((model(x,a,b)-y)**2)/np.sum((y-y.mean())**2):.12g}',
             f'- Positive-angle asymptote A/√B = {a/np.sqrt(b):.12g}' if b > 0 else '- No positive-B saturation.',
             '- Four initial guesses converged to the same objective value.', '',
             '| Set (converged cases only) | n | Model RMSE | Cubic RMSE |',
             '|---|---:|---:|---:|']
    for name, group in df[df.status == 'CONVERGED'].groupby('group', sort=False):
        lines.append(f'| {name} | {len(group)} | {np.sqrt(np.mean(group.residual**2)):.9g} | {np.sqrt(np.mean(group.cubic_residual**2)):.9g} |')
    lines += ['', 'The comparison cubic includes an intercept and is independently fitted to the same 17 training points.',
              'The 22° and 24° saved SIMPLE traces are unconverged; their final-iteration values are provisional and excluded from RMSE summaries.',
              'Saturation is a model assumption, not a demonstrated large-angle physical limit.', '',
              '| Angle | Observed | Predicted | Observed − predicted | Status |',
              '|---|---:|---:|---:|---|']
    for row in df[df.angle_deg >= 8].itertuples():
        lines.append(f'| {row.angle_deg:g} | {row.observed:.8f} | {row.predicted:.8f} | {row.residual:+.8f} | {row.status} |')
    lines += ['', 'Run: `conda run -n research python bizarre_experiments/fit_saturating.py`', '',
              '[Diagnostics through 12°](saturating_diagnostics_12deg.png) · [Extrapolation through 24°](saturating_extrapolation_24deg.png) · [CSV](saturating_predictions.csv)']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    for limit, filename in [(12, 'saturating_diagnostics_12deg.png'), (24, 'saturating_extrapolation_24deg.png')]:
        fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
        grid = np.linspace(0, limit, 600)
        t = np.deg2rad(grid)
        axes[0].plot(grid, model(t,a,b), label=r'$Ax/\sqrt{1+Bx^2}$', color='C0')
        axes[0].plot(grid, cubic[0]+cubic[1]*t+cubic[2]*t**3, '--', color='C3', label='Intercept + x + x³')
        subset = df[df.angle_deg <= limit]
        for name, marker, color in [('training','o','black'), ('interpolation','x','C2'),
                                    ('extrapolation_8_12','s','C1'), ('extrapolation_12_24','s','C1')]:
            group = subset[(subset.group == name) & (subset.status == 'CONVERGED')]
            if len(group):
                axes[0].scatter(group.angle_deg, group.observed, s=25, marker=marker, color=color, label=name.replace('_',' '))
                axes[1].scatter(group.angle_deg, group.residual, s=25, marker=marker, color=color)
        axes[1].plot(subset.angle_deg, subset.cubic_residual, '--', color='C3', label='Cubic residual')
        provisional = subset[subset.status != 'CONVERGED']
        if len(provisional):
            axes[0].scatter(provisional.angle_deg, provisional.observed, facecolors='none', edgecolors='C1', label='Unconverged final iterate')
            axes[1].scatter(provisional.angle_deg, provisional.residual, facecolors='none', edgecolors='C1')
        for ax in axes:
            ax.axvspan(8, limit, color='gray', alpha=0.07)
            ax.axvline(8, ls=':', color='gray')
            ax.grid(alpha=0.25)
        axes[0].set(ylabel='$C_L$', title=f'Steady lift: trained 0–8°, evaluated through {limit}°')
        axes[0].legend(fontsize=8)
        axes[1].axhline(0, lw=0.7, color='black')
        axes[1].set(xlabel='Angle of attack (degrees)', ylabel='Observed − predicted $C_L$')
        axes[1].legend()
        fig.tight_layout()
        fig.savefig(OUT/filename, dpi=180)
        plt.close(fig)


if __name__ == '__main__':
    main()
