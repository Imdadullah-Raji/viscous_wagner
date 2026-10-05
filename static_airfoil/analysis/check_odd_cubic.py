"""Validate the fixed 0–12 degree OLS model against additional steady cases."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from steady_summary import load, verdict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis' / 'figs'
data = pd.read_csv(ROOT / 'lift_coeffs.csv', sep=r'\s+')
angles = data.iloc[:, 0].to_numpy()
y = data.iloc[:, 1].to_numpy()

def design(degrees):
    t = np.deg2rad(np.asarray(degrees))
    return np.column_stack((t, t**3))

model = LinearRegression().fit(design(angles), y)
rows = []
for case in (ROOT / 'steady_re500').glob('steady_aoa*'):
    if not case.is_dir():
        continue
    angle = float(case.name.removeprefix('steady_aoa'))
    trace = load(str(case))
    if trace is None:
        continue
    it, cd, cl = trace
    status, ranges = verdict(it, cl, 1e-4)
    tail = cl[it >= 0.9 * it[-1] + 0.1 * it[0]]
    group = ('training' if np.any(np.isclose(angle, angles)) else
             'interpolation' if angles.min() < angle < angles.max() else 'extrapolation')
    prediction = model.predict(design([angle]))[0]
    rows.append(dict(angle_deg=angle, observed=cl[-1], predicted=prediction,
                     residual=cl[-1]-prediction, group=group, status=status,
                     tail_min=tail.min(), tail_max=tail.max()))
df = pd.DataFrame(rows).sort_values('angle_deg')
df.to_csv(OUT / 'odd_cubic_validation.csv', index=False)
report = []
for group in ['training', 'interpolation', 'extrapolation']:
    subset = df[(df.group == group) & (df.status == 'CONVERGED')]
    report.append(f'{group}: n={len(subset)}, RMSE={np.sqrt(np.mean(subset.residual**2)):.8g}, '
                  f'max abs error={subset.residual.abs().max():.8g}')
peak = np.rad2deg(np.sqrt(-model.coef_[0] / (3 * model.coef_[1])))
report.append(f'Polynomial maximum at {peak:.8g} degrees')
report.append(df[df.angle_deg >= 11.5].to_string(index=False))
text = '\n'.join(report)
print(text)
(OUT / 'odd_cubic_validation.txt').write_text(text + '\n')

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex='col')
for column, limit in enumerate([12, 24]):
    ax, residual_ax = axes[:, column]
    grid = np.linspace(0, limit, 500)
    ax.plot(grid, model.predict(design(grid)), color='C0', label='Fixed OLS fit (trained 0–12°)')
    for group, marker, color in [('training', 'o', 'black'),
                                  ('interpolation', 'x', 'C2'),
                                  ('extrapolation', 's', 'C1')]:
        subset = df[(df.group == group) & (df.status == 'CONVERGED') & (df.angle_deg <= limit)]
        ax.scatter(subset.angle_deg, subset.observed, marker=marker, color=color, s=24, label=group)
        residual_ax.scatter(subset.angle_deg, subset.residual, marker=marker, color=color, s=24)
    provisional = df[(df.status != 'CONVERGED') & (df.angle_deg <= limit)]
    if len(provisional):
        midpoint = (provisional.tail_min + provisional.tail_max)/2
        half = (provisional.tail_max-provisional.tail_min)/2
        ax.errorbar(provisional.angle_deg, midpoint, yerr=half, fmt='D', mfc='none',
                    color='C3', capsize=4, label='Unconverged: last 10% iteration range')
        residual_ax.errorbar(provisional.angle_deg, midpoint-provisional.predicted,
                             yerr=half, fmt='D', mfc='none', color='C3', capsize=4)
    for panel in [ax, residual_ax]:
        panel.axvline(11.5, color='gray', ls=':', label='Shedding onset ≈11.5°' if panel is ax else None)
        if limit > 12:
            panel.axvspan(12, limit, color='C1', alpha=0.06)
        panel.grid(alpha=0.2)
    ax.set_title('Interpolation: independent omitted cases' if limit == 12 else 'Extrapolation along steady branch')
    ax.set_ylabel('$C_L$')
    ax.legend(fontsize=7)
    residual_ax.axhline(0, color='black', lw=0.7)
    residual_ax.set(xlabel='Angle of attack (degrees)', ylabel='Observed − predicted $C_L$')
fig.tight_layout()
fig.savefig(OUT / 'odd_cubic_validation.png', dpi=200)
plt.close(fig)
