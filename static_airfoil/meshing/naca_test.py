import matplotlib.pyplot as plt
from airfoil import naca4
af = naca4("2412", n_points=200, x_trunc=0.98)
plt.plot(*af.upper.T, '.-', ms=2); plt.plot(*af.lower.T, '.-', ms=2)
plt.plot([af.te_upper[0], af.te_lower[0]], [af.te_upper[1], af.te_lower[1]], 'k-')
plt.axis('equal'); plt.show()