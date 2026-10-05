import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from pathlib import Path 
from scipy import stats

PATH = Path('/home/raji/Research/Thesis/static_airfoil/lift_coeffs.csv')
print(PATH)
df = pd.read_csv(PATH, delimiter=" ")

aoa = df.iloc[:,0].to_numpy()
mask = (aoa<=2)

aoa_deg = aoa[mask]
aoa_rad = aoa_deg*np.pi/180
cl = df.iloc[:,1].to_numpy()
cl = cl[mask]
# print(df.head())
# print(aoa_deg)
# print(aoa_rad)
print(cl)


result = stats.linregress((aoa_rad), (cl))

print(f"Slope: {result.slope:.4f}")
print(f"Slope/2$\pi$: {result.slope/(2*np.pi):.4f}")
print(f"Intercept:{result.intercept:.4f}")
print(f"R-squared: {result.rvalue**2:.4f}")


cl_pred = result.slope*aoa_rad+result.intercept
# cl_pred = np.exp(result.intercept)*(aoa_rad**result.slope)
print(cl_pred)

residuals = cl - cl_pred

fig, ax = plt.subplots(1,2,figsize=(10,4))

ax[0].scatter(aoa_deg, residuals)
ax[0].set_title('Residual plot')
ax[0].set_xlabel('Angle of Attack(deg)')
ax[0].set_ylabel('cl-cl_pred')

stats.probplot(residuals, dist="norm", plot=ax[1])
ax[1].set_title('Normal Q-Q Plot')

plt.tight_layout()
plt.show()


