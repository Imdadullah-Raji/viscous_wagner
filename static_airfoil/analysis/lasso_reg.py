import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import Lasso
from sklearn.pipeline import make_pipeline
import pandas as pd

PATH = "/home/raji/Research/Thesis/static_airfoil/lift_coeffs.csv"

df = pd.read_csv(PATH, delimiter=" ")

alpha = df.iloc[:,0].to_numpy()*(np.pi/180)
cl = df.iloc[:, 1].to_numpy()
X = alpha.reshape(-1,1)

model = make_pipeline(
    PolynomialFeatures(degree=9, include_bias=False),
    StandardScaler(),
    Lasso(alpha=1e-4, max_iter=100000)
)

model.fit(X,cl)

cl_pred = model.predict(X)

poly = model.named_steps["polynomialfeatures"]
lasso = model.named_steps["lasso"]
scaler = model.named_steps["standardscaler"]

means = scaler.mean_
scales = scaler.scale_


names = poly.get_feature_names_out(["alpha"])

for name, coef in zip(names, (lasso.coef_)/scales):
    print(f"{name:10s}: {coef:+.8e}")

print("intercept:", lasso.intercept_-np.sum(lasso.coef_*means/scales))


plt.plot(alpha, cl_pred)
plt.plot(alpha, cl)
plt.show()