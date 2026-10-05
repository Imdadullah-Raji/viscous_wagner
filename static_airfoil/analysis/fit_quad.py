import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import pandas as pd

PATH = '/home/raji/Research/Thesis/static_airfoil/lift_coeffs.csv'

df = pd.read_csv(PATH, delimiter=" ")

alpha = df.iloc[:,0].to_numpy()
cl = df.iloc[:,1].to_numpy()
mask = alpha<=10

alpha_full = alpha*(np.pi/180)
alpha = alpha[mask]*(np.pi/180)
cl_full = cl
cl = cl[mask]


X = alpha.reshape(-1,1)
y = cl
poly = PolynomialFeatures(degree=3, include_bias=False)
X_poly = poly.fit_transform(X)

model = LinearRegression()
model.fit(X_poly,y)

X_full_poly = poly.fit_transform(alpha_full.reshape(-1,1))

cl_pred = model.predict(X_poly)
cl_pred_full = model.predict(X_full_poly)

print("intercept:", model.intercept_)
print("coefficients:", model.coef_)



residual = cl - cl_pred

plt.plot(alpha, residual)
plt.show()
plt.plot(alpha,cl )
plt.plot(alpha, cl_pred)
plt.show()

score = model.score(X_poly, y)
print(score)

print("Cl:" ,cl)
print("Cl_pred: ",cl_pred)