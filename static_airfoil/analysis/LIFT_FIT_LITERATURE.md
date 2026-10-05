# Lift-curve fitting at ultra-low Reynolds numbers

Literature search and modelling notes — 5 October 2026

## Purpose and scope

This report examines polynomial and alternative representations of lift coefficient versus angle of attack, motivated by the present two-dimensional NACA0012 campaign at chord Reynolds number **Re = 500**.

The quantity of interest is the lift of the **steady Navier–Stokes solution**, including its unstable continuation beyond the reported vortex-shedding onset near 11.5°. This is distinct from the time-averaged lift of a shedding flow, a transient force after startup, or a quasi-steady coefficient inferred from a revolving wing.

The search found relevant polynomial-fitting precedents and nonlinear low-Re lift models. **It did not establish a published polynomial correlation specifically for the unstable steady NACA0012 branch at Re = 500.** This is a limitation of the present search, not proof that no such work exists.

## Relevant literature

### 1. Nonlinear viscous model for a steady flat plate

**Liu, T., Wang, S., and He, G. (2017). “Explicit Role of Viscosity in Generating Lift.” AIAA Journal, 55(11), 3990–3994.**

- DOI: [10.2514/1.J055907](https://doi.org/10.2514/1.J055907)
- Geometry and conditions: a stationary two-dimensional flat plate; comparison with direct numerical simulations at **Re = 200**.
- Model: a viscous-flow derivation gives a nonlinear lift relation of the form

  $$C_L = 2\pi\theta F(\theta)\cos\theta.$$

  Here, the correction factor depends on quantities including Reynolds number, a velocity ratio, and trailing-edge pressure loading. It is not simply a universal constant fitted to the lift slope.
- Relevance: this is the closest identified source to the present **steady-flow** modelling question. It provides a physical route to nonlinear lift behaviour at Reynolds numbers of a few hundred.
- Limitation: its flat-plate geometry and modelling assumptions differ from NACA0012; it does not supply a validated fit for the present unstable branch.

The formula and its Re = 200 comparison are also presented in Section 5, Eqs. (12)–(13), of [*Evolutionary understanding of airfoil lift*](https://link.springer.com/article/10.1186/s42774-021-00089-4). That exposition is useful for reading the correction factor, but the 2017 paper is the original reference above.

### 2. Approximately sinusoidal lift at Re = 192

**Wang, Z. J. (2000). “Vortex shedding and frequency selection in flapping flight.” Journal of Fluid Mechanics, 410, 323–341.**

- DOI: [10.1017/S0022112099008071](https://doi.org/10.1017/S0022112099008071)
- Open author copy: [PDF](https://dragonfly.tam.cornell.edu/publications/2000_JFM_Wang_FlpS.pdf)
- Geometry and conditions: translating two-dimensional elliptical sections; Section 3.2.2 and Fig. 6 compare lift versus angle at **Re = 192**, measured at a specified time after startup.
- Model observation: the lift curve is approximately proportional to

  $$C_L \sim A\sin(2\theta).$$

- Relevance: a direct low-Re precedent for a compact trigonometric representation of nonlinear lift.
- Limitation: the force depends on elapsed time and vortex development. These data are not the unstable steady Navier–Stokes branch. The paper's discussion of a lift maximum near 45° should not be transferred to the present airfoil.

### 3. Explicit sinusoidal fit for a revolving wing

**Santhanakrishnan, A., Jones, S. K., Dickson, W. B., Peek, M., Kasoju, V. T., Dickinson, M. H., and Miller, L. A. (2018). “Flow Structure and Force Generation on Flapping Wings at Low Reynolds Numbers Relevant to the Flight of Tiny Insects.” Fluids, 3(3), 45.**

- DOI and article: [10.3390/fluids3030045](https://doi.org/10.3390/fluids3030045)
- Geometry and conditions: three-dimensional revolving wings; an explicit force-coefficient fit is given at **Re = 130**.
- Model: Eq. (3) gives

  $$C_L = 0.18 + 1.36\sin(2.16\alpha-9.5).$$

  In this empirical expression, angles use the degree convention. For implementation with a radians-based sine function and angle of attack $\theta$ in radians, write the argument as $2.16\theta-9.5\pi/180$.
- Relevance: an explicit example of fitting a non-polynomial lift law at ultra-low Re.
- Limitation: revolving-wing forces and three-dimensional leading-edge-vortex behaviour differ from a stationary two-dimensional airfoil. The offset and phase shift are not appropriate symmetry constraints for NACA0012 without independent justification.

### 4. Piecewise polynomial airfoil fits

**Xie, C., Tao, G., and Wu, Z. (2021). “Performance Prediction and Design of Stratospheric Propeller.” Applied Sciences, 11(10), 4698.**

- DOI and article: [10.3390/app11104698](https://www.mdpi.com/2076-3417/11/10/4698)
- Geometry and conditions: FX63 airfoil database spanning Reynolds numbers beginning at **1,000**, with example fitted curves at **8,000** and **200,000**.
- Method: piecewise cubic and quartic polynomials represent lift and drag over specified angle intervals; see Eqs. (27)–(30).
- Relevance: a direct precedent for polynomial representations of airfoil coefficient data, especially **local or piecewise** representations rather than unrestricted global extrapolation.
- Limitation: the study uses RANS with a Spalart–Allmaras model for the airfoil database. It is a fitting-method precedent, not a validation of laminar Re = 500 aerodynamics or an unstable steady branch.

## Connection to the present fits

The following results are from local computations in this workspace, not from the cited papers. Angles in all polynomial equations below are in **radians**.

### Cubic fitted over 0–8°

Using the 17 original CSV observations through 8°, ordinary least squares gives

$$C_L = 0.000086595934 + 3.267395529\theta - 28.235887546\theta^3.$$

| Diagnostic set | Number of points | RMSE |
|---|---:|---:|
| Training observations, 0–8° | 17 | 0.000102739 |
| Additional cases within the fitted range | 14 | 0.0000826863 |
| Extrapolation cases, above 8° through 12° | 12 | 0.0117707 |

Training $R^2$ is 0.999999248. At 12°, the prediction is 0.425004 versus computed steady lift 0.452294, an underprediction of approximately 6.03%. The polynomial has a maximum at 11.253°.

The interpolation errors are small but systematically structured. Extrapolation errors grow before the reported shedding onset. The polynomial maximum is therefore not, by itself, evidence of a stability transition.

Artifacts:

- [Fit and residual plots](figs/odd_cubic_8deg_diagnostics.png)
- [Diagnostic data](figs/odd_cubic_8deg_diagnostics.csv)
- [Fit report](figs/odd_cubic_8deg_fit.txt)
- [Reproducible script](fit_odd_cubic_8deg.py)

### Cubic fitted over 0–12°

The earlier fit to all 23 CSV observations was

$$C_L = 0.002632829 + 3.184219881\theta - 24.279871653\theta^3.$$

Its training RMSE was 0.00218016. It peaks near 11.98° and extrapolates poorly: at 20° it predicts approximately 0.08145 compared with the computed steady value 0.48208. Extending the fitting interval changes the coefficients appreciably but does not resolve the cubic's inability to describe the wider branch.

See the [validation plot](figs/odd_cubic_validation.png) and [validation data](figs/odd_cubic_validation.csv).

The saved SIMPLE traces at 22° and 24° remain oscillatory in iteration space. Their final values were shown only as provisional comparisons, with the last 10% iteration ranges plotted. This numerical convergence issue does not establish nonexistence of the unstable steady solution.

## What a sinusoidal fit would imply

This section is a mathematical inference from the local fit, not a result claimed by the cited authors.

For a free-frequency sinusoid,

$$A\sin(B\theta) = AB\theta - \frac{AB^3}{6}\theta^3 + O(\theta^5).$$

Matching the linear and cubic coefficients of the 0–8° fit, while setting aside its small intercept, requires

$$B = \sqrt{-\frac{6b_3}{b_1}} \approx 7.20.$$

The commonly discussed $A\sin(2\theta)$ form instead fixes $B=2$, so it imposes a cubic-to-linear coefficient ratio of $-2/3$. The present fitted ratio is approximately $-8.64$.

Thus, fitting only the amplitude of $\sin(2\theta)$ cannot match both the local slope and curvature in these data. A free $B$ allows that local match but does not guarantee correct extrapolation or identify the shedding onset.

## Candidate models and validation plan

The literature supports exploring nonlinear representations; it does not select a universal law for this dataset. The following are proposed empirical candidates, not established Re = 500 correlations from the papers above.

| Candidate | Form | Main issue to test |
|---|---|---|
| Odd polynomial | $a_1\theta+a_3\theta^3+a_5\theta^5$ | Whether an extra term improves held-out accuracy without unstable extrapolation |
| Odd rational model | $(a_1\theta+a_3\theta^3)/(1+b_2\theta^2)$ | Denominator must remain nonzero over the validation interval; far-field behaviour is not automatically bounded |
| Saturating model | $A\tanh(B\theta)$ | Whether saturation represents the observed branch rather than merely improving visual agreement |
| Trigonometric model | $A\sin(B\theta)$ | Whether its imposed turning behaviour agrees with independent steady data |
| Piecewise polynomial or shape-preserving interpolant | Local representations between computed angles | Useful for interpolation; does not establish predictive extrapolation |

For a symmetric airfoil on a reflection-symmetric steady branch, an odd representation enforces $C_L(-\theta)=-C_L(\theta)$ and $C_L(0)=0$. The fits requested so far include an intercept; a constrained zero-intercept fit should be compared explicitly rather than substituted silently.

Recommended comparison procedure:

1. Preserve the existing 0–8° training set for a controlled comparison between model families.
2. Evaluate interpolation on omitted steady cases inside that interval.
3. Evaluate extrapolation separately on 8–12° and on the converged saved cases through 20°.
4. Report RMSE, maximum absolute error, residual trends, and predicted turning points. Avoid relative errors near zero lift.
5. Assess sensitivity to the upper training angle. A model whose extrapolation changes strongly with that choice should be described as a local surrogate.
6. Keep steady-branch data separate from shedding means. Treat 22° and 24° as provisional validation points until converged steady solutions are available.

**Working conclusion:** the present odd cubic is an accurate local approximation through 8°, but its extrapolation is inadequate. Liu et al. provide the most relevant physical modelling reference; the trigonometric and piecewise-polynomial papers provide alternative functional-form precedents with substantially different flow conditions. Any wider-range correlation should be selected using independent steady-branch validation rather than training $R^2$ alone.
