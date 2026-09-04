# Static Airfoil Simulation Plan

## Modeling Assumptions 

The simulation is 2D, hence any 3D spanwise instability, vortex stretching cannot appear. The highest Re =1000 and  highest AoA is 40 degrees. 


## Sampling 

- *limit cycle shedding* :32 snapshots per shedding cycle, there are more frequencies to capture other than the shedding frequencies. 

- *Decision Concering the Transient*: The transient depends on Initial conditions. The transient decays fast in the beginning and slowly converges towards the limit cycle. The transient can be useful for operator 
learning frameworks, but that is a study parameterized 
by initial conditions. My current plan is not a study 
parameterized by initial conditions, rather a study parameterized by actuation. **So I decide to discard the transient**

## Storage

- Original simulation data in NetCDF format 
- TBD: Interpolation onto cartesian grids and POD of that data vs cell-volume-weighted POD on original mesh.
- After all the simulations are done a global POD basis
will be constructed. 

## Post-processing 

- Fields u, v, p, . Needs to cropped around the domain of interest later; sampled at least 32 times per shedding cycle 
- $C_p$, sampling rate is each simulation step
- Integrated loads $C_L$ and $C_D$, *calculated at every simulation step, including the transient.*
- wall shear stress $\tau_w$ for tracking seperation, sampling rate: each simulation step


## Experiment Matrix

- Re = 100, 200, 300, 400, 500
- AoA = 20:1:40 





