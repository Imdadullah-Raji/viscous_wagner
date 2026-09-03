# Campaign Plan 

As of September 3rd, we have less than 4 weeks left before our thesis predefence. Without a solid plan we would be in a danger of squandering our limited time, a sizeable portion of which will be occupied by simulation time itself. Here I layout a plan of our entire campaign of about 20 days: the simulations to run, variables to be measured and analyses to be run. 



## Types of Simulations 

All simulations will be run on the laminar solver in OpenFOAM. No
turbulence closure model is needed. Divergence scheme second-order upwind, time-stepping second-order implicit. 
### Unactuated

#### Static Airfoil

**Purpose**: 
- Establishing a baseline, validation against literature. 
- Ensuring that no meshing problems are present that may carry over to dynamic mesh(AMI) simulations.
- Capturing static stall behaviour.

Essentially two categories: pre-stall and post-stall. The prestall simulations are essentially steady state, so simulations are cheap. The post-stall has periodic vortex shedding, transient simulations are necessary.

| variables |max    |min        |steps                      |
|-----------|-------|-----------|---------------------------|
| Re        |100    |5000       | 100, 500, 1000, 2500, 5000|
| AoA       |0      |40         |`num_steps=40`             |


*total number of simulations*: 40*5 =200

The prestall AoAs do not need transient simulations. They need negligible time compared to the post-stall ones. The post-stall 
AoAs generally start after 20 degrees. Hence the compute costly 
simulations number around 100. 


#### Pitching Airfoil

**Purpose**
- Studying dynamic stall

5 angle of attacks and 4 reduced frequencies, and 3 mean+amp combinations. The experiment matrix is sized $5\times4\times3=60$
 

| variables |max    |min        |steps                      |
|-----------|-------|-----------|---------------------------|
| Re        |100    |5000       | 100, 500, 1000, 2500, 5000|
| k         |0.25   |2          | 0.25, 0.5, 1,2            |
|AoA        |-      |-          | $15\pm 10$, $15\pm 15$, $15\pm20$|



#### Step change in AoA 

**Purpose**
- Comparing it to Wagner function
- Comparing it to 2021 Dawson, Brunton paper, and fit their model to 
see if it holds. 

Two broad regimes. 
1. When the step does not cross into static stall region.
2. When it does

*Modeling Assumption*: The step can be modeled by a tanh function, 
keeping it smooth and differentiable. 

### Actuated
TBD

## Post-Processing 

Things to be stored. 
- Fields( $u$, $v$, $p$)
- Pressure Coefficient $C_p$
- Lift and drag coefficients $C_L, C_D$
-  Wall-shear stress $\tau_w$

### Storage Issue and Proper Orthogonal Decomposition

Field data needs to be compressed before they can be stored efficiently 
within our limited storage capacity. We will do proper orthogonal
decompositions of data for compression. 


