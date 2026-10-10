This simulation workspace has got very large. It is a heterogeneous mixture of 
- few originial openfoam case files, 
- case generation bash scripts 
- generated OpenFOAM cases from those bash scripts 
- python mesher files using gmsh API
- python analysis files
- log files 
- markdown files: instructions, notes, results summary
- Simulation results, large fields and small force history files. 

We need to put the scripts under version control. The scripts include. 

- mesher files
- python analysis files
- case generation bash scripts 
- sample OpenFOAM case for steady and unsteady run. 
- Result markdowns 


This folder is too messy to put under version contrlol. I want a clean repository at the grandparent directory `/home/raji/Research/Repositories/`. 


## Goal of this project

In that `Repositories` a new folder should be created. The research program that I am currently conducting is the step response of a airfoil at ultra low Re: which I call **The Viscous Wagner Problem**. The research is inspired by 2021 paper by SL Bruton and STM Dawson: *Improved Approximation to Wagner Function by Sparse Identification of Nonlinear Dynamics*. Currently I am setting the steady baseline, so far I have found that, from the very onset the lift curve is described by a odd cubic( higher terms are small near the origin) from the very onset. This is the first significant difference from the potential flow linear slope of $2\pi\alpha$.
>Most of the git commit messages should reflect the research intent rather than software implementation or bug fix messages. 

## Types of Simulations

### Pre-October Simualtions
This folder contains several types of simulation. The initial simulations are older than the intent of attacking the viscous wagner problem. Initially there was no intent at all beyond seeing the NACA0012 shed. Then a claude agent reported period doubling, searching the literature I hit Kurtulus(2016) and related papers. - Rossi, Durante et al
- Gupta et al, *Two- and three-dimensional wake transitions of a NACA0012 airfoil*, DOI: 10.1017/jfm.2022.958 
These simulations save fields over 20 complete shedding cycles, with period doubled flows being twice as long. 

### September Last Week

Large angle pitching simulations were done. These can be found at the parent directory. These saves fields too, the goal was to validate against literature, some paper. Details to be found at `/home/raji/Research/Thesis/reports/predefense`(a git repo)


### October Simulations
Then the simulations targeting the first portion of viscous wagner, avoiding large angle steps: setting the steady baseline.

- First, LASSO regression fits were done on Re=1000 cases, preshedding. This was prbably a misstep, but I intend to document it. 
-I had fit all powers. Then I decided to fit odd powers from the symmetry of the problem.  
- A problem arose from powers of $\alpha$ being higly correleated in the interval of fit, resulting in very high condition numbers[ CITATION: *Check the docs*] Decision was taken that Legendre polynomials will be fit. It gave stable results. Its condition number was also quite small, around 6-7.  
-  A power law type scaling was observed in the cubic coefficient. The fits were done accordingly for NACA0012

- The entire process had been repeated for NACA0006 and NACA0012. Some plate cases might be in `/home/raji/Research/Thesis/plate_case` fetch them too. 


## Machines and CFD Code Version Issue

Simulations have been conducted on 4 different machines, all foundation OpenFOAMs
1. This one, 6 core Intel core i5 12th Gen 12400, OpenFOAM-dev
2. A 8 core AMD ryzen, OpenFOAMv14.
3. M2 MacBook Air running ubuntu on a 6 core multipass VM, OpenFOAM13
4. A google cloud c3-standard-8, 4 cores, running OpenFOAMv13 


Bottom line: Make a folder at `/home/raji/Research/repositories`, copy codes there, not data. The commit messages should reflect sequential research decisions, mostly. To determine the sequence look at the last modifications dates.

**In the commit messages, also flag explicitly if a code is claude generated. Like, most of analysis scripts**


## Hypothesis Source: where each idea came from

- LASSO regression: myself,  LASSO path inspection: ChatGPT ; Q-Q plot, residual plots, ChatGPT again, I did not knew the concepts.

- orthogonal fit idea and saturating, suggest by Claude Opus 5.5 in a conversation. The Legendre choice and the specific saturating fit in `bizarre_experiments` is mine. 

- Cubic is sufficient to describe nonlinearity at origin: myself

- The Re scaling law of a3 for NACA0012, suggested by claude code mid-session.(claude opus 5.5) **It calculated the fit by itself, without prompting**

- **UNTESTED**: I added one entry you didn't list: the explanation that the plate's a3 peak comes from its sharp leading edge: Claude Opus 5.5 in a claude code session. 
