# EASC 711: Numerical Methods Final Project   


Numerical experiments with the advection diffusion equation using finite
difference methods. This is the final project repository for EASC 711: Numerical
Methods. Repository structure is loosely based on the [`CFD Python`](https://github.com/barbagroup/CFDPython)  repository.


## Notebooks
>  To run these notebooks interactively launch a binder session with: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andrewdnolan/AdvDiff/master)

- [1-Diffusion](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Diffusion_1D.ipynb): Crank-Nicolson solution to the diffusion equation in one dimension.
- [2-Advection](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Advection_1D.ipynb): Lax-Wendroff, Beam-Splitting, Slope Limiting methods advection equation in one dimension.


## To Do:   
- Binder:
  - [x] export environment
  - [x] set up binder link
  - [ ] fix `requirments.txt` file

- Diffusion Equation:  
   - [ ] stability analysis
   - [x] python implementation
   - [x] Boundary conditions are wonky
   - [x] Periodic Boundary Conditions ?
   - [ ] references section in notebook
   - [ ] experiments to determine most efficient solver

- Avection Equation:
  - [x] port info from [MATH 709 Final Proj](https://github.com/andrewdnolan/MATH-709-Final-Project)
  - [x] state mathematical problem and numerical methods to solve  
  - [x] Beam-Splitting Method
  - [ ] Analytical Solution
  - [ ] TVD method??
  - [ ] stability analysis, CFL condition
  - [ ] Reference section

- Fractional Step Methods:
  - [ ] scale parameters to glaciological scales
  - [ ] Peclet Number based experiments
  - [ ] Stability and convergence experiments
