# EASC 711: Numerical Methods Final Project   


Numerical experiments with the advection diffusion equation using finite
difference methods. This is the final project repository for EASC 711: Numerical
Methods. Repository structure is loosely based on the [`CFD Python`](https://github.com/barbagroup/CFDPython)  repository.  


## Notebooks
>  To run these notebooks interactively launch a binder session with: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andrewdnolan/AdvDiff/master)

- [1-Diffusion](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Diffusion_1D.ipynb): Crank-Nicolson solution to the diffusion equation in one dimension.
- [2-Advection](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Advection_1D.ipynb): UpWind, Lax-Wendroff, and Beam-Splitting methods advection equation in one dimension.
- [3-AdvDiff](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/AdvDiff_1D.ipynb): First and Second order operator splitting for solving the Advection-Diffusion equation in one dimension.


## To Do:   
- Binder:
  - [x] set up binder link
  - [ ] fix `requirments.txt` file
  - [ ] export environment  

- Diffusion Equation:  
   - [x] python implementation
   - [x] Boundary conditions are wonky
   - [x] Periodic Boundary Conditions ?
   - [ ] stability analysis
   - [ ] references section in notebook

- Avection Equation:
  - [x] port info from [MATH 709 Final Proj](https://github.com/andrewdnolan/MATH-709-Final-Project)
  - [x] state mathematical problem and numerical methods to solve  
  - [x] Beam-Splitting Method
  - [x] Analytical Solution
  - [x] stability analysis, CFL condition
  - [ ] Reference section

- Fractional Step Methods:
  - [ ] scale parameters to glaciological scales
  - [ ] Peclet Number based experiments
  - [ ] Stability and convergence experiments
