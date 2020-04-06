# EASC 711: Numerical Methods Final Project   

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andrewdnolan/AdvDiff/master)

Numerical experiments with the advection diffusion equation using finite
difference methods.

## Notebooks
- [1-Diffusion](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Diffusion_1D.ipynb): Crank-Nicolson solution to the diffusion equation in one dimension.
- [2-Advection](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Advection_1D.ipynb)): Lax-Wendroff, Beam-Splitting, Slope Limiting methods advection equation in one dimension.


## To Do:   
- Binder:
  - [x] export environment
  - [ ] set up binder link
  - [ ] make `setup.py` file for `advdiff` package
- Iterative Solvers:
    - [ ] add GMRES
    - [ ] add conjugate gradient
    - [ ] sparse matrix implementation?

- Diffusion Equation:  
   - [ ] stability analysis
   - [x] python implementation
   - [ ] Boundary conditions are wonky
   - [ ] Periodic Boundary Conditions ?
   - [ ] references section in notebook
   - [ ] experiments to determine most efficient solver

- Avection Equation:
  - [ ] port info from [MATH 709 Final Proj](https://github.com/andrewdnolan/MATH-709-Final-Project)
  - [ ] state mathematical problem and numerical methods to solve  
  - [ ] stability analysis, CFL condition


- Fractional Step Methods:
  - [ ] scale parameters to glaciological scales
  - [ ] Peclet Number based experiments
  - [ ] Stability and convergence experiments
