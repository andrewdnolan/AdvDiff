# EASC 711: Numerical Methods Final Project   


Numerical experiments with the advection diffusion equation using finite
difference methods. This is the final project repository for EASC 711: Numerical
Methods. Repository structure is loosely based on the [`CFD Python`](https://github.com/barbagroup/CFDPython)  repository.  


## Notebooks
>  To run these notebooks interactively launch a binder session with: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/andrewdnolan/AdvDiff/master)

- [1-Diffusion](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Diffusion_1D.ipynb): Crank-Nicolson solution to the diffusion equation in one dimension.
- [2-Advection](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/Advection_1D.ipynb): UpWind, Lax-Wendroff, and Beam-Splitting methods advection equation in one dimension.
- [3-AdvDiff](https://nbviewer.jupyter.org/github/andrewdnolan/AdvDiff/blob/master/notebooks/AdvDiff_1D.ipynb): First and Second order operator splitting for solving the Advection-Diffusion equation in one dimension.

## References  

- Barba et al., (2018). _CFD Python: the 12 steps to Navier-Stokes equations_. Journal of Open Source Education, 1(9), 21, https://doi.org/10.21105/jose.00021

- Langtangen and  Linge. Finite Difference Computing with PDEs - A Modern Software Approach, Texts in Computational Science and Engineering, Springer, 2016, https://doi.org/10.1007/978-3-319-55456-3  

- LeVeque, Randall J. 2002. \textit{Finite Volume Methods for Hyperbolic Problems}. Cambridge Texts in Applied Mathematics. Cambridge: Cambridge University Press. doi:10.1017/CBO9780511791253.

- LeVeque, Randall J. 2007. _Finite Difference Methods for Ordinary and Partial Differential Equations: Steady-State and Time-dependent Problems_. Classics in Applied Mathematics. Society of Industrial and Applied Mathematics (SIAM). doi:10.1137/1.9780898717839

## To Do:  

- Fractional Step Methods:
  - [x] Peclet Number based experiments
  - [x] Stability and convergence experiments
