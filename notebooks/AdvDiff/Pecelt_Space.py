#!/usr/bin/env python3

import sympy
import numpy as np
import scipy.linalg as LA
from cycler import cycler
import matplotlib.pyplot as plt
from advdiff.model import AdvDiff
from sympy.utilities.lambdify import lambdify

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14

x, kappa, t = sympy.symbols('x kappa t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * kappa * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * kappa * (t + 1))))

phiprime = phi.diff(x)

u = -2 * kappa * (phiprime / phi) + 4
ufunc = lambdify((t, x, kappa), u)

params = {'L':2*np.pi,'nx':50,'nt':5,'linear':False}
coeffs = {'κ':1e-1, 'a':10., 'σ':1.0}

nxs    = np.arange(50,1550,300)
kappas = np.linspace(1e-2,1e0,10)
error_kappas  = np.zeros((kappas.shape[0],nxs.shape[0]))


params = {'L':2*np.pi,'nx':50,'nt':5,'linear':False}
coeffs = {'κ':1e-1, 'a':10., 'σ':1.0}

for i, nx in enumerate(nxs):
    for j, kappa in enumerate(kappas):

        params['nx'] = nx
        coeffs['κ'] = kappa
        model = AdvDiff(params,coeffs)
        model.U[0,:] = ufunc(0,model.x,model.κ)

        Strang_Beam_kap = model.Strang('BeamWarming','w')
        # Goundov_Beam_kap = model.Goundov('BeamWarming','w')
        # Strang_Lax_kap = model.Strang('LaxWendroff','w')
        # Goundov_Lax_kap = model.Goundov('LaxWendroff','w')

        exact  =  np.zeros((model.nt,model.nx))
        for k in range(0,model.nt):
            exact[k,:] = ufunc(model.dt*k,model.x,model.κ)

        # error_kappas[0,j,i] = LA.norm(Goundov_Beam_kap- exact,np.inf)
        # error_kappas[1,j,i] = LA.norm(Goundov_Lax_kap- exact,np.inf)

        error_kappas[j,i] = LA.norm(Strang_Beam_kap- exact,np.inf)
        # error_kappas[3,j,i] = LA.norm(Strang_Lax_kap- exact,np.inf)

new_colors = [plt.get_cmap('viridis')(1. * i/len(kappas)) for i in range(len(kappas))]
prop_cycle=cycler('color', new_colors)

fig, ax = plt.subplots(1,1,sharey=True,figsize=(8,6))
ax.set_prop_cycle(prop_cycle)
for i, kappa in enumerate(kappas):
    ax.loglog(1/nxs,error_kappas[i,:],label="PE={}".format(kappa))

ax.grid()
plt.legend()
ax.set_ylabel(r'$||U - u||_\infty$',fontsize='large')
ax.set_xlabel(r'$\Delta x$')

plt.tight_layout()
plt.show()
