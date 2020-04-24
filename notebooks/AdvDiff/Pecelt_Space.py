#!/usr/bin/env python3

import sympy
import numpy as np
import scipy.linalg as LA
from cycler import cycler
import matplotlib.pyplot as plt
from advdiff.model import AdvDiff
from sympy.utilities.lambdify import lambdify

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.titlesize'] = 16

x, kappa, t = sympy.symbols('x kappa t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * kappa * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * kappa * (t + 1))))

phiprime = phi.diff(x)

u = -2 * kappa * (phiprime / phi) + 4
ufunc = lambdify((t, x, kappa), u)

params = {'L':2*np.pi,'nx':50,'nt':5,'linear':False}
coeffs = {'κ':1e-1, 'a':10., 'σ':1.0}

nxs    = np.arange(50,1550,100)
kappas = np.linspace(1e-2,1e0,50)
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

        exact  =  np.zeros((model.nt,model.nx))
        for k in range(0,model.nt):
            exact[k,:] = ufunc(model.dt*k,model.x,model.κ)


        error_kappas[j,i] = LA.norm(Strang_Beam_kap- exact,np.inf)


xx, yy = np.meshgrid(1/nxs,kappas)
Pe = 4*xx/yy

from matplotlib import ticker, cm

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(17, 8))

pes = ax[0].contourf(1/nxs,kappas,Pe,levels = np.logspace(-3,1,9), locator=ticker.LogLocator(), cmap=cm.PuBu_r)
cbar0 = fig.colorbar(pes, ax = ax[0])
cbar0.ax.set_title('$Pe$',fontsize='large', pad = 20)
ax[0].set_ylabel(r'$\kappa$',fontsize='xx-large')
ax[0].set_xlabel(r'$\Delta x$',fontsize='xx-large')
ax[0].set_title('Peclet Number',fontsize='large')

cs = ax[1].contourf(1/nxs,kappas,error_kappas, levels = np.logspace(-1,2,19), locator=ticker.LogLocator(), cmap=cm.viridis)
cbar1 = fig.colorbar(cs, ax = ax[1])
cbar1.ax.set_title('$||U - u||_\infty$',fontsize='large', pad = 20)
ax[1].set_xlabel(r'$\Delta x$',fontsize='xx-large')
ax[1].set_title('Error',fontsize='large')

plt.tight_layout()
