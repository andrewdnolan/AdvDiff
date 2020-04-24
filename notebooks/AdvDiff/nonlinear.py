#!/usr/bin/env python3

import sympy
import numpy as np
import scipy.linalg as LA
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

nxs    = np.arange(50,1550,300)
kappas = np.linspace(1e-2,1e0,10)
error_loose  = np.zeros((4,nxs.shape[0]))
error_strict  = np.zeros((4,nxs.shape[0]))
error_kappas  = np.zeros((4,kappas.shape[0],nxs.shape[0]))


for i, nx in enumerate(nxs):
    params['nx'] = nx
    coeffs['dt'] = 10*(params['L']/nx) * coeffs['κ']
    model = AdvDiff(params,coeffs)
    model.U[0,:] = ufunc(0,model.x,model.κ)

    Strang_Beam = model.Strang('BeamWarming','w')
    Goundov_Beam = model.Goundov('BeamWarming','w')
    Strang_Lax = model.Strang('LaxWendroff','w')
    Goundov_Lax = model.Goundov('LaxWendroff','w')

    exact  =  np.zeros((model.nt,model.nx))
    for j in range(0,model.nt):
        exact[j,:] = ufunc(model.dt*j,model.x,model.κ)

    error_loose[0,i] = LA.norm(Goundov_Beam- exact,np.inf)
    error_loose[1,i] = LA.norm(Goundov_Lax- exact,np.inf)

    error_loose[2,i] = LA.norm(Strang_Beam- exact,np.inf)
    error_loose[3,i] = LA.norm(Strang_Lax- exact,np.inf)

params = {'L':2*np.pi,'nx':50,'nt':5,'linear':False}
coeffs = {'κ':1e-1, 'a':10., 'σ':1.0}

for i, nx in enumerate(nxs):

    params['nx'] = nx
    model = AdvDiff(params,coeffs)
    model.U[0,:] = ufunc(0,model.x,model.κ)

    Strang_Beam_Str = model.Strang('BeamWarming','w')
    Goundov_Beam_Str = model.Goundov('BeamWarming','w')
    Strang_Lax_Str = model.Strang('LaxWendroff','w')
    Goundov_Lax_Str = model.Goundov('LaxWendroff','w')

    exact  =  np.zeros((model.nt,model.nx))
    for j in range(0,model.nt):
        exact[j,:] = ufunc(model.dt*j,model.x,model.κ)

    error_strict[0,i] = LA.norm(Goundov_Beam_Str- exact,np.inf)
    error_strict[1,i] = LA.norm(Goundov_Lax_Str- exact,np.inf)

    error_strict[2,i] = LA.norm(Strang_Beam_Str- exact,np.inf)
    error_strict[3,i] = LA.norm(Strang_Lax_Str- exact,np.inf)

params = {'L':2*np.pi,'nx':50,'nt':5,'linear':False}
coeffs = {'κ':1e-1, 'a':10., 'σ':1.0}

fig, ax = plt.subplots(1,2,sharey=True,figsize=(12,6))

ax[0].loglog(1/nxs, error_loose[0],'D:',color='olivedrab',label='Goundov: CN-BeamWarming', markersize=8)
ax[0].loglog(1/nxs, error_loose[1],'D:',color='olivedrab',label='Goundov: CN-LaxWendroff', markersize=8)
ax[0].loglog(1/nxs, error_loose[2],'D-',color='slategrey',label='Strang: CN-BeamWarming', markersize=8)
ax[0].loglog(1/nxs, error_loose[3],'D-',color='slategrey',label='Strang: CN-LaxWendroff', markersize=8)

ax[1].loglog(1/nxs, error_strict[0],'D:',color='olivedrab',label='Goundov: CN-BeamWarming', markersize=8)
ax[1].loglog(1/nxs, error_strict[1],'D:',color='olivedrab',label='Goundov: CN-LaxWendroff', markersize=8)
ax[1].loglog(1/nxs, error_strict[2],'D-',color='slategrey',label='Strang: CN-BeamWarming', markersize=8)
ax[1].loglog(1/nxs, error_strict[3],'D-',color='slategrey',label='Strang: CN-LaxWendroff', markersize=8)

ax[0].grid()
plt.legend()
ax[0].set_ylabel(r'$||U - u||_\infty$',fontsize='large')
ax[0].set_xlabel(r'$\Delta x$')
ax[0].set_title(r'$\Delta t = 10 \; \Delta x \kappa$')

ax[1].grid()
ax[1].set_xlabel(r'$\Delta x$')
ax[1].set_title(r'$\Delta t = \Delta x \kappa$')

plt.tight_layout()
plt.show()
