#!/usr/bin/env python3

import sys
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt

sys.path.append('../')
from advdiff.model import advection

plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['text.usetex'] = True

def η(x,t,a=3.):
    return 5*np.sin(2*((x-(a*t))-np.pi))+5

params     = {'L':2*np.pi,'nx':50,'nt':100}
adv_params_nx0 = {'a':3, 'σ':0.5}

nxs = np.arange(10,4100,100)
err_nx0 = np.zeros((4,nxs.shape[0]))

for i, nx in enumerate(nxs):
    params['nx'] = nx
    model = advection(params,adv_params_nx0)

    model.U[0,:] = η(model.x,0)

    UpWind = model.run('UpWind','w')
    FTCS   = model.run('FTCS','w')
    LaxWendroff = model.run('LaxWendroff','w')
    BeamWarming = model.run('BeamWarming','w')

    exact = np.zeros_like(model.U)
    for j in range(0,model.nt):
        exact[j,:] = η(model.x,(j*model.dt))

    err_nx0[0,i] = LA.norm(FTCS - exact,np.inf)
    err_nx0[1,i] = LA.norm(UpWind - exact,np.inf)
    err_nx0[2,i] = LA.norm(LaxWendroff - exact,np.inf)
    err_nx0[3,i] = LA.norm(BeamWarming - exact,np.inf)

params     = {'L':2*np.pi,'nx':50,'nt':100}
adv_params_nx1 = {'a':3, 'σ':0.01}

nxs = np.arange(10,4100,100)
err_nx1 = np.zeros((4,nxs.shape[0]))

for i, nx in enumerate(nxs):
    params['nx'] = nx
    model = advection(params,adv_params_nx1)

    model.U[0,:] = η(model.x,0)

    UpWind = model.run('UpWind','w')
    FTCS   = model.run('FTCS','w')
    LaxWendroff = model.run('LaxWendroff','w')
    BeamWarming = model.run('BeamWarming','w')

    exact = np.zeros_like(model.U)
    for j in range(0,model.nt):
        exact[j,:] = η(model.x,(j*model.dt))

    err_nx1[0,i] = LA.norm(FTCS - exact,np.inf)
    err_nx1[1,i] = LA.norm(UpWind - exact,np.inf)
    err_nx1[2,i] = LA.norm(LaxWendroff - exact,np.inf)
    err_nx1[3,i] = LA.norm(BeamWarming - exact,np.inf)

params     = {'L':2*np.pi,'nx':100,'nt':100}
adv_params_sig = {'a':3, 'σ':0.01}

sigmas = np.arange(0.0001,1,0.01)
err_sig = np.zeros((4,sigmas.shape[0]))

for i, sigma in enumerate(sigmas):
    adv_params_sig['σ'] = sigma
    model = advection(params,adv_params_sig)
    model.U[0,:] = η(model.x,0)

    UpWind = model.run('UpWind','w')
    FTCS   = model.run('FTCS','w')
    LaxWendroff = model.run('LaxWendroff','w')
    BeamWarming = model.run('BeamWarming','w')

    exact = np.zeros_like(model.U)
    for j in range(0,model.nt):
        exact[j,:] = η(model.x,(j*model.dt))

    err_sig[0,i] = LA.norm(FTCS - exact,np.inf)
    err_sig[1,i] = LA.norm(UpWind - exact,np.inf)
    err_sig[2,i] = LA.norm(LaxWendroff - exact,np.inf)
    err_sig[3,i] = LA.norm(BeamWarming - exact,np.inf)

fig, ax = plt.subplots(1,3,sharey=True,figsize=(18,6))

ax[0].loglog(1/nxs, err_nx0[0],'x-',label='FTCS',   markersize=8)
ax[0].loglog(1/nxs, err_nx0[1],'o-',label='UpWind', markersize=8)
ax[0].loglog(1/nxs, err_nx0[2],'d-.',label='LaxWendroff',   markersize=8)
ax[0].loglog(1/nxs, err_nx0[3],'^:', label='BeamWarming',   markersize=8)
ax[0].grid()
ax[0].set_ylabel(r'$||U - u||_\infty$',fontsize='large')
ax[0].set_xlabel(r'$\Delta x$')
ax[0].set_title(r'nt=100, $\sigma$={}'.format(adv_params_nx0['σ']))

ax[1].loglog(1/nxs, err_nx1[0],'x-',label='FTCS',   markersize=8)
ax[1].loglog(1/nxs, err_nx1[1],'o-',label='UpWind', markersize=8)
ax[1].loglog(1/nxs, err_nx1[2],'d-.',label='LaxWendroff',   markersize=8)
ax[1].loglog(1/nxs, err_nx1[3],'^:', label='BeamWarming',   markersize=8)
ax[1].grid()
ax[1].set_xlabel(r'$\Delta x$')
ax[1].set_title(r'nt=100, $\sigma$={}'.format(adv_params_nx1['σ']))

ax[2].loglog(sigmas, err_sig[0],'x-',label='FTCS',   markersize=8)
ax[2].loglog(sigmas, err_sig[1],'o-',label='UpWind', markersize=8)
ax[2].loglog(sigmas, err_sig[2],'d-.',label='LaxWendroff',   markersize=8)
ax[2].loglog(sigmas, err_sig[3],'^:' ,label='BeamWarming',   markersize=8)
ax[2].grid()
ax[2].set_ylim(None,10**4)
ax[2].set_xlabel(r'$\sigma$')
ax[2].set_title(r'nt={}, nx={}'.format(params['nt'],params['nx']))

plt.legend()
plt.tight_layout()
