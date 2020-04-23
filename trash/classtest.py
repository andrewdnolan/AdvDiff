#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.model import Diffusion, AdvDiff
#from advdiff.heuristic import AdvDiff

params = {'L':4.*np.pi,'nx':100,'nt':100}
coeffs = {'κ':0.3, 'σ':0.75, 'a':3.}

def initcond(test, M = 1,P = 1,x_0 = np.pi,t=0.1):
    return (M / (P *(4*np.pi*test.κ*t)**0.5)) * np.exp(-(((test.x-(x_0 + test.a*t))**2)/(4*test.κ*t)))

nxs = np.arange(100,2500,400)
err = np.zeros((4,nxs.shape[0]))

for i, nx in enumerate(nxs):
    params['nx'] = nx
    test = AdvDiff(params,coeffs)
    test.U[0,:] = initcond(test)
    FTCS  = test.Goundov('FTCS','w')
    UpWind  = test.Goundov('UpWind','w')
    LaxWendroff = test.Goundov('LaxWendroff','w')
    BeamWarming  = test.Goundov('BeamWarming','w')

    analytical = np.zeros((params['nt'],params['nx']))
    for t in range(0,params['nt']):
        analytical[t,:] = initcond(test,t=(t*test.dt+0.1))
        
    err[0,i] = LA.norm(analytical - FTCS, np.inf)
    err[1,i] = LA.norm(analytical - UpWind, np.inf)
    err[2,i] = LA.norm(analytical - LaxWendroff, np.inf)
    err[3,i] = LA.norm(analytical - BeamWarming, np.inf)

# for t in range(params['nt']):
#     plt.plot(test.x, second[t,:])
# plt.show()

fig, ax = plt.subplots(1,1,figsize=(6,4))
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.loglog(1/nxs,err[0],'kx-',label='CN & FTCS',markersize=4,basey=10)
ax.loglog(1/nxs,err[1],'ko-',label='CN & UpWind',markersize=4,basey=10)
ax.loglog(1/nxs,err[2],'k-.', label='CN & LaxWendroff',markersize=4,basey=10)
ax.loglog(1/nxs,err[3],'k:',label='CN & BeamWarming',markersize=4,basey=10)

ax.grid()
#
#ax.axhline(1e-16,nxs[0],nxs[-1])
ax.set_ylabel(r'$|U_{num} - U_{anl}|_\infty$')
ax.set_xlabel('$\Delta x$')
ax.set_title(r'Goundov Splitting (nt=100)')
plt.legend()

# fig, ax = plt.subplots(1,1,figsize=(6,4))
# ax.plot(test.x,UpWind[-1,:],'x')
# ax.plot(test.x,LaxWendroff[-1,:],'x')
# ax.plot(test.x,BeamWarming[-1,:],'x')
# ax.plot(test.x,analytical[-1,:])
plt.tight_layout()
plt.show()
