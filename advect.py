#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.model import advection

plt.rcParams['text.usetex'] = True

params     = {'L':2*np.pi,'nx':100,'nt':100}
adv_params = {'a':3, 'σ':0.7, 'dt':0.0025}

nxs = np.arange(10,210,10)
err = np.zeros((4,nxs.shape[0]))

for i, nx in enumerate(nxs):
    params['nx'] = nx
    model = advection(params,adv_params)
    model.U[0,:] = 5*np.sin((model.x-np.pi/2))+5
    model.run('UpWind')

    analytical = np.zeros((params['nt'],params['nx']))
    for t in range(0,params['nt']):
        analytical[t,:] = 5*np.sin(((model.x-(model.a*model.dt*t))-np.pi/2))+5

    err[1,i] = LA.norm(model.U - analytical,2)

    model = advection(params,adv_params)
    model.U[0,:] = 5*np.sin((model.x-np.pi/2))+5
    model.run('FTCS')
    err[0,i] = LA.norm(model.U - analytical,2)

    model = advection(params,adv_params)
    model.U[0,:] = 5*np.sin((model.x-np.pi/2))+5
    model.run('LaxWendroff')
    err[2,i] = LA.norm(model.U - analytical,2)

    model = advection(params,adv_params)
    model.U[0,:] = 5*np.sin((model.x-np.pi/2))+5
    model.run('BeamWarming')
    err[3,i] = LA.norm(model.U - analytical,2)

fig, ax = plt.subplots(1,1,figsize=(6,4))
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(nxs,err[0],'kx-',label='FTCS',   markersize=4)
ax.plot(nxs,err[1],'ko-',label='UpWind', markersize=3)
ax.plot(nxs,err[2],'k-.',label='LaxWendroff')
ax.plot(nxs,err[3],'k:',label='BeamWarming')

#ax.axhline(1e-16,nxs[0],nxs[-1])
ax.set_ylabel(r'$|U_{num} - U_{anl}|_2$')
ax.set_xlabel('N')
ax.set_title(r'nt=100, dt=0.001')

plt.legend()
plt.tight_layout()
plt.show()
