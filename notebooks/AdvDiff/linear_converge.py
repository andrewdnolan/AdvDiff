#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.model import AdvDiff

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14

params = {'L':100.,'nx':1000,'nt':25}
coeffs = {'κ':1e-4, 'a':3., 'σ':1.0, 'tol':1e-4 }

#######################################
############   High Pe   ##############
#######################################
def η(x, t=0.1, κ=3e-6, a=3e-6, x_0=50):
    return (4 / (1*(4*np.pi*κ*t)**0.5)) * np.exp(-(((x-(x_0 + a*t))**2)/(4*κ*t)))

#nxs    = np.logspace(2,3.5,3)
kappas = np.logspace(-6,4,20)*3
error  = np.zeros(kappas.shape[0])


for i,kappa in enumerate(kappas):
    #print('('+str(i+1)+'/' +str(len(kappas))+') Kappa ----->'+str(kappa))
    # for j, nx in enumerate(nxs):
    #     print('\t'+'('+str(i+1)+'/' +str(len(kappas))+') nx ----->'+str(nx)
    coeffs['κ'] = kappa
    #params['nx'] = int(nx)
    model = AdvDiff(params,coeffs)
    model.U[0,:] = η(model.x,model.dt,model.κ,model.a)
    Beam = model.Goundov('BeamWarming','w')
    anal = np.zeros_like(Beam)
    for t in range(model.nt):
        anal[t,:] = η(model.x,(t*model.dt)+model.dt,model.κ,model.a)

    error[i] = LA.norm(Beam-anal,np.inf)

# params = {'L':100.,'nx':1000,'nt':25}
# coeffs = {'κ':3e-1, 'a':3., 'σ':1.0, 'tol':1e-8 }
#
# nxs = np.logspace(2,3.5,10)
# error_nx  = np.zeros(nxs.shape[0])
#
# for j, nx in enumerate(nxs):
#     #print('('+str(i+1)+'/' +str(len(kappas))+') Kappa ----->'+str(kappa))
#     # for j, nx in enumerate(nxs):
#     # print('\t'+'('+str(j+1)+'/' +str(len(nxs))+') nx ----->'+str(nx))
#     params['nx'] = int(nx)
#     #params['nx'] = int(nx)
#     model = AdvDiff(params,coeffs)
#     model.U[0,:] = η(model.x,model.dt,model.κ,model.a)
#     Beam = model.Goundov('BeamWarming','w')
#     anal_nx = np.zeros_like(Beam)
#     print((model.a * model.dx)/model.κ)
#     for t in range(model.nt):
#         anal_nx[t,:] = η(model.x,(t*model.dt)+model.dt,model.κ,model.a)
#
#     error_nx[j] = LA.norm(Beam-anal_nx,np.inf)


fig, ax = plt.subplots(1,1,sharey=True,figsize=(8,6))

ax.loglog((model.a * model.dx)/kappas, error,'x-',label='CN-BeamWarming', markersize=8)

ax.grid()
plt.legend()
ax.set_title(r'$\Delta x =$ {:.3f}'.format(model.dx) )
ax.set_ylabel(r'$||U - u||_\infty$',fontsize='large')
ax.set_xlabel(r'$Pe$')
plt.show()
