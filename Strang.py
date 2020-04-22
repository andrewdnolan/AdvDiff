#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.solvers import TDMA
from advdiff.model import advection


########################################################
#################   Init. Constant   ###################
########################################################
κ  = 3                   # Diffusivity
a  = 3                    # wave speed
L  = np.pi*24.             # Domain Length
nx = 1000                  # Num. grid cells
dx = L/(nx-1)              # grid spacing

nt = 8000                   # Num time steps
σ  = .75                   # courant number
dt = (σ*dx)**2/κ           # time step

M = 1
P = 1
x_0 = 2*np.pi

########################################################
##################   Init. Domain   ####################
########################################################
r  = (κ*dt)/(2*dx**2) /2   # matrix const.
x  = np.linspace(dx,L,nx)  # spatial grid
u  = np.zeros((nt+1,nx))   # solution array
exact  =  np.zeros((nt+1,nx))
t = 0.1
u[0,:] = (M / (P *(4*np.pi*κ*t)**0.5)) * np.exp(-(((x-(x_0 + a*t))**2)/(4*κ*t)))
exact[0,:] = (M / (P *(4*np.pi*κ*t)**0.5)) * np.exp(-(((x-(x_0 + a*t))**2)/(4*κ*t)))
# init. condition

# Mat. A
A = np.diagflat([[(1+2*r) for __ in range(nx)]]) + \
    np.diagflat([[  (-r)  for __ in range(nx-1)]],1) +\
    np.diagflat([[  (-r)  for __ in range(nx-1)]],-1)

# Mat. B
B = np.diagflat([[(1-2*r) for __ in range(nx)]]) + \
    np.diagflat([[  (r)   for __ in range(nx-1)]],1) +\
    np.diagflat([[  (r)   for __ in range(nx-1)]],-1)

u_star  = np.zeros((nt+1,nx))
u_star2 = np.zeros((nt+1,nx))

err = np.zeros(nt)
for t in range(0,nt):
    # U^* &= \mathcal{N}_{\mathcal{A}}(U^n,k)
    b     = B.dot(u[t,:])                 # left vect.
    b[0]  = b[1]                          # boundary cond.
    b[-1] = b[-2]                         # boundary cond.
    u_star[t,:]  = TDMA(A,b).flatten()    # itterative solv.
    u_star[t,0]  = u_star[t,1]            # boundary cond.
    u_star[t,-1] = u_star[t,-2]           # boundary cond.

    #U^{n+1} &= \mathcal{N}_{\mathcal{B}}(U^*,k)
    u_star2[t+1,1:-1] = u_star[t,1:-1] - ((dt*a)/dx) * (u_star[t,1:-1] - u_star[t,0:-2])

    b     = B.dot(u_star2[t+1,:])                 # left vect.
    b[0]  = b[1]                          # boundary cond.
    b[-1] = b[-2]                         # boundary cond.
    u[t+1,:]  = TDMA(A,b).flatten()    # itterative solv.
    u[t+1,0]  = u[t+1,1]            # boundary cond.
    u[t+1,-1] = u[t+1,-2]           # boundary cond.

    exact[t+1,:] = (M / (P *(4*np.pi*κ*(0.1 + t*dt))**0.5)) * np.exp(-(((x-(x_0 + a*(0.1 + t*dt)))**2)/(4*κ*(0.1 + t*dt))))
    #plt.plot(x,u[t,:])
    err[t] = LA.norm(u[t,:] - exact[t,:],2)

plt.plot(x,exact[-1,:])
plt.plot(x,u[t-1,:],'D')
plt.show()
