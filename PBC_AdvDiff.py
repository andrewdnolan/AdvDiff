#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.solvers import TDMA
from advdiff.model import advection

def f(x,t):
    D = 3e-3
    k = 0.5
    U = 3e-2
    A = 1.
    '''Analytical Solution to the Advec. Diff. Eqn.
    '''
    return np.exp(-D*k**2*t)*np.sin(2*np.pi*k*(x-U*t))

########################################################
#################   Init. Constant   ###################
########################################################
κ  = 3e-3                  # Diffusivity
a  = 3e-2                  # wave speed
L  = np.pi*2.              # Domain Length
nx = 500                   # Num. grid cells
dx = L/(nx-1)              # grid spacing

nt = 2                     # Num time steps
σ  = .75                   # courant number
dt = (σ*dx)**2/κ           # time step


########################################################
##################   Init. Domain   ####################
########################################################
r  = (κ*dt)/(2*dx**2)      # matrix const.
x  = np.linspace(dx,L,nx)  # spatial grid
u  = np.zeros((nt,nx))   # solution array
exact  =  np.zeros((nt,nx))

u[0,:] = f(x,0)
exact[0,:] = f(x,0)
# init. condition

# Mat. A
A = np.diagflat([[(1+2*r) for __ in range(nx)]]) + \
    np.diagflat([[  (-r)  for __ in range(nx-1)]],1) +\
    np.diagflat([[  (-r)  for __ in range(nx-1)]],-1)
    
A[0,-1] = -r
A[-1,0] = -r

# Mat. B
B = np.diagflat([[(1-2*r) for __ in range(nx)]]) + \
    np.diagflat([[  (r)   for __ in range(nx-1)]],1) +\
    np.diagflat([[  (r)   for __ in range(nx-1)]],-1)
# B[0,-1] = r
# B[-1,0] = r

u_star = np.zeros((nt,nx))

for t in range(0,nt-1):
    # U^* &= \mathcal{N}_{\mathcal{A}}(U^n,k)
    b     = B.dot(u[t,:])                 # left vect.
    u_star[t,:]  = LA.solve(A,b).flatten()    # itterative solv.

    u[t+1,1:] = u_star[t,1:] - ((dt*a)/dx)\
                        * (u_star[t,1:] - u_star[t,0:-1])

    u[t+1,0] = u_star[t,0] - ((dt*a)/dx)\
                        * (u_star[t,0] - u_star[t,-1])

    exact[t+1,:] = f(x,t*dt+dt)

plt.plot(x,exact[-1,:])
plt.plot(x,u[-1,:],'x')
plt.show()
