#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.solvers import TDMA
from advdiff.model import advection

########################################################
#################   Init. Constant   ###################
########################################################

a  = 3                     # wave speed
L  = np.pi*2.              # Domain Length
nx = 100                   # Num. grid cells
dx = L/(nx-1)              # grid spacing

nt = 100                   # Num time steps
σ  = .75                   # courant number
dt = 3e-3 #(σ*dx)**2/a     # time step


def η(x,t):
    return 5*np.sin(((x-(a*t))-np.pi))+5

########################################################
##################   Init. Domain   ####################
########################################################
x  = np.linspace(dx,L,nx)  # spatial grid
u  = np.zeros((nt,nx))   # solution array
exact  =  np.zeros((nt,nx))

u[0,:] = η(x,0)
exact[0,:] = η(x,0)


for t in range(0,nt-1):
    u[t+1,1:] = u[t,1:] - ((dt*a)/dx)\
                        * (u[t,1:] - u[t,0:-1])

    u[t+1,0] = u[t,0] - ((dt*a)/dx)\
                        * (u[t,0] - u[t,-1])

    exact[t+1,:] = η(x,t*dt+dt)


    plt.plot(x,u[t,:],'x')
    plt.plot(x,exact[t,:])

plt.show()
