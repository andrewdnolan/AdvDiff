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


#################   Init. Constant   ###################
κ  = 0.07                  # Diffusivity
L  = 2*np.pi               # Domain Length
nx = 201                   # Num. grid cells
dx = L/(nx-1)              # grid spacing

nt = 113                   # Num time steps
dt = dx * κ                # time step

##################   Init. Domain   ####################

r  = (κ*dt)/(2*dx**2)      # matrix const.
x  = np.linspace(dx,L,nx)  # spatial grid

numer  = np.zeros((nt,nx))   # solution array
exact  =  np.zeros((nt,nx))

numer[0,:] = ufunc(0,x,κ)
exact[0,:] = ufunc(0,x,κ)

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
B[0,-1] = r
B[-1,0] = r

u_star = np.zeros((nt,nx))
u_star2 = np.zeros((nt,nx))

error = np.zeros((nt,nx))
Pe    = np.zeros((nt,nx))
for j in range(0,nt-1):
    b     = B.dot(numer[j,:])                 # left vect.
    u_star[j,:]  = LA.solve(A,b)    # itterative solv.
    u_star[j,-1] = u_star[j,-1]

    numer[j+1,1:] = u_star[j,1:] - ((dt*u_star[j,1:])/dx) * (u_star[j,1:] - u_star[j,0:-1])
    numer[j+1,0] = u_star[j,0] - ((dt*u_star[j,0])/dx)* (u_star[j,0] - u_star[j,-1])
    exact[j+1,:] = ufunc(dt*j,x,κ)

    error[j,:] = np.abs(exact[j,:]-numer[j,:])
    Pe[j,:] = (exact[j,:] * dx) / κ

from matplotlib import animation, rc
rc('animation', html='jshtml')

fig, ax = plt.subplots(1,figsize=(15,8))

ax.set_xlim(0, L)
ax.set_ylim(-0.2, 6.25)
ax.set_ylabel('$Pe$',fontsize='x-large')
ax.set_xlabel(' x ')

ax2 = ax.twinx()
ax2.set_xlim(0, L)
ax2.set_ylim(0, 5.75)
ax2.set_ylabel('$|U-u|$',rotation=270,fontsize='x-large',labelpad=25)
ax2.set_xlabel(' x ')

line0, = ax.plot([], [], 'x-', lw=3, color = 'darkslategrey',label='$Pe$')
line1, = ax2.plot([], [], 'D-',lw=3, color = 'darkmagenta',label='Error')


line = [line0,line1]

# added these three lines
labs = [l.get_label() for l in line]
ax.legend(line, labs, loc=0)

def animate(i):
    global error,Pe
    line[0].set_data(x, Pe[i,:])
    line[1].set_data(x, error[i,:])
    return line

anim = animation.FuncAnimation(fig, animate,
                               frames=range(0,nt), interval=100, blit=True)

plt.close()
anim
