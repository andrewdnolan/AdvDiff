# Global
import sys
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt

# Local
sys.path.append('../')
from advdiff.solvers import TDMA
from advdiff.plot import animation


########################################################
#################   Init. Constant   ###################
########################################################
a  = 3                     # wave speed
L  = 2*np.pi               # Domain Length
nx = 60                    # Num. grid cells
dx = L/(nx-1)              # grid spacing

nt = 1000                  # Num time steps
σ  = 0.5                   # courant number
dt = (σ*dx)**2/a           # time step

########################################################
#################      Init. Cond      #################
########################################################

def η(x,t,a=3.):
    return 5*np.sin(2*((x-(a*t))-np.pi))+5

########################################################
##################   Init. Domain   ####################
########################################################
x  = np.linspace(dx,L,nx)  # spatial grid
u  = np.zeros((5,nt,nx))   # (num methods) X (nx) X (nt)
u[:,0,:] = η(x,0)          # init. condition

for j in range(0,nt-1):
    # FTCS (0)
    u[0,j+1,1:-1]  = u[0,j,1:-1] - ((dt*a)/(dx*2)) * (u[0,j,2:] - u[0,j,0:-2])
    u[0,j+1, 0] = u[0,j,0]  - ((dt*a)/(dx*2)) * (u[0,j,1] - u[0,j,-1])
    u[0,j+1,-1] = u[0,j,-1] - ((dt*a)/(dx*2)) * (u[0,j,0] - u[0,j,-2])

    # Up-wind (1)
    u[1,j+1,1:] = u[1,j,1:] - ((dt*a)/dx) * (u[1,j,1:] - u[1,j,0:-1])
    u[1,j+1,0]  = u[1,j,0]  - ((dt*a)/dx) * (u[1,j,0]  - u[1, j,-1])

    # Lax-Wendroff (2)
    u[2,j+1,1:-1] = u[2,j,1:-1] - ((dt*a)/(dx*2))*(u[2,j,2:] - u[2,j,0:-2]) + ((dt*a)**2/(2*dx**2)) *(u[2,j,2:] - 2*u[2,j,1:-1] + u[2,j,0:-2])
    u[2,j+1,0] =  u[2,j,0]  - ((dt*a)/(dx*2)) *(u[2,j,1] - u[2,j,-1])  + ((dt*a)**2/(2*dx**2)) *(u[2,j,1] - 2*u[2,j,0]  + u[2,j,-1])
    u[2,j+1,-1] = u[2,j,-1] - ((dt*a)/(dx*2)) *(u[2,j,0] - u[2,j,-2]) + ((dt*a)**2/(2*dx**2)) *(u[2,j,0] - 2*u[2,j,-1] + u[2,j,-2])

    # Beam Warming (3)
    u[3,j+1,2:] = u[3,j,2:] - ((dt*a)/(dx*2))* (3*u[3,j,2:] - 4*u[3,j,1:-1] + u[3,j,0:-2])+ ((dt*a)**2/(2*dx**2))*(u[3,j,2:]- 2*u[3,j,1:-1] + u[3,j,0:-2])
    u[3,j+1, 0] = u[3,j,0]  - ((dt*a)/(dx*2))* (3*u[3,j,0]  - 4*u[3,j,-1] + u[3,j,-2]) + ((dt*a)**2/(2*dx**2)) * (u[3,j,0] - 2*u[3,j,-1] + u[3,j,-2])
    u[3,j+1, 1] = u[3,j,1]  - ((dt*a)/(dx*2))* (3*u[3,j,1]  - 4*u[3,j,0]  + u[3,j,-1]) + ((dt*a)**2/(2*dx**2)) * (u[3,j,1] - 2*u[3,j,0]  + u[3,j,-1])
    u[3,j+1, 2] = u[3,j,2]  - ((dt*a)/(dx*2))* (3*u[3,j,2]  - 4*u[3,j,1]  + u[3,j, 0]) + ((dt*a)**2/(2*dx**2)) * (u[3,j,2] - 2*u[3,j,1]  + u[3,j, 0])

    #Anlytical (4)
    u[4,j+1,:] = η(x,(j*dt)+dt)

from matplotlib import animation, rc
rc('animation', html='jshtml')

fig, ax = plt.subplots(1,figsize=(15,8))

ax.set_xlim(0, L)
ax.set_ylim(-0.2, 10.2)

ax.set_ylabel('Amp.')
ax.set_xlabel(' x ')

line0, = ax.plot([], [], 'x', lw=2, color='lightblue',label='FTCS')
line1, = ax.plot([], [], 'o', lw=2, color='#707BA7',label='Upwind')
line2, = ax.plot([], [], '3', lw=2, color='g',label='Lax-Wefford')
line3, = ax.plot([], [], '^', lw=2, color='b',label='BeanWarming')
line4, = ax.plot([], [], lw=2, color='k',label='Analtical',zorder=0)

line = [line0,line1, line2,line3,line4]
ax.legend(loc=3)

def animate(i):
    global x,u
    line[0].set_data(x, u[0,i,:])
    line[1].set_data(x, u[1,i,:])
    line[2].set_data(x, u[2,i,:])
    line[3].set_data(x, u[3,i,:])
    line[4].set_data(x, u[4,i,:])
    return line

anim = animation.FuncAnimation(fig, animate,
                               frames=range(0,nt,10), interval=100, blit=True)
plt.close()
anim
