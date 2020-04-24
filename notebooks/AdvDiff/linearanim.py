# Global
import sys
import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt

# Local
sys.path.append('../')
from advdiff.plot import animation
from advdiff.model import advection,Diffusion
plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['text.usetex'] = True

def η(x, t=0.1, κ=3e-6, a=3e-6, x_0=np.pi):
    return (1 / ((4*np.pi*κ*t)**0.5)) * np.exp(-(((x-(x_0 + a*t))**2)/(4*κ*t)))

params = {'L':6.*np.pi,'nx':500,'nt':1191}
adv_coef  = {'a':3, 'σ':0.1}
diff_coef = {'κ':3e-1, 'σ':0.1, 'tol':1e-6}

Adv  = advection(params,adv_coef)
Adv.U[0,:] = η(Adv.x,0.1,diff_coef['κ'],adv_coef['a'])

Diff = Diffusion(params,diff_coef)
Diff.U[0,:] = η(Adv.x,0.1,diff_coef['κ'],adv_coef['a'])

UStar = np.zeros_like(Diff.U)
U = np.copy(Diff.U)  # numerical sol. array
u = np.copy(Diff.U)  # analytical sol. array

for t in range(params['nt']-1):
    UStar[t,:] = Diff.crank_nicolson(t,U,Adv.dt)
    U[t+1,:]   = Adv.BeamWarming(t,UStar)

    u[t+1,:] = η(Adv.x,(t*Adv.dt)+0.1, diff_coef['κ'],adv_coef['a'])

from matplotlib import animation, rc
rc('animation', html='jshtml')

fig, ax = plt.subplots(1,figsize=(15,8))

ax.set_xlim(0, params['L'])
ax.set_ylim(-0.2, 1.75)

ax.set_ylabel('Amp.')
ax.set_xlabel(' x ')

line0, = ax.plot([], [], lw=2, color = 'darkslategrey',label='Analytical')
line1, = ax.plot([], [], 'D',lw=2, color = 'darkseagreen',label='CN-Beam')

line = [line0,line1]
ax.legend(loc=2)

def animate(i):
    global Diff, u, U
    line[0].set_data(Diff.x, u[i,:])
    line[1].set_data(Diff.x, U[i,:])
    return line

anim = animation.FuncAnimation(fig, animate,
                               frames=range(0,Diff.nt,10), interval=60, blit=True)

plt.close()
