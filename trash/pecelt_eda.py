
# Global
import sys
import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titlesize'] = 14

# Local
sys.path.append('../')
from advdiff.model import advection,Diffusion

def η(x, t=0.1, κ=3e-6, a=3e-6, x_0=50):
    return (4 / (1*(4*np.pi*κ*t)**0.5)) * np.exp(-(((x-(x_0 + a*t))**2)/(4*κ*t)))

# def η(x, t=0.1, κ=3e-6, a=3e-6, x_0=50):
#     return ( np.exp(-(((x - a*t))/(4*κ*t))))

params = {'L':100,'nx':1000,'nt':25}
adv_coef  = {'a':3e0, 'σ':0.5}
diff_coef = {'κ':2e-2, 'σ':0.5}

dx = params['L']/params['nx']
dt = (1*dx)**2/adv_coef['a']
x  = np.linspace(0,params['L'],params['nx'])

y0  =  η(x, dt, diff_coef['κ'],adv_coef['a'])
y50 = η(x,10*dt+dt, diff_coef['κ'],adv_coef['a'])

plt.plot(x,y0)
plt.plot(x,y50)

plt.show()
print((adv_coef['a']* dx)/diff_coef['κ'])
