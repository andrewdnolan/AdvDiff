#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA
import matplotlib.pyplot as plt
from advdiff.model import AdvDiff

params = {'L':12.*np.pi,'nx':100,'nt':100}
coeffs = {'κ':0.3, 'σ':0.5, 'a':3.}

def η(test, M = 1,P = 1,x_0 = np.pi,t=0.1):
    return (M / (P *(4*np.pi*test.κ*t)**0.5)) * np.exp(-(((test.x-(x_0 + test.a*t))**2)/(4*test.κ*t)))

test1 = AdvDiff(params,coeffs)
test2 = AdvDiff(params,coeffs)
test1.U[0,:] = η(test1,t=0.1)
test2.U[0,:] = η(test2,t=0.1)

LaxWendroff_G = test1.Goundov('LaxWendroff','w')
LaxWendroff_S = test2.Strang('LaxWendroff','w')

anal = np.zeros_like(LaxWendroff_S)
for t in range(test1.nt-1):
    anal[t+1,:] = η(test1,t=(t*test1.dt)+0.1)

plt.plot(test1.x,anal[-1,:])
plt.plot(test1.x,LaxWendroff_G[-1,:],'x')
plt.plot(test1.x,LaxWendroff_S[-1,:],'+')
plt.show()

error =np.zeros((2,test1.nt,test1.nx))
error[0,:,:] = LaxWendroff_G - anal # FTCS (0)
error[1,:,:] = LaxWendroff_S - anal # Up-wind (1)

names = ['Goundov','Strang']

print('Error (infty norm)')
for i,method in zip(range(error.shape[0]), names):
    print('\t{}: \n\t\t{}'.format(method,LA.norm(error[i],np.inf)))
