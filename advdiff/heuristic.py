#!/usr/bin/env python3

import numpy as np
import scipy.linalg as LA

class model:

    def parameters(self,L,nx,nt):
        self.L  = L
        self.nx = nx
        self.dx = L/(nx-1)
        self.nt = nt

    def domain(self):
        self.tol = 1e-6
        self.x   = np.linspace(self.dx,self.L,self.nx)
        self.U   = np.zeros((self.nt+1,self.nx))

    def __init__(self,params):
        self.parameters(**params)
        self.domain()

class implicit:

    def init_matrix(self):
        pass

class explicit:
    pass

class AdvDiff(model):

    def coefficients(self,κ,σ,a,dt=None):
        self.κ  = κ
        self.a  = a
        self.σ  = σ
        if dt != None:
            self.dt = dt
        else:
            self.dt = (σ*self.dx)**2/a

    def __init__(self,params,adv_params):
        super().__init__(params)
        self.coefficients(**adv_params)
        #self.implicit = self.implicit()

    def init_CN(self,dt=0.001):

        # exception to deal with fractional time step for Strang splitting
        if dt == None:
            dt = self.dt
        else:
            pass

        self.r  = (self.κ*dt)/(2*self.dx**2)      # matrix const.

        self.A = np.diagflat([[(1+2*self.r) for __ in range(self.nx)]]) + \
            np.diagflat([[  (-self.r)  for __ in range(self.nx-1)]],1) +\
            np.diagflat([[  (-self.r)  for __ in range(self.nx-1)]],-1)

        self.B = np.diagflat([[(1-2*self.r) for __ in range(self.nx)]]) + \
            np.diagflat([[  (self.r)   for __ in range(self.nx-1)]],1) +\
            np.diagflat([[  (self.r)   for __ in range(self.nx-1)]],-1)

    def crank_nicolson(self,t,w=False,dt=None):
        """ Crank Nicolson

        Keyword arguments:
        -----------------
            w -- {'w' or True}
                 write solution to output array. Else if you leave it blank or
                 set to equal False will be written to self.U
        """
        from advdiff.solvers import TDMA

        if w == 'w' or w == True:

            U = self.U.copy()

            b     = self.B.dot(U[t,:])            # left vect.
            b[0]  = b[1]                     # boundary cond.
            b[-1] = b[-2]                    # boundary cond.
            U[t+1,:]  = TDMA(self.A,b).flatten()  # itterative solv.
            U[t+1,0]  = U[t+1,1]             # boundary cond.
            U[t+1,-1] = U[t+1,-2]            # boundary cond.

            return U[t+1,:]

        else:
            b     = self.B.dot(self.U[t,:])            # left vect.
            b[0]  = b[1]                          # boundary cond.
            b[-1] = b[-2]                         # boundary cond.
            self.U[t+1,:]  = TDMA(self.A,b).flatten()  # itterative solv.
            self.U[t+1,0]  = self.U[t+1,1]        # boundary cond.
            self.U[t+1,-1] = self.U[t+1,-2]       # boundary cond.

    def Goundov(self):

        self.init_CN()
        u_star = np.zeros((nt+1,nx))
               =
        for t in range(0,nt):
        pass
