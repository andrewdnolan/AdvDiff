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
        self.U   = np.zeros((self.nt,self.nx))

    def __init__(self,params):
        self.parameters(**params)
        self.domain()

    def FTCS(self,j,ustar=None):

        if isinstance(ustar, np.ndarray):
            ustar[j+1,1:-1]  = ustar[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
                                * (ustar[j,2:] - ustar[j,0:-2])

            ustar[j+1, 0] = ustar[j,0] - ((self.dt*self.a)/(self.dx*2))\
                                * (ustar[j,1] - ustar[j,-1])

            ustar[j+1,-1] = ustar[j, -1 ] - ((self.dt*self.a)/(self.dx*2))\
                                * (ustar[j,0] - ustar[j,-2])

            return ustar[j+1,:]

        else:
            self.U[j+1,1:-1]  = self.U[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
                                * (self.U[j,2:] - self.U[j,0:-2])
            self.U[j+1, 0] = self.U[j,0] - ((self.dt*self.a)/(self.dx*2))\
                                * (self.U[j,1] - self.U[j,-1])
            self.U[j+1,-1] = self.U[j, -1 ] - ((self.dt*self.a)/(self.dx*2))\
                                * (self.U[j,0] - self.U[j,-2])

    def UpWind(self,j,ustar=None):
        """ For only one time step!!!!
        """
        if isinstance(ustar, np.ndarray):
            ustar[j+1,1:] = ustar[j,1:] - ((self.dt*self.a)/self.dx)\
                                * (ustar[j,1:] - ustar[j,0:-1])

            ustar[j+1,0] = ustar[j,0] - ((self.dt*self.a)/self.dx)\
                                * (ustar[j,0] - ustar[j,-1])
            return ustar[j+1,:]

        else:
            self.U[j+1,1:] = self.U[j,1:] - ((self.dt*self.a)/self.dx)\
                                * (self.U[j,1:] - self.U[j,0:-1])

            self.U[j+1,0] = self.U[j,0] - ((self.dt*self.a)/self.dx)\
                                * (self.U[j,0] - self.U[j,-1])

    def LaxWendroff(self,j,ustar=None):
        if isinstance(ustar, np.ndarray):

            ustar[j+1,1:-1] = ustar[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
            *(ustar[j,2:] - ustar[j,0:-2]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(ustar[j,2:] - 2*ustar[j,1:-1] + ustar[j,0:-2])

            ustar[j+1,0] = ustar[j,0] - ((self.dt*self.a)/(self.dx*2))\
            *(ustar[j,1] - ustar[j,-1]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(ustar[j,1] - 2*ustar[j,0] + ustar[j,-1])

            ustar[j+1,-1] = ustar[j,-1] - ((self.dt*self.a)/(self.dx*2))\
            *(ustar[j,0] - ustar[j,-2]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(ustar[j,0] - 2*ustar[j,-1] + ustar[j,-2])

            return ustar[j+1,:]

        else:
            self.U[j+1,1:-1] = self.U[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
            *(self.U[j,2:] - self.U[j,0:-2]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(self.U[j,2:] - 2*self.U[j,1:-1] + self.U[j,0:-2])

            self.U[j+1,0] = self.U[j,0] - ((self.dt*self.a)/(self.dx*2))\
            *(self.U[j,1] - self.U[j,-1]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(self.U[j,1] - 2*self.U[j,0] + self.U[j,-1])

            self.U[j+1,-1] = self.U[j,-1] - ((self.dt*self.a)/(self.dx*2))\
            *(self.U[j,0] - self.U[j,-2]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(self.U[j,0] - 2*self.U[j,-1] + self.U[j,-2])

    def BeamWarming(self,j,ustar=None):
        if isinstance(ustar, np.ndarray):
            ustar[j+1,2:] = ustar[j,2:] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*ustar[j,2:] - 4*ustar[j,1:-1] + ustar[j,0:-2]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(ustar[j,2:] \
                        - 2*ustar[j,1:-1] + ustar[j,0:-2])

            ustar[j+1, 0] = ustar[j,0] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*ustar[j,0] - 4*ustar[j,-1] + ustar[j,-2]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(ustar[j,0] \
                        - 2*ustar[j,-1] + ustar[j,-2])

            ustar[j+1, 1] = ustar[j,1] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*ustar[j,1] - 4*ustar[j,0] + ustar[j,-1]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(ustar[j,1] \
                        - 2*ustar[j,0] + ustar[j,-1])

            ustar[j+1, 2] = ustar[j,2] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*ustar[j,2] - 4*ustar[j,1] + ustar[j, 0]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(ustar[j,2] \
                        - 2*ustar[j,1] + ustar[j, 0])

            return ustar[j+1,:]

        else:
            self.U[j+1,2:] = self.U[j,2:] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*self.U[j,2:] - 4*self.U[j,1:-1] + self.U[j,0:-2]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(self.U[j,2:] \
                        - 2*self.U[j,1:-1] + self.U[j,0:-2])

            self.U[j+1, 0] = self.U[j,0] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*self.U[j,0] - 4*self.U[j,-1] + self.U[j,-2]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(self.U[j,0] \
                        - 2*self.U[j,-1] + self.U[j,-2])

            self.U[j+1, 1] = self.U[j,1] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*self.U[j,1] - 4*self.U[j,0] + self.U[j,-1]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(self.U[j,1] \
                        - 2*self.U[j,0] + self.U[j,-1])

            self.U[j+1, 2] = self.U[j,2] - ((self.dt*self.a)/(self.dx*2))\
                        * (3*self.U[j,2] - 4*self.U[j,1] + self.U[j, 0]) \
                        + ((self.dt*self.a)**2/(2*self.dx**2))*(self.U[j,2] \
                        - 2*self.U[j,1] + self.U[j, 0])

class advection(model):
    '''1-D Advection Class

    Keyword arguments:
    -----------------
        params -- dictionary of numerical parameters
            L  -- length of domain
            nx -- number of grid cell
            nt -- number of time steps

        adv_params -- dictionary of parameters Adv. Eqn.
            a  -- advection velocity (scalar)
            σ  -- courant numner
            dt -- (optitonal) time step. If the time step is not
                  specificed will calculate based on the CFL cond.
    \n
    '''
    def adv_parameters(self,a,σ,dt=None):
        self.a  = a
        self.σ  = σ
        if dt != None:
            self.dt = dt
        else:
            self.dt = (σ*self.dx)/a

    def __init__(self,params,adv_params):
        super().__init__(params)
        self.adv_parameters(**adv_params)

    def boundray_conditions(self,type):
        if   type == 'Periodic':
            self.BC = 'PBC'
        elif type == 'Dirichlet':
            self.BC = 'DBC'
        elif type == 'Neumann':
            self.BC = 'NBC'
        else:
            print("\n\tBoundary Conditions Must be:\n\t\
            Periodic, Dirichlet, or Neumann\n")

    def solver(self,solver):
        self.solver = solver

    def run(self,solver='UpWind',w=None):
        """Run the Integration Scheme

        This function numerically solves the advection equation and writes
        the solution to self.U.

        NOTE: you must set the initial condition self.U[0,:] for schemes to
        run.

        Keyword arguments:
        -----------------
            solver -- {'FTCS', 'UpWind', 'LaxWendroff', 'BeamWarming'}
                FTCS   -- 1st order symetirc (Foward Time Centered Space)
                Upwind -- 1st order non-symetric
                LaxWendroff -- 2nd order symetirc
                BeamWarming -- 2nd order non-symetric

        Example:
        --------
        >>>> import numpy as np
        >>>> from advdiff.model import advection

        >>>> params     = {'L':2*np.pi,'nx':100,'nt':100}
        >>>> adv_params = {'a':3, 'σ':0.7, 'dt':0.0025}

        >>>> model = advection(params,adv_params)
        >>>> model.U[0,:] = 5*np.sin((model.x-np.pi/2))+5
        >>>> model.run('UpWind')

        """

        if LA.norm(self.U[0,:]) == 0:
            raise Exception('Initial condition not set')
        else:
            #self.solver(solver)
            if isinstance(w,str):
                U = self.U.copy()
                for t in range(0,self.nt-1):
                    U[t+1,:] = self.__getattribute__(solver)(t,U)
                return U
            else:
                for t in range(0,self.nt-1):
                    self.__getattribute__(solver)(t)

class Diffusion(model):
    '''1-D Diffusion Class

    Keyword arguments:
    -----------------
        params -- dictionary of numerical parameters
            L  -- length of domain
            nx -- number of grid cell
            nt -- number of time steps

        adv_params -- dictionary of parameters Adv. Eqn.
            κ  -- diffusivity (scalar)
            σ  -- courant numner
            dt -- (optitonal) time step. If the time step is not
                  specificed will calculate based on the CFL cond.
    '''

    def diff_parameters(self,κ,σ,dt=None):
        self.κ  = κ
        self.σ  = σ
        if dt != None:
            self.dt = dt
        else:
            self.dt = (σ*self.dx)/κ

    def __init__(self,params,diff_params):
        super().__init__(params)
        self.diff_parameters(**diff_params)

    def crank_nicolson(self,t,U=None,dt=None):
        """ Crank Nicolson

        Keyword arguments:
        -----------------
            w -- {'w' or True}
                 write solution to output array. Else if you leave it blank or
                 set to equal False will be written to self.U
        """
        from advdiff.solvers import TDMA

        # exception to deal with fractional time step for Strang splitting
        if dt == None:
            dt = self.dt
        else:
            pass

        r  = (self.κ*dt)/(2*self.dx**2)      # matrix const.

        if hasattr(self, 'A') == False:
            self.A = np.diagflat([[(1+2*r) for __ in range(self.nx)]]) + \
                np.diagflat([[  (-r)  for __ in range(self.nx-1)]],1) +\
                np.diagflat([[  (-r)  for __ in range(self.nx-1)]],-1)
            self.A[0,-1] = -r
            self.A[-1,0] = -r
            self.B = np.diagflat([[(1-2*r) for __ in range(self.nx)]]) + \
                np.diagflat([[  (r)   for __ in range(self.nx-1)]],1) +\
                np.diagflat([[  (r)   for __ in range(self.nx-1)]],-1)
            self.B[0,-1] = r
            self.B[-1,0] = r

        if isinstance(U, np.ndarray):
            b  = self.B.dot(U[t,:])            # left vect.
            U[t+1,:] = LA.solve(self.A,b).flatten()  # itterative solv.
            return U[t+1,:]
        else:
            b = self.B.dot(self.U[t,:])            # left vect.
            self.U[t+1,:] = LA.solve(self.A,b).flatten()  # itterative solv.

    def run(self,solver='crank_nicolson',w=None):
        if LA.norm(self.U[0,:]) == 0:
            raise Exception('Initial condition not set')
        else:
            #self.solver(solver)
            if isinstance(w,str):
                U = self.U.copy()
                for t in range(0,self.nt-1):
                    U[t+1,:] = self.__getattribute__(solver)(t,U)
                return U
            else:
                for t in range(0,self.nt-1):
                    self.__getattribute__(solver)(t)
class AdvDiff(model):
    '''1-D Advection Diffusion Class

    Keyword arguments:
    -----------------
        params -- dictionary of numerical parameters
            L  -- length of domain
            nx -- number of grid cell
            nt -- number of time steps

        coeffs -- dictionary of parameters Adv. Eqn.
            κ  -- diffusivity (scalar)
            σ  -- courant numner
            dt -- (optitonal) time step. If the time step is not
                  specificed will calculate based on the CFL cond.
    '''

    def coefficients(self,κ,σ,a,dt=None):
        self.κ  = κ
        self.a  = a
        self.σ  = σ

        if dt != None:
            self.dt = dt
        else:
            self.dt = (σ*self.dx)/a

    def __init__(self,params,adv_params):
        super().__init__(params)
        self.coefficients(**adv_params)
        self.set_solvers()

    def set_solvers(self, explicit='LaxWendroff'):
        self.explict = explicit

    def init_CN(self,dt=None):

        # exception to deal with fractional time step for Strang splitting
        if dt == None:
            dt = self.dt
        else:
            pass

        self.r  = (self.κ*dt)/(2*self.dx**2)      # matrix const.

        self.A = np.diagflat([[(1+2*self.r) for __ in range(self.nx)]]) + \
            np.diagflat([[  (-self.r)  for __ in range(self.nx-1)]],1) +\
            np.diagflat([[  (-self.r)  for __ in range(self.nx-1)]],-1)
        self.A[0,-1] = -self.r
        self.A[-1,0] = -self.r

        self.B = np.diagflat([[(1-2*self.r) for __ in range(self.nx)]]) + \
            np.diagflat([[  (self.r)   for __ in range(self.nx-1)]],1) +\
            np.diagflat([[  (self.r)   for __ in range(self.nx-1)]],-1)
        self.B[0,-1] = self.r
        self.B[-1,0] = self.r

    def crank_nicolson(self,t,U=None,dt=None):
        """ Crank Nicolson

        Keyword arguments:
        -----------------
            w -- {'w' or True}
                 write solution to output array. Else if you leave it blank or
                 set to equal False will be written to self.U
        """
        from advdiff.solvers import TDMA

        if isinstance(U, np.ndarray):
            b  = self.B.dot(U[t,:])            # left vect.
            U[t+1,:] = LA.solve(self.A,b).flatten()  # itterative solv.
            return U[t+1,:]

        else:
            b = self.B.dot(self.U[t,:])            # left vect.
            self.U[t+1,:] = LA.solve(self.A,b).flatten()  # itterative solv.

    def Goundov(self,explict='UpWind',w=True):
        self.init_CN()
        u_star = np.zeros((self.nt,self.nx))

        if w == 'w' or w == True:
            U = self.U.copy()
            for t in range(0,self.nt-1):
                u_star[t,:] = self.crank_nicolson(t,U)
                U[t+1,:] = self.__getattribute__(explict)(t,u_star)
            return U

        else:
            for t in range(0,self.nt-1):
                u_star[t,:] = self.crank_nicolson(t)
                self.__getattribute__(explict)(t,u_star)

    def Strang(self,explict='UpWind',w=True):
        self.init_CN(self.dt/2.)
        u_star  = np.zeros((self.nt,self.nx))
        u_star2 = np.zeros((self.nt,self.nx))

        if w == 'w' or w == True:
            U = self.U.copy()
            for t in range(0,self.nt-1):
                u_star[t,:]  = self.crank_nicolson(t,U)
                u_star2[t,:] = self.__getattribute__(explict)(t,u_star)
                U[t+1,:] = self.crank_nicolson(t,u_star2)
            return U

        else:
            for t in range(0,self.nt-1):
                u_star[t,:]  = self.crank_nicolson(t)
                u_star2[t,:] = self.__getattribute__(explict)(t,u_star)
                self.crank_nicolson(t,u_star2)
