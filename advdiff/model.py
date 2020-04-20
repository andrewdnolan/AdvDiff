import os
import json
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

    def init_cond(self):
        pass

    def analytical(self,):
        pass

    def dump(self):
        '''Dump params. to HDF5 file
        '''
        with open(os.path.join(self.out_dir,'params.json'),'w') as fp:
            json.dump()
        pass

    def log(self):
        '''Store info for log file
        '''
        pass

class advection(model):

    def adv_parameters(self,a,σ,dt=None):
        self.a  = a
        self.σ  = σ
        if dt != None:
            self.dt = dt
        else:
            self.dt = (σ*self.dx)**2/a

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

    def FTCS(self):

        for j in range(0,self.nt-1):
            self.U[j+1,1:-1]  = self.U[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
                                * (self.U[j,2:] - self.U[j,0:-2])

            self.U[j+1, 0] = self.U[j+1,-2]
            self.U[j+1,-1] = self.U[j+1,1]

    def UpWind(self):
        for j in range(0,self.nt-1):
            self.U[j+1,1:-1] = self.U[j,1:-1] - ((self.dt*self.a)/self.dx)\
                                * (self.U[j,1:-1] - self.U[j,0:-2])

            self.U[j+1, 0] = self.U[j+1,-2]
            self.U[j+1,-1] = self.U[j+1,1]

    def LaxWendroff(self):
        for j in range(0,self.nt-1):
            self.U[j+1,1:-1] = self.U[j,1:-1] - ((self.dt*self.a)/(self.dx*2))\
            *(self.U[j,2:] - self.U[j,0:-2]) + ((self.dt*self.a)**2/(2*self.dx**2))\
            *(self.U[j,2:] - 2*self.U[j,1:-1] + self.U[j,0:-2])

            self.U[j+1, 0] = self.U[j+1,-2]
            self.U[j+1,-1] = self.U[j+1,1]

    def BeamWarming(self):
        for j in range(0,self.nt-1):
            for i in range(1,self.nx-1):
                self.U[j+1,i] = self.U[j,i] - ((self.dt*self.a)/(self.dx*2))\
                            * (3*self.U[j,i] - 4*self.U[j,i-1] + self.U[j,i-2]) \
                            + ((self.dt*self.a)**2/(2*self.dx**2))*(self.U[j,i] \
                            - 2*self.U[j,i-1] + self.U[j,i-2])

            self.U[j+1, 0] = self.U[j+1,-2]
            self.U[j+1,-1] = self.U[j+1,1]

    def run(self,solver='UpWind'):
        if LA.norm(self.U[:,0]) == 0:
            raise Exception('Initial condition not set')
        else:
            self.__getattribute__(solver)()
