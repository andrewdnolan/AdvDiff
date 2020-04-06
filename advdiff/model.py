import os
import json

class model(object):

    def __init__(self,out_dir='./output'):
        self.out_dir = out_dir

    def parameters(self,paramfile):
        '''Initalize Model Parameters
        '''
        with open(paramfile,'r') as fp:
            params = json.load(fp)
        self.parameters
        pass

    def numerical(self):
        '''Initalize Numerical Parameters
        '''
        self.tol = 1e-6
        pass

    def dump(self):
        '''Dump params. to HDF5 file
        '''
        self.parameters =
        with open(os.path.join(self.out_dir,'params.json'),'w') as fp:
            json.dump()
        pass

    def log(self):
        '''Store info for log file
        '''
        pass
