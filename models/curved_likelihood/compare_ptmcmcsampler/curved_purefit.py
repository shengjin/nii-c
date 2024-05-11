#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pylab as plt
import time
import corner
import glob
import scipy.optimize as so
import scipy.linalg as sl
from PTMCMCSampler import PTMCMCSampler

from mpi4py import MPI



class CurvedLikelihood(object):
    
    def __init__(self):
        
        self.pmin = np.array([-10, -10])
        self.pmax = np.array([10, 10])
    
    def lnlikefn(self, x):
        
        # x = x[0],  y = x[1]
        
        ll = np.exp(-x[0]**2 - (9+4*x[0]**2 + 9*x[1])**2) + 0.5 * np.exp(-8*x[0]**2-8*(x[1]-2)**2)
        
        return np.log(ll)
    
    
    def lnpriorfn(self, x):
        
        if np.all(self.pmin < x) and np.all(self.pmax > x):
            return 0.0
        else:
            return -np.inf  
        return 0.0

    
cl = CurvedLikelihood()


p0 = np.array([-0.1, -0.5])

cov = np.diag([1.0, 1.0])


start = time.time()      
sampler = PTMCMCSampler.PTSampler(2, cl.lnlikefn, cl.lnpriorfn, np.copy(cov), outDir='./chains')

sampler.sample(p0, 1000000, ladder=np.array([1.0,0.3,0.1,0.03]), burn=100000, thin=1, SCAMweight=10, AMweight=10, DEweight=10, NUTSweight=10, HMCweight=10, MALAweight=0, HMCsteps=50, HMCstepsize=0.08)

end = time.time()
multi_time = end - start 
print(multi_time) 










