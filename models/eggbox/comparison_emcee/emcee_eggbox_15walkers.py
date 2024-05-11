###########################
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import time
import math

from multiprocessing import Pool

from multiprocessing import cpu_count

##############
import os
os.environ["OMP_NUM_THREADS"] = "1"
# Some builds of NumPy (including the version included with Anaconda) will automatically parallelize some operations using something like the MKL linear algebra. This can cause problems when used with the parallelization methods described here so it can be good to turn that off (by setting the environment variable OMP_NUM_THREADS=1, for example).
##############





############################################
############################################
############################################
# prior
a_min=              0
a_max=              50
b_min=              0
b_max=              50

def log_prior_a(a):
    if a_min < a < a_max:
        return np.log(1/(b_max-b_min))
    return -np.inf

def log_prior_b(b):
    if b_min < b < b_max:
        return np.log(1/(b_max-b_min))
    return -np.inf


def log_prior(theta):
    a = theta[0]
    b = theta[1]
    log_a = log_prior_a(a)
    log_b = log_prior_b(b)
    log_prior = log_a + log_b
    #print("a,b,d,la,lb,ld:", a, b, d, log_a, log_b, log_d)
    return log_prior


def log_likelihood(theta):
    a = theta[0]
    b = theta[1]
    #
    logll = math.sin(b/2.0)*math.cos(a/2.0)+1
    return logll


def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)


theta = np.random.rand(2)
print(log_prior(theta),log_likelihood(theta),log_probability(theta))




pos = np.random.randn(15, 2)
nwalkers, ndim = pos.shape

N_threads = 15
with Pool(N_threads) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(), pool=pool)
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    start = time.time()
    sampler.run_mcmc(pos, 405000, progress=True);
    end = time.time()
    multi_time = end - start
    print(multi_time)
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))
    print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

ncpu = cpu_count()
print("{0} CPUs".format(ncpu))


flat_samples = sampler.get_chain(discard=5000, thin=1, flat=True)
#print(flat_samples.shape)
#np.savetxt('emcee.out', np.transpose([theta[:,0], theta[:,1]]))
np.savetxt('emcee.out', flat_samples) #np.transpose([a,b,c])







