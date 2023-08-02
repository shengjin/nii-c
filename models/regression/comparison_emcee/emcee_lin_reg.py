###########################
import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import time

from multiprocessing import Pool

from multiprocessing import cpu_count

##############
import os
os.environ["OMP_NUM_THREADS"] = "1"
# Some builds of NumPy (including the version included with Anaconda) will automatically parallelize some operations using something like the MKL linear algebra. This can cause problems when used with the parallelization methods described here so it can be good to turn that off (by setting the environment variable OMP_NUM_THREADS=1, for example).
##############



# data_file
datafile = "xy.dat"
data_delim = " "
header = 0


# true model
a = 3.0
b = 5
d = 2

data = np.genfromtxt(datafile, delimiter=data_delim, skip_header=header)
x = data[:,0]
y = data[:,1]

plt.scatter(x, y)
plt.plot(x, a*x+b, "k", alpha=0.3, lw=3)
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("xy.png")



############################################
############################################
############################################
# prior
a_min=              -20
a_max=              20
b_min=              -50
b_max=              50
d_min=              0.01
d_max=              30

def log_prior_a(a):
    if a_min < a < a_max:
        return np.log(1/(a_max-a_min)) 
    return -np.inf

def log_prior_b(b):
    if b_min < b < b_max:
        return np.log(1/(b_max-b_min))
    return -np.inf

def log_prior_d(d):
    if d_min < d < d_max:
        return np.log(1/(d*np.log(d_max/d_min)))
    return -np.inf

def log_prior(theta):
    a = theta[0]
    b = theta[1]
    d = theta[2]
    log_a = log_prior_a(a)
    log_b = log_prior_b(b)
    log_d = log_prior_d(d)
    log_prior = log_a + log_b +log_d
    #print("a,b,d,la,lb,ld:", a, b, d, log_a, log_b, log_d)
    return log_prior


def log_likelihood(theta, x, y):
    a = theta[0]
    b = theta[1]
    d = theta[2]
    #
    # model: y = ax+b+N(0,d) 
    # 
    PI = 3.141592653589793
    #
    len_xy = len(x)
    #
    #
    logll = 0;
    for i in range(len(x)):
        sig_power = d**2.0
        logll_one = -0.5*np.log(2*PI) -0.5*np.log(sig_power) - (x[i]*a+b-y[i])**2.0/2.0/sig_power
        logll = logll + logll_one
    return logll


def log_probability(theta, x, y):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y)


theta = np.random.rand(3)
#why random theta return positive log_p :::::????
#theta[0] = 3
#theta[1] = 5
#theta[2] = 2
print(log_prior(theta),log_likelihood(theta,x,y),log_probability(theta,x,y))




pos = np.random.randn(8, 3)
nwalkers, ndim = pos.shape

N_threads = 8
with Pool(N_threads) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y), pool=pool)
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    start = time.time()
    sampler.run_mcmc(pos, 100000, progress=True);
    end = time.time()
    multi_time = end - start
    print(multi_time)
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))



ncpu = cpu_count()
print("{0} CPUs".format(ncpu))


flat_samples = sampler.get_chain(discard=10000, thin=15, flat=True)
print(flat_samples.shape)
np.savetxt('emcee.out', flat_samples) #np.transpose([a,b,c])
#np.savetxt('emcee.out', flat_samples[:,0,:]) #np.transpose([a,b,c])




