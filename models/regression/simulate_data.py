##########
import numpy as np

a = 3.5
b = 4
d = 2

x = np.random.rand(1000)*10
y = x*a+b
noise = np.random.randn(1000)*d

np.savetxt("xy.dat", np.transpose([x, y+noise]))


