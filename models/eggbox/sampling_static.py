
import numpy as np
import corner
import math

font = {'family' : 'serif', #monospace
'weight' : 'bold', #bold
'size'   : 30,
}


x_dim = 100
y_dim = 100
x_max = 50.0
y_max = 50.0
x_min = 0.0
y_min = 0.0

dx = (x_max-x_min)/x_dim
dy = (y_max-y_min)/y_dim

x_grid = np.linspace(x_min, x_max, num=x_dim, endpoint=False) + dx/2
y_grid = np.linspace(y_min, y_max, num=y_dim, endpoint=False) + dy/2

n_para = 2
ndim = n_para
skipnum = 000 #00000

fname_ch = "%s%s" % ("./", "chain_all.dat")
print(fname_ch)
# read chain 
i_ch_tup = np.genfromtxt(fname_ch, skip_header=skipnum, usecols=(0, 1))
#print(i_ch_tup.shape)
#### change to array
i_ch = np.asarray(i_ch_tup)

m = i_ch.shape[0]
n = i_ch.shape[1]
print(m,n)

x_y_NumInBin = np.zeros((x_dim*y_dim, 4), float)
for i in range(x_dim):
    for j in range(y_dim):
        x_y_NumInBin[i*y_dim+j,0] = x_grid[i]
        x_y_NumInBin[i*y_dim+j,1] = y_grid[j]
        x_y_NumInBin[i*y_dim+j,3] = math.cos(x_grid[i]/2.0)*math.sin(y_grid[j]/2.0) + 1
#np.savetxt("xy0.dat", np.transpose([x_y_NumInBin[:,0],x_y_NumInBin[:,1],x_y_NumInBin[:,2]]))


for i in range(m):
    x_ind = int(i_ch[i,0]/dx)
    y_ind = int(i_ch[i,1]/dy)
    x_y_NumInBin[x_ind*y_dim+y_ind,2] = x_y_NumInBin[x_ind*y_dim+y_ind,2] + 1 

np.savetxt("xy_fit_mod.dat", np.transpose([x_y_NumInBin[:,0],x_y_NumInBin[:,1],x_y_NumInBin[:,2],x_y_NumInBin[:,3]]))
print(x_y_NumInBin[:,2].max())
