
import numpy as np

n_parm = 2
comp_arr = np.zeros((n_parm,2), float)

#filename = "gaussian_prop.chain0"
filename = "chain3.dat"
new_filename = "chain3.d"

file = open(filename,'r')
old_line = file.readline()
file.close()

file_new = open(new_filename,'a')


file = open(filename,'r')
iline = 0
while True:
    #
    print(iline)
    #
    next_line = file.readline()
    if not next_line:
        break;
    #
    all_diff = 0
    #
    str_split_new = next_line.split()
    str_split_old = old_line.split()
    #
    for i in range(n_parm):
        #print(str_split_new[i], str_split_old[i])
        comp_arr[i,0] = str_split_new[i]
        comp_arr[i,1] = str_split_old[i]
    #
    all_diff = np.sum(np.abs(comp_arr[:,0]-comp_arr[:,1]))
    if all_diff > 0:
        #print(comp_arr[:,0])
        #print(next_line.strip())
        file_new.write(next_line)
    #
    old_line = next_line
    iline = iline + 1

file.close()
file_new.close()


