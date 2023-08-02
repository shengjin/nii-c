
import numpy as np
import corner


font = {'family' : 'serif', #monospace
'weight' : 'bold', #bold
'size'   : 30,
}


n_para = 6
ndim = n_para
skipnum = 1000000 #00000

fname_ch = "%s%s" % ("./results082900_debug/", "chain3.dat")
print(fname_ch)
# read chain 
# number of n_para
i_ch_tup = np.genfromtxt(fname_ch, skip_header=skipnum, usecols=(0, 1, 2, 3, 4, 5))
#print(i_ch_tup.shape)
#### change to array
i_ch = np.asarray(i_ch_tup)

m = i_ch.shape[0]
n = i_ch.shape[1]
print(m,n)


figure = corner.corner(i_ch,
        labels = [
    r"$\gamma$",
    r"$log(C)$",
    r"$log(\sigma_{R_{0}})$",
    r"$log(\sigma_{R_{obs}})$",
    r"$log(\sigma_{F_{obs}})$",
    r"${\varphi}$"
            ],
        label_kwargs={"fontsize": 18},
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 13}, #smooth=True
        ) 

figure.savefig("pic_corner.png", dpi=300)


