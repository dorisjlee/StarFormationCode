import matplotlib
matplotlib.use("Agg")
import pylab
import numpy as np
import matplotlib.pyplot as plt
import yt 
yt.funcs.mylog.setLevel(50) #coerce output null
import matplotlib.cm as cm 
from matplotlib.colors import LogNorm
import os
os.chdir("../../project/FLASH4.3_2/object/")
G = 6.67e-8 #cgs
a = 28730.5 #cm/s

i =0 
ds = yt.load("sod_hdf5_chk_{}".format(str(i).zfill(4)))
lev = 5
dim = 2**(lev+3)
all_data = ds.covering_grid(level=lev, left_edge=[0,0.0,0.0],dims=[dim,dim,dim])
dens_arr =  np.array(all_data["density"])

margin = dim/3
start  = margin
end = dim-margin
plt.figure()
plt.imshow(dens_arr[start:end,dim/2,start:end],cmap=cm.jet,norm=LogNorm())
cell_size = int((ds.domain_width/dim)[0].in_cgs())
print "cell size: ", cell_size
print end-start
print "looping through: ", (end-start)**3

# xi_range = np.logspace(-3,5)
xi_range = np.logspace(0,50,num=10)
r_range = xi_range/1.057E-17
#let dr = cell_size 
rcloud =1.59886e18 #xi = 16.90 (/1.057E-17 conv factor) 
xyzrange = np.arange(start,end)
sum_args_list = []
for ri in r_range:
    print "Looking at radius: ", ri
    sum_args = 0
    for i in xyzrange:
        for j in xyzrange: 
            for k in xyzrange: 
                r = np.sqrt(i**2+j**2+k**2)*cell_size
                #print r
		if np.isclose(r,ri,atol=cell_size):
                    print "inside:" , r
                    sum_args+=r**2*dens_arr[i][j][k]*cell_size
    sum_args_list.append(sum_args)
sum_args_list = np.array(sum_args_list)

print sum_args_list
plt.loglog(xi_range,4*pi*G*sum_args_list/rcloud/a)
plt.savefig('mass.png')
