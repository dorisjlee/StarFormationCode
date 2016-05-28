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
import sys
lev = int(sys.argv[1])
fat_fname=100
print "mass plots for fat 100"
os.chdir("../../project/FLASH4.3_3/object/")
G = 6.67e-8 #cgs
a = 28730.5 #cm/s
timestep= 30
ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
dim = 2**(lev+3)
cell_size = int((ds.domain_width/dim)[0].in_cgs())
print "cell size: ", cell_size
ctr =dim/2
dr=cell_size
boxlen=int(ds.domain_width[0])
ratio = boxlen/dim
#xi_range = np.logspace(0,1.04,num=20)
xi_range = np.logspace(-1.5,1.5,num=30)
r_range = xi_range/1.057E-17

def plot_MR(timestep):
    ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
    all_data = ds.covering_grid(level=lev, left_edge=[0,0.0,0.0],dims=[dim,dim,dim])
    dens_arr =  np.array(all_data["density"])
    x,y,z =  np.indices((dens_arr.shape))
    r = np.sqrt((x - ctr)**2 + (y -ctr)**2+(z -ctr)**2)*cell_size
    print "Working on t = ", timestep
    sum_args_list = []
    confident_blockcount_lst = []
    for ri in r_range[::-1]:
        ix,iy,iz =  np.where(np.isclose(r,ri,atol=dr))
        print "At radius xi= ",ri*1.057E-17,", number of blocks within dr: ",len(ix)
	confident_blockcount_list.append(len(ix))
        #val =np.sum(r[ix][iy][iz]*dens_arr[ix][iy][iz]*dr)
	sum_args_list.append(np.sum(r[ix,iy,iz]*dens_arr[ix,iy,iz]*dr))
    sum_args_list = np.array(sum_args_list)
    print "confidence_blockcount_list: ", confident_blockcount_lst
    np.savetxt("fast_sum_args_list{0}_lev{1}.txt".format(timestep,lev),sum_args_list)
    print "sum_args_list: ",sum_args_list
    plt.loglog(xi_range,4*np.pi*G*sum_args_list/a,label= "t={}".format(timestep))
    plt.savefig('fast_mass{0}_lev{1}.png'.format(timestep,lev))

plt.figure()
tlst = [0,10,20,30,40,50]
for t in tlst :
    plot_MR(t)
plt.legend(loc='upper left')
plt.savefig("fast_MRplot.png")
