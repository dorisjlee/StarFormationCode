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
#os.chdir("../../project/FLASH4.3_3/object/")
#os.chdir("../../project/FLASH4.3_3/lev8sink/")
#os.chdir("../../project/FLASH4.3_3/lev8nosink/")
#os.chdir("../../project/FLASH4.3_3/low_res_no_sink/")
#os.chdir("../../project/FLASH4.3_2/object/fat{}/".format(fat_fname))
os.chdir("../../project/FLASH4.3_3/hi_res_beta_0.01sink_CFL0.3_Jeans/")
G = 6.67e-8 #cgs
a = 28730.5 #cm/s
timestep= 30
ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
dim = 2**(lev+3)
cell_size = int((ds.domain_width/dim)[0].in_cgs())
#cell_size = 15625000000000000
#cell_size = 39062500000000000
#cell_size = 3906250000000000
print "cell size: ", cell_size
ctr =dim/2
dr=cell_size
boxlen=int(ds.domain_width[0])
ratio = boxlen/dim
xi_range = np.logspace(-1.5,1.5,num=30)
#xi_range = np.logspace(0,1.04,num=20)
r_range = xi_range/1.057E-17
dv =  dr**3 #volume element 
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
	confident_blockcount_lst.append(len(ix))
        #val =np.sum(r[ix][iy][iz]*dens_arr[ix][iy][iz]*dr)
	#sum_args_list.append(np.sum(r[ix,iy,iz]*dens_arr[ix,iy,iz]*dr))
	sum_args_list.append(np.sum(dens_arr[ix,iy,iz]*dv))
    sum_args_list = np.array(sum_args_list)
    print "confidence_blockcount_list: ", confident_blockcount_lst
    np.savetxt("compare_new_fast_sum_args_list{0}_lev{1}.txt".format(timestep,lev),sum_args_list)
    print "sum_args_list: ",sum_args_list
    plt.loglog(xi_range,4*np.pi*G*sum_args_list/a,label= "t={}".format(timestep))
    plt.savefig('mass{0}_lev{1}.png'.format(timestep,lev))

plt.figure()
#tlst = [0,10,20,30,40,50]
#tlst = [0,28,30,32,40]
#tlst = [34,36,38]
#tlst = [0,28,30,32,40]
#tlst = [0,28,30,32,40]
#tlst = [3034]
#for t in tlst :
END_TIME = 348
for t in np.arange(210,END_TIME):
    if t%30 ==0:
	print "time: ", t
    	plot_MR(t)
plt.legend(loc='upper left')
plt.savefig("compare_fast_MRplot.png")
