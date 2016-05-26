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
timestep= 0 
ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
dim = 2**(lev+3)
cell_size = int((ds.domain_width/dim)[0].in_cgs())
print "cell size: ", cell_size
ctr =dim/2
dr=cell_size
boxlen=int(ds.domain_width[0])
ratio = boxlen/dim
#xi_range = np.logspace(0,1.04,num=20)
xi_range = np.logspace(0,1.04,num=20)
r_range = xi_range/1.057E-17
end  = dim/2
#let dr = cell_size
#xyzrange = np.arange(start,end)
def plot_MR(timestep):
	ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
	all_data = ds.covering_grid(level=lev, left_edge=[0,0.0,0.0],dims=[dim,dim,dim])
	dens_arr =  np.array(all_data["density"])
	print "Working on t = ", timestep
	sum_args_list = []
	confident_blockcount_lst = []
	confident_blockcount=0
	for ri in r_range:
#	    print "Looking at radius: ", ri
		start = int(ri/ratio)
		margin = [int(ri/ratio)+70 if int(ri/ratio)+70 <127 else 0][0]
		sum_args = 0
		xyzrange = np.arange(start,end)
		print "ri : ",ri
		print "Looping over: ", len(xyzrange)**3

		for i in xyzrange:
			for j in xyzrange:
				for k in xyzrange:
					r = np.sqrt((i-ctr)**2+(j-ctr)**2+(k-ctr)**2)*cell_size
					if np.isclose(r,ri,atol=dr):#atol is +/-
						sum_args+=r*dens_arr[i][j][k]*dr
						confident_blockcount+=1
	    	print sum_args
	    	print "At radius ",ri*1.057E-17,", number of blocks within dr: ",confident_blockcount
	    	confident_blockcount_lst.append(confident_blockcount)
	    	sum_args_list.append(sum_args)
	sum_args_list = np.array(sum_args_list)
	print "confidence_blockcount_list: ", confident_blockcount_lst
	np.savetxt("sum_args_list{0}_lev{1}.txt".format(timestep,lev),sum_args_list)
	print "sum_args_list: ",sum_args_list
	plt.loglog(xi_range,4*np.pi*G*sum_args_list/a,label= "t={}".format(timestep))
	plt.savefig('mass{0}_lev{1}.png'.format(timestep,lev))

plt.figure()
tlst = [30]
for t in tlst :
	plot_MR(t)
plt.legend(loc='upper left')
plt.savefig("MRplot.png")
