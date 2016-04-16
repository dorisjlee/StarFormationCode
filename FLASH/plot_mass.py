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
os.chdir("../../project/FLASH4.3_2/object/fat1000/")
#print "not fat 1000 one" 
#os.chdir("../../project/FLASH4.3_2/object/")
G = 6.67e-8 #cgs
a = 28730.5 #cm/s
timestep= 0 
ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
lev = 5
dim = 2**(lev+3)
#all_data = ds.covering_grid(level=lev, left_edge=[0,0.0,0.0],dims=[dim,dim,dim])
#dens_arr =  np.array(all_data["density"])
margin = dim/3
start  = margin
end = dim-margin
#start = 0
#end = dim
cell_size = int((ds.domain_width/dim)[0].in_cgs())
print "cell size: ", cell_size
print end-start
print "looping through: ", (end-start)**3
boxlen = 1e19
ctr =dim/2
dr=cell_size
xi_range = np.logspace(-2,np.log10(17),num=20)
r_range = xi_range/1.057E-17
	#let dr = cell_size
xyzrange = np.arange(start,end)
def plot_MR(timestep):
	ds = yt.load("sod_hdf5_chk_{}".format(str(timestep).zfill(4)))
	all_data = ds.covering_grid(level=lev, left_edge=[0,0.0,0.0],dims=[dim,dim,dim])
	dens_arr =  np.array(all_data["density"])
	print "Working on t = ", timestep
	sum_args_list = []
	for ri in r_range:
#	    print "Looking at radius: ", ri
	    sum_args = 0
	    for i in xyzrange:
		for j in xyzrange:
		    for k in xyzrange:
	#		print "(i,j,k)",i,j,k
			r = np.sqrt((i-ctr)**2+(j-ctr)**2+(k-ctr)**2)*cell_size
	#		print "(r,ri): ",r,ri
			if np.isclose(r,ri,atol=dr):#atol is +/-
	#                     if r<ri:
	  #                  print "hit r :" , r
			    sum_args+=r*dens_arr[i][j][k]*dr
	    print sum_args
	    sum_args_list.append(sum_args)
	sum_args_list = np.array(sum_args_list)

	np.savetxt("sum_args_list{}.txt".format(timestep),sum_args_list)
	print "sum_args_list: ",sum_args_list
	plt.loglog(xi_range,4*np.pi*G*sum_args_list/a,label= "t={}".format(timestep))
	#plt.loglog(xi_range,4*np.pi*G*sum_args_list/a)
	plt.savefig('mass{}.png'.format(timestep))

plt.figure()
#tlst = [0,10,20,30]
tlst = [22, 24, 26, 28, 30, 31]
for t in tlst :
	plot_MR(t)
plt.legend(loc='upper left')
plt.savefig("MRplot.png")
