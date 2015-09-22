import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.display import display
from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *
import yt
yt.funcs.mylog.setLevel(50) #coerce output null

def plot_time_slice(physical_quantity,timestep,text="",title=""):
    ds= yt.load("output_{0}/info_{0}.txt".format(str(timestep).zfill(5)))
    slc = yt.SlicePlot(ds, "z",physical_quantity)
    slc.set_axes_unit('pc')
    slc.set_cmap(physical_quantity,"rainbow")
    slc.set_font_size(20)
    if title!="":
    	slc.annotate_title(title)
    if text!="":
	slc.annotate_text((0.1, 0.1),text, coord_system='axis')	
    slc.annotate_text((0.05, 0.05),"timestep: {}".format(timestep), coord_system='axis')
    slc.annotate_text((0.05, 0.02),"time: {} Myrs".format(timestep*61793.091/1000000.), coord_system='axis')
    slc.annotate_velocity()
    slc.annotate_grids()
    slc.show()

def density_radial_profile(timestep):
    ds= yt.load("output_0000{0}/info_0000{0}.txt".format(timestep))
    c = ds.find_max("density")[1]
    ax = 0 # Cut through x axis
    # cutting through the y0,z0 such that we hit the max density
    ray = ds.ortho_ray(ax, (c[1], c[2]))
    srt = np.argsort(ray['x'])
    plt.figure()
    plt.subplot(211)
    plt.loglog(np.array(ray['x'][srt]), np.array(ray['density'][srt]))
    plt.title("Timestep {}".format(timestep),fontsize=13)
    plt.xlabel("log Radius",fontsize=13)
    plt.ylabel('log Density',fontsize=13)
    # plt.subplot(212)
