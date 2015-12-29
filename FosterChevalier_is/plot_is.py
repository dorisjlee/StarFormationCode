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

def plot_time_slice(physical_quantity,timestep,zmin="",zmax="" ,text="",title="",zoom_factor="",velocity=True,grid=False,save=False,log=True,plot_size=5):
    ds= yt.load("output_{0}/info_{0}.txt".format(str(timestep).zfill(5)))
    slc = yt.SlicePlot(ds, "z",physical_quantity)
    slc.set_axes_unit('pc')
    slc.set_figure_size(plot_size)
    slc.set_font_size(plot_size*2.5)
    if (log==False):
	slc.set_log(physical_quantity, False)
    if zoom_factor!="":
        slc.zoom(zoom_factor)
    if zmin!="" and zmax!="":
        slc.set_zlim(physical_quantity, zmin,zmax)
    slc.set_cmap(physical_quantity,"rainbow")
    if title!="":
    	slc.annotate_title(title)
    if text!="":
	slc.annotate_text((0.1, 0.1),text, coord_system='axis')	
    slc.annotate_text((0.05, 0.05),"timestep: {}".format(timestep), coord_system='axis')
    slc.annotate_text((0.05, 0.02),"time: {} Myrs".format(timestep*61793.091/1000000.), coord_system='axis')
    if (velocity):
        slc.annotate_velocity()
    if(grid): 
	slc.annotate_grids()
    if (save):
       name  =str(timestep)#physical_quantity[:3]+str(timestep)
       slc.save(name)
    else:
       slc.show()

def check_IC_profiles(timestep=1):
    from mpl_toolkits.axes_grid1 import AxesGrid
    ds= yt.load("output_{0}/info_{0}.txt".format(str(timestep).zfill(5)))
    fig = plt.figure()
    grid = AxesGrid(fig, ( (0.1, 0.1, 0.8, 0.8)),
                    nrows_ncols = (1, 3),
                    axes_pad = 1.0,
                    label_mode = "1",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")
    fields = ['density','pressure','temperature']
    slc = yt.SlicePlot(ds, 'z', fields)
    slc.set_log('pressure', False)
    slc.set_axes_unit('pc')
    slc.annotate_text((0.05, 0.02),"time: {} Myrs".format(timestep*61793.091/1000000.), coord_system='axis')
    slc.annotate_velocity()
    slc.set_font_size(12)
    for i, field in enumerate(fields):
        plot = slc.plots[field]
        plot.figure = fig
        slc.set_cmap(field,"rainbow")
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        slc._setup_plots()



def density_radial_profile(timestep):
    ds= yt.load("output_{0}/info_{0}.txt".format(str(timestep).zfill(5)))
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
