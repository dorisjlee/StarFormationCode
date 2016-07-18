import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
from yt.units import second, g, cm ,dyne
from yt.visualization.fixed_resolution import FixedResolutionBuffer
G = 6.674e-8*cm**3/second**2/g
ctr = 5e18*cm
yt.mylog.setLevel(50)
import numpy as np
import os
# Particle Clean File  
# cp ../source/Simulation/SimulationMain/unitTest/SinkMomTest/utils/clean_flashdat.py . 
# python clean_sinks_evol.py
os.chdir("../object")
SAVE_PATH = "/global/homes/d/dorislee/StarFormationCode/FLASH/"
def plot_dens(i,fname="sod",plane="z", velocity=False,grid=False,zmin ="",zmax="",magnetic=False, particle=False,zoom=""):
    ds = yt.load("{0}_hdf5_chk_{1}".format(fname,str(i).zfill(4)))
    physical_quantity="density"
    slc = yt.SlicePlot(ds, plane,physical_quantity)#,center=(0.5,0.5,0.5))
    slc.set_figure_size(5)
    if zoom!="": slc.zoom(zoom)
    if grid: slc.annotate_grids()
    if velocity: slc.annotate_velocity(normalize=True)
    if magnetic: slc.annotate_magnetic_field()
    slc.set_cmap("all","rainbow")
    if zmin!="" and zmax!="": slc.set_zlim(physical_quantity,zmin,zmax)
    if particle : 
        #os.system("cp ../source/Simulation/SimulationMain/unitTest/SinkMomTest/utils/clean_sinks_evol.py .")
        os.system("python clean_sinks_evol.py")
        data =np.loadtxt("sinks_evol.dat_cleaned",skiprows=1)
        pcl_indx_at_t = np.where(np.isclose(int(ds.current_time.in_cgs()),data[:,1]))[0]
        print "Number of sink particles: " , len(pcl_indx_at_t)
        pcl_pos_at_t = data[pcl_indx_at_t,2:5]
        for pos in pcl_pos_at_t:
            slc.annotate_marker(pos, coord_system='data',marker='.',plot_args={'color':'black','s':3})
    slc.save(SAVE_PATH+"{0}_{1}_zoom_4.png".format(ds,physical_quantity))
def plot_var(i,physical_quantity,fname="sod",cut="z",velocity=False,grid=False,zmin ="",zmax="",particle=False):
    ds = yt.load("{0}_hdf5_chk_{1}".format(fname,str(i).zfill(4)))
    slc = yt.SlicePlot(ds, cut,physical_quantity)#,center=(0.5,0.5,0.5))
    slc.set_figure_size(5)
    if grid: slc.annotate_grids()
    if velocity: slc.annotate_velocity()
    slc.set_cmap("all","rainbow")
    if zmin!="" and zmax!="": slc.set_zlim(physical_quantity,zmin,zmax)
    #slc.show()
    slc.save(SAVE_PATH+"{0}_{1}.png".format(ds,physical_quantity))

def all_direction_slices(i,fname="sod",physical_quantity="density",zmin="",zmax="",zoom=""):
    from mpl_toolkits.axes_grid1 import AxesGrid
    ds = yt.load("{0}_hdf5_chk_{1}".format(fname,str(i).zfill(4)))
    fig = plt.figure()
    grid = AxesGrid(fig, ( (0, 0, 0.8, 0.8)),
                    nrows_ncols = (1, 3),
                    axes_pad = 1.0,
                    label_mode = "1",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")
    direction = ['x','y','z']
    for i, direc in enumerate(direction):
        slc = yt.SlicePlot(ds,direc, physical_quantity)
        slc.set_axes_unit('pc')
        slc.set_font_size(12)
	if zmin!="" and zmax!="": slc.set_zlim(physical_quantity,zmin,zmax)
	if zoom!="": slc.zoom(zoom)
        slc.annotate_magnetic_field()
        plot = slc.plots[physical_quantity]
        plot.figure = fig
        slc.set_cmap(physical_quantity,"rainbow")
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        slc._setup_plots()
    plt.savefig(SAVE_PATH+"{0}_alldir_{1}_zoom_4.png".format(ds,physical_quantity))
START_TIME = 0
END_TIME = 494
#plot_dens(288,velocity=True,particle=True)
#plot_dens(0,velocity=True,particle=True,zmin=5e-21,zmax = 1e-17,zoom=4)
#for i in np.arange(START_TIME,END_TIME,20):
#    plot_dens(i,velocity=True,particle=True,zmin=5e-21,zmax = 1e-17,zoom=4)
#plot_dens(END_TIME,velocity=True,particle=True,zmin=5e-21,zmax = 1e-17,zoom=4)

all_direction_slices(0,zmin=5e-21,zmax = 1e-17,zoom=4)
for i in np.arange(START_TIME,END_TIME,20):
    all_direction_slices(i,zmin=5e-21,zmax = 1e-17,zoom=4)
all_direction_slices(END_TIME,zmin=5e-21,zmax = 1e-17,zoom=4)

##plot_dens(340,grid=True,zoom=4)
