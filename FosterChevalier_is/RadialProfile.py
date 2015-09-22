import matplotlib
matplotlib.use('Agg')

import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

from yt.mods import *

import numpy as np
import sys

# What file are we looking for?
fNum = sys.argv[1]

# In file name
fName = 'data.' + fNum.rjust(4,'0') + '.3d.hdf5'
print("Opening file " + fName)

# Out file name
sName = 'RadialVel_' + fNum.rjust(4,'0') + '.png'
print("Will save to file " + sName)

# Plotting what?
plotMe = 'radial-velocity' #density'

# Max and Min
plotMax = 8
plotMin = 0

# Load file
pf = yt.load(fName)
fig = plt.figure()
time = pf.current_time.convert_to_units('yr')


# Sphere
sp = pf.sphere(pf.domain_center, (0.1,"pc"))
rp = yt.create_profile(sp,'radius','radial_velocity',units={'radius':'pc'},logs={'radius':False})

ax = fig.add_subplot(111)
ax.set_xlim(0,0.1)
ax.set_ylim(plotMin,plotMax)
ax.plot(rp.x.value,rp["radial_velocity"].in_units("km/s").value)
ax.set_xlabel(r"Radius   (pc)")
ax.set_ylabel(r"Radial Velocity    (km/s)")




plt.suptitle("Time = " + str(time))
fig.savefig(sName)


