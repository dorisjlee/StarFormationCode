import plot_is
import os
import yt 
import numpy as np
yt.funcs.mylog.setLevel(50)
os.chdir("/project/projectdirs/astro250/doris/ramses3/trunk/ramses/bin/")
for i in np.arange(1,10):
#    if i%10==0:
#        plot_time_slice("density", i)
    plot_is.plot_time_slice("density", i,save=True)
