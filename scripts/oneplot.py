from yt.mods import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import sys 

# Order: 
# fname
# res
# bdy
# size
# slice
# plot
# log
# title (no spaces!)

# data.0001.3d.hdf5 128 per 0.2pc x Density 1 Log(density)

#Give me a file name
#pf=load('data.0.2pc.128per.hdf5')
pf=load(sys.argv[1])

# Give me a resolution
res = int(sys.argv[2])

# Boundary condition
bc = sys.argv[3]

# Size
sz = sys.argv[4]

# Slice
sd = sys.argv[5]
if(sd=='x'): slicedirection=0
elif(sd=='y'): slicedirection=1
else: slicedirection=2

# Plot
plotme = sys.argv[6]

# Log?
logit=0
if( int(sys.argv[7]) != 0 ): logit=1

# Title
plttitle = sys.argv[8]

#Output name
#foutname = '128per.0.2pc.x.dens.pdf'
if(logit): pme = 'log'+ plotme
else: pme = plotme

foutname = str(res) + bc + "." + sz + "." + sd + "." + pme + ".pdf"



if slicedirection==0:
	xdir = 1
	ydir = 2
elif slicedirection==1:
	xdir = 0
	ydir = 2
else:
	xdir=0
	ydir=1
	slicedirection=2
	

#dd = pf.h.all_data()
#dd2 = pf.h.all_data()



##Choose some weird data spot say x-y slice at z=0.75*domain_size
c = 0.0*pf.domain_right_edge[slicedirection]
sl = pf.h.slice(slicedirection,c) ##Get the Slice
w = [pf.domain_left_edge[xdir],pf.domain_right_edge[xdir],pf.domain_left_edge[ydir],pf.domain_right_edge[ydir]] #Choose the width, here I choose the entire domain
frb1 = FixedResolutionBuffer(sl,w,(res,res))  #Create FixedResolution Buffer

if(plotme=='gravitational-potential'):
	frb1[plotme] = frb1[plotme] - frb1[plotme].min()

##THIS IS NOT ANOTHER FRB. IT IS ACTUALLY A NUMPY ARRAY WITH THE DENSITY DIFFERENCES
#frb1[plotme] = frb1[plotme] - frb1[plotme].min()
#frb2[plotme] = frb2[plotme] - frb2[plotme].min()
#frb3 = np.fabs(frb1[plotme]-frb2[plotme])
#frb3 = 2.0*np.fabs(frb1[plotme]-frb2[plotme])/np.fabs(frb1[plotme]+frb2[plotme])

#Try ticks
ticks = np.arange(0,res+1,res/5.0)
rad  = 0.5*(pf.domain_right_edge[ydir]-pf.domain_left_edge[ydir])*pf['pc']
radt = np.linspace(-1,1,num=len(ticks))
radt = radt * (rad)
radt = np.array(radt,dtype='float')
for i in range(len(radt)):
    radt[i] = '%.2f'%(radt[i])
xticks = []
yticks = []

for i in range(len(radt)):
    xticks.append(radt[i])
    yticks.append(radt[len(radt)-i-1])
    
##Plot each FRB and the difference, include a colorbar
fig=plt.figure()
plt.clf()
ax=plt.gca()
if(logit): image = plt.imshow(np.log10(frb1[plotme]))
else: image = plt.imshow(frb1[plotme])
image.set_cmap('Paired')
plt.colorbar(image)
#mpl.colors.Normalize(vmin=-0.5, vmax=2.5)

ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel('pc')
ax.set_ylabel('pc')
ax.set(xticklabels=xticks)
ax.set(yticklabels=yticks)


plt.title(plttitle)
plt.savefig(foutname)

