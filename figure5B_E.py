#####################################################################################
#%% Getting chunks of tACS spikes Figures3B-E
#####################################################################################

import numpy as np
import pandas as pd
from listcells import *
from plot_utils import mm, color_by_cellname, colors, get_continuous_cmap

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


################
### READING DATA
df = pd.read_pickle(f"results/allCell_w10Hz.pkl")
amp = 10
cn, freq, dt = L5[0], 10, 1e3 #
# cn, freq, dt = PV[4], 40 , 0.4e3 

cell = df[(df["cell_name"]==cn) & (df["tacs_freq"]==freq)]
spiketACS = cell[cell["tacs_amp"]==amp]["t_spikes_tacs"].values[0]
spikeSHAM = cell[cell["tacs_amp"]==0]["t_spikes_tacs"].values[0]


###############
### spike Data
times = np.arange(120e3,360e3+dt,dt)
tacs_events, sham_events = [], []
for t in times[1:]:
    loc = (spiketACS > t-dt)*(spiketACS < t)
    tacs_events.append(spiketACS[loc]-t+dt)

    loc = (spikeSHAM > t-dt)*(spikeSHAM < t)
    sham_events.append(spikeSHAM[loc]-t+dt)
    

### COLORMAP CREATION
top = mpl.colormaps['Oranges'].resampled(128)
bottom = mpl.colormaps['Blues_r'].resampled(128)
gray = np.ones((128,4))
gray[:,:3] *=0.1
newcolors = np.vstack( ( np.vstack((1-bottom(np.linspace(0, 0.5, 64)), gray)),
                         1-top(np.linspace(0.5, 1, 64)) ) )
newcmp = ListedColormap(newcolors, name='OrangeGrayBlue')
newcmp = get_continuous_cmap([[0.0, 0.06, 0.9],
                              [0.6,0.6,0.6], 
                              [0.9, 0.06, 0.0]
                              ])

###############
###FIGURE   
figwidth = 90 
figheight = 70*figwidth/110
fig, ax = plt.subplots(3, 1, figsize=(figwidth*mm,figheight*mm),dpi=250,  sharex=True, 
                       gridspec_kw={'height_ratios': [2,3,3]})
fig.subplots_adjust(hspace=0.02)

t = np.arange(0,dt, 0.05)
y = np.sin(2*np.pi*freq*t/1e3)

###############
## Signal tACS
points = np.array([t, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=newcmp, norm =  plt.Normalize(-1, 1), antialiaseds =True)
lc.set_array(y)
lc.set_linewidth(2)
lc.set(rasterized=True)
## TACS waveform
ax[0].add_collection(lc)
ax[0].plot(t,y, path_effects = [patheffects.SimpleLineShadow(offset=(0,0), linewidth=5, alpha=0.1)])
# ax[0].plot(t, y , rasterized=True, color="k", lw=0.8)

###############
## Raster plot
# ax[1].eventplot(tacs_events[40:50], color="orange", lineoffsets=1, linelengths=1)#len(tacs_events[10:20])/200)
# [ax[1].axvline(x =   1e3/(4*freq) + i*1000/freq, color = 'tab:red', alpha=0.2, linewidth=3) for i in range(round(freq * dt/1e3))]
# ax[2].eventplot(sham_events[40:50], color="tab:blue", lineoffsets=1, linelengths=1)#len(tacs_events[10:20])/200)
# [ax[2].axvline(x =   1e3/(4*freq) + i*1000/freq, color = 'tab:red', alpha=0.2, linewidth=3) for i in range(round(freq * dt/1e3))]
# col = color_by_cellname(cn)
ms, alf = 1.6*figwidth/110 , 0.2
for i,e in enumerate(tacs_events):
    loc = np.random.uniform(0, 1, len(e))
    ax[1].plot(e, loc, ls=" ", marker="o", markersize=ms,  
               alpha=alf, rasterized=True, color = "orange" )#color=col,)
    
for i,e in enumerate(sham_events):
    loc = np.random.uniform(0, 1, len(e))
    ax[2].plot(e, loc, ls=" ", marker="o", markersize=ms, 
                alpha=alf, rasterized=True, color = "tab:blue")#color="gray",)
ax[1].margins(x=0, y=0.02)
ax[2].margins(x=0, y=0.02)

#############################################
### Adding usefull legend and text
# ax[0].set_ylabel("E (V/m)", fontsize=13)
tscale = 50
# ax[0].plot([20,20+tscale], [-0.8,-0.8], "k")
# ax[0].text(25+tscale, -0.8, f"{tscale} ms", va="center")

ax[1].text(-4,0.5, "tACS", ha='right', va='center', rotation="vertical")
ax[2].text(-4,0.5, "control", ha='right', va='center', rotation="vertical")
ax[1].set_xticks([0, dt/2, dt], ["0", f"{round(dt/2e3,1)}", f"{dt/1e3}"])

#############################################
#### SET axes visibility and labels/limits
ax[0].spines['left'].set_bounds(-1, 1)
ax[0].spines['left'].set( lw=0.7)
ax[0].tick_params(axis="y",direction="out", width=0.7)

ax[1].spines['bottom'].set_bounds(0, dt)
ax[1].spines['bottom'].set(linestyle=":", lw=0.7)
ax[1].tick_params(axis="x",direction="in", width=0.7)
ax[1].set_xlim((-10,dt+2))

ax[2].spines['bottom'].set_bounds(0, dt)
ax[2].spines['bottom'].set( lw=0.7)
ax[2].tick_params(axis="x",direction="out", width=0.7)

[ax[i].spines['right'].set_visible(False) for i in [0,1,2]]
[ax[i].spines['top'].set_visible(False) for i in [0,1]]
[ax[i].spines['bottom'].set_visible(False) for i in [0]]
[ax[i].spines['left'].set_visible(False) for i in [1,2]]
[ax[i].xaxis.set_visible(False) for i in [0,1]]
[ax[i].yaxis.set_visible(False) for i in [1,2]]
[ax[2].spines[loc].set_visible(False) for loc in ["left", "right", "top"]]

ax[2].set_xlabel("time (s)")
plt.savefig(f"figures/{cn}_raster_{freq}Hz_V3.jpg", bbox_inches='tight', dpi=250)
plt.savefig(f"figures/{cn}_raster_{freq}Hz_V3.svg",bbox_inches='tight', dpi=250)
plt.show()

# %%
