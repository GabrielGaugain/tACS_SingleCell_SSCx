#%%
import numpy as np
import scipy.signal as signal
import os
import pickle as pkl
import efel

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from plot_utils import get_continuous_cmap

#%%sinusoidal plot with gradient color for extrema
###################################################


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

fig, ax = plt.subplots(figsize=(40,7))
T = 0.5e3
t = np.arange(0,T, 0.05)
y = (t>T/2)*np.sin(2*np.pi*10*(t-T/2)/1e3)
points = np.array([t, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=newcmp, norm =  plt.Normalize(-1, 1), antialiaseds =True)
lc.set_array(y)
lc.set_linewidth(10)
ax.add_collection(lc)
ax.plot(t,y, path_effects = [patheffects.SimpleLineShadow(offset=(0,0), linewidth=20, alpha=0.1)])
# ax.plot(t, (t>T/2)*np.sin(2*np.pi*10*(t-T/2)/1e3), lw=3, cmap = plt.get_cmap("jet"),)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_xlim((0,T))
ax.set_ylim((-2,2))
plt.show()


#%% Reading data 
########################
import pandas as pd
from listcells import *
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_pickle(f"results/allCell_w10Hz.pkl")
cn = L5[1]
freq, amp = 10,10
cell = df[(df["cell_name"]==cn) & (df["tacs_freq"]==freq)]
spiketACS = cell[cell["tacs_amp"]==amp]["t_spikes_tacs"].values[0]
spikeSHAM = cell[cell["tacs_amp"]==0]["t_spikes_tacs"].values[0]

dt = 3. / freq *1e3
times = np.arange(120e3,360e3+dt,dt)
tacs_events, sham_events = [], []
for t in times[1:]:
    loc = (spiketACS > t-dt)*(spiketACS < t)
    tacs_events.append(spiketACS[loc]-t+2*dt )

    loc = (spikeSHAM > t-dt)*(spikeSHAM < t)
    sham_events.append(spikeSHAM[loc]-t+dt)


# #%% Figure 1D
fig, ax = plt.subplots(3,1,sharex=True, gridspec_kw={'height_ratios': [1,2, 1]}, dpi=1000)
fig.subplots_adjust(hspace=-0.)

### Waveform
t = np.arange(0,2*dt, 0.05)
y = np.sin(2*np.pi*10*(t-dt)/1e3)*(t>dt)
points = np.array([t, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=newcmp, norm =  plt.Normalize(-1, 1), antialiaseds =True)
lc.set_array(y)
lc.set_linewidth(2)

## TACS waveform
ax[0].add_collection(lc)
ax[0].plot(t,y, path_effects = [patheffects.SimpleLineShadow(offset=(0,0), linewidth=5, alpha=0.1)])
# ax[0].plot(t, np.sin(2*np.pi*freq*t/1e3) )

## Raster plot ####################################
# ax[1].eventplot(tacs_events, color="k", lineoffsets=0.2, linelengths=len(tacs_events)/250)
ms, alf = 2.5, 0.2
for i,e in enumerate(tacs_events):
    loc = np.random.uniform(0, 1, len(e))
    ax[1].plot(e, loc, ls=" ", marker="o", markersize=ms, color="gray", alpha=alf)
for i,e in enumerate(sham_events):
    loc = np.random.uniform(0, 1, len(e))
    ax[1].plot(e, loc, ls=" ", marker="o", markersize=ms, color="gray", alpha=alf)
ax[1].margins(x=0, y=0.02)

## PTSH / bin count ####################################
dbin = 10
hist_tacs,_,_ = ax[2].hist(np.concatenate(tacs_events), 
                           bins=np.arange(dt,2*dt+0.1,dbin), density=True, color="gray")
hist_sham,_,_ = ax[2].hist(np.concatenate(sham_events), 
                           bins=np.arange(0,dt+0.1,dbin), density=True, color="gray")

# ax[2].plot([0,2*dt],[np.mean(hist_sham)]*2, "k",
#            [0,2*dt],[np.mean(hist_sham)-2*np.std(hist_sham)]*2, "k--",
#            [0,2*dt],[np.mean(hist_sham)+2*np.std(hist_sham)]*2, "k--", lw=0.8,)
ax[2].invert_yaxis()
# ax[2].set_ylabel("PSTH")


#### Adding usefull legend and text
# ax[0].set_ylabel("E (V/m)", fontsize=13)
tscale = 20
ax[0].plot([20,20+tscale], [-0.8,-0.8], "k")

ax[0].text(25+tscale, -0.8, f"{tscale} ms", va="center")

ax[1].text(-20,0.5, "Neural spikes",fontsize=13, ha='center', va='center', rotation="vertical")

ax[2].set_xticks([0, dt, 2*dt],[" ", "", " "])

ax[2].text(dt/2, 1.3*np.max(hist_tacs), "SHAM",fontsize=13, ha='center')
ax[2].text(3*dt/2, 1.3*np.max(hist_tacs) , "tACS",fontsize=13, ha='center')
ax[2].text(-20,np.mean(hist_sham), "bin count",
           fontsize=13,rotation="vertical", ha='center',va='center', )

#### SET axes visibility and labels/limits
ax[0].spines['left'].set_bounds(-1, 1)
ax[0].tick_params(axis="y",direction="out")
ax[2].spines['bottom'].set_bounds(0, 2*dt)
ax[2].tick_params(axis="x",direction="out")
ax[2].set_xlim((-10,2*dt+2))
[ax[i].spines['right'].set_visible(False) for i in [0,1,2]]
[ax[i].spines['top'].set_visible(False) for i in [0,1,2]]
[ax[i].spines['bottom'].set_visible(False) for i in [0,1]]
[ax[i].spines['left'].set_visible(False) for i in [1,2]]
[ax[i].xaxis.set_visible(False) for i in [0,1]]
[ax[i].yaxis.set_visible(False) for i in [1,2]]


plt.savefig("figures/figure1/raster/test_rasterplot.svg")
plt.savefig("figures/figure1/raster/test_rasterplot.jpg", dpi=1000,bbox_inches='tight')
plt.show()

# # #%% try to smooth psth with gaussian kernel
# t = np.arange(0,400.001,0.05)
# def gaussian(x, mean, sigma = 10):
#     return np.exp(-(x - mean)**2/sigma**2 )/sigma

# smt = np.zeros(t.shape)
# sig = 8
# for spike in np.concatenate(tacs_events):
#     smt += gaussian(t, spike,sigma=sig )
# for spike in np.concatenate(sham_events):
#     smt += gaussian(t, spike,sigma=sig )
# plt.plot(t,smt/(np.concatenate(tacs_events).size + np.concatenate(sham_events).size))
# plt.ylabel("PSTH")
# plt.show()