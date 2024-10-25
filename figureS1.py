#%%
import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
import pickle as pkl
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

from astropy.stats import rayleightest
from listcells import *
from BBPresults import load_postres, calc_PLV, PPC
from plot_utils import color_by_cellname, mm, colors
plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025
# tstart = time.time()

markers = ["x", "o", "v", "s", "d"]
import neuron


#%% PLV by cell for one freq and linear regression for each cell to get slopes
###############################################################################
## Getting 2D dataframe of PLV[ cell_name, tacs_amp]
from AberraCell import AberraCell

Params = {"temperature": 36, "tstop":50,"cell_id": 15  , 
          "specie_type" : 1, "myelinate_ax": 1,
                "save_res" : False, "rec_axons":True,
                "syn_weight" :0. , "syn_freq" :50,
        "iclamp":{"amp":1.5, "dur":1e3, "del":1e3}
        }


# Params = {"temperature": 36, "tstop": 10, "weight": 0.2  , 
#           "cell_id": 15 , "specie_type": 2, "myelinated": 1,
#           "tacs_amp": 0., "tacs_freq": 10, "tacs_phi": 0, "tacs_dur": None, "tacs_del": None,
#           "clamp_amp": 0.5, "clamp_dur":100, "clamp_del":0,
#           "save_rec": True
#           }

cell = AberraCell(**Params)


# cell.initialize()
# cell.run()
# ax = cell.plot_results(plot_ax=True)
# ax.set_xlim((18,22))
# # cell.plot2D(show_legend =True )
# plt.show()

#%% AP initiation location

savefig = False
cell.tstop = 50

cell.initialize()
spike_timer=0
while (neuron.h.t < cell.tstop-neuron.h.dt/2) and (spike_timer==0): 
    neuron.h.fadvance()
    for sec in cell.cell.axonal:
        if (sec(0.5).v > 0):
            spike_timer = 1
            print(f"AP initiated at time {neuron.h.t} ms")
            AP_site = [sec(0.5).x_xtra,sec(0.5).y_xtra,sec(0.5).z_xtra]
            break


fig,ax = cell.plot2D(show_legend =False, show_scale=False )
ax.plot(AP_site[0], AP_site[1], "ro", fillstyle='none',markersize=10, markeredgewidth=2)

if savefig:
    if cell.myelinate_ax:
        plt.savefig("figures/FiguresSupp/figureS1/APloc_withmyelin.svg")
        plt.savefig("figures/FiguresSupp/figureS1/APloc_withmyelin.png", dpi=200,transparent=True)
    else:
        plt.savefig("figures/FiguresSupp/figureS1/APloc_withoutmyelin.svg")
        plt.savefig("figures/FiguresSupp/figureS1/APloc_withoutmyelin.png", dpi=200,transparent=True)

plt.show()



#%% ## Custom recording in the myelinated tree
recordings = cell.recordings
dt = neuron.h.dt
if cell.myelinate_ax:
    secs = ["axon[0]","Node[0]","Node[2]", "Unmyelin[0]", "Unmyelin[1]", "Unmyelin[2]"]
else:
    secs = ["axon[0]","axon[2]","axon[4]","axon[6]"]

for sec in secs:

    if ("Node" in sec) or ("Unm" in sec):
        recordings[sec] =  neuron.h.Vector().record(eval(f"neuron.h.{sec}(0.5)._ref_v"), dt)
    else:
        recordings[sec] =  neuron.h.Vector().record(eval(f"neuron.h.cell.{sec}(0.5)._ref_v"), dt)
  
cell.tstop = 2e3

cell.initialize()
cell.run()

#%%
savefig =False
distance = neuron.h.distance
distance(sec = cell.cell.soma[0])

from matplotlib import colormaps
cmap = colormaps["viridis"]
from plot_utils import mm

d= []
%matplotlib widget
fig, ax = plt.subplots(figsize=(120*mm,80*mm)) #nrows = 2)
ax.plot(cell.recordings["t"], cell.recordings["soma"], label="soma", color = cmap(0))
for sec in secs:
    if ("Node" in sec) or ("Unm" in sec):
        d.append(distance(eval(f"neuron.h.{sec}(0.5)")))
    else:
        d.append(distance(eval(f"neuron.h.cell.{sec}(0.5)"))) 

d= np.array(d)
# d /=d.max()
dmax = 400
d /=dmax
for i,sec in enumerate(secs):
    ax.plot(cell.recordings["t"], cell.recordings[sec], label=sec, color =cmap(d[i]) )

ax.legend(frameon=False)
ax.set_xlabel("time (ms)")
ax.set_ylabel(r"$V_m$ (mV)")
ax.set_yticks([-80,-40,0,40])
# ax.set_xticks([[-80,-40,0,40]])


# ax.set_xlim((19,22))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=dmax))
plt.colorbar(sm,ax=ax, label = "distance to soma (µm)", ticks=[0,200,400])

# ax.set_xlim((14.5,17))
# ax.set_xlim((19,21))

# if not cell.myelinate_ax:
#     ax.set_xlim((45,47))
# if cell.myelinate_ax and (cell.specie_type==2):
#     ax.set_xlim((19,21))
# if cell.myelinate_ax and (cell.specie_type==1):
#    ax.set_xlim((14.5,16.5))

if savefig:
    if cell.myelinate_ax:    
        plt.savefig(f"figures/PB_Aberra/aberra_propa_{'human' if cell.specie_type==2 else 'rat'}.svg")
        plt.savefig(f"figures/PB_Aberra/aberra_propa_{'human' if cell.specie_type==2 else 'rat'}.jpg", dpi=250)
    else:
        plt.savefig("figures/FiguresSupp/figureS1/BBP_propa.svg")
        plt.savefig("figures/FiguresSupp/figureS1/BBP_propa.jpg", dpi=250)
