#%%
import numpy as np
import pandas as pd
import scipy.signal as signal
import matplotlib.pyplot as plt
import pickle as pkl
import matplotlib as mpl
from matplotlib.collections import LineCollection

from astropy.stats import rayleightest
from listcells import *
from BBPresults import load_postres, calc_PLV, PPC
from plot_utils import color_by_cellname, mm
plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025
# tstart = time.time()

markers = ["x", "o", "v", "s", "d"]

#%% PLV by cell for one freq and linear regression for each cell to get slopes
###############################################################################
## Getting 2D dataframe of PLV[ cell_name, tacs_amp]
import pandas as pd
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
amps = [0,1,2,3,4,5,10]

dfs = pd.read_pickle(f"results/allCell_w10Hz.pkl")
dfs.reset_index(inplace=True)

# #%% Calculating FR during tACS
FR = []
for i in range(len(dfs)):
    t = dfs["t_spikes_tacs"][i]
    FR.append( t.size/240 )
dfs["FR_tACS"] = FR


#%% select one freq
freq = 5
cl = "PV"
for freq in [5,10,20,40]:
    df = dfs[dfs["tacs_freq"] == freq]

    fig, axs = plt.subplots(2,2,figsize=(150*mm, 80*mm), dpi=250, sharex=True, sharey=True)
    for i, cl in enumerate(["PC", "PV", "SST", "VIP"]):
        
        cell_list = eval(cl) #np.unique(df["cell_name"])
        ax = axs[i//2, i%2]
        d = {}

        for c in cell_list:
            d[c] = df[df["cell_name"]==c].sort_values("tacs_amp")["FR_tACS"].values
            d = pd.DataFrame(d, columns=cell_list, index=amps)
        
        dFR =  d if i==0 else  pd.concat([d,dFR], axis=1)
            
        
        
        #### PLOT
        # fig, ax = plt.subplots(figsize=(90*mm, 55*mm), dpi=250)
        for i,c in enumerate(cell_list):
            ax.plot(amps, d[c].values/d[c].values[0],  color = color_by_cellname(c),
                    ls = " ", 
                    marker = markers[i%5],  ms=3, markerfacecolor='none', markeredgewidth=0.5)#clip_on=False,

        ax.set_ylim((0.95,1.05))
        ax.set_xlim((0,10.5))
        ax.spines["bottom"].set_bounds(0,10)
        ax.set_xticks([0,5,10])
        ax.set_yticks([0.95,1,1.05])
        [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]
    # [axs[1,i].set_xlabel("Electric field (V/m)") for i in range(2)]
    # [axs[i,0].set_ylabel("Relative firing rate") for i in range(2)]

    fig.text(0.5, 0.0, "Electric field (V/m)", ha='center')
    fig.text(0.02, 0.5, "Relative firing rate", va='center', rotation='vertical')

    plt.subplots_adjust(wspace=0.03, hspace=0.1)

    # plt.savefig(f"results/figures/FR_S1/allCell_{freq}HztACS.jpg",dpi=250, bbox_inches='tight')
    # plt.savefig(f"results/figures/FR_S1/allCell_{freq}HztACS.svg")

    # plt.savefig(f"results/figures/FR_S1/{cl}_{freq}HztACS.jpg",dpi=250, bbox_inches='tight')
    # plt.savefig(f"results/figures/FR_S1/{cl}_{freq}HztACS.svg")
plt.show()




#%%
freq = 5
for freq in [5,10,20,40]:
    df = dfs[dfs["tacs_freq"] == freq]
    d = {}
    cell_list = PC + VIP + SST + PV 
    for c in cell_list:
        d[c] = df[df["cell_name"]==c].sort_values("tacs_amp")["FR_tACS"].values.T
    d = pd.DataFrame(d, columns=cell_list, index=amps)
    d = d.transpose()
    d.to_csv(f"results/Tables/FR_{freq}Hz_tACS.csv")