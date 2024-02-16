#%%
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import pickle as pkl
from astropy.stats import rayleightest
import pandas as pd
from listcells import *
from BBPresults import load_postres, calc_PLV, PPC
from plot_utils import color_by_cellname, interp, mm
import matplotlib as mpl
from matplotlib.collections import LineCollection
plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025


##### MAIN script #########
amps = np.array([0,1,2,3,4,5,10])

#%% MEAN PLV for all cell types #### sub figure in FIGURE 5

df = pd.read_pickle(f"results/AllCell_natural.pkl")
freq= 10

fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)

cell_lists = ["L5","L6_TPC" , "L23" , "SST" , "VIP" , "PV"]
# cell_lists = ["L23"]
meanPLVs = np.zeros((len(cell_lists),7))
for j, cl in enumerate(cell_lists):
    cell_list = eval(cl)
    # sort_values("tacs_amp")
    PLVs = []
    for c in cell_list:
        PLVs.append(df[df["cell_name"]==c]["PLV"].fillna(0).values)
    PLVs = np.array(PLVs)
    ## for additionnal individual plots
    # for i,c in enumerate(cell_list):
    #     ax.plot(amps, PLVs[i,:], color = color_by_cellname(c), lw = 0.6, ls=':', 
    #             ms = 3., markerfacecolor="none" ,marker = markers[i%5], clip_on=False)
        
    meanPLVs[j] = np.nanmean(PLVs, axis = 0)
    ax.plot(amps, meanPLVs[j] , 
                label=cl, color = color_by_cellname(cell_list[0]), 
                marker= "x", markerfacecolor='none', markersize = 4,
                linestyle = '-', lw=0.8 , clip_on=False)


ax.set_xlabel("amplitude (V/m)")
ax.set_xlim((0, amps[-1]))
ax.set_ylabel("PLV")
# ax.set_ylim((0,0.2))
ax.spines["bottom"].set_bounds(0,10)
[ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
plt.legend(frameon = False)
plt.savefig(f"results/figures/PLVplots/Natural_meanPLVs_{freq}Hz_tACS.jpg", dpi=250, bbox_inches='tight')
plt.savefig(f"results/figures/PLVplots/Natural_meanPLVs_{freq}Hz_tACS.svg")
plt.show() 



#%%

sign = []
for amp in [0,1,2,3,4,5,10]:
    temp = df[df["tacs_amp"]==amp][["cell_name", "pval"]]
    sign.append(temp[temp["pval"]<.05]["cell_name"].values)