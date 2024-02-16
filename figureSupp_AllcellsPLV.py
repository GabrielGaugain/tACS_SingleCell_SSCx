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


amps = np.array([0,1,2,3,4,5,10])
## Select from 
markers = ["x", "o", "v", "s", "d"]
markers_PV = ["x", "d", "D", "v", "<", ">", "s", "p", "*"]


#%%  PLV for all cell types #### sub figure in FIGURE supp ??, all subplot combined

dfs =pd.read_pickle(f"results/allCell_w10Hz.pkl")
freq= 40
df = dfs[dfs["tacs_freq"]==freq]

cell_lists = ["L5","L6_TPC" , "L23" , "SST" , "VIP" , "PV"]
labels = ["L5 PC", "L6 TPC", "L2/3 PC",  "SST" , "VIP" , "PV"]

###
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

fig, axs = plt.subplots(3,2,figsize=(145*mm, 160*mm), dpi=250, sharex=True, sharey=True)
fig.tight_layout()

for i, cl in enumerate(cell_lists):
    ax = axs[i//2][i%2]
    df = dfs[dfs["tacs_freq"]==freq]

    # meanPLVs = np.zeros((len(cell_lists),7))
    PLVs = np.zeros((len(cl),7))
    cell_list = eval(cl)
    PLVs = []
    for c in cell_list:
        PLVs.append(df[df["cell_name"]==c]["PLV"].values)
    PLVs = np.array(PLVs)
    ## for additionnal individual plots
    for j,c in enumerate(cell_list):
        ax.plot(amps, PLVs[j,:], color = color_by_cellname(c), lw = 0.6, ls=':', 
                ms = 3., markerfacecolor="none" ,marker = markers[j%5], clip_on=False)
        
    # meanPLVs[j] = PLVs.mean(axis = 0)
    # ax.plot(amps, meanPLVs[j] , 
    #             label=labels[j], color = color_by_cellname(cell_list[0]), 
    #             marker= "x", markersize = 4, markerfacecolor='none',
    #             linestyle = ':', lw=1. , clip_on=False)

    # ax.set_yticks([0,0.1,0.2,0.3], ["0","0.1", "0.2", "0.3"])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        
    if i>=2:
        ax.set_xlabel("amplitude (V/m)")
    else:
        [xtlab.set_visible(False) for xtlab in ax.get_xticklabels()]
    if not i%2:
        ax.set_ylabel("PLV")
    else:
        [ytlab.set_visible(False) for ytlab in ax.get_yticklabels()]
    
    # ax.set_ylim((0,round_to_half(meanPLVs.max())))
    ax.set_xlim((0, amps[-1]))
    # ax.set_ylim((0,0.3))
    
    ax.spines["bottom"].set_bounds(0,10)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
    if i==0:
        ax.legend(frameon = False,ncol=2 ,columnspacing=0.8, loc =(0.02,0.5) )#"upper left")

    ax.annotate(labels[i], xy=(0.3,0.8),xycoords='axes fraction',
                 fontsize=13,  ha="center", color=color_by_cellname(cell_list[0]))
    # ax.tick_params( length=4, width=1.5,)
    # ax.tick_params(which="minor",  length=2.5, width=1.2,)
plt.savefig(f"results/figures/PLVplots/allPLVs_{freq}Hz_tACS.jpg", dpi=250, bbox_inches='tight')
plt.savefig(f"results/figures/PLVplots/allPLVs_{freq}Hz_tACS.svg")
# plt.savefig(f"results/figures/testlegend.svg")
plt.show() 
