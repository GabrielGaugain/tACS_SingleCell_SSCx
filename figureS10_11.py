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



#%% PLV by cell for one freq and linear regression for each cell to get slopes
###############################################################################
## Getting 2D dataframe of PLV[ cell_name, tacs_amp]
import pandas as pd
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

custom_lines = {"PC":[Line2D([0], [0], color=colors[c], ls=":", lw=2) for c in ["L5", "L6_TPC", "L23"]] ,
                "VIP":[Line2D([0], [0], color=colors["VIP"][i], ls=":", lw=2) for i in [0,5]] ,
                "SST" : [Line2D([0], [0], color=colors["SST"], ls=":", lw=2),] ,
                "PV":[Line2D([0], [0], color=colors["PV"][i], ls=":", lw=2) for i in [0,1,3,-1] ] ,
                }
dleg = {"PC":['L5 PC','L6 PC','L2/3 PC'],   "PV": ['LBC', 'NBC', 'SBC', 'ChC'], 
        "VIP":['BP bNAC','BP cAC',],        "SST":["MC",]}

#%% select one freq

df =pd.read_pickle(f"results/tACS/allfreq_{10}amp.pkl")

metric = "PLV" #"PLV"
freqs = np.arange(5,51,5)

fig, axs = plt.subplots(2,2,figsize=(250*mm, 190*mm), dpi=250, sharex=True, sharey=True)
for i, cl in enumerate(["PC", "PV", "SST", "VIP"]):
    
    cell_list = eval(cl) #np.unique(df["cell_name"])
    cell_names = eval(cl)
    ax = axs[i//2, i%2]        
    
    #### PLOT
    # fig, ax = plt.subplots(figsize=(90*mm, 55*mm), dpi=250)
    for i,cell_name in enumerate(cell_names):
        plv = df[df["cell_name"]==cell_name]["PLV"].values 
        # Lps.append(Lp)
        ax.plot(freqs, plv, color = color_by_cellname(cell_name), label=cell_name,
                marker= markers[i%5], markerfacecolor='none', markersize= 4 ,
                        linestyle = ':', lw=0.8,clip_on=False )

    # for i,c in enumerate(cell_list):
    #     ax.plot(amps, d[c].values/d[c].values[0],  color = color_by_cellname(c),
    #             ls = " ", 
    #             marker = markers[i%5],  ms=3, markerfacecolor='none', markeredgewidth=0.5)#clip_on=False,

    # ax.set_ylim((0.95,1.05))
    ax.set_xlim((0,50))
    ax.spines["bottom"].set_bounds(0,50)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]
    # [axs[1,i].set_xlabel("Electric field (V/m)") for i in range(2)]
    # [axs[i,0].set_ylabel("Relative firing rate") for i in range(2)]
    if cl == "PC":
        ax.legend(custom_lines[cl], dleg[cl], frameon=False,fontsize=12, ncol=2, bbox_to_anchor=(1.1, 1.05))
    else:
        ax.legend(custom_lines[cl], dleg[cl], frameon=False,fontsize=12, ncol=2)


fig.text(0.5, 0.05, "frequency (Hz)",fontsize=14,  ha='center')
fig.text(0.05, 0.5, "PLV",fontsize=14, va='center', rotation='vertical')

plt.subplots_adjust(wspace=0.1, hspace=0.1)

plt.savefig(f"figures/Figures Supp/figureS10_11/PLV_x_freq_amp10.jpg",dpi=250, bbox_inches='tight')
plt.savefig(f"figures/Figures Supp/figureS10_11/PLV_x_freq_amp10.svg")
plt.show()




#%%
fig, axs = plt.subplots(2,2,figsize=(250*mm, 190*mm), dpi=250, sharex=True, sharey=False)
for i, cl in enumerate(["PC", "PV", "SST", "VIP"]):
    
    cell_list = eval(cl) #np.unique(df["cell_name"])
    cell_names = eval(cl)
    ax = axs[i//2, i%2]        
    
    #### PLOT
    # fig, ax = plt.subplots(figsize=(90*mm, 55*mm), dpi=250)
    for i,cell_name in enumerate(cell_names):
        plv = df[df["cell_name"]==cell_name]["PPh"].values 
        # Lps.append(Lp)
        ax.plot(freqs, plv, color = color_by_cellname(cell_name), label=cell_name,
                marker= markers[i%5], markerfacecolor='none', markersize= 4 ,
                        linestyle = ':', lw=0.8,clip_on=False )

    ax.set_xlim((0,50))
    ax.spines["bottom"].set_bounds(0,50)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]

    if cl == "PC":
        ax.legend(custom_lines[cl], dleg[cl], frameon=False,fontsize=12, ncol=2)
        ax.set_yticks([ np.pi/4,np.pi/2, 3*np.pi/4], [r"$\pi/4$",r"$\pi/2$", r"$3\pi/4$"])
    else:
        ax.legend(custom_lines[cl], dleg[cl], frameon=False,fontsize=12, ncol=2)
        ax.set_yticks([-3*np.pi/4,-np.pi/2, -np.pi/4, 0,  np.pi/4], [r"$-3\pi/4$",r"$-\pi/2$",r"$-\pi/4$",r"0",r"$\pi/4$"])

fig.text(0.5, 0.05, "frequency (Hz)",fontsize=14,  ha='center')
fig.text(0.05, 0.5, "PPh",fontsize=14, va='center', rotation='vertical')

# ax[1].spines["left"].set_bounds(ax[1].get_ylim()[0],3*np.pi/4 )
# [ ax[i].spines["bottom"].set_bounds(0,50) for i in range(2)]
# [ ax[0].spines[pos].set_visible(False) for pos in ['right', 'top']]
# [ ax[1].spines[pos].set_visible(False) for pos in ['right', 'top']]

# ax[0].grid(which="major", linestyle='--')
# ax[0].grid(which="minor", linestyle=':')
# plt.minorticks_on()

plt.subplots_adjust(wspace=0.2, hspace=0.1)
plt.savefig(f"figures/Figures Supp/figureS10_11/PPh_x_freq_amp10.jpg",dpi=250, bbox_inches='tight')
plt.savefig(f"figures/Figures Supp/figureS10_11/PPh_x_freq_amp10.svg")
plt.show()
