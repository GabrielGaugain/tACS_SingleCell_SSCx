# %% PLV plot for all freq 
import os
import numpy as np
import matplotlib.pyplot as plt
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

l = []
amp = 10

# cells = L5+L23+L6_TPC+VIP+SST+PV
# for cell in cells:
#     try:
#         df = load_postres(cell, pattern=dict(amp=amp), resdir =f"results/{cell}/post_processing/w10Hz", w10=True)
#         msg = [f"data missing for cell {cell} and freq {freq}" for freq in np.arange(5,51,5) if freq not in df["tacs_freq"].values ]
#         for m in msg: print(m)
        
#         df.drop(df[df["tacs_amp"]==0].index, inplace=True)
#         df.sort_values("tacs_freq", inplace=True, ignore_index = True) ## sort by tacs amp 0->10 V/m
#         df = calc_PLV(df, drop_data=False, PPC=False)
#         if len(msg)==0: l.append(df)
#     except:
#         print(f"didnt work for cell {cell}")
# dfs = pd.concat(l)
# dfs.to_pickle(f"results/allfreq_{amp}amp.pkl")
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

dfs =pd.read_pickle(f"results/allfreq_{amp}amp.pkl")

# #%%
metric = "PLV" #"PLV"

dec = np.arange(-1,1,0.3)
fig, ax = plt.subplots(2,1,figsize=(90*mm, 100*mm), dpi=250, sharex=True)
plt.subplots_adjust(wspace=0., hspace=0.02)
savefig = True #False

cl = ["L5", "L6_TPC", "L23"]#
labels = ["L5 PC", "L6 TPC", "L2/3 PC"]
cln = "PC"

# cl = ["VIP", "SST", "PV"] #
# labels = cl
# cln = "Inhib"
data = {}

for i,cells in enumerate(cl):
    cell_list = eval(cells)
    PLVs, PPhs = [], []
    for c in cell_list:
        PLVs.append(dfs[dfs["cell_name"]==c][metric].values)
        PPhs.append(dfs[dfs["cell_name"]==c]["PPh"].values)
    PLVs, PPhs = np.array(PLVs), np.array(PPhs)

    freqs = np.arange(5,51,5)

    # for i in range(PLVs.shape[0]):
    #     plt.plot(freqs, PLVs[i], ":", lw=1.5)
        # plt.scatter(freqs, PLVs[i])
    # plt.plot(freqs, PLVs.mean(axis=0))
    # plt.errorbar(freqs,PLVs.mean(axis=0), [PLVs.mean(axis=0)- PLVs.min(axis=0),PLVs.max(axis=0)-PLVs.mean(axis=0)],
    #              capsize=3 )
    ax[0].errorbar(freqs+dec[i], PLVs.mean(axis=0), PLVs.std(axis=0)/np.sqrt(PLVs.shape[0]) ,color=color_by_cellname(cell_list[0]),
                capsize=3, marker = "o", markersize=3.,
                lw=.7,   label =  labels[i] )
    ax[1].errorbar(freqs+dec[i], PPhs.mean(axis=0), PPhs.std(axis=0)/np.sqrt(PPhs.shape[0]) ,color=color_by_cellname(cell_list[0]),
                capsize=3, marker = "o", markersize=3.,
                lw=.7,  label = labels[i] )
ax[0].set_ylabel("PLV")
ax[1].set_ylabel("PPh")
ax[1].set_xlabel("frequency (Hz)")
# 
ax[0].legend(frameon= False, ncol=3)

ax[1].set_xlim(0,52)
ax[0].set_ylim(0,0.2)
major, minor = [ 0, 0.1, 0.2,0.3], [0.05,0.15,0.25]
ax[0].set_yticks(major )
ax[0].set_yticks(minor, minor = True )


ax[1].set_yticks([ np.pi/4,np.pi/2, 3*np.pi/4], [r"$\pi/4$",r"$\pi/2$", r"$3\pi/4$"])
# ax[1].set_yticks([-3*np.pi/4,-np.pi/2, -np.pi/4, 0,  np.pi/4], [r"$-3\pi/4$",r"$-\pi/2$",r"$-\pi/4$",r"0",r"$\pi/4$"])
# ax[1].set_yticks([-np.pi/2, -np.pi/4, 0, np.pi/4], [r"$-\pi/2$",r"$-\pi/4$",r"$0$",r"$\pi/4$",])
ax[1].spines["left"].set_bounds(ax[1].get_ylim()[0],3*np.pi/4 )
[ ax[i].spines["bottom"].set_bounds(0,50) for i in range(2)]
[ ax[0].spines[pos].set_visible(False) for pos in ['right', 'top']]
[ ax[1].spines[pos].set_visible(False) for pos in ['right', 'top']]


# ax[0].grid(which="major", linestyle='--')
# ax[0].grid(which="minor", linestyle=':')
# plt.minorticks_on()

if savefig:
    plt.savefig(f"results/figures/spectrumPLV/{cln}_meanSpec_V2.jpg", dpi=250,bbox_inches='tight')
    plt.savefig(f"results/figures/spectrumPLV/{cln}_meanSpec_V2.svg")

plt.show()



#%% FIGURE

fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)

new_f = np.arange(5,50.1,0.5)
mi,ma = interp(freqs, np.min(PLVs, axis=0), new_f, type_interpolation="quadratic"), interp(freqs, np.max(PLVs, axis=0), new_f,type_interpolation="quadratic")
m, sd =interp(freqs, np.mean(PLVs, axis=0), new_f, type_interpolation="quadratic"), interp(freqs, np.std(PLVs, axis=0), new_f, type_interpolation="quadratic") 

for i in range(PLVs.shape[0]):
    ax.plot(new_f, interp(freqs, PLVs[i], new_f, type_interpolation="quadratic") ,
             color = color_by_cellname(cell_list[i]), ls=":", lw=0.8)
ax.plot(new_f, m ,color = color_by_cellname(cell_list[i]))
# ax.fill_between(new_f,mi , ma, alpha=0.5 )
# plt.fill_between(new_f,m - sd , m+sd, alpha=0.5 )
# plt.yticks([np.pi/4, np.pi/2, 3*np.pi/4], [r"$\pi/4$",r"$\pi/2$",r"$3\pi/4$"])
plt.xlim((5,50))
plt.xlabel("frequency (HZ)")
plt.ylabel("PLV ")
if not os.path.exists("figures/spectrumPLV") : os.mkdir("figures/spectrumPLV")
plt.savefig(f"results/figures/spectrumPLV/{cells}_PLVspec.jpg", dpi=250,bbox_inches='tight')
plt.savefig(f"results/figures/spectrumPLV/{cells}_PLVspec.svg")
plt.show()

# # %% PPH

fig, ax = plt.subplots(figsize=(90*mm,70*mm), dpi=250)

for i in range(PPhs.shape[0]):
    ax.plot(new_f, interp(freqs, PPhs[i], new_f, type_interpolation="quadratic"), 
                    color = color_by_cellname(cell_list[i]), ls=":", lw=1)

ax.plot(new_f, interp(freqs, np.mean(PPhs, axis=0), new_f, type_interpolation="quadratic"), color =color_by_cellname(cell_list[i]) )
# plt.fill_between(new_f, interp(freqs[0], np.min(PPh, axis=1), new_f, type_interpolation="quadratic"),
                        # interp(freqs[0], np.max(PPh, axis=1), new_f,type_interpolation="quadratic"), alpha=0.5 )
# plt.yticks([np.pi/4, np.pi/2, 3*np.pi/4], [r"$\pi/4$",r"$\pi/2$",r"$3\pi/4$"])
plt.yticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], [r"$-\pi$",r"$-\pi/2$","0", r"$\pi/2$", r"$\pi$"])
plt.xlim((5,50))
plt.xlabel("frequency (HZ)")
plt.ylabel("PPh (rad)")
plt.savefig(f"results/figures/spectrumPLV/{cells}_PPhspec.jpg", dpi=250,bbox_inches='tight')
plt.savefig(f"results/figures/spectrumPLV/{cells}_PPhspec.svg")
plt.show()
