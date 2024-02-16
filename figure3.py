#%%
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import preprocessing
from listcells import *
from plot_utils import mm, color_by_cellname
from matplotlib.lines import Line2D


plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
#%% FIGURE  3A ##############################################

d = pd.read_pickle("results/AllCell_polarisationLength.pkl")
freqs = d["freqs"]

cell_list = [PC, SST, PV,VIP]

fig,axs = plt.subplots(nrows=4, ncols=1, figsize=(90*mm,4*70*mm))

for k,cell_names in enumerate(cell_list):
    
    ax = axs[k]
    for i,cell_name in enumerate(cell_names):

        color = color_by_cellname(cell_name)
        Lp = d[cell_name]
        # Lps.append(Lp)
        ax.plot(freqs, Lp, color = color, label=cell_name)

    ax.set_xlim((0, 50))
    if k==3:
        ax.set_xlabel(r"frequency (Hz)")
    if k!=3:
        ax.set_xticklabels([])

    ax.set_ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
    
    plt.tight_layout()
    ax.set_xlim((0, 50))
    if cell_names ==PC:
        custom_lines = [Line2D([0], [0], color=color_by_cellname(c), ls="-", lw=2) for c in [L5[0], L6_TPC[0], L23[0]]]
        ax.legend(custom_lines, ['L5 PC','L6 PC','L2/3 PC'], frameon=False)
    elif cell_names == VIP:
        custom_lines =[Line2D([0], [0], color=color_by_cellname(VIP[i]), ls="-", lw=2) for i in [0,5]]
        ax.legend(custom_lines, ['BP bNAC','BP cAC',], frameon=False)
    elif cell_names == SST:
        custom_lines =[Line2D([0], [0], color=color_by_cellname(SST[0]), ls="-", lw=2),] 
        ax.legend(custom_lines, ['MC',], frameon=False)
    elif cell_names ==PV:
        custom_lines =[Line2D([0], [0], color=color_by_cellname(PV[i]), ls="-", lw=2) for i in [0,1,3,-1] ]
        ax.legend(custom_lines, ['LBC', 'NBC', 'SBC', 'ChC'], frameon=False, fontsize = 9)
  
# plt.title("L5 PC polarisation lengths")
# plt.savefig(f"results/figures/polarLength_{cell_list}_V2.jpg", dpi = 250)
# plt.savefig(f"results/figures/polarLength_{cell_list}_V2.svg")
plt.savefig(f"results/figures/polarLength_ALL_V1.jpg", dpi = 250)
plt.savefig(f"results/figures/polarLength_ALL_V1.svg")
plt.show()    


#%% FIGURE  3B ##############################################
x = pd.read_pickle('results/AllCell_polarisationLength.pkl')
df = pd.DataFrame(x)
keys = ['L5_TTPC2_cADpyr232_1', 'L23_PC_cADpyr229_2','L6_TPC_L4_cADpyr231_3','L4_BP_bNAC219_2','L4_MC_cACint209_4','L4_SBC_cACint209_1']
df1 = df[keys].copy()
l23=df1.iloc[:,1].values
l23=preprocessing.normalize([l23],norm='max')
l5=preprocessing.normalize([df1.iloc[:,0].values],norm='max')
l6=preprocessing.normalize([df1.iloc[:,2].values],norm='max')
vip=preprocessing.normalize([df1.iloc[:,3].values],norm='max')
sst=preprocessing.normalize([df1.iloc[:,4].values],norm='max')
pv=preprocessing.normalize([df1.iloc[:,5].values],norm='max')

freqs = df["freqs"].values


df_lam=pd.DataFrame({'L5 PC_s': 0.137*l5[0,:],'L5 PC_d': 0.071*l5[0,:],'L6 PC_s': 0.117*l6[0,:],'L6 PC_d': 0.046*l6[0,:],'L23 PC_s': 0.032*l23[0,:],'L23 PC_d': 0.120*l23[0,:], })
df_VIP=pd.DataFrame({'L23 VIP': 0.2*vip[0,:],'L4 VIP': 0.025*vip[0,:], 'L5 VIP': 0.085*vip[0,:],'L6 VIP': 0.1*vip[0,:]})
df_SST=pd.DataFrame({'L23 SST': 0.08*sst[0,:],'L4 SST': 0.07*sst[0,:], 'L5 SST': 0.08*sst[0,:],'L6 SST': 0.06*sst[0,:]})
df_PV=pd.DataFrame({'L23 PV': 0.09*pv[0,:],'L6 PV': 0.07*pv[0,:]})

pola = [[df]]

from plot_utils import mm
# fig, axes = plt.subplots(4,1, sharex=True, figsize=(155*mm,250*mm))

fig,axs = plt.subplots(nrows=4, ncols=1, figsize=(90*mm,4*70*mm))#, sharex=True)
cell_list = ["PC", "SST", "PV", "VIP"]

lam_list = [df_lam, df_SST, df_PV, df_VIP]
cmaps = ["BuPu",'YlGn','Reds','autumn_r']

for k,cell_names in enumerate(cell_list):
    lam = lam_list[k]
    ax = axs[k]
    for i in range(lam.shape[1]):
        ls="-"

        l = lam.iloc[:,i]
        if "PC" in l.name:
            if "L23" in l.name:
                color = color_by_cellname(L23[0]) 
            elif "L5" in l.name:
                color = color_by_cellname(L5[0]) 
            else:
                color = color_by_cellname(L6_TPC[0]) 
            if "d" in l.name:
                ls=":"
            else:
                ls="-"
        else:
            color = mpl.colormaps[cmaps[k]](0.5 + 0.5/lam.shape[1]*(i+1))  
            if "SST" in l.name:
                    color = mpl.colormaps[cmaps[k]](0.5 + 0.5/lam.shape[1]*i)  
                    if "L5" in l.name:
                        ls = (0, (3, 5, 1, 5))
                        ls=":"
            if "VIP" in l.name:
                    color = mpl.colormaps[cmaps[k]](0.5/lam.shape[1]*(i+1))  


        ax.plot(freqs,l.values, label = l.name, color = color,ls=ls)

    ax.legend(frameon=False)

    ax.set_xlim((0.1,50))
    if k!=3:
        ax.set_xticklabels([])
    if k==3:
        ax.set_xlabel(r"frequency (Hz)")
    ax.set_ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
    plt.tight_layout()

  
plt.savefig(f"results/figures/polarLengthSimp_ALL_V1.jpg", dpi = 250)
plt.savefig(f"results/figures/polarLengthSimp_ALL_V1.svg")
plt.show()    

