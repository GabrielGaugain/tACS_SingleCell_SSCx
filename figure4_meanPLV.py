#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from listcells import *
import pickle as pkl
from plot_utils import color_by_cellname, interp, mm
plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025


##### MAIN script #########
amps = np.array([0,1,2,3,4,5,10])
def round_to_half(a):
    return (a*100//5 +1)*5/100
#%% MEAN PLV for all cell types #### sub figure in FIGURE 5

dfs =pd.read_pickle(f"results/allCell_w10Hz.pkl")
freq= 40
df = dfs[dfs["tacs_freq"]==freq]

fig, ax = plt.subplots(figsize=(90*mm, 70*mm), dpi=250)

cell_lists = ["L5","L6_TPC" , "L23" , "SST" , "VIP" , "PV"]
labels = ["L5 PC", "L6 TPC", "L2/3 PC",  "SST" , "VIP" , "PV"]
data= {"amps":amps}
meanPLVs = np.zeros((len(cell_lists),7))
for j, cl in enumerate(cell_lists):
    cell_list = eval(cl)
    # sort_values("tacs_amp")
    PLVs = []
    for c in cell_list:
        PLVs.append(df[df["cell_name"]==c]["PLV"].values)
    PLVs = np.array(PLVs)
    ## for additionnal individual plots
    # for i,c in enumerate(cell_list):
    #     ax.plot(amps, PLVs[i,:], color = color_by_cellname(c), lw = 0.6, ls=':', 
    #             ms = 3., markerfacecolor="none" ,marker = markers[i%5], clip_on=False)
    meanPLVs[j] = PLVs.mean(axis = 0)
    data[labels[j]] = meanPLVs[j]

    ax.plot(amps, meanPLVs[j] , 
                label=labels[j], color = color_by_cellname(cell_list[0]), 
                marker= "x", markersize = 4, markerfacecolor='none',
                linestyle = ':', lw=0.8 , clip_on=False)


ax.set_xlabel("amplitude (V/m)")
ax.set_xlim((0, amps[-1]))
ax.set_ylabel("PLV")
# ax.set_ylim((0,round_to_half(meanPLVs.max())))
ax.set_ylim((0,0.3))
ax.spines["bottom"].set_bounds(0,10)
[ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
plt.legend(frameon = False,ncol=2 ,columnspacing=0.8)
# plt.savefig(f"results/figures/PLVplots/meanPLVs_{freq}Hz_tACS.jpg", dpi=250, bbox_inches='tight')
# plt.savefig(f"results/figures/PLVplots/meanPLVs_{freq}Hz_tACS.svg")
# plt.savefig(f"results/figures/testlegend.svg")
plt.show() 

with open(f"results/meanPLV_{freq}Hz.pkl", "wb") as f:
    pkl.dump(data, f)


#%% MEAN PLV for all cell types #### sub figure in FIGURE 5, all subplot combined

dfs =pd.read_pickle(f"results/allCell_w10Hz.pkl")
freqs = [5,10,20,40]

cell_lists = ["L5","L6_TPC" , "L23" , "SST" , "VIP" , "PV"]
labels = ["L5 PC", "L6 TPC", "L2/3 PC",  "SST" , "VIP" , "PV"]

###
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

fig, axs = plt.subplots(2,2,figsize=(145*mm, 100*mm), dpi=250, sharex=True, sharey=True)
fig.tight_layout()

for i, freq in enumerate(freqs):
    ax = axs[i//2][i%2]
    df = dfs[dfs["tacs_freq"]==freq]

    meanPLVs = np.zeros((len(cell_lists),7))

    for j, cl in enumerate(cell_lists):
        cell_list = eval(cl)
        PLVs = []
        for c in cell_list:
            PLVs.append(df[df["cell_name"]==c]["PLV"].values)
        PLVs = np.array(PLVs)
        ## for additionnal individual plots
        # for i,c in enumerate(cell_list):
        #     ax.plot(amps, PLVs[i,:], color = color_by_cellname(c), lw = 0.6, ls=':', 
        #             ms = 3., markerfacecolor="none" ,marker = markers[i%5], clip_on=False)
            
        meanPLVs[j] = PLVs.mean(axis = 0)
        ax.plot(amps, meanPLVs[j] , 
                    label=labels[j], color = color_by_cellname(cell_list[0]), 
                    marker= "x", markersize = 4, markerfacecolor='none',
                    linestyle = ':', lw=1. , clip_on=False)

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
    ax.set_ylim((0,0.3))
    
    ax.spines["bottom"].set_bounds(0,10)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
    if i==0:
        ax.legend(frameon = False,ncol=2 ,columnspacing=0.8, loc =(0.02,0.5) )#"upper left")

    ax.annotate(f"{freq} Hz",xy=(0.5,1), fontsize=16, xycoords='axes fraction', ha="center")
    # ax.tick_params( length=4, width=1.5,)
    # ax.tick_params(which="minor",  length=2.5, width=1.2,)
# plt.savefig(f"results/figures/PLVplots/meanPLVs_4freqs_tACS.jpg", dpi=250, bbox_inches='tight')
# plt.savefig(f"results/figures/PLVplots/meanPLVs_4freqs_tACS.svg")
# plt.savefig(f"results/figures/testlegend.svg")
plt.show() 




#%% MEAN PLV for all cell types #### sub figure in FIGURE 5, all subplot combined with also simplified aside

dfs =pd.read_pickle(f"results/allCell_w10Hz.pkl")
freqs = [5,10,20,40]

cell_lists = ["L5","L6_TPC" , "L23" , "SST" , "VIP" , "PV"]
labels = ["L5 PC", "L6 TPC", "L2/3 PC",  "SST" , "VIP" , "PV"]

colors = [ color_by_cellname(eval(cl)[0]) for cl in cell_lists]

###
import matplotlib
matplotlib.rcParams.update({'font.size': 10})
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

fig, axs = plt.subplots(4,2,figsize=(145*mm, 200*mm), dpi=250, sharex=True, sharey=True)
fig.tight_layout()

for i, freq in enumerate(freqs):
    ax = axs[i][0]
    df = dfs[dfs["tacs_freq"]==freq]

    meanPLVs = np.zeros((len(cell_lists),7))

    for j, cl in enumerate(cell_lists):
        cell_list = eval(cl)
        PLVs = []
        for c in cell_list:
            PLVs.append(df[df["cell_name"]==c]["PLV"].values)
        PLVs = np.array(PLVs)
        ## for additionnal individual plots
        # for i,c in enumerate(cell_list):
        #     ax.plot(amps, PLVs[i,:], color = color_by_cellname(c), lw = 0.6, ls=':', 
        #             ms = 3., markerfacecolor="none" ,marker = markers[i%5], clip_on=False)
            
        meanPLVs[j] = PLVs.mean(axis = 0)
        ax.plot(amps, meanPLVs[j] , 
                    label=labels[j], color = color_by_cellname(cell_list[0]), 
                    marker= "x", markersize = 4, markerfacecolor='none',
                    linestyle = ':', lw=1. , clip_on=False)

    # ax.set_yticks([0,0.1,0.2,0.3], ["0","0.1", "0.2", "0.3"])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        
    ax.set_xlabel("amplitude (V/m)")
    ax.set_ylabel("PLV")
    
    # ax.set_ylim((0,round_to_half(meanPLVs.max())))
    ax.set_xlim((0, amps[-1]))
    ax.set_ylim((0,0.3))
    
    ax.spines["bottom"].set_bounds(0,10)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
    if i==0:
        ax.legend(frameon = False,ncol=2 ,columnspacing=0.8, loc =(0.02,0.5) )#"upper left")
    yloc = 1 if i==0 else 0.8
    ax.annotate(f"{freq} Hz",xy=(0.5,yloc), fontsize=16, xycoords='axes fraction', ha="center")


### simplified model results

# test = pd.read_csv("results/simplified/Type_3_layer_5_F.csv")
# test.drop(test.columns[0], axis=1, inplace=True)
# plv = test[test["freq"]==freq]["PLV"]
# test = pd.read_pickle("results/simplified/allPLV_simplified.pkl")
test = pd.read_pickle("results/simplified/meanPLV_simplified.pkl")

simp_list = ["L5_PC","L6_PC", "L23_PC", "SST", "VIP", "PV"]


for i, freq in enumerate(freqs):
    ax = axs[i][1]

    temp = test[test["freq"]==freq]
    PLVs = []
    for j,c in enumerate(simp_list):
        PLVs.append(temp[c].values)
        ax.plot(np.arange(0,10.1,1), temp[c].values , 
                    label=labels[j], color = colors[j], 
                    marker= "x", markersize = 4, markerfacecolor='none',
                    linestyle = ':', lw=1. , clip_on=False)

    # ax.set_yticks([0,0.1,0.2,0.3], ["0","0.1", "0.2", "0.3"])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        
    ax.set_xlabel("amplitude (V/m)")
    ax.set_ylabel("PLV")
    
    # ax.set_ylim((0,round_to_half(meanPLVs.max())))
    ax.set_xlim((0, amps[-1]))
    ax.set_ylim((0,0.3))
    
    ax.spines["bottom"].set_bounds(0,10)
    [ ax.spines[pos].set_visible(False) for pos in ['right', 'top']]    
    # if i==0:
    #     ax.legend(frameon = False,ncol=2 ,columnspacing=0.8, loc =(0.02,0.5) )#"upper left")
    yloc = 1 if i==0 else 0.8
    ax.annotate(f"{freq} Hz",xy=(0.5,yloc), fontsize=16, xycoords='axes fraction', ha="center")
    # ax.tick_params( length=4, width=1.5,)
    # ax.tick_params(which="minor",  length=2.5, width=1.2,)
plt.savefig(f"figures/Figure4/meanPLV_detailed_simplified.jpg", dpi=250, bbox_inches='tight')
plt.savefig(f"figures/Figure4/meanPLV_detailed_simplified.svg")
# plt.savefig(f"results/figures/testlegend.svg")
plt.show() 




#%%
type_mapping = {"3":"PC", "10":"SST", "20": "PV", "30": "VIP"}
import os

test = pd.read_csv("results/simplified/Type_3_layer_5_F.csv")
test.drop(test.columns[0], axis=1, inplace=True)
test.drop("PLV", axis=1, inplace=True)
test.drop("Pph", axis=1, inplace=True)
testmean = test.copy()



for typ in [3,10,20,30]:
    celln = type_mapping[f"{typ}"]
    for layer in [1,23,4,5,6]:
        layn = f"L{layer}"
        fn = f"results/simplified/Type_{typ}_layer_{layer}_F.csv"
        if os.path.exists(fn):
            temp = pd.read_csv(fn)
            plv = temp["PLV"]

            if np.all((temp["Amp"] == test["Amp"]) * (temp["freq"] == test["freq"])):
                test[f"{layn}_{celln}"] = plv
        
test.drop(0, inplace=True)
test.to_pickle("results/simplified/allPLV_simplified.pkl")




testmean["L23_PC"] = test["L23_PC"]
testmean["L5_PC"] = test["L5_PC"]
testmean["L6_PC"] = test["L6_PC"]
for celln in ["PV", "SST", "VIP"]:
    col = test[test.columns[np.where([celln in cn for cn in  test.columns])]]
    testmean[celln] = col.mean(axis = 1)


testmean.to_pickle("results/simplified/meanPLV_simplified.pkl")