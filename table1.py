


#%%
from listcells import *
import pandas as pd
import numpy as np

cl = L5+L6_TPC+L23 + VIP+SST+PV
cell_list = ["L5 PC"]*5 + ["L6 TPC"]*5 + ["L2/3 PC"]*5 + ["VIP"]*10 +["SST"]*5 + ["PV"]*len(PV) 
d= {"cell list":cell_list, "cell name":cl }


df = pd.read_pickle("results/allCell_w10Hz.pkl")
freq10 = []
for c in cl:
    temp = df[(df["cell_name"]==c)*(df["tacs_amp"]==0)]["mean_frequency"].values[0]
    if type(temp) is np.ndarray:
        temp = temp[0]
    freq10.append(round(temp,2))
d["alf frequency"] = freq10

df = pd.read_pickle(f"results/AllCell_natural_10Hz_tACS.pkl")
freqnat = []
for c in cl:
    temp = df[(df["cell_name"]==c)*(df["tacs_amp"]==0)]["mean_frequency"].values[0]
    if type(temp) is np.ndarray:
        temp = temp[0]

    freqnat.append(round(temp,2))
d["natural frequency"] = freqnat

# pd.DataFrame(d).to_csv("results/Tables/tableS1_frompython.csv")
#%%

d["synaptic weight"] = []
for cn in d["cell name"]:
    if cn in alpha10Hz.keys():
        d["synaptic weight"].append(alpha10Hz[cn])
    else:
        d["synaptic weight"].append(np.nan)

