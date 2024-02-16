#%%
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import os, time,efel
import pickle as pkl
import pandas as pd
from listcells import *
from astropy.stats import rayleightest

df = pd.read_pickle("results/allCell_w10Hz.pkl")

#%% Gett significativity at freq Hz for each amp
# ###############################################################################
freq=40
for amp in [0,1,2,]:
    temp = df[(df["tacs_freq"]==freq)&(df["tacs_amp"]==amp)]
    print(f"significant entrainment at {amp}V/m:\n",temp[temp["pval"]<0.01]["cell_name"].values, "\n\n")
