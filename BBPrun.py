import os, sys
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
from BBPcell import BBPcell
from listcells import *

from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()

amps = np.array([0,1,2,3,4,5,10])
freq = 10
amp = amps[RANK]
tacs_params = { "tacs_amp": amp, "tacs_freq": freq, 
                "tacs_phi":0, "tacs_dur": 120e3, "tacs_del": 1e3}
Params = {"temperature": 36, "tstop":121e3,"cell_name": L5[0]  , 
        "use_BBPsynapses":True,  "fullAxon":False, "save_res":False,
        "tacs_params": tacs_params}

cell = BBPcell(**Params)
cell.set_tACS()
# cell.plot2D(linewidth=1)#,savefig=True, figname=f"figures/morpho/{cell.cell_name}"  )
# plt.show()

#%%

# alf = 0.25 if RANK == 0 else alpha10Hz[cell.cell_name]
alf=0.24
for isyn in range(cell.synapses.nsyn):
        cell.synapses.synapses[isyn].Fac = 0
        cell.synapses.synapses[isyn].Dep = 0
        if  cell.synapses.df["synapse_type"][isyn] >= 100:
                cell.synapses.weights[isyn]  = alf*cell.synapses.df["weight"][isyn]
        else:
                cell.synapses.weights[isyn]  = 1/alf*cell.synapses.df["weight"][isyn]
cell.synapses.set_weights()

sys.stdout.flush()
cell.initialize(-65)
cell.run()
cell.get_spike_timing()
print(cell.mean_frequency, " Hz activity")

sys.stdout.flush()

# # ## SAVING RESULTS  
fn =f"{cell.resdir}/{freq}Hz_tACS_amp{amp}_withoutPlast.pkl" #if RANK == 0 else f"{cell.resdir}/test_withPlast.pkl"
cell.save_results(fn=fn, recordings=True)
