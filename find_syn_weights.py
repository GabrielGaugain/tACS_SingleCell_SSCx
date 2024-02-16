#%%
import time
from neuron import h, nrn
import os, sys
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
import efel

from BBPcell import BBPcell

## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results"
if not os.path.exists(resdir): os.mkdir(resdir)
h.celsius = 36

from listcells import *
from BBPcell import BBPcell


alf_L5 = np.array([[ 0.35, 0.9],
                  [ 0.35, 0.9],
                  [ 0.35, 0.9],
                  [ 1, 1.4],
                  [ 0.35, 0.9]])

# alf_chc = np.array([[ 0.15, 1],
#                   [ 0.15, 1],
#                   [ 0.15, 1]])
# L4_BP_cACint209_1
# "L4_NBC_cNAC187_1"
# L23_PC_cADpyr229_1
# chc = ["L4_ChC_dNAC222_1","L4_ChC_cACint209_1" ,"L4_ChC_cNAC187_1"]
Params = {"temperature": 36, "tstop":10e3,
        "cell_name": L5[0],
        "use_BBPsynapses":True, "syn_freq": None, "syn_events":None,
        "save_res" : False, 
        "tacs_params":None
}
cell = BBPcell(**Params)
for isyn in range(cell.synapses.nsyn):
        cell.synapses.synapses[isyn].Fac = 0
        cell.synapses.synapses[isyn].Dep = 0

def run_step(cell, alf, plot_v = False):

    # ind_exc = [i for i in range(cell.synapses.nsyn) if cell.synapses.synapses_type[i]==0]
    # ind_netstim_exc = np.unique( np.array(cell.synapses.synapses_netstim)[ind_exc])
    # for i in ind_netstim_exc:
    #     cell.synapses.netstim_list[i].interval =  1000/1000
    # cell.synapses.netstim_list[0].interval = cell.synapses.netstim_list[1].interval = \
    # cell.synapses.netstim_list[3].interval = 1000/50
        
    for isyn in range(cell.synapses.nsyn):
        if  cell.synapses.df["synapse_type"][isyn] >= 100:
            cell.synapses.weights[isyn]  = alf*cell.synapses.df["weight"][isyn]
        else:
            cell.synapses.weights[isyn]  = 1/alf*cell.synapses.df["weight"][isyn]
    cell.synapses.set_weights()

    sys.stdout.flush()

    cell.initialize(-65)
    cell.run()
    cell.get_spike_timing()
    freq = cell.mean_frequency

    if freq is None: freq = 0

    if plot_v:
        plt.plot(cell.recordings["t"], cell.recordings["axon[0]"] )
        plt.title(f" alpha ={alpha10Hz} ; FR = {freq}")
        plt.show()

    return float(freq)

#%% INIT the binary search
# eps, iter = 0.5 , 0
# alf =np.array([0.5,2] )
# # alf = alf_L5[1]
# freq = np.zeros((2,))
# freq[0] =  run_step(cell, alf[0]) 
# freq[1] =  run_step(cell, alf[1]) 

eps, iter = 0.5 , 0
alf =np.array([0.25, 0.27,  3] )

# alf = alf_L5[1]
freq = np.zeros((alf.size,))

## MPI init
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()
f=  run_step(cell, alf[RANK], plot_v=False)

COMM.Barrier()
freqs = COMM.gather(f, root=0) 
COMM.Barrier()
if RANK == 0:
    print(freqs, "Hz for alpha = ", alf)

# print(freq)
# while freq[1]< 10-eps or freq[0]>10+eps:

#     if freq[1]< 9.5:
#         freq[0], alf[0] = freq[1], alf[1] #alf_min = alf_max
#         alf[1] *= 10 

#         freq[1] =  run_step(cell, alf[1]) 

#     elif freq[0]>10.5:
#         freq[1], alf[1] = freq[0], alf[0] # alf_max = alf_min
#         alf[0] /=10

#         freq[0] =  run_step(cell, alf[0]) 
# ## END while()

# print("init done; alf = ", alf,"\nfreq = ", freq, "for cell ", cell.cell_name)
#%%

# while np.all(np.abs(10-freq)>0.5 ): 

#     iter +=1
#     w = (10 - freq[0])/(freq[1]-freq[0])
#     alpha = w*alf[1] + (1-w)*alf[0]  

#     f = run_step(cell, alpha) 
#     print(f"\t {iter}th iteration; freq = {f}; alpha = {alpha}")

#     # f = run_step(cell, np.mean(alf)) 
#     # print(f"\t {iter}th iteration; freq = {f}; alpha = {np.mean(alf)}")

#     if f >10+eps:
#         freq[1] = f
#         alf[1] = alpha
#     elif f<10-eps:
#         alf[0] = alpha
#         freq[0] = f

#     else:
#         break

#     if iter>15:
#         break
# ## END while()

# print(f"Took {iter} iteration to find weights to generate {f} Hz activity with alpha = {alpha} for cell {cell.cell_name} \n")
# np.save(f"{cell.cell_name}_10Hz_alpha.npy", np.array(alpha))

MPI.Finalize()

# %%
