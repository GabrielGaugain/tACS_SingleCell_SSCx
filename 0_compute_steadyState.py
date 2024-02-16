
#%%
import argparse, sys, os
from neuron import h
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from BBPcell import BBPcell
from listcells import *


## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results"
if not os.path.exists(resdir): os.mkdir(resdir)
h.celsius = 36
# # %%
## Parser for calling the script through cmd line or script to run in loops 
# parser = argparse.ArgumentParser()
# parser.add_argument("-n","--cell_name", help="name of the cell to simulate", required=True)
# parser.add_argument("--fullAxon", help="Whether considering or not the full axonal tree", default=0)
# args = parser.parse_args()
# cell_name = args.cell_name
# fullAxon = bool(int(args.fullAxon))

## MPI init
# from mpi4py import MPI
# COMM = MPI.COMM_WORLD
# SIZE = COMM.Get_size()
# RANK = COMM.Get_rank()

# cell_list = L5 + L23+VIP+SST + PV
# Ncells = len(cell_list)
# for i in range(Ncells/SIZE):
#     cell_name = cell_list[i*SIZE + RANK]

# cell_name = L5[0]
fullAxon = False
cell_name = PV[5]#"L4_NBC_dNAC222_1"
Params = {"temperature": 36, "tstop":60e3,
        "cell_name": cell_name,
        "use_BBPsynapses":False, "save_res" : False, "fullAxon":fullAxon}

cell_name = Params["cell_name"]
cell = BBPcell(**Params)

#%%
# cell.recordings = cell.set_recordings(dt = 0.5)
if cell.fullAxon:
    print(f"Computing Steady state for {cell.cell_name} with full axonal tree")
else:
    print(f"Computing Steady state for {cell.cell_name} with cutted axonal tree (2 compartment)")

sys.stdout.flush()
# cell.compute_steady_state(savestate=False)
cell.initialize()
cell.run()

cell.save_state(f"{cell.resdir}/States/steadyState_{cell.cell_name}.pkl")
#%% Testing it
cell.tacs_params = { "tacs_amp": 1, "tacs_freq": 10, "tacs_phi":0,  
                            "tacs_dur": 10e3, "tacs_del": 10}
cell.set_tACS()
cell.initialize()
cell.get_steady_state()
       
cell.tstop = 250
cell.run()

plt.plot(cell.recordings["t"],cell.recordings["axon[0]"])
plt.show()