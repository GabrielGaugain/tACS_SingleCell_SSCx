#%%
# import os
from listcells import *
import subprocess, os
import numpy as np
import pickle as pkl
#%%
## MPI init
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()
# print(SIZE, " cores")

##############################################################
### For tACS script lauching with multiple amp, freq and cells

amps = np.array([0, 1, 2, 3, 4, 5, 10]) # 1 to 5 V/m
# freqs = np.array([15, 25, 30, 35, 45, 50])# freqs
# # Nsim = amps.size * freqs.size
# # freq = freqs[ RANK%freqs.size ]
amp = amps[RANK] #V/m

cell_name = PV[4]
cmd = f"python run_tACS.py -n {cell_name} -amp {amp} -freq {10}"
process = subprocess.run(cmd, shell=True)


# freq = 140
# # for amp in [0,2,5,10]:  
# for cell_name in  SST + PV:
#     cmd = f"python run_tACS.py -n {cell_name} -amp {amp} -freq {freq}"
#     process = subprocess.run(cmd, shell=True)


# if RANK==0:
#     cmd = f"python run_tACS.py -n {cell_name} -amp {0} -freq {10} --use_alf 0"
#     process = subprocess.run(cmd, shell=True)
# elif RANK==1:
#     cmd = f"python run_tACS.py -n {cell_name} -amp {0} -freq {10}"
#     process = subprocess.run(cmd, shell=True)

##################################################
### For polarisationLength script or steady states
# import time
# time.sleep(60*10)
# cell_list = PV[:3] + PV[5:]
# Ncells = len(cell_list)

# for i in range(round(Ncells/SIZE)):
#     i_cell = i*SIZE + RANK
#     if i_cell < len(cell_list):    
#         cell_name = cell_list[i_cell]
#         # cmd = f"python 0_compute_steadyState.py -n {cell_name} --fullAxon 0"
#         cmd = f"mpiexec -n 8 python 3_polarization_length.py -n {cell_name} "
#         process = subprocess.run(cmd, shell=True)
# for cell_name in L23+L6_TPC:
#     cmd = f"mpiexec -n 5 python 3_polarization_length.py -n {cell_name} --df 1 --fmax 100 --fmin 51"
#     process = subprocess.run(cmd, shell=True)
##################################################
## Generate all shape plots
# cell_list = L5 + L23 + L4_LBC + SST + VIP + PV
# cell_list = L6_TPC
# N = len(cell_list)//SIZE +1
# for i in range(N):
#     icell = RANK + i*SIZE 
#     if icell<len(cell_list):
#         cell_name = cell_list[icell]
#         cmd = f"python 5_plot_morpho.py -n {cell_name} "
#         process = subprocess.run(cmd, shell=True)
# for c in L5+L23+L6_TPC + L1_NGC+ L4_LBC + SST + VIP + PV:
#     subprocess.run(f"python 5_plot_morpho.py -n {c} ")



# %%
MPI.Finalize()