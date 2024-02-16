
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

## Parser for calling the script through cmd line or script to run in loops 
parser = argparse.ArgumentParser()
parser.add_argument("-n","--cell_name", help="name of the cell to simulate", required=True)
parser.add_argument("--fullAxon", help="Whether considering or not the full axonal tree", default=0)
args = parser.parse_args()
cell_name = args.cell_name
fullAxon = bool(int(args.fullAxon))
#%%

Params = {"temperature": 36, "tstop":10, "cell_name": cell_name,
        "use_BBPsynapses":False, "save_res" : False, "fullAxon":fullAxon}


cell = BBPcell(**Params)


fn = "results/cell_bounds.pkl"
if os.path.exists(fn):
    
        with open(fn,"rb") as f:
                d = pkl.load(f)
else:
        d = {}

bounds = np.array(cell._get_boundaries())

bounds = np.vstack((bounds,np.array([cell.cell.soma[0](0.5).x_xtra,cell.cell.soma[0](0.5).y_xtra,cell.cell.soma[0](0.5).z_xtra])))
d[cell.cell_name] = np.array(bounds)

with open(fn,"wb") as f:
    pkl.dump(d, f)
sys.stdout.flush()
# cell.plot2D(savefig=True, figname=f"figures/morpho/{cell.cell_name}" )
