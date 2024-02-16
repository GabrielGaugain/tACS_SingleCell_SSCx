import os, sys, argparse
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np

from BBPcell import BBPcell
from listcells import  alpha10Hz


## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results"
if not os.path.exists(resdir): os.mkdir(resdir)


## Parser for calling the script through cmd line or script to run in loops 
parser = argparse.ArgumentParser()

parser.add_argument("-n","--cell_name", help="name of the cell to simulate", required=True)
parser.add_argument("-amp","--tacs_amp", help="Amplitude of the applied tACS", required= True)
parser.add_argument("-freq","--tacs_freq", 
                    help="Frequency of the applied tACS", default=10 ) ## optionnal but can be changed

## optionnal args
parser.add_argument("--tstop", help="duration of the simulation", default=360e3 )
parser.add_argument("-dur","--tacs_dur", help="duration of the applied tACS", default=240e3 )
parser.add_argument("-del","--tacs_del", help="offeset of the applied tACS", default=120e3 )
parser.add_argument("--temperature", help="temperature at which the simulation is runned", default=36 )
parser.add_argument("--use_alf", help="either using synaptic param to generate 10Hz activity or default", default=1 )
parser.add_argument("--overwrite", help="Overwrite existing results if any", default=0 )

args = parser.parse_args()

overwrite = bool(int(args.overwrite))
use_alf = int(args.use_alf)

tacs_params = { "tacs_amp": float(args.tacs_amp), "tacs_freq": int(args.tacs_freq), 
                "tacs_phi":0, "tacs_dur": float(args.tacs_dur), "tacs_del": float(args.tacs_del)}

# ## Main part of the script to run the simulation giver the cell's name and tacs parameters
Params = {"temperature": float(args.temperature), "tstop":float(args.tstop),
        "cell_name": args.cell_name,
        "use_BBPsynapses":True, "syn_freq": None, "syn_events":None, "save_res" : False, 
        "tacs_params":tacs_params
}

cell = BBPcell(**Params)
cell.set_tACS()


# resdir = f"{cell.resdir}/tACS/{cell.tacs_params['tacs_freq']}Hz"
# if not os.path.exists(resdir) : os.mkdir(resdir)   

if use_alf:
    # resdir = f"{cell.resdir}/post_processing/w10Hz"
    resdir = f"{cell.resdir}/tACS/{round(cell.tacs_params['tacs_freq'])}Hz"
    if not os.path.exists(resdir) : os.mkdir(resdir)   

    fn = f"{resdir}/" +\
            f"amp{cell.tacs_params['tacs_amp']}_" +\
            f"freq{cell.tacs_params['tacs_freq']}_" +\
            f"dur{round(cell.tstop/1e3)}s_without_ax_w10Hz.pkl" 
else:
    resdir = f"{cell.resdir}/post_processing/natural"
    if not os.path.exists(resdir) : os.mkdir(resdir)   

    fn = f"{resdir}/" +\
            f"amp{cell.tacs_params['tacs_amp']}_" +\
            f"freq{cell.tacs_params['tacs_freq']}_" +\
            f"dur{round(cell.tstop/1e3)}s_without_ax_natural.pkl"   
print(fn)
print(args.cell_name in alpha10Hz.keys())
##Simulate only if results does not exist
# if (not os.path.exists(fn) )and (not os.path.exists(fn[:-4]+"_postproc.pkl")) or overwrite:

### Seting synaptic weights to enforce 10Hz cellular activity
if args.cell_name in alpha10Hz.keys() and use_alf:
    alpha = alpha10Hz[args.cell_name]

    for isyn in range(cell.synapses.nsyn):
        if  cell.synapses.df["synapse_type"][isyn] >= 100:
            cell.synapses.weights[isyn]  = alpha*cell.synapses.df["weight"][isyn]
        else:
            cell.synapses.weights[isyn]  = 1/alpha*cell.synapses.df["weight"][isyn]
    cell.synapses.set_weights()

elif not use_alf : 
    print("Entered parameter use_alf was 0 => unaffect synaptic weight to generate default/natural activity ")  

else:
    print("use_alpha=", use_alf)
    print("No alpha entered for setting syn weight to generate 10Hz cellular activity and --use_alf =1 .... Exiting ")  
    exit()

sys.stdout.flush()

cell.initialize(-65)
cell.run()
cell.get_spike_timing()
print(cell.mean_frequency, " Hz acitivity")
# with open("test_freq.txt", "a") as f:
#     f.write(f"{cell.cell_name} \t {alpha} \t {cell.mean_frequency} Hz\n")
sys.stdout.flush()
# # ## SAVING RESULTS  
cell.save_results(fn=fn, recordings=True)

# else:
#     print("existing results, aborting simulation ...")