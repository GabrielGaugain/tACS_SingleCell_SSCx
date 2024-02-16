#%%
import argparse, sys, os, time
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




#%%

def results_list2dict(d):
    assert type(d)==list and type(d[0])==dict

    freqs, Lp, amps = d[0]["freqs"], d[0]["Lp"], [d[0]["amp"]]
    rec = d[0]["rec"]   
    for i in range(1, len(d)):
        freqs = np.concatenate( (freqs, d[i]["freqs"]))
        amps += [d[i]["amp"]] 
        Lp = np.concatenate( (Lp, d[i]["Lp"]))
        for key in rec.keys():
            rec[key] =  rec[key] + d[i]["rec"][key]

    amps = np.unique(np.array(amps))
    return {"amp": amps, "freqs":freqs, "Lp": Lp, "rec":rec}


def calc_Lp(Params, freqs = np.arange(0, 50, 1), amp = 1, plot_LP = False):

    ## Instanciate the neuron
    cell_name = Params["cell_name"]
    cell = BBPcell(**Params)
    cell.recordings = cell.set_recordings(dt = 0.5)
    print(cell.fullAxon)
    sys.stdout.flush()
    # cell.compute_steady_state()
    
    t, axon, soma = [],[],[]
    Lp = np.zeros(freqs.shape)

    for i,freq in enumerate(freqs):

        cell.tstop = 4e3 if freq<5 else 2e3
        cell.tacs_params = { "tacs_amp": amp, "tacs_freq": freq, "tacs_phi":0,  
                            "tacs_dur": 10e3, "tacs_del": 50}
        cell.set_tACS()
        tACS_period = 1/freq*1e3 if freq >0 else 0

        # #%%
        ### RUNNING
        cell.initialize()
        cell.get_steady_state()
        sys.stdout.flush()
        cell.run()

        # #%%S
        ti, v_axon, v_soma = np.array(cell.recordings["t"].to_python()), np.array(cell.recordings["axon[0]"].to_python()), np.array(cell.recordings["soma"].to_python())
        t.append(ti)
        axon.append(v_axon)
        soma.append(v_soma)
        if freq !=0:
            Lp[i] =  np.abs(np.max(axon[i][t[i]>cell.tstop-tACS_period]) - 
                            np.min(axon[i][t[i]>cell.tstop-tACS_period])
                            )/2
        elif freq ==0:
            Lp[i] = np.abs(axon[i][-1] -axon[i][0] )
        
        else:
            print(f"Error : freq = {freq} is an invalide value for frequency")
        # plt.plot(t[i], axon_v[i])

        print(f"done calculating polarisation length for freq = {freq}")
        sys.stdout.flush()


    # #%% Post Processing
    if plot_LP:
        plt.plot(freqs, Lp)
        plt.title(cell_name)
        plt.xlabel(r"frequency (Hz)")
        plt.ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
        plt.xlim((freqs[0], freqs[-1]))
        plt.savefig(f"results/{cell_name}/figures/polarisationLength/polarisation_length.jpg", dpi = 250)
        plt.savefig(f"results/{cell_name}/figures/polarisationLength/polarisation_length.svg")
        # plt.show()  

    d = {"amp" :amp  ,"freqs":freqs, "Lp":Lp, "rec":{"t":t, "axon": axon, "soma":soma}}

    # plt.show()
    return Lp, d


def compute_pola_spectrum(Params, freqs, save_Lp = True):


    suffix = "_fullAxon" if Params["fullAxon"] else ""

    if SIZE >1:
        ## scatter frequencies over available RANK | order for easier gathering
        print(f"mpi called, freqs array spread over {SIZE} cores")

        freqbycore = freqs.size//SIZE +1
        if (RANK+1)* freqbycore <= freqs.size:
            sub_freq = freqs[RANK*freqbycore : (RANK+1)*freqbycore]
        else:
            sub_freq = freqs[RANK*freqbycore:]

        # sub_freq = freqs[np.arange(RANK, freqs.size , SIZE )]
        Lp, d = calc_Lp(Params, freqs=sub_freq, plot_LP=False)

        data = COMM.gather(d,root=0)

        if save_Lp and RANK==0:
            print("saving results.")
            data = results_list2dict(data)
            with open(f"results/{args.cell_name}/polarisationLength/PolarisationLength_{freqs.min()}_{freqs.max()}HzSpectrum{suffix}.pkl", "wb") as f:
                pkl.dump(data, f)    
        print(f"Lp calculations finished for {args.cell_name} cell ")


    else:
        ## Looping on freq array
        print("Looping over all frequency (one core compuation)")
        print(save_Lp)
        Lp, d = calc_Lp(Params, freqs=freqs, plot_LP=False)

        if save_Lp:
            with open(f"results/{args.cell_name}/polarisationLength/PolarisationLengthSpectrum_{freqs.min()}_{freqs.max()}Hz{suffix}.pkl", "wb") as f:
                pkl.dump(d, f)

        print(f"Lp calculations finished for {args.cell_name} cell ")    

    return Lp


def compute_pola_amp(Params, amps, freq = 10, save_Lp = True):


    suffix = "_fullAxon" if Params["fullAxon"] else ""
    freqs = np.array([freq])

    if SIZE >1:
        ## scatter frequencies over available RANK | order for easier gathering
        print(f"mpi called, freqs array spread over {SIZE} cores")

        amp = amps[RANK]
        Lp, d = calc_Lp(Params, freqs=freqs, amp=amp)

        data = COMM.gather(d,root=0)

        if save_Lp and RANK==0:
            print("saving results.")
            data = results_list2dict(data)
            with open(f"results/{args.cell_name}/polarisationLength/PolarisationLength_{freq}Hz_multAmp{suffix}.pkl", "wb") as f:
                pkl.dump(data, f)    
        print(f"Lp calculations finished for {args.cell_name} cell ")


    else:
        ## Looping on freq array
        print("Looping over all frequency (one core compuation)")
        print(save_Lp)
        Lp, d  = [], []
        for i, amp in enumerate(amps):
            a, b = calc_Lp(Params, freqs=freqs, amp=amp)
            Lp.append(a)
            d.append(b)

        d = results_list2dict(d)        
        if save_Lp:
            with open(f"results/{args.cell_name}/polarisationLength/PolarisationLength_{freq}Hz_multAmp{suffix}.pkl", "wb") as f:
                pkl.dump(d, f)

        print(f"Lp calculations finished for {args.cell_name} cell ")    

    return Lp

#%%
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()
#%%
# Params = {"temperature": 36, "tstop":10e3,
#         "cell_name": VIP[0],
#         "use_BBPsynapses":False, "syn_freq": None, "syn_events":None, "save_res" : False, 
#         "tacs_params":{ "tacs_amp":1, "tacs_freq": 0, 
#         "tacs_phi":0, "tacs_dur": 10e3, "tacs_del": 50}
# }
# freqs = np.array([0,1,10])
# Lp, d = calc_Lp(Params, freqs=freqs, plot_LP=True, save =True)
# %
# cell_name= L5[0] #VIP[0]
# with open(f"results/{cell_name}/PolarisationLength.pkl", "rb") as f:
#     d = pkl.load(f)
# freqs, Lp = d["freqs"], d["Lp"]
# plt.plot(freqs, Lp)
# # plt.title(cell_name)
# plt.xlabel(r"frequency (Hz)")
# plt.ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
# plt.xlim((freqs[0], freqs[-1]))
# plt.savefig(f"results/{cell_name}/figures/polarisation_length.jpg", dpi = 250)
# plt.savefig(f"results/{cell_name}/figures/polarisation_length.svg")
# plt.show()  
#% CHECKING TIME COURSE
# i = 1
# rec = d["rec"]
# plt.plot(rec["t"][i], rec["axon"][i])
# plt.show()


#%%
## Parser for calling the script through cmd line or script to run in loops 
parser = argparse.ArgumentParser()

parser.add_argument("-n","--cell_name", help="name of the cell to simulate", required=True)
# parser.add_argument("--freq", help="frequency of tACS to compute amplitude dependancy", default=10)
parser.add_argument("--fmax", help="maximum frequency to consider (from 0 to fmax)", default=50)
parser.add_argument("--fmin", help="maximum frequency to consider (from 0 to fmax)", default=0)
parser.add_argument("--df", help="maximum frequency to consider (from 0 to fmax)", default=0.25)
parser.add_argument("--fullAxon", help="Whether considering or not the full axonal tree", default=0)
parser.add_argument("--save", help="SAving the results into a pickle file (default 1 for True)", default=1 )

args = parser.parse_args()
save_Lp = bool(int(args.save))

Params = {"temperature": 36, "tstop":4e3,
        "cell_name": args.cell_name,
        "use_BBPsynapses":False, "syn_freq": None, "syn_events":None, "save_res" : False, "fullAxon":bool(int(args.fullAxon)),
        "tacs_params":{ "tacs_amp":1, "tacs_freq": 10, 
                "tacs_phi":0, "tacs_dur": 4e3, "tacs_del": 50}
}

resdir = f"results/{args.cell_name}/polarisationLength"
if not os.path.exists(resdir): 
    if RANK==0:
        os.mkdir(resdir)
    time.sleep(1)

##################################################
df, fmax, fmin = float(args.df), float(args.fmax), float(args.fmin)
print(df, "Hz sampling ")
freqs= np.arange(fmin, fmax+df, df) # tACS frequency range over which to calculate Lp


# %%
if __name__ == "__main__":

    # amps = np.array([0.2, 0.5, 1, 2, 3, 5, 10 ])
    # compute_pola_amp(Params, amps=amps, freq=freq, save_Lp=save_Lp)
    compute_pola_spectrum(Params, freqs)