#%%
import time, sys, os, efel 
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle as pkl
import numpy as np
import pandas as pd 
from listcells import *
mm = 1/25.4


def plot_Lp_spectrum(cell_names, save_fig = False):

    freq, Lps = [], []
    if type(cell_names)== str : cell_names = [cell_names]

    for cell_name in cell_names:
        try:
        
            with open(f"results/{cell_name}/PolarisationLength/PolarisationLengthSpectrum.pkl", "rb") as f:
                d = pkl.load(f)

            freqs, Lp = d["freqs"], d["Lp"]
            Lps.append(Lp)
            freq.append(freqs)

            plt.plot(freqs, Lp, label=cell_name)

    
        except:
            print(f"Failed to read and plot polarisation length data for {cell_name}")

    plt.xlabel(r"frequency (Hz)")
    plt.ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
    plt.xlim((0, 50))
    if save_fig:
        plt.savefig(f"results/{cell_name}/figures/polarisation_length.jpg", dpi = 250)
        plt.savefig(f"results/{cell_name}/figures/polarisation_length.svg")
    plt.show()     

    return Lps, freq
# Lps, f = plot_Lp_spectrum(VIP)

def results_list2dict(d):
    assert type(d)==list and type(d[0])==dict

    freqs, Lp = d[0]["freqs"], d[0]["Lp"]
    rec = d[0]["rec"]   
    for i in range(1, len(d)):
        freqs = np.concatenate( (freqs, d[i]["freqs"]))
        Lp = np.concatenate( (Lp, d[i]["Lp"]))
        for key in rec.keys():
            rec[key] =  rec[key] + d[i]["rec"][key]

    return {"freqs":freqs, "Lp": Lp, "rec":rec}


def Lp_from_recordings(d, part = "axon"):
    
    rec = d["rec"]
    freqs = d["freqs"]
    Lp = np.zeros(freqs.shape)
    
    for i,freq in enumerate(freqs):
        tstop = 4e3 if freq<5 else 2e3
        tACS_period = 1/freq*1e3 if freq >0 else 0

        t = rec["t"]
        v = rec[part]

        if freq !=0:
            Lp[i] =  np.abs(np.max(v[i][t[i]>tstop-tACS_period]) - 
                            np.min(v[i][t[i]>tstop-tACS_period])
                            )/2
        elif freq ==0:
            Lp[i] = np.abs(v[i][-1] -v[i][0] )
        
        else:
            print(f"Error : freq = {freq} is an invalide value for frequency")

    return Lp





from matplotlib.lines import Line2D
# colors_PV = ["red"] + 2*["xkcd:deep red"] + 3*["xkcd:salmon"] +3*[ "xkcd:tomato"]
rmap = mpl.colormaps["Reds"]
colors = {"L23":"xkcd:carolina blue", "L5":"blue", "L6":"purple", 
          "VIP":["orange"]*5 + ["gold"]*5, 
          "SST":"xkcd:dark pastel green",
          "PV":[rmap(0.99)] + 2*[rmap(0.8)] + 3*[rmap(0.6)] + 3*[rmap(0.3)]}

custom_lines = {"PC":[Line2D([0], [0], color=colors[c], ls=":", lw=2) for c in ["L5", "L6", "L23"]] ,
                "VIP":[Line2D([0], [0], color=colors["VIP"][i], ls=":", lw=2) for i in [0,5]] ,
                "SST" : [Line2D([0], [0], color=colors["SST"], ls=":", lw=2),] ,
                "PV":[Line2D([0], [0], color=colors["PV"][i], ls=":", lw=2) for i in [0,1,3,-1] ] ,
                }

# %% Polarisation length with frequency

d = pd.read_pickle("results/AllCell_polarisationLength.pkl")
freqs = d["freqs"]

# cell_names = [L5[0]]
color = colors["L23"]
cell_names = PC
cell_list = "PC" if cell_names==PC else "VIP" if cell_names==VIP else "SST"if cell_names==SST else PV is cell_names==PV
suffix = "" #fullAxon

fig,axs = plt.subplots(nrows=2, ncols=2, figsize=(150*mm,2*60*mm))
# fig, ax = plt.subplots(figsize=(90*mm,70*mm))
for k,cell_names in enumerate([PC, VIP, SST, PV]):
    ax = axs.flatten()[k]
    for i,cell_name in enumerate(cell_names):
        # with open(f"results/{cell_name}/polarisationLength/PolarisationLengthSpectrum.pkl", "rb") as f:
        #     d = pkl.load(f)
        #     d = results_list2dict(d) if type(d) == list
                
        if cell_name in PC:
            color = colors[cell_name.split("_")[0]]
        else:
            if cell_name in SST:  
                color = colors["SST"]  
            elif cell_name in PV:
                color =colors["PV"][PV.index(cell_name)]
            else:
                color =colors["VIP"][VIP.index(cell_name)]

        Lp = d[cell_name]
        # Lps.append(Lp)
        ax.plot(freqs, Lp, color = color, label=cell_name)
        ax.set_xlim((0,50))

    ax.set_xlabel(r"frequency (Hz)")
    ax.set_ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
    plt.tight_layout()
    plt.xlim((0, 50))

    if cell_names ==PC:
        plt.legend(custom_lines["PC"], ['L5 PC','L6 PC','L2/3 PC'], frameon=False)
    elif cell_names == VIP:
        plt.legend(custom_lines["VIP"], ['BP bNAC','BP cAC',], frameon=False)
    elif cell_names == SST:
        plt.legend(custom_lines["SST"], ['MC',], frameon=False)
    elif cell_names ==PV:
        plt.legend(custom_lines["PV"], ['LBC', 'NBC', 'SBC', 'ChC'], frameon=False, fontsize = 9)
  
# plt.title("L5 PC polarisation lengths")
# plt.savefig(f"results/figures/polarLength_{cell_list}_V2.jpg", dpi = 250)
# plt.savefig(f"results/figures/polarLength_{cell_list}_V2.svg")
plt.show()    




#%% CHECKING TIME COURSE 
cell_name = PV[0]
# with open(f"results/{cell_name}/polarisationLength/PolarisationLengthSpectrum.pkl", "rb") as f:
# with open(f"results/{cell_name}/PolarisationLength.pkl", "rb") as f:
list_d= {}
for cell_name in L5:
    with open(f"results/{cell_name}/PolarisationLength/PolarisationLengthSpectrum.pkl", "rb") as f:
        d = pkl.load(f)

    list_d[cell_name] = d

with open(f"results/L5PC_PolarisationLengthSpectrums.pkl", "wb") as f:
    pkl.dump(list_d, f)


#%% Phase plot of Lp #################################

    plot_sig = False
    phi = []
    # freq = d["freqs"][i]
    # for freq in [10,20,40]:
    freq = 20
    i= np.where((d["freqs"]==freq))[0][0]

    for i,freq in enumerate(d["freqs"]):

        t, v_axon = d["rec"]["t"][i] , d["rec"]["axon"][i]
        tacs = d["Lp"][i]*np.sin(2*np.pi*freq*(t-50)*1e-3)*(t>50)/4 +v_axon[0]
        imax = np.argmax(v_axon[t> t[-1]-1e3/freq]) + np.where(t> t[-1]-1e3/freq)[0][0]
        phi.append( (2*np.pi*freq*(t[imax]-50)*1e-3 )%(2*np.pi) -np.pi/2 )

        if plot_sig:
            # plt.plot(t, v_axon-v_axon[0], label = "vm")
            plt.plot(t, v_axon, label = "vm")
            # tacs= tacs +v_axon[0]
            plt.plot(t, tacs, label="tACS")


            # plt.vlines(t[imax], -74.9, -74.8, "k")

            plt.xlim((t[0], t[-1]))
            plt.xlabel("t (ms)")
            plt.ylabel("Vm (mV)")
            plt.legend()
            plt.xlim(1800,2000)
            plt.show()


    plt.plot(d["freqs"][1:], np.array(phi)[1:])
plt.yticks([-np.pi, -np.pi/2,  0, np.pi/4, np.pi/2, np.pi], [r"$-\pi$",r"$-\pi/2$",r"$0$",r"$\pi/4$",r"$\pi/2$",r"$\pi$",])
plt.ylabel("phase")
plt.xlabel("freq. (Hz)")
plt.show()
# print(phi[0]*180/np.pi, " degree for phase")


#%% Test to extract phase of polarisation (dephasing due to membrane capacitance)
# phi = []
# for i in range(1,len(d["freqs"])):
#     t, v_axon = d["rec"]["t"][i] , d["rec"]["axon"][i]
#     freq = d["freqs"][i]
#     imax = np.argmin(v_axon[t> t[-1]-1e3/freq]) + np.where(t> t[-1]-1e3/freq)[0][0]
#     phi.append( (2*np.pi*freq*(t[imax]-50)*1e-3 -np.pi/2)% (2*np.pi) )
# plt.plot(freqs[1:], 180/np.pi*np.array(phi))
# plt.xlim((0,50))
# plt.show()
# from numpy import fft
# f, sf =fft.fftfreq(v_axon.size, d=t[1]*1e-3 ) , fft.fft(v_axon)/v_axon.size 
# phi, a = np.angle(sf[np.where(f==10)]), np.abs(sf[np.where(f==10)])
# plt.plot(f[1:sf.size//2],np.abs(2*sf[1:sf.size//2]) )
# plt.xlim((0,20))
# plt.show()

# %% plot of pola length with amplitude for 7 freqs

amps = np.array([0.2, 0.5, 1, 2, 3, 5, 10 ])
# cell_names = [L5[0]]
color = colors["L5"]
cell_names = L5
suffix = "" #fullAxon
freq, Lps = [], []
if type(cell_names)== str : cell_names = [cell_names]
markers = ["x", "o", "v", "s", "d","<"]
# cell_name = L5[0]
from scipy.stats import pearsonr, linregress
r= []
for i,cell_name in enumerate(cell_names):

    for j,f in enumerate(np.array([5, 10, 20, 40, 50, 100])):
        fn = f"results/{cell_name}/polarisationLength/PolarisationLength_{float(f)}Hz_multAmp{suffix}.pkl"
        
        try:
            with open(fn, "rb") as f:
                d = pkl.load(f)
            if type(d) == list:
                d = results_list2dict(d)
                with open(fn, "wb") as f:
                    pkl.dump(d, f)

            freqs, Lp = d["freqs"], d["Lp"]
            Lps.append(Lp)
            freq.append(freqs)
            
            plt.plot(amps, Lp, color = color[0], marker=markers[j] , label=cell_name)

            r.append(linregress(amps, Lp ))

        except:
            print(f"Failed to read {fn} and plot polarisation length data for {cell_name}")

    plt.xlabel(r"frequency (Hz)")
    plt.ylabel(r"Polarisation length $\lambda_{p}$  (mm)")
    # plt.xlim((2, 50))

    # from matplotlib.lines import Line2D
    # custom_lines = [Line2D([0], [0], color=color[0], ls=":", lw=2),
    #                 Line2D([0], [0], color=color[1], ls=":", lw=2),
    #                 Line2D([0], [0], color=color[3], ls=":", lw=2),
    #                 Line2D([0], [0], color=color[-1], ls=":", lw=2),]
    # plt.legend(custom_lines, ['BP bNAC', 'BP cAC'], frameon=False)
    # plt.legend(custom_lines, ['MC cAC'], frameon=False)
    # plt.legend(custom_lines, ['L2/3 PC'], frameon=False)
    # plt.legend(custom_lines, ['LBC', 'NBC', 'SBC', 'ChC'], frameon=False)
    # plt.title("L5 PC polarisation lengths")
    # plt.savefig(f"results/figures/polarLength_PV.jpg", dpi = 250)
    # plt.savefig(f"results/figures/polarLength_PV.svg")
    plt.show()    