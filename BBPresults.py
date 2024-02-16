#%%
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import os, time,efel
import pickle as pkl
import pandas as pd
from listcells import *
from astropy.stats import rayleightest
from scipy.stats import pearsonr, linregress

plt.rcParams.update({        
    "xtick.direction": 'in',
    "ytick.direction": 'in',
})
dt = 0.025
tstart = time.time()

## Another approx from Mardia (1972)
def Rayprob(r, n, approx = "Mardia72"):
    R = n*r
    z = R**2/n
    assert approx in ["Mardia72", "Zar99"], "approx should be <<Mardia72>> or <<zar99>>" 
    if approx == "Mardia72":
        pval = np.exp(-z)*(1 + (2*z-z**2)/(4*n) -(24*z -132*z**2 + 76*z**3-9*z**4)/(288*n**2))
    elif approx == "Zar99":
        pval = np.exp(np.sqrt(1 + 4 * n + 4 * (n ** 2 - R ** 2)) - (1 + 2 * n))
    return pval


### hand made Rayleight test approx
n, r = 928,  0.06
# compute p value using approxation in Zar 1999, p. 617
print(Rayprob(r,n,approx="Zar99"))


def PPC(phi):
    ## calc PPC in addition
    N = len(phi)
    norm_factor = 2/(N*(N-1))
    ppc = 0
    for i in range(N-1):
        # ppc += np.sum( np.cos(phi[i] - phi[i+1:]) )
        ppc += np.sum( np.cos(phi[i])*np.cos(phi[i+1:]) + np.sin(phi[i])*np.sin(phi[i+1:]) )

    return ppc*norm_factor

def listdf2tupledf(dfs):
    new = {}
    for j in range(len(dfs)):
        for i in range(len(dfs[j])):
            temp = dfs[j][i].copy()
            cell_name = temp.pop("cell_name")
            amp = temp.pop("tacs_amp")
            new[(cell_name, amp)] = temp
    return pd.DataFrame(new)

def pre_process_results(cell_name, fn, save_postproc = True):
    if type(fn)==str:
        fn = [fn]
    p ,ds = [], []
    resdir = f"results/{cell_name}"
    print(f"reading results files...")
    for i in range(len(fn)):
        with open(f"{fn[i]}","rb") as f:
            p.append(pkl.load(f))
            ds.append(p[i]["recordings"])


    print("extracting spike timing")
    traces = []
    for i in range(len(ds)):
        trace = {}
        if "axon0" in ds[i].keys():
            ds[i]["axon[0]"] =ds[i]["axon0"] 
        trace["T"] = ds[i]["t"]
        trace["V"] = ds[i]["axon[0]"]
        trace["stim_start"] = [0]
        trace["stim_end"] = [ds[i]["t"][-1]]
        traces.append(trace)

    feats = ["peak_time", "mean_frequency"]
    tr_res = efel.getFeatureValues(traces, feats)

    t_spikes = [t["peak_time"] for t in tr_res]
    print([t["mean_frequency"] for t in tr_res], " Hz")


    fn = [os.path.split(filename)[-1]  for filename in fn]        
    for i, param in enumerate(p):
        param.pop("recordings")
        param["t_spikes"] = t_spikes[i]
        param["mean_frequency"] = tr_res[i]["mean_frequency"]
            
        # Saving the initial files with tspike and without recordings and saving recordings in another file for portability
        if save_postproc:
                with open(f"{resdir}/post_processing/{fn[i][:-4]}_postproc.pkl","wb") as f:
                    pkl.dump(param, f)
    
    return p
## END pre_process_results()

def calc_phase_spikes(p):
    tacs_p = [ i["tacs_params"] for i in p]
    t_spikes = [ param["t_spikes"] for param in p]
    print([ f"{param['mean_frequency']} Hz\n" for param in p ])

    # #%% IF one wants to compare sham to tACS with same phase etc
    for i in range(len(tacs_p)):
        if tacs_p[i] is None:
            tacs_p[i] = {}
            tacs_p[i]["tacs_del"] =tacs_p[i+1]["tacs_del"] 
            tacs_p[i]["tacs_dur"] =tacs_p[i+1]["tacs_dur"] 
            tacs_p[i]["tacs_phi"] =tacs_p[i+1]["tacs_phi"] 
            tacs_p[i]["tacs_freq"] =tacs_p[i+1]["tacs_freq"] 
            tacs_p[i]["tacs_amp"] =0 

    # # %% Phase analysis
    t_spikes_tacs = []
    for i in range(len(t_spikes)):
        itacs =  (t_spikes[i]> tacs_p[i]["tacs_del"])*(t_spikes[i]< tacs_p[i]["tacs_del"] + tacs_p[i]["tacs_dur"])
        t_spikes_tacs.append( t_spikes[i][itacs] )

    twopi = 2*np.pi
    phase_spikes =  []
    for i in range(len(t_spikes_tacs)):
        phase_spikes.append( (twopi*tacs_p[i]["tacs_freq"]* 1e-3*(t_spikes_tacs[i]-tacs_p[i]["tacs_del"]) + tacs_p[i]["tacs_phi"]) % twopi )

    return phase_spikes
## END calc_phase_spikes


def load_postres(cell_name,pattern=dict(freq=10), w10 = True, resdir=""):

    data = []
    if resdir == "":
        resdir = f"results/{cell_name}/post_processing"
    try:
        all_fn =  [filename for filename in os.listdir(resdir) if os.path.isfile(f"{resdir}/{filename}") ]
        # print(all_fn)
        ## select only w10Hz files
        if w10:
            all_fn = [ f for f in all_fn if "w10Hz" in f]
        
        ## select files of the interest frequency
        for k in pattern.keys():
            print(f"{k}{pattern[k]}_")
            fn = [f for f in all_fn if (f"{k}{pattern[k]}" in f)or(f"{k}{float(pattern[k])}" in f)  ]
        # for k in pattern:
        #     fn = [f for f in all_fn if (k in f)  ]

        # fn = [f for f in all_fn if (f"freq{freq}_" in f)or(f"freq{float(freq)}_" in f)  ]
        ## adding control if not file found with amp 0
        if not np.any(["amp0" in fname for fname in fn]) and (np.any(["control" in fname for fname in all_fn]) or np.any(["amp0" in fname for fname in all_fn])):
            fn.append(np.array(all_fn)[np.array(["control" in fname for fname in all_fn] + np.array(["amp0" in fname for fname in all_fn]))][0])

        print(fn)
        ## read files and append data to listl
        for filename in fn:
            # print(filename)
            with open(f"{resdir}/{filename}","rb") as f:
                d = pkl.load(f) 

            ## removing useless infos
            for key in ["use_BBPsynapses", "syn_freq","save_synapses", "syn_events", "save_rec", "save_res", "resdir"]:
                if key in d.keys():
                    d.pop(key)
            for key in d["tacs_params"].keys():
                d[key] = d["tacs_params"][key]
            d.pop("tacs_params")
            data.append(d)

    except:
        print("didnt find any file")
    return pd.DataFrame(data)
    # return data
## END load_postres()


def calc_PLV(df, drop_data= True, PPC=False):

    t_spikes_tacs = []
    for i in range(len(df)):
        if df["t_spikes"][i] is not np.NaN:
            itacs =  (df["t_spikes"][i] > df["tacs_del"][i])*(df["t_spikes"][i] < df["tacs_del"][i] + df["tacs_dur"][i])
            t_spikes_tacs.append( df["t_spikes"][i][itacs] )
        else:
            t_spikes_tacs.append(np.NaN)
    df["t_spikes_tacs"] = t_spikes_tacs

    twopi = 2*np.pi
    phase_spikes = (twopi*df["tacs_freq"]* 1e-3*(df["t_spikes_tacs"] - df["tacs_del"]) + df["tacs_phi"]) % twopi 

    meanVecs, pvals, PPCs = [], [], []
    for i in range(len(df)):
        meanVecs.append(np.mean(np.exp(1j*phase_spikes[i] ) ) )        
        pvals.append(rayleightest(phase_spikes[i])) if not np.all(np.isnan(phase_spikes[i])) else pvals.append(np.NaN)
        if PPC : PPCs.append(PPC(phase_spikes[i]))
    # print(phase_spikes)

    df["meanVec"] = meanVecs
    df["PLV"] = np.abs(meanVecs)
    df["PPh"] = np.angle(meanVecs)
    df["pval"] =  pvals

    N = np.array([ i.size if not np.all(np.isnan(i)) else np.NaN  for i in phase_spikes])
    df["ubPLV"] = np.sqrt(  (df["PLV"]**2 * N -1 )/(N-1))# if N != np.NaN else np.NaN 
    df["Coh"] = [np.abs( np.sum(np.sin(i), axis=-1) /\
                 ( np.sqrt( i.size * np.sum(np.sin(i)**2,  axis=-1) ) ) ) 
                 if not np.all(np.isnan(i)) else np.NaN
                 for i in phase_spikes]## Coherence
    if PPC: df["PPC"] =  PPCs

    if drop_data:
        df.drop(["tsim", "tstop", "t_spikes_tacs", "t_spikes"], axis=1)

    return df
## END calc_PLV


def allcell_toDF(freqs =[5,10,20,40] ,two_levels=False, natural = False, savedf = False):
   
    cell_list = L5  + L23+ L4_LBC  + L1_NGC + L6_TPC + VIP + SST+PV #+ L23[:3]+[L23[4]]
    r, params =[], []
    # freq=40

    for freq in freqs:
        for cell_name in cell_list:
            # print(cell_name)
            if natural:
                df = load_postres(cell_name, pattern=dict(freq=freq), resdir =f"results/{cell_name}/post_processing/natural", w10=False)
            else:
                df = load_postres(cell_name, pattern=dict(freq=freq), resdir =f"results/{cell_name}/post_processing/w10Hz")
                if len(df)== 0:
                    print(cell_name, " do not have file")
                    break
            for amp in [0,1,2,3,4,5,10]:
                if amp not in df["tacs_amp"].values:
                    print(cell_name, f" missed {freq}Hz ", amp, " V/m")
                    nanline = { ind:np.NaN for ind in df.columns  }
                    nanline["cell_name"] = cell_name
                    nanline["tacs_amp"] = amp
                    nanline["tacs_freq"] = freq
                    df.loc[len(df)] = nanline
            df.sort_values("tacs_amp", inplace=True, ignore_index = True) ## sort by tacs amp 0->10 V/m
            df.at[0,"tacs_freq"]= freq
            df = calc_PLV(df, drop_data=False)
            # print(df["mean_frequency"])
            if two_levels:
                df = df.to_dict(orient='records')
            params.append(df)
            if len(df)<7:
                print(cell_name)
            # r.append( np.corrcoef(df["PLV"] ,df["tacs_amp"]) )
            # r.append( pearsonr(df["PLV"] ,df["tacs_amp"]) )
            # r.append(linregress(df["tacs_amp"],df["PLV"] ))

        ########## OLD version to get multiindex dataframe (cell_name, tacs_amp)
        if two_levels and savedf:
            all_dfs = listdf2tupledf(params)
            all_dfs.to_pickle(f"results/allCell_natural_{freq}Hz_tACS.pkl") if natural else all_dfs.to_pickle(f"results/allCell_w10Hz_{freq}Hz_tACS.pkl")
            return all_dfs
    dfs = pd.concat(params)
    if not two_levels and savedf:
        dfs.to_pickle(f"results/AllCell_w10Hz.pkl" ) if not natural else dfs.to_pickle(f"results/AllCell_natural.pkl" )

    return dfs



#################### MAIN playground ########################################################################

### cheat sheet for multiindex dataframe selections
# ###############################################################################
# idx = pd.IndexSlice
# df[L5].loc[:,idx[:,10]]    OR    df[L5].iloc[:, dfs.columns.get_level_values(1)==10]

# # %% Preprocessing of voltages files
# def main():
# t_spikes, phase_spikes, ds, p = [], [], [], []
# cell_name = L23[0]
# print(f"post processing of {cell_name}")
# resdir = f"results/{cell_name}"
# suffix = "_without_ax_natural"
# ### Pre-process results to get spike timing and remove recordings (heavy)
# fn = [f"{resdir}/tACS/amp{amp}_freq10_dur360s{suffix}.pkl" for amp in np.array([0,1,2,3,5,10])]
# # fn = [f"{resdir}/tACS/{filename}" for filename in os.listdir(f"{resdir}/tACS")]
# p = pre_process_results(cell_name, fn, save_postproc = True)


#%% LOADING POST PROCESSING RESULTS (SPIKE EVENTS AND METRICS)
# ###############################################################################
# if __name__ == "__main__":
# df = allcell_toDF(freqs=[5,10,20,40], two_levels=False, natural=False) #, savedf=True)
# df.to_pickle("results/allCell_w10Hz.pkl")

# df = allcell_toDF(freqs=[1], two_levels=False, natural=False) 
# df = allcell_toDF(freqs=[10], two_levels=False, natural=True)
# df.to_pickle(f"results/AllCell_natural_10Hz_tACS.pkl" )

# df = allcell_toDF(freqs=[0], two_levels=False, natural=False)
# df.to_pickle("results/allCell_w10Hz_tDCS.pkl")

#%% importing data from another dataset to combine datasets
# df = pd.read_pickle("results/allCell_w10Hz.pkl")
# amp=2
# dfreq = pd.read_pickle(f"results/allfreq_{amp}amp.pkl")
# freqs = np.arange(5,51,5)
# for c in PC + VIP + SST:
#     for f in freqs:
#         if f not in dfreq[dfreq["cell_name"]==c]["tacs_freq"].values:
#             # print(c, f)
#             if f in df[(df["cell_name"]==c)*(df["tacs_amp"]==amp)]["tacs_freq"].values :
#                 dfreq = pd.concat( (dfreq,df[(df["cell_name"]==c)*(df["tacs_amp"]==amp)*(df["tacs_freq"]==f)])  )              

# dfreq.to_pickle(f"results/allfreq_{amp}amp_V2.pkl")

#%% Loading pickled data => which cell are sign entrained at each amp
# ###############################################################################
# freq = 40
# dfs = pd.read_pickle(f"results/allCell_w10Hz.pkl")
# df = dfs[dfs["tacs_freq"]==freq]
# sign = []
# for amp in [0,1,2,3,4,5,10]:
#     temp = df[df["tacs_amp"]==amp][["cell_name", "pval"]]
#     sign.append(temp[temp["pval"]<.05]["cell_name"].values)

# # plt.plot([0,1,2,3,4,5,10],df[df["cell_name"]=="L4_NBC_dNAC222_1"]["PLV"].values)
# # plt.show()

# # %% Perform permutation non-parametric test to test if PLV_tacs & PLV_sham distribution differs ( stat diff, stat higher) 
# ###############################################################################
# from scipy.stats import permutation_test
# def phi(df, twocol=True):
    
#     if twocol:
#         return (2*np.pi*df.loc["tacs_freq"] * 1E-3* (df.loc["t_spikes_tacs"] -df.loc["tacs_del"]))%(2*np.pi)
#     else:
#         return (2*np.pi*df["tacs_freq"].values * 1E-3* (df["t_spikes_tacs"].values -df["tacs_del"].values))%(2*np.pi)

# def dPLV(a,b, axis):
#     return np.abs(np.mean(np.exp(1j*a), axis=axis)) -np.abs(np.mean(np.exp(1j*b), axis=axis)) 

# # df = pd.read_pickle(f"results/allCell_w10Hz.pkl")
# df = pd.read_pickle(f"results/allCell_w10Hz_{freq}Hz_tACS.pkl")
# cn = "L4_SBC_cACint209_1"#"L4_NBC_dNAC222_1"
# # temp = df[df["cell_name"]==cn]
# # x = phi(temp[ temp["tacs_amp"]==1], twocol=False ) 
# # y = phi(temp[ temp["tacs_amp"]==0] , twocol=False) 
# x = phi(df[(cn,1)])  
# y = phi(df[(cn,0)] ) 

# res = permutation_test((x, y), dPLV, vectorized=True,
#                        n_resamples=1000, alternative='greater')

# a = plt.hist(res.null_distribution,bins =30 )
# plt.plot([res.statistic, res.statistic], [0, a[0].max()], color="orange")
# plt.text(res.statistic+1e-3,a[0].max(), f"pval ={round(res.pvalue,2)}", ha ="left", color="orange" )
# plt.show()

#%% ISI calculation
# ###############################################################################
# df = pd.read_pickle(f"results/allCell_w10Hz.pkl")
# df.reset_index(inplace=True)
# ISI= []
# for i in range(len(df)):
#     t = df["t_spikes_tacs"][i]
#     ISI = t[1:] - t[:-1]
#     ISIs.append(ISI)
# df["ISI"] =ISIs

# freq = 5
# for c in L4_LBC:
#     isimax, disi = 500, 10
#     a,b = np.histogram(np.array(df[(df["cell_name"]==c) & (df["tacs_amp"]==10)& (df["tacs_freq"]==freq) ]["ISI"])[0], bins = np.arange(0,isimax,disi))
#     c,d = np.histogram(np.array(df[(df["cell_name"]==c) & (df["tacs_amp"]==0)& (df["tacs_freq"]==freq) ]["ISI"])[0], bins = np.arange(0,isimax,disi))
#     plt.bar((b[:-1] + b[1:])/2, a, width= disi)
#     plt.bar((b[:-1] + b[1:])/2, c, width= disi, alpha = 0.5)
#     plt.show()
# from scipy.linalg import toeplitz
# toep = toeplitz(a[1:], a[:-1])
# plt.imshow(toep)



#%%
# d= {}
# for freq in [5,10,20,40]:
#     plvs = pd.read_csv(f"results/PLV_w10Hz_{freq}Hz_tACS.csv")

#%% PLV by cell for one freq and linear regression for each cell to get slopes
###############################################################################
## Getting 2D dataframe of PLV[ cell_name, tacs_amp]
# import pandas as pd
# import pickle as pkl
# import numpy as np
# import matplotlib.pyplot as plt
# freq=10
# df = pd.read_pickle(f"results/allCell_w10Hz.pkl")
# df_freq = df[df["tacs_freq"]==freq]
# cell_names =  L5#np.unique(df["cell_name"])
# d={}
# dpval = {}
# amps = np.array([0,1,2,3,4,5,10])
# for cell_name in cell_names:
#     d[cell_name] = df_freq[df_freq["cell_name"]==cell_name].sort_values("tacs_amp")["PLV"].values
#     dpval[cell_name] = df_freq[df_freq["cell_name"]==cell_name].sort_values("tacs_amp")["pval"].values
# d = pd.DataFrame(d, columns=cell_names, index=amps)
# dpval = pd.DataFrame(dpval, columns=cell_names, index=amps)

#### WITH OLD DF with multiple indexes
# PLVs = dfs.loc["PLV"]
# amps = np.array([0,1,2,3,4,5,10])
# temp = []
# cell_names = []
# PLVs = dfs.loc["PLV"]
# for cell_name in np.unique(PLVs.reset_index(level=0).level_0 ):
#     temp.append(PLVs[cell_name].values)
#     cell_names.append(cell_name)
# temp = np.vstack(temp)
# d = pd.DataFrame(d, columns=cell_names, index=amps)
# d.to_csv(f"results/PLV_w10Hz_{freq}Hz_tACS.csv")
# # df = pd.read_csv(f"results/PLV_w10Hz_{freq}Hz_tACS.csv")



