#%%
from neuron import h
import os
import numpy as np
import pickle as pkl
import pandas as pd

## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results/"
if not os.path.exists(resdir): os.mkdir(resdir)
h.celsius = 36

ProbGABAAB_params = [   'tau_r_GABAA', 'tau_d_GABAA', 'tau_r_GABAB', 'tau_d_GABAB', 
                        'Use', 'Dep', 'Fac', 'e_GABAA', 'e_GABAB', 'u0', 
                        'GABAB_ratio', 'i', 'i_GABAA', 'i_GABAB', 'g_GABAA', 'g_GABAB', 'g',
                        'A_GABAA_step', 'B_GABAA_step', 'A_GABAB_step', 'B_GABAB_step', 
                        'Rstate', 'tsyn_fac', 'u', 'rng',
                        'A_GABAA', 'B_GABAA', 'A_GABAB', 'B_GABAB']

ProbAMPANMDA_params = [ 'tau_r_AMPA', 'tau_d_AMPA', 'tau_r_NMDA', 'tau_d_NMDA', 
                        'Use', 'Dep', 'Fac', 'e', 'mg', 'u0', 'NMDA_ratio', 
                        'i', 'i_AMPA', 'i_NMDA', 'g_AMPA', 'g_NMDA', 'g', 
                        'A_AMPA_step', 'B_AMPA_step', 'A_NMDA_step', 'B_NMDA_step', 
                        'Rstate', 'tsyn_fac', 'u', 'rng',
                        'A_AMPA', 'B_AMPA', 'A_NMDA', 'B_NMDA']

def read_tsv(filename):

    col = ["synapse_id" , "pre_cell_id" ,  "pre_mtype" , "sectionlist_id" , 
          "sectionlist_index" , "seg_x" , "synapse_type" ,  
          "dep" , "fac" , "use" , "tau_d" , "delay" , "weight" ]

    return pd.read_csv(filename, sep='\t', names = col,  skiprows=[0])  


class BBPsynapses:

    def __init__(self,cell_name, post_cell, list_synevents = None, 
                 base_seed = 0, activated = True, freq = None) -> None:

        self.post_cell = post_cell
        self.post_cell_name = cell_name
        self.synapses_folder = f"./cells/{self.post_cell_name}/synapses"
        self.synapses = []
        self.input_events = []
        self.rng_list = []
        self.weights = []
        self.delays = []

        ## Creation of all synapses (ex and in) from tsv data file of the cell
        self.init_synapses(base_seed)

        ## Attributing synaptic event to each synapses with either precomputed spike events
        ## or netstim artificial spiking cell with a proper netcon for each synapse
        if list_synevents is not None:
            self.set_vecstim(list_synevents, activated)
        else:
            self.create_netstims(activated, freq)

   ## END __init__()

    def init_synapses(self, base_seed) -> None:
        print(f"Synapses creation from BBP data for cell {self.post_cell_name}")
        df = read_tsv(f'{self.synapses_folder}/synapses.tsv')

        self.nsyn = df.shape[0]
        self.df = df
        self.load_synconf()
        self.synapses_type = []
        self.get_gid()
        # for index, row in df.iterrows():
        for isyn in range(self.nsyn):

            ## creating the synapse from tsv file info
            sec = self.selec_section(df, isyn)
            seg_x = df["seg_x"][isyn] +1e-3 if df["seg_x"][isyn]==0 else df["seg_x"][isyn] - 1e-3 if df["seg_x"][isyn]==1 else df["seg_x"][isyn]

            if df["synapse_type"][isyn] <100:
                synapse = h.ProbGABAAB_EMS( sec( seg_x ) )
                synapse.tau_d_GABAA  = df["tau_d"][isyn]
                rng = h.Random()                                                      
                rng.MCellRan4( isyn*100000+100, self.gid+250+base_seed )                
                rng.lognormal(0.2, 0.1)                                                 
                synapse.tau_r_GABAA = rng.repick()
                self.synapses_type.append(1)
            else:
                synapse = h.ProbAMPANMDA_EMS( sec( seg_x ) )
                synapse.tau_d_AMPA  = df["tau_d"][isyn]
                self.synapses_type.append(0)
            synapse.Use = np.abs(df["use"][isyn])
            synapse.Dep = np.abs(df["dep"][isyn])
            synapse.Fac = np.abs(df["fac"][isyn])
            self.weights.append(df["weight"][isyn])
            self.delays.append(df["delay"][isyn])

            ## Synapse additionnal configuration
            for s in self.synconf[isyn].split("%s")[1:]:
                s = "synapse"+s
                exec(s)

            ## Create the random number generator for the synapse
            rng = h.Random()                                                          
            rng.MCellRan4( isyn*100000+100, self.gid+250+base_seed )                    
            rng.uniform(0,1)                                                            
            synapse.setRNG( rng )                                                       
            self.rng_list.append(rng) 
            self.synapses.append(synapse)

        print(f"Initialized {self.nsyn} synapses")
    ## END init_synapses()

    def selec_section(self, dataframe, index ): 
        sec_id = dataframe["sectionlist_id"][index] 
        sec_index = dataframe["sectionlist_index"][index]
        # id, index):
        if sec_id ==0:
            sec = self.post_cell.soma[ sec_index ]
        elif sec_id ==1:
            sec = self.post_cell.dend[ sec_index ]
        elif sec_id ==2:
            sec = self.post_cell.apic[ sec_index ]
        elif sec_id ==3:
            sec = self.post_cell.axon[ sec_index ]
        else:
            print(f"error : sectionlist_id {sec_id} not supported ")
        
        return sec
    ## END selec_section()

    def set_vecstim(self, stim_events, activated = True ):
        assert len(stim_events) == self.df["pre_cell_id"].unique().size, "stim event list should have the same elements number as synapses"
        print("Attributing given synaptic events to all synapses")

        self.vecstim_list = []
        self.netcon_list = []
        self.pre_cell_ids = []

        for i in range(self.nsyn):
            pre_cell_id = self.df["pre_cell_id"][i]

            if pre_cell_id in self.pre_cell_ids:
                j = self.pre_cell_ids.index(pre_cell_id) 
                stim = self.vecstim_list[ j ]

            else:
                ## equivalent to the netstim for the ith synapse
                self.pre_cell_ids.append(pre_cell_id)
                j = self.pre_cell_ids.index(pre_cell_id) 

                stim = h.VecStim()
                stim.play(stim_events[j])
                self.vecstim_list.append(stim)


            ## connecting the stim to the synapse
            nc = h.NetCon(stim, self.synapses[i] )
            nc.delay = self.delays[i]
            if activated:
                nc.weight[0] = self.weights[i]	
            else:
                nc.weight[0] = 0
            self.netcon_list.append(nc)
    ## END set_vecstim(self)
    
    def create_netstims(self, activated, freq = None):
        print("Creating netstim for all synapses as artificial spiking cells")        
        self.netstim_list = []
        self.netcon_list = []
        self.pre_cell_ids = []
        self.synapses_netstim = []
        for i in range(self.nsyn):  
            pre_cell_id = self.df["pre_cell_id"][i]

            if pre_cell_id in self.pre_cell_ids:
                stim = self.netstim_list[ self.pre_cell_ids.index(pre_cell_id) ]
            else:
                stim = h.NetStim()
                stim.interval = 100 if freq is None else 1000/freq # ms (mean) time between spikes
                stim.number = 1e20  # (average) number of spikes
                stim.start = 0 # ms (mean) start time of first spike
                stim.noise = 1 # range 0 to 1. Fractional randomness
                self.netstim_list.append(stim)
                self.pre_cell_ids.append(pre_cell_id)

            self.synapses_netstim.append(self.pre_cell_ids.index(pre_cell_id))
            ## connecting the stim to the synapse
            nc = h.NetCon(stim, self.synapses[i] )
            nc.delay = self.delays[i]
            if activated:
                nc.weight[0] = self.weights[i]	
            else:
                nc.weight[0] = 0

            self.netcon_list.append(nc)
    ## END create_netstims()

    def record_synapses(self):
        self.syn_records = []
        ## need to record only synapses with different netcon (diff presyn cell)
        ## cause netcon.record() only record one instance of a shared netstim...
        pre_cell_ids = self.df["pre_cell_id"].drop_duplicates() # first shared netstim occurence with index
        
        for i in range(len(self.pre_cell_ids)):
            self.syn_records.append( h.Vector())
            self.netcon_list[pre_cell_ids.index[i]].record( self.syn_records[i] ) 
    ## END record_synapses()

    def activate_synapses(self):
        for i in range(self.nsyn):
            self.netcon_list[i].weight[0] = self.weights[i]
    ## END activate_synapses()

    def desactivate_synapses(self):
        for i in range(self.nsyn):
            self.netcon_list[i].weight[0] = 0
    ## END desactivate_synapses()

    def set_weights(self):
        for i in range(self.nsyn):
            self.netcon_list[i].weight[0] = self.weights[i]
          
    ## END set_weights

    def change_weight_by_alpha(self, alpha = 0.05):
        for isyn in range(self.nsyn):
            if self.df["synapse_type"][isyn] <100:
                self.weights[isyn]  = (1-alpha)*self.df["weight"][isyn]
            else:
                self.weights[isyn]  = (1+alpha)*self.df["weight"][isyn]

        self.set_weights()
    ## END change_weight_by_alpha()

    def get_gid(self):

        with open(f"{self.synapses_folder}/synapses.hoc","r") as f:
            lines = f.readlines()
        for i,l in enumerate(lines):
            if "gid = " in l:
                line = l
        import re 
        self.gid = int(re.findall(r'\d+', line)[0])
    ## END get_gid()

    def load_synconf(self):
        
        conf, gids = [],[]
        with open(f"{self.synapses_folder}/synconf.txt", "r") as f:
            lines = f.readlines()
        for i in range(len(lines)//2):
            conf.append(lines[2*i][:-1])
            gids.append(lines[2*i+1])
        
        gids = [ np.array(gids[i].split()[:-1]).astype(int) for i in range(len(gids))]

        self.synconf = ["" for x in range(self.nsyn)]
        for i in range(self.nsyn):
            for j in range(len(conf)):
                if i in gids[j]:
                    self.synconf[i] += conf[j]
    ## END load_synconf()

    def get_state(self):
        state = {"synapses":[], "weights":np.array([i.weight[0] for i in self.netcon_list])}

        for i, syn in enumerate(self.synapses):
            d={}

            if "ProbGABAAB" in syn.hname():
                for p in ProbGABAAB_params:
                    exec(f"d['{p}'] = syn.{p}")
                
            elif "ProbAMPANMDA" in syn.hname():
                for p in ProbAMPANMDA_params:
                    exec(f"d['{p}'] = syn.{p}")
            else:
                d["nosyn"] = "Error: neither ProbGABAAB nor ProbAMPANMDA synapse instance ..."

            state["synapses"].append(d)

        return state
    ## END get_state()
    
    def save_state(self, fn):
        s = self.get_state()
        with open(fn, "wb") as f:
            pkl.dump(s,f)
    ## END save_state()

    def set_state(self, state):
        
        if type(state) == str:
            with open(state,"rb") as f:
                state = pkl.load(f)

        self.weights  = state["weights"]
        for i, syn in enumerate(self.synapses):
            s = state['synapses'][i]
            try:
                if "ProbGABAAB" in syn.hname():
                    for p in ProbGABAAB_params:
                        exec(f"syn.{p} = s['{p}']")
                    
                elif "ProbAMPANMDA" in syn.hname():
                    for p in ProbAMPANMDA_params:
                        exec(f"syn.{p} = s['{p}']")
                else:
                    print(f"carreful : {p} is neither ProbAMPANMDA nor ProbGABAAB")

            except:
                print(f"WARNING : couldn't set the {i}th synapse to the state")
                pass

## END BBPsynapses






def main():
    return

if __name__=="__main__":
    main()
# %%
