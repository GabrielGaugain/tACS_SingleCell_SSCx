#%%
import time
import neuron
# from neuron import h
from neuron.units import mV, ms
import os
import numpy as np
import matplotlib.pyplot as plt
# from listcells import  L5, L23
import pickle as pkl
from matplotlib.patches import Circle

## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results"
if not os.path.exists(resdir): os.mkdir(resdir)
neuron.h.celsius = 36

## LOADING USEFULL HOC SCRIPTS
neuron.h.load_file("nrngui.hoc")
neuron.h.load_file("aberra_scripts/mySettings.hoc")
neuron.h.load_file("aberra_scripts/interpCoordinates.hoc")
neuron.h.load_file("aberra_scripts/setPointers.hoc")
neuron.h.load_file("Estim.hoc")
neuron.h.load_file("aberra_scripts/cellChooser.hoc")
neuron.h.load_file("aberra_scripts/editMorphology.hoc")
neuron.h.load_file("aberra_scripts/OneSynapticInput.hoc")


cell_names = [f"L1_NGC-DA_bNAC219_{i}" for i in range(1,6)]+ [f"L23_PC_cADpyr229_{i}" for i in range(1,6)] + [f"L4_LBC_cACint209_{i}" for i in range(1,6)]+  [f"L5_TTPC2_cADpyr232_{i}" for i in range(1,6)]+ [f"L6_TPC_L4_cADpyr231_{i}" for i in range(1,6)]

class AberraCell:
    """
    Class to deal with BBP template cells in python. 
    Usefull hoc script for cell instantiation
    are used as well as for extracellular pilotting. 
    Params :
    --------
            termperature (float) : temperature at which the simulation is run in NEURON.
            cell_id (str) : name of the cell ( BBP nomenclature).
            tstop (float) : Duration of the simulation. 
            save_res (bool) :  Whether or not the recordings are saved
            steady_state (bool) : whether or not calculating/ setting the steady state
            use_BBPsynapses (bool) : Whether or not BBP synapses are loaded
            save_synapses (bool) : whether to save synaptic events or not
            syn_freq (float) : mean frequency of all presynaptic cells
            syn_event (list) : list of neuron Vector where synaptic events are stored to drive synapses
            tacs_params (dict) : parameters for tACS stimulation. Default: no tACS (or SHAM)
            iclamp (dict) : parameter for a soma current clamp. Default: no current clamp inserted 
    """

    def __init__(self, temperature = 36,tstop = 200,
                 cell_id = None , specie_type = 1, myelinate_ax = 1,
                save_res = True, rec_axons=False,
                syn_weight = None, syn_freq = 50, syn_events = None, save_synapses = False,
                tacs_params = None, iclamp = None): #, fullAxon = False, use_BBPsynapses = False, steady_state = None,
        time_init = time.time()
        ## setting NEURON temperature (affecting channel kinetics)
        neuron.h.celsius = temperature
        neuron.h.tstop = tstop
        self.cell_id = cell_id
        self.cell_name = cell_names[cell_id]
        self.specie_type = specie_type
        self.myelinate_ax = myelinate_ax
        neuron.h.myelinate_ax =myelinate_ax

        self.tstop = tstop
        # self.use_BBPsynapses = use_BBPsynapses
        # self.syn_events= syn_events
        self.syn_freq = syn_freq
        self.save_synapses = save_synapses
        self.save_res = save_res
        self.tacs_params = tacs_params
        # self.fullAxon = fullAxon
        
        for sec in neuron.h.allsec():
            neuron.h.delete_section(sec=sec)

        # self.comment_axon()

        # self.cell = None
        self.instanciate_AberraCell() ## loading BBP cell


        ## getting steady state | not very usefull coz u have to init after and then call back this function
        # if steady_state is True:
        #     self.get_steady_state()
        ## Instanciate class containing BBP synapses informations 
        # if self.use_BBPsynapses:
        #     self.synapses = BBPsynapses(self.cell_name, self.cell, 
        #                                 activated = True, list_synevents=syn_events, freq = syn_freq)
        if syn_weight is not None:
            self.syn_weight = syn_weight
            self.syn = self.set_single_synapse_pyr(weight = syn_weight, freq = syn_freq)
            self.syn_events = self.syn[0]
        # self.set_results_directories() ## Creating folders to save data in

        ## set recordings
        self.recordings = self.set_recordings(record_axons= rec_axons)
        # self.recordings["cell_name"] =  self.cell_name  ## Already contained in the new dict saved

        # if save_synapses:
        #     self.synapses.record_synapses()

        ## Extracellular stimulation with sinusoidal wave (for tACS in our case)
        if (tacs_params is not None) and (tacs_params["tacs_amp"] !=0) : 
            self.set_tACS() ## to be called again if something is changed

        ## Soma current clamp to generate activity
        if iclamp is not None:
            self.add_soma_clamp(iclamp["amp"], iclamp["dur"], iclamp["del"])

        print(f"Took {round(time.time() -time_init,2)} seconds to init {self.cell_name}.")
    ## END __init__()


    ############################
    ##   Basic running tools  ##
    ############################

    def instanciate_AberraCell(self):

        neuron.h.setParams(self.specie_type) ## set param for either AdultRats(1) of AdultHuman (2)
        neuron.h.cell_chooser(self.cell_id) 
        os.chdir("..") ## if cell_choser and how script in a subfolder (here aberra_script folder)
        self.cell = neuron.h.cell

        # ## Inserting extracellular and xtra mechanism 
        # for sec in neuron.h.cell.all:
        #     if consider_ax or "axon" not in sec.name():
        #         sec.insert("extracellular")
        #         sec.insert("xtra")
        ## Assigning sections coordinates and call setPointers to enable xtra stimulation mech
        # neuron.h.getcoords()
        # neuron.h.setpointers()
    ## END instanciate_BBPcell()
    
    def initialize(self, vinit = None):
        neuron.h.finitialize() if vinit is None else neuron.h.finitialize(vinit)
    ## END init()

    def run(self):
        ## Lauching simulation
        tsim = time.time() 
        neuron.h.frecord_init()
        print(f"Running simulation of {self.tstop} ms of {self.cell_name} activity")
        neuron.h.continuerun(self.tstop)
        self.tsim = round((time.time() - tsim)/60,2)
        print(f"Simulation took {self.tsim} min to simulate")
        if self.save_res:
            self.save_results()
    ## END run()

    def run_save_eachdt(self, Dt = 10000):

        t_inter = np.arange(Dt,neuron.h.tstop+1, Dt)
        for t_i in t_inter:
            neuron.h.frecord_init() 

            print(t_i," ms")
            neuron.h.continuerun(t_i * ms)

            with open(f"{resdir}/{self.cell_name}/freq{self.tacs_params['tacs_freq']}_amp{self.tacs_params['tacs_amp']}_dur{self.tstop}_t{t_i}.pkl", "wb") as f:
                pkl.dump(self.recordings, f)

        return
    ## END run_save_eachdt()
    
    ############################
    ##   Recordings and data  ##
    ############################

    def set_recordings(self, dt = neuron.h.dt, record_axons=False):

        recordings = {  
            "t": neuron.h.Vector().record(neuron.h._ref_t, dt),
            "soma": neuron.h.Vector().record(neuron.h.cell.soma[0](0.5)._ref_v, dt),
            "axon[0]": neuron.h.Vector().record(neuron.h.cell.axon[0](0.5)._ref_v, dt),
        }
        if self.cell_id>=15 :
            recordings["apic[10]"] =  neuron.h.Vector().record(neuron.h.cell.apic[10](0.5)._ref_v, dt)

        if record_axons:
            if self.myelinate_ax ==1:
                recordings["Node[0]"] =  neuron.h.Vector().record(neuron.h.Node[0](0.5)._ref_v, dt)
                recordings["Node[1]"] =  neuron.h.Vector().record(neuron.h.Node[1](0.5)._ref_v, dt)
                recordings["Unmyelin[0]"] =  neuron.h.Vector().record(neuron.h.Unmyelin[0](0.5)._ref_v, dt)
        
            else: 
                recordings["axon[1]"] =  neuron.h.Vector().record(neuron.h.cell.axon[1](0.5)._ref_v, dt)
                recordings["axon[2]"] =  neuron.h.Vector().record(neuron.h.cell.axon[2](0.5)._ref_v, dt)
        ## END if
        
        return recordings
    ## END set_recordings()

    def get_all_params(self, recordings = True):
        param = {}
        for k in self.__dict__.keys():
            if type(self.__dict__[k]) is not neuron.hoc.HocObject:
                param[k] = self.__dict__[k]
        ## removing hoc obj ref that cannot be pickled
        if "cell" in param.keys(): param.pop("cell") 
        if "synapses" in param.keys(): param.pop("synapses")

        if not recordings:
            param.pop("recordings")

        return param

    def save_results(self, fn = None, recordings = True):
        t_presave = time.time()
        
        ## saving Input parameters as object attr + recordings
        res = self.get_all_params(recordings)

        if self.save_synapses:
            res["syn_events"] = self.synapses.syn_records

        ## Naming the results file according to tACS params
        if fn is None:
            if (self.tacs_params is None) or (self.tacs_params['tacs_amp'] == 0):
                fn = f"{self.resdir}/control_dur{round(neuron.h.tstop/1e3)}s.pkl"
            else:
                fn = f"{self.resdir}/tACS/amp{self.tacs_params['tacs_amp']}_freq{self.tacs_params['tacs_freq']}_dur{round(neuron.h.tstop/1e3,2)}s.pkl"

        ## Saving results with cell param
        with open(fn, "wb") as f:
            pkl.dump(res, f)

        print(f'Took {round(time.time() - t_presave)} seconds to save the results')
    ## END save_results()

    def save_recordings(self, fn = None):
        
        if fn is None:
            if (self.tacs_params is None) or (self.tacs_params['tacs_amp'] == 0):
                fn = f"{self.resdir}/control_dur{round(neuron.h.tstop/1e3)}s_recordings.pkl"
            else:
                fn = f"{self.resdir}/tACS/amp{self.tacs_params['tacs_amp']}_freq{self.tacs_params['tacs_freq']}_dur{round(neuron.h.tstop/1e3,2)}s_recordings.pkl"

        ## Saving recordings
        with open(fn, "wb") as f:
            pkl.dump(cell.recordings, f)        
    ## END save_recordings()

    ##########################
    ##   Synaptic inputs    ##
    ##########################

    def set_single_synapse_pyr(self, weight = 0.4, n_comp = 10, x_comp = 0.5, freq=50.,  vec = None):

        if vec is None:
            ## Spike generator
            stim = neuron.h.NetStim()
            stim.interval = 1000/freq  # ms (mean) time between spikes
            stim.number = 100000000000  # (average) number of spikes
            stim.start = 0 # ms (mean) start time of first spike
            stim.noise = 1 # range 0 to 1. Fractional randomness
            stim.seed(1)
        else:
            stim = None ## To be coded to input a pre computed pre-spiketrain

        ## Synapse param
        synapse = neuron.h.Exp2Syn(neuron.h.cell.apic[n_comp](x_comp)) ## in the apical tree for pyramidal cell 
        synapse.e = 0 * mV 
        synapse.tau1 = 2*ms
        synapse.tau2 = 10*ms

        ## Synaptic connection with spike generator
        nc = neuron.h.NetCon(stim, synapse )
        nc.threshold = 10
        nc.delay = 0
        nc.weight[0] = weight	

        ## recording spike synaptic events
        syn_events = neuron.h.Vector()
        nc.record(syn_events)

        return syn_events, stim, synapse, nc, 
    ## END set_single_synapse()

    ##########################
    ##   Point processes    ##
    ##########################

    def add_soma_clamp(self, amp = 0.7, duration = 300, delay = 0):

        iclamp = neuron.h.IClamp(0.5, sec = neuron.h.cell.soma[0])
        iclamp.delay = delay
        iclamp.dur = duration
        iclamp.amp = amp
        self.iclamp = iclamp
    ## END add_clamp()

    def setstim_python(self, amp=1., freq = 10., phi = 0., dur = 1000., offset =0. ):
        """
        code any time waveform for the stimulation into a vector here
        """
        # n = int(neuron.h.tstop/neuron.h.dt +1)
        ti =np.arange(0,self.tstop+neuron.h.dt,neuron.h.dt)
        if freq == 0.:
            # print(f"tDCS stimulation set from {offset}ms to {offset+dur}ms ")
            s = amp * ( (ti>offset)&(ti<= offset+dur) )
        else:
            print(f"{freq}Hz tACS stimulation set from {offset}ms to {offset+dur}ms with amp {amp}V/m ")
            s = amp * np.sin( 2*np.pi* freq * (ti - offset) *1e-3 + phi ) * ( (ti>offset)&(ti<= offset+dur) )

        ## Assigning stim_time and stim_amp and attach the stim (link to pointer stim of xtra mech)
        neuron.h.stim_time = neuron.h.Vector(ti)
        neuron.h.stim_amp = neuron.h.Vector(s)
        neuron.h.attach_stim()
        return s
    ## END setstim_python()

    def set_tACS(self):
        ## Extracellular stimulation (tACS or tDCS if f == 0 )
        if ( self.tacs_params["tacs_dur"] is None) & ( self.tacs_params["tacs_del"] is None):
            self.tacs_params["tacs_dur"] = 3*neuron.h.tstop/4 
            self.tacs_params["tacs_del"] = neuron.h.tstop/4
            
        s = self.setstim_python(amp = self.tacs_params["tacs_amp"], 
                                freq = self.tacs_params["tacs_freq"], 
                                phi = self.tacs_params["tacs_phi"], 
                                dur = self.tacs_params["tacs_dur"], 
                                offset = self.tacs_params["tacs_del"])
        # neuron.h.setstim(neuron.h.DEL,neuron.h.DUR,20.0,neuron.h.FREQ)
        neuron.h.calcesE(0,-1,0) ## 1 V/m, normalized with amp param set into the waveform vector
        self.recordings["tacs_waveform"] = neuron.h.stim_amp

        return s
    ## END set_tACS()

    ##########################
    ##   Usefull methods    ##
    ########################## 

    def get_xyz(self):
        self.xyz = []
        self.diam = []
        for sec in neuron.h.allsec():
            for i in range(sec.n3d()):
                self.xyz.append([sec.x3d(i) , sec.y3d(i) , sec.z3d(i) ])
                self.diam.append( sec.diam3d(i) )

        return self.xyz, self.diam
    ## END get_xyz()

    def save_state(self, fn = None, savesyn = False):
        l = []

        for sec in self.cell.all:
            dic = sec.psection()
            v = np.array([ seg.v for seg in sec])
            l.append( {"v":v, 
                       "density_mechs": dic["density_mechs"], 
                       "ions": dic["ions"]} )

        if fn == None:
            fn = f"{self.resdir}/States/state_t{neuron.h.t}ms.pkl"

        with open(fn, "wb") as f:
            pkl.dump(l, f)

        if savesyn:
            self.synapses.save_state(f"{fn[:-4]}_synapses.pkl")
    ## END save_state()

    def restore_state(self, d):
        if type(d) == str:
            with open(d, "rb") as f:
                d = pkl.load(f)
        try:
            for i,sec in enumerate(self.cell.all):
                for j, seg in enumerate(sec):
                    seg.v = d[i]["v"][j]

                    ## setting density mechanisms with taking car to extracellular as it is calling differently
                    mechanisms = d[i]["density_mechs"]
                    for m in mechanisms:
                        if m != "extracellular" and m!="xtra":
                            for sub_m in mechanisms[m]:
                                exec( f"seg.{m}.{sub_m} = {mechanisms[m][sub_m][j]}"  )
                        elif m== "extracellular":
                            exec( f"seg.vext[0] = {mechanisms[m]['vext'][j][0]}"  )
                            exec( f"seg.vext[1] = {mechanisms[m]['vext'][j][1]}"  )
                            exec( f"seg.i_membrane = {mechanisms[m]['i_membrane'][j]}"  )

                    ## setting ions concentrations
                    ions = d[i]["ions"]
                    for ion in ions:
                        for ci in ions[ion]:
                            exec( f"seg.{ion}_ion.{ci} = {ions[ion][ci][j]}"  )

        except:
            print("Could not affect density mechs and ions to each seg")
            pass
    ## END restore_state()
    
    def compute_steady_state(self, tol = 1e-6, savestate = True):

        print(f"Computing steady state, it may take a while...")
        tss = time.time()
        neuron.h.init(-70)
        ## init t-1 voltages list
        previous_v = []    
        for i,sec in enumerate(neuron.h.cell.all):
            previous_v.append(sec(0.5).v)
        ## simulating till convergence ( dV < 1e-6 mV= 1 nV )
        cond = True
        while cond and neuron.h.t<20e3:
            neuron.h.fadvance()

            cond=False
            ## setting current voltage to compare with previous voltage
            curent_v =[]
            for i,sec in enumerate(neuron.h.cell.all):
                curent_v.append(sec(0.5).v)
                if (not cond) & ( np.abs(curent_v[i]-previous_v[i]) > tol ):
                    cond=True
            previous_v = curent_v
            
        if cond and neuron.h.t>60e3:
            print("failed to have a steady state in less than 20 s of activity..")
            exit("SteadyState not reached")

        print(f"Steady state reached at t = {neuron.h.t} ms and took {round((time.time() -tss)/60, 2)} min to simulate")
        neuron.h.continuerun(neuron.h.t + 5000)
        ## Saving the corresponding steady state in cell folder "steadyState"
        suffix = "_fullAxon" if self.fullAxon else ""
        if savestate:
            self.save_state(f"{highdir}/results/{self.cell_name}/States/steadyState_{self.cell_name}{suffix}.pkl")
    ## END compute_steady_state()

    def get_steady_state(self):
        # print("getting steady state")
        suffix = "_fullAxon" if self.fullAxon else ""
        steady_fn = f"{highdir}/results/{self.cell_name}/States/steadyState_{self.cell_name}{suffix}.pkl"

        if os.path.exists(steady_fn):
            # print(f"getting steady state : {steady_fn}")
            self.restore_state(steady_fn)
        else:
            self.compute_steady_state()
    ## END get_steady_state()


    def plot2D(self, plan = "xy", show_scale = True, linewidth = None,
               colors = ["lightgray", "firebrick","coral", "royalblue", ], ## soma, dendrites, apic, axon
               savefig= False, figname=None, show_legend=False):
        
        (ix, iy) = (0,1) if (plan == "xy")or(plan == "yx") else (2,1) if (plan == "zy")or(plan == "yz") else (0,2)
    
        def plot_seclist(ax, seclist, color, name="", 
                    circle = False, plan="xy", lw = 1):

            for sec in seclist:
                x,y,z,d = [],[],[],[]
                for i in range(sec.n3d()):
                    x.append(sec.x3d(i))
                    y.append(sec.y3d(i))
                    z.append(sec.z3d(i))
                    d.append(sec.diam3d(i))
                X = np.array([x,y,z,d])
                if circle:
                    [ ax.add_patch( Circle((X[ix,i], X[iy,i]), 0.75*X[3,i], 
                                                        fc = color, lw =lw, ec=color, 
                                                        zorder=0) ) 
                                for i in range(sec.n3d()) ]  
                else:
                    ax.plot(X[ix], X[iy], color=color,  linewidth=lw)
            return ax
        ### END plt_seclist()

        Xmin, Xmax = self._get_boundaries()
        xmin, ymin = Xmin[ix], Xmin[iy]
        l = [Xmax[i] - Xmin[i] for i in range(3)]
        max_length = max([l[ix],l[iy]])

        if self.myelinate_ax:
            parts = ["soma", "dend", "apic", "unmyelin", "myelin", "nodes" ]
            if len(colors)==4:
                colors += ["dimgrey","fuchsia"]

        fig, ax = plt.subplots(figsize=(150/27.4, 150/27.4), dpi=150)
        # colors = ["red", "royalblue","navy", "coral", ] ## soma, dendrites, apical, axon
        if type(colors) == str:
            print("colors entered as a string")
            colors = [ colors for i in range(len(parts))]

        
        ax = plot_seclist(ax, neuron.h.Myelin, color=colors[4], name="myelin", lw=1)
        ax = plot_seclist(ax, neuron.h.Unmyelin, color=colors[3], name="Unmyelin", lw=1.5)
        ax = plot_seclist(ax, neuron.h.Node, color=colors[5], name="Nodes", 
                          circle = True, lw=1.6)

        ax = plot_seclist(ax, self.cell.soma, color=colors[0], name="soma", 
                          circle=True, lw= 0.1,)
        ax = plot_seclist(ax, self.cell.dend, color=colors[1], name="dend", lw=0.5)
        ax = plot_seclist(ax, self.cell.apic, color=colors[2], name="apical", lw=0.5)
        
        # ax.axis(ax.axis('equal')) 
        ax.set_aspect('equal', "box")
        [ax.spines[side].set_visible(False) for side in ['right','left', 'top', 'bottom']]
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


        ## Scale
        if show_scale:
            lbar = 100
            if max_length<1000 : lbar = 50
            ax.plot([xmin-150-lbar, xmin-150], [ymin-10, ymin-10],'k')
            ax.text( xmin-140, ymin-10,f'{lbar} µm',va='center')

            ax.plot([xmin-150-lbar, xmin-150-lbar],[ymin-10, ymin-10+lbar],'k')
            ax.text( xmin-150-lbar, ymin+5+lbar,f'{lbar} µm',ha='center')
        ##legend
        if show_legend:
            from matplotlib.lines import Line2D
            custom_lines=[Line2D([0], [0], color=colors[i], lw=1) if parts[i] != "soma" else Line2D([0], [0],  marker='o',color=colors[i], lw=0) for i in range(len(colors)) ]
            plt.legend(custom_lines, parts,
                    bbox_to_anchor=(0.95,0.95), borderaxespad=0,
                        frameon=False, fontsize=9)

        if savefig :
            # morpho = os.path.split(self.morphology)[-1]
            if figname is None:
                plt.savefig(f"{self.resdir}/figures/morpho{plan}_{self.cell_name}.jpg", dpi = 500)
                plt.savefig(f"{self.resdir}/figures/morpho_{plan}_{self.cell_name}.eps", dpi = 500)
                plt.savefig(f"{self.resdir}/figures/morpho_{plan}_{self.cell_name}.svg")
            else: 
                plt.savefig(f"{figname}.jpg", dpi = 250)
                # plt.savefig(f"{figname}.eps", dpi = 500)
                plt.savefig(f"{figname}.svg")

        # plt.show()

        return 
    ## END plot_2D()
    
    def set_results_directories(self):
        self.resdir = f"results/{self.cell_name}"
        if not os.path.exists(self.resdir): os.mkdir(self.resdir)

        dirs = ["States", "post_processing", "tACS","figures"]
        for d in dirs:
            if not os.path.exists(f"{self.resdir}/{d}"): os.mkdir(f"{self.resdir}/{d}")
    ## END set_results_directories()

    def _get_boundaries(self):
        self.x_min = min([ min([sec.x3d(i) for i in range(sec.n3d())]) for sec in self.cell.all if sec.n3d()!=0 ])
        self.y_min = min([ min([sec.y3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.z_min = min([ min([sec.z3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])

        self.x_max = max([ max([sec.x3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.y_max = max([ max([sec.y3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.z_max = max([ max([sec.z3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        xyz_max = np.array([self.x_max, self.y_max, self.z_max])
        xyz_min = np.array([self.x_min, self.y_min, self.z_min])

        return xyz_min, xyz_max
    ## END _get_boundaries()

    ##############################
    ##   Specific running tools ##
    ##############################


    def get_spike_timing(self, sec="axon[0]"):
        
        print("getting spike timing")
        import efel
        traces, trace = [] , {}
        trace["T"] = self.recordings["t"]
        trace["V"] = self.recordings[sec]
        trace["stim_start"] = [0]
        trace["stim_end"] = [self.recordings["t"][-1]]
        traces.append(trace)

        feats = ["peak_time", "mean_frequency"]
        tr_res = efel.getFeatureValues(traces, feats)
        
        self.mean_frequency = tr_res[0]["mean_frequency"][0] if tr_res[0]["mean_frequency"] is not None else 0
        self.t_spikes = tr_res[0]["peak_time"]
    ## END get_spike_timing()

    def get_steady_state(self):
        # print("getting steady state")
        suffix = "_fullAxon" if self.fullAxon else ""
        steady_fn = f"{highdir}/results/{self.cell_name}/States/steadyState_{self.cell_name}{suffix}.pkl"

        if os.path.exists(steady_fn):
            # print(f"getting steady state : {steady_fn}")
            self.restore_state(steady_fn)
        else:
            self.compute_steady_state()
    ## END get_steady_state()


    def plot2D(self, plan = "xy", show_scale = True, linewidth = None,
               colors = ["lightgray", "firebrick","coral", "royalblue", ], ## soma, dendrites, apic, axon
               savefig= False, figname=None, show_legend=False):
        
        (ix, iy) = (0,1) if (plan == "xy")or(plan == "yx") else (2,1) if (plan == "zy")or(plan == "yz") else (0,2)
    
        def plot_seclist(ax, seclist, color, name="", 
                    circle = False, plan="xy", lw = 1):

            for sec in seclist:
                x,y,z,d = [],[],[],[]
                for i in range(sec.n3d()):
                    x.append(sec.x3d(i))
                    y.append(sec.y3d(i))
                    z.append(sec.z3d(i))
                    d.append(sec.diam3d(i))
                X = np.array([x,y,z,d])
                if circle:
                    [ ax.add_patch( Circle((X[ix,i], X[iy,i]), 0.75*X[3,i], 
                                                        fc = color, lw =lw, ec=color, 
                                                        zorder=0) ) 
                                for i in range(sec.n3d()) ]  
                else:
                    ax.plot(X[ix], X[iy], color=color,  linewidth=lw)
            return ax
        ### END plt_seclist()

        Xmin, Xmax = self._get_boundaries()
        xmin, ymin = Xmin[ix], Xmin[iy]
        l = [Xmax[i] - Xmin[i] for i in range(3)]
        max_length = max([l[ix],l[iy]])

        if self.myelinate_ax:
            parts = ["soma", "dend", "apic", "unmyelin", "myelin", "nodes" ]
            if len(colors)==4:
                colors += ["dimgrey","fuchsia"]

        fig, ax = plt.subplots(figsize=(150/27.4, 150/27.4), dpi=150)
        # colors = ["red", "royalblue","navy", "coral", ] ## soma, dendrites, apical, axon
        if type(colors) == str:
            print("colors entered as a string")
            colors = [ colors for i in range(len(parts))]

        if self.myelinate_ax:
            ax = plot_seclist(ax, neuron.h.Myelin, color=colors[4], name="myelin", lw=1)
            ax = plot_seclist(ax, neuron.h.Unmyelin, color=colors[3], name="Unmyelin", lw=1.5)
            ax = plot_seclist(ax, neuron.h.Node, color=colors[5], name="Nodes", 
                          circle = True, lw=1.6)
        else:
            ax = plot_seclist(ax, cell.cell.axonal, color=colors[3], name="axon", lw=1)

        ax = plot_seclist(ax, self.cell.soma, color=colors[0], name="soma", 
                          circle=True, lw= 0.1,)
        ax = plot_seclist(ax, self.cell.dend, color=colors[1], name="dend", lw=0.5)
        ax = plot_seclist(ax, self.cell.apic, color=colors[2], name="apical", lw=0.5)
        
        # ax.axis(ax.axis('equal')) 
        ax.set_aspect('equal', "box")
        [ax.spines[side].set_visible(False) for side in ['right','left', 'top', 'bottom']]
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)


        ## Scale
        if show_scale:
            lbar = 100
            if max_length<1000 : lbar = 50
            ax.plot([xmin-150-lbar, xmin-150], [ymin-10, ymin-10],'k')
            ax.text( xmin-140, ymin-10,f'{lbar} µm',va='center')

            ax.plot([xmin-150-lbar, xmin-150-lbar],[ymin-10, ymin-10+lbar],'k')
            ax.text( xmin-150-lbar, ymin+5+lbar,f'{lbar} µm',ha='center')
        ##legend
        if show_legend:
            from matplotlib.lines import Line2D
            custom_lines=[Line2D([0], [0], color=colors[i], lw=1) if parts[i] != "soma" else Line2D([0], [0],  marker='o',color=colors[i], lw=0) for i in range(len(colors)) ]
            plt.legend(custom_lines, parts,
                    bbox_to_anchor=(0.95,0.95), borderaxespad=0,
                        frameon=False, fontsize=9)

        if savefig :
            # morpho = os.path.split(self.morphology)[-1]
            if figname is None:
                plt.savefig(f"{self.resdir}/figures/morpho{plan}_{self.cell_name}.jpg", dpi = 500)
                plt.savefig(f"{self.resdir}/figures/morpho_{plan}_{self.cell_name}.eps", dpi = 500)
                plt.savefig(f"{self.resdir}/figures/morpho_{plan}_{self.cell_name}.svg")
            else: 
                plt.savefig(f"{figname}.jpg", dpi = 250)
                # plt.savefig(f"{figname}.eps", dpi = 500)
                plt.savefig(f"{figname}.svg")

        # plt.show()

        return fig,ax
    ## END plot_2D()
    
    def set_results_directories(self):
        self.resdir = f"results/{self.cell_name}"
        if not os.path.exists(self.resdir): os.mkdir(self.resdir)

        dirs = ["States", "post_processing", "tACS","figures"]
        for d in dirs:
            if not os.path.exists(f"{self.resdir}/{d}"): os.mkdir(f"{self.resdir}/{d}")
    ## END set_results_directories()

    def _get_boundaries(self):
        self.x_min = min([ min([sec.x3d(i) for i in range(sec.n3d())]) for sec in self.cell.all if sec.n3d()!=0 ])
        self.y_min = min([ min([sec.y3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.z_min = min([ min([sec.z3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])

        self.x_max = max([ max([sec.x3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.y_max = max([ max([sec.y3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        self.z_max = max([ max([sec.z3d(i) for i in range(sec.n3d())]) for sec in self.cell.all  if sec.n3d()!=0])
        xyz_max = np.array([self.x_max, self.y_max, self.z_max])
        xyz_min = np.array([self.x_min, self.y_min, self.z_min])

        return xyz_min, xyz_max
    ## END _get_boundaries()




# from listcells import *
Params = {"temperature": 36, "tstop":50,"cell_id": 15  , 
          "specie_type" : 2, "myelinate_ax": 0,
                "save_res" : False, "rec_axons":True,
                "syn_weight" :0.0 , "syn_freq" :50,
        "iclamp":{"amp":1.5, "dur":80, "del":10}
        }


# Params = {"temperature": 36, "tstop": 10, "weight": 0.2  , 
#           "cell_id": 15 , "specie_type": 2, "myelinated": 1,
#           "tacs_amp": 0., "tacs_freq": 10, "tacs_phi": 0, "tacs_dur": None, "tacs_del": None,
#           "clamp_amp": 0.5, "clamp_dur":100, "clamp_del":0,
#           "save_rec": True
#           }

cell = AberraCell(**Params)

# #%%
# cell.initialize()
# cell.run()
# ax = cell.plot_results(plot_ax=True)
# ax.set_xlim((18,22))
# # cell.plot2D(show_legend =True )
# plt.show()
#%% AP initiation location
cell.tstop = 50

cell.initialize()
spike_timer=0
while (neuron.h.t < cell.tstop-neuron.h.dt/2) and (spike_timer==0): 
    neuron.h.fadvance()
    for sec in cell.cell.axonal:
        if (sec(0.5).v > 0):
            spike_timer = 1
            print(f"AP initiated at time {neuron.h.t} ms")
            AP_site = [sec(0.5).x_xtra,sec(0.5).y_xtra,sec(0.5).z_xtra]
            break
fig,ax = cell.plot2D(show_legend =False, show_scale=False )
ax.plot(AP_site[0], AP_site[1], "ro", fillstyle='none',markersize=10, markeredgewidth=2)
if cell.myelinate_ax:
    plt.savefig("figures/PB_Aberra/APloc_withmyelin.svg")
    plt.savefig("figures/PB_Aberra/APloc_withmyelin.png", dpi=200,transparent=True)
else:
    plt.savefig("figures/PB_Aberra/APloc_withoutmyelin.svg")
    plt.savefig("figures/PB_Aberra/APloc_withoutmyelin.png", dpi=200,transparent=True)

plt.show()

#%%
## Custom recording in the myelinated tree
recordings = cell.recordings
dt = neuron.h.dt

# secs = ["axon[0]","Node[0]","Node[2]", "Unmyelin[0]", "Unmyelin[1]", "Unmyelin[2]"]
secs = ["axon[0]","axon[2]","axon[4]","axon[6]"]

for sec in secs:

    if ("Node" in sec) or ("Unm" in sec):
        recordings[sec] =  neuron.h.Vector().record(eval(f"neuron.h.{sec}(0.5)._ref_v"), dt)
    else:
        recordings[sec] =  neuron.h.Vector().record(eval(f"neuron.h.cell.{sec}(0.5)._ref_v"), dt)
  
cell.tstop = 50

cell.initialize()
cell.run()


#%%
distance = neuron.h.distance
distance(sec = cell.cell.soma[0])
from matplotlib import colormaps
cmap = colormaps["viridis"]
from plot_utils import mm

d= []
fig, ax = plt.subplots(figsize=(120*mm,80*mm)) #nrows = 2)
ax.plot(cell.recordings["t"], cell.recordings["soma"], label="soma", color = cmap(0))
for sec in secs:
    if ("Node" in sec) or ("Unm" in sec):
        d.append(distance(eval(f"neuron.h.{sec}(0.5)")))
    else:
        d.append(distance(eval(f"neuron.h.cell.{sec}(0.5)"))) 
d= np.array(d)
# d /=d.max()
dmax = 400
d /=dmax
for i,sec in enumerate(secs):
    ax.plot(cell.recordings["t"], cell.recordings[sec], label=sec, color =cmap(d[i]) )

ax.legend(frameon=False)
ax.set_xlabel("time (ms)")
ax.set_ylabel(r"$V_m$ (mV)")
ax.set_yticks([-80,-40,0,40])
# ax.set_xticks([[-80,-40,0,40]])


# ax.set_xlim((19,22))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=dmax))
plt.colorbar(sm, label = "distance to soma (µm)", ticks=[0,200,400])

# ax.set_xlim((14.5,17))
# ax.set_xlim((19,21))
# plt.savefig("figures/PB_Aberra/aberra_propa_human.svg")
# plt.savefig("figures/PB_Aberra/aberra_propa_human.jpg", dpi=250)

ax.set_xlim((45,47))
plt.savefig("figures/PB_Aberra/BBP_propa_V2.svg")
plt.savefig("figures/PB_Aberra/BBP_propa_V2.jpg", dpi=250)

#%%
# if __name__== "__main__":
   
    # Params = {"temperature": 36, "tstop":100,"cell_id": 15  , 
    #         "specie_type" : 2, "myelinate_ax": 1,
    #                 "save_res" : False, 
    #                 "syn_weight" :0.02 , "syn_freq" :50
    #         }
#             # "iclamp":{"amp":0.2, "dur":80, "del":10}
#             # }
#     cell = BBPcell(**Params)
#     # cell.run()




