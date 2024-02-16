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

## GLOBAL CONSTANTS
highdir = os.getcwd()
resdir = "./results"
if not os.path.exists(resdir): os.mkdir(resdir)
neuron.h.celsius = 36

## LOADING USEFULL HOC SCRIPTS
neuron.h.load_file("nrngui.hoc")
neuron.h.load_file("interpCoordinates.hoc")
neuron.h.load_file("setPointers.hoc")
neuron.h.load_file("Estim.hoc")
from BBPsynapses import BBPsynapses


class BBPcell:
    """
    Class to deal with BBP template cells in python. 
    Usefull hoc script for cell instantiation
    are used as well as for extracellular pilotting. 
    Params :
    --------
            termperature (float) : temperature at which the simulation is run in NEURON.
            cell_name (str) : name of the cell ( BBP nomenclature).
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

    def __init__(self, temperature = 36, cell_name = None , 
                tstop = 200, save_res = True, steady_state = None,
                use_BBPsynapses = False, syn_freq = None, syn_events = None, save_synapses = False,
                tacs_params = None, iclamp = None, fullAxon = False ):
        time_init = time.time()
        ## setting NEURON temperature (affecting channel kinetics)
        neuron.h.celsius = temperature
        neuron.h.tstop = tstop
        self.cell_name = cell_name
        self.tstop = tstop
        self.use_BBPsynapses = use_BBPsynapses
        self.syn_freq = syn_freq
        self.syn_events= syn_events
        self.save_synapses = save_synapses
        self.save_res = save_res
        self.tacs_params = tacs_params
        self.fullAxon = fullAxon
        
        for sec in neuron.h.allsec():
            neuron.h.delete_section(sec=sec)

    
        # self.comment_axon()

        # self.cell = None
        self.instanciate_BBPcell(consider_ax=False) ## loading BBP cell


        ## getting steady state | not very usefull coz u have to init after and then call back this function
        if steady_state is True:
            self.get_steady_state()

        ## Instanciate class containing BBP synapses informations 
        if self.use_BBPsynapses:
            self.synapses = BBPsynapses(self.cell_name, self.cell, 
                                        activated = True, list_synevents=syn_events, freq = syn_freq)

        self.set_results_directories() ## Creating folders to save data in

        ## set recordings
        self.recordings = self.set_recordings()
        # self.recordings["cell_name"] =  self.cell_name  ## Already contained in the new dict saved

        if save_synapses:
            self.synapses.record_synapses()

        ## Extracellular stimulation with sinusoidal wave (for tACS in our case)
        # if (tacs_params is not None) and (tacs_params["tacs_amp"] !=0) : 
        #     self.set_tACS() ## to be called again if something is changed

        ## Soma current clamp to generate activity
        if iclamp is not None:
            self.add_soma_clamp(iclamp["amp"], iclamp["dur"], iclamp["del"])


        print(f"Took {round(time.time() -time_init,2)} seconds to init {self.cell_name}.")
    ## END __init__()


    ############################
    ##   Basic running tools  ##
    ############################

    def instanciate_BBPcell(self, consider_ax =False):
        ## CREATING BBP CELL
        os.chdir(f"cells/{self.cell_name}") ## moving to cell dir to use hoc scripts
        neuron.h.load_file("createsimulation.hoc")
        neuron.h.create_cell(False)
        os.chdir(highdir)## moving back to the top dir
        print(f"{self.cell_name} loaded.")
 
        self.cell = neuron.h.cell
        # self.set_axon_xyz()
        # ## Inserting extracellular and xtra mechanism 
        for sec in neuron.h.cell.all:
            if consider_ax or "axon" not in sec.name():
                sec.insert("extracellular")
                sec.insert("xtra")
        ## Assigning sections coordinates and call setPointers to enable xtra stimulation mech
        neuron.h.getcoords()
        neuron.h.setpointers()
    ## END instanciate_BBPcell()

    def comment_axon(self):
        from comment_script import replaceAll
        fn = f"cells/{self.cell_name}/template.hoc"

        if self.fullAxon:
            replaceAll(fn, "    replace_axon()", "    //replace_axon()")
        else:
            replaceAll(fn, "    //replace_axon()", "    replace_axon()")
        pass
    ## END comment_axon()
    
    def set_axon_xyz(self):
        n = 0#self.cell.soma[0].n3d()-1
        x0, y0, z0 = self.cell.soma[0].x3d(n), self.cell.soma[0].y3d(n), self.cell.soma[0].z3d(n)
        L= self.cell.axon[0].L
        self.cell.axon[0].pt3dadd(x0, y0-0.5, z0, self.cell.axon[0].diam)
        self.cell.axon[0].pt3dadd(x0, y0-0.5-self.cell.axon[0].L, z0, self.cell.axon[0].diam)


        self.cell.axon[1].pt3dadd(x0, y0-0.5-self.cell.axon[0].L, z0, self.cell.axon[1].diam)
        self.cell.axon[1].pt3dadd(x0, y0-0.5-self.cell.axon[0].L-self.cell.axon[1].L, z0, self.cell.axon[1].diam)
    ## END set_axon_xyz()
    
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

    def set_recordings(self, dt = neuron.h.dt):

        recordings = {  
            "t": neuron.h.Vector().record(neuron.h._ref_t, dt),
            "soma": neuron.h.Vector().record(neuron.h.cell.soma[0](0.5)._ref_v, dt),
            "axon[0]": neuron.h.Vector().record(neuron.h.cell.axon[0](0.5)._ref_v, dt),
        }
        ## specific section for pyramidal cells
        # if "TTPC" in self.cell_name:
        #     recordings["apic[10]"] =  neuron.h.Vector().record(neuron.h.cell.apic[10](0.5)._ref_v)
        return recordings
    ## END set_recordings()

    def get_all_params(self, recordings = True):
        param = {}
        for k in self.__dict__.keys():
            if type(self.__dict__[k]) is not neuron.hoc.HocObject:
                param[k] = self.__dict__[k]
            elif k == "iclamp":
                param["iclamp"] = {"amp":self.iclamp.amp, "del": self.iclamp.delay, "dur": self.iclamp.dur}

        
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

        return stim, synapse, nc, syn_events
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
               parts = ["soma", "dend", "apic", "axon" ], colors = ["lightgray", "firebrick","coral", "royalblue", ], ## soma, dendrites, apic, axon
               savefig= False, figname=None, show_legend=False):

        (ix, iy) = (0,1) if (plan == "xy")or(plan == "yx") else (2,1) if (plan == "zy")or(plan == "yz") else (0,2)
        Xmin, Xmax = self._get_boundaries()
        xmin, ymin = Xmin[ix], Xmin[iy]
        l = [Xmax[i] - Xmin[i] for i in range(3)]
        max_length = max([l[ix],l[iy]])

        fig, ax = plt.subplots()
        # colors = ["red", "royalblue","navy", "coral", ] ## soma, dendrites, apical, axon
        if type(colors) == str:
            print("colors entered as a string")
            colors = [ colors for i in range(len(parts))]

        from matplotlib.patches import Circle
        for sec in self.cell.all:

            for i,part in enumerate(parts):
                if part in sec.name(): 
                    c = colors[i]
                    p = part  

            x,y,z,d = [],[],[],[]
            for i in range(sec.n3d()):
                x.append(sec.x3d(i))
                y.append(sec.y3d(i))
                z.append(sec.z3d(i))
                d.append(sec.diam3d(i))
            
            X = np.array([x,y,z,d])
            if p=="soma":
                fac = 0.75 if linewidth is None else 1.*linewidth
                [ ax.add_patch( Circle((X[ix,i], X[iy,i]), fac*X[3,i], 
                                        fc = c,lw= 0.2 , ec="gray", zorder=5) ) 
                 for i in range(sec.n3d()) ] 
                # i=0
                # ax.add_patch( Circle((X[ix,i], X[iy,i]), 0.75*sec.diam, 
                #              fc = c,lw= 0.2 , ec="gray", zorder=5) ) 

            elif p == "axon" and len(self.cell.axon)==2:
                lw = 2 if linewidth is None else 4*linewidth
                n = 0#self.cell.soma[0].n3d()-1
                x0, y0, z0 = self.cell.soma[0].x3d(n), self.cell.soma[0].y3d(n), self.cell.soma[0].z3d(n)
                L, d=self.cell.axon[0].L+self.cell.axon[1].L , self.cell.axon[0].diam
                X = np.array([[x0,y0-0.5,z0,d],[x0,y0-0.5-L,z0,d]]).T
                ax.plot(X[ix], X[iy], color=c, linewidth=lw)
            else:
                lw =0.5 if linewidth is None else linewidth
                ax.plot(X[ix], X[iy], color=c, linewidth=lw)


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
            plt.legend(custom_lines, ["Soma", "Dendrites","Apical", "Axon"],
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

    def run_control_simulation(self, dt = 30e3):
        ## Control or SHAM simulation running with segmented output for state restoring
        ts = time.time()
        dt = dt ## each 30s = 30e3 ms
        time2save = np.arange(dt, self.tstop+1, dt)
        self.initialize(-65)
        # cell.get_steady_state()
        fn_dir = f"{self.resdir}/States/"
        for t in time2save:
            neuron.h.continuerun(t)
            cell.save_state(f"{fn_dir}/SHAM_state_{t//1e3}s_cutaxon.pkl", savesyn =True)
            print(f"Done simulating {t//1e3} s and saving current voltages...")

        self.save_results(f"{self.resdir}/control_dur{t//1e3}s_cutaxon")     

        pass
    ## END run_control_simulation(self)

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

## END BBPcell


# from listcells import *
# Params = {"temperature": 36, "tstop":50,"cell_name": L23[0]  , 
#         "use_BBPsynapses":False,  "fullAxon":True,
#         "iclamp":{"amp":0.7, "dur":100, "del":1}
#         }
# cell = BBPcell(**Params)

# # # %%
# import blenderspike_py
# # cell.set_axon_xyz()
# cell_list = cell.cell.soma[0].wholetree()
# cell_list.pop()
# cell_list.pop()

# recorder = blenderspike_py.CellRecorder(cell_list)
# from blenderneuron import neuronstart

# cell.initialize()
# cell.get_steady_state()
# cell.run()

# recorder.save_pickle(filename=f"test_blenderspike_{cell.cell_name}_withaxon.pkl", FRAME_NUM=300)

# plt.figure()
# plt.plot(cell.recordings["t"].as_numpy(),cell.recordings["soma"].as_numpy())
# plt.show()


# cell.plot2D(savefig=True, figname="L23[0]_black", linewidth = 1.5,
#              show_scale=False ,show_legend=False,  colors = "k", 
#             )

# #%%
# with open("results/cell_15/synaptic_inputs480s_V2.pkl", "rb") as f:
#     events_list = pkl.load(f)
# e = [i.as_numpy() for i in events_list]
# e = [neuron.h.Vector(i[i>50]-50) for i in e]
# events_list = e
# neuron.h.dt = 0.002
# "L4_BP_bNAC219_1"
# Params = {"temperature": 36, "tstop":500,
#         "cell_name":  "L4_BP_bNAC219_1",  
#         "steady_state": False,"save_res": False,
#         "use_BBPsynapses":True, "syn_freq": None, "syn_events":None,
#         "tacs_params":{"tacs_amp": 10, "tacs_freq": 10, "tacs_phi":0,  
#         "tacs_dur": 500, "tacs_del":50 }
#         }
#         # "iclamp":{"amp":0.2, "dur":80, "del":10}
#         # }
# cell = BBPcell(**Params)
# # # cell.plot2D()
# # # cell.synapses.desactivate_synapses()
# cell.set_tACS()
# cell.init(-68.5)
# cell.run()
# cell.get_spike_timing()
# cell.save_results(recordings=False)
# cell.save_recordings()
# plt.plot(cell.recordings["t"].as_numpy(),cell.recordings["soma"].as_numpy())
# plt.show()

# plt.plot(cell.recordings["t"].as_numpy(),cell.recordings["tacs_waveform"])
# plt.show()
# plt.switch_backend("qtagg")
# cell.plot2D()


#%%
# if __name__== "__main__":
   
#     Params = {"temperature": 36, "tstop":360e3,
#             "cell_name": "L5_TTPC2_cADpyr232_1" , 
#             "use_BBPsynapses":True, "syn_freq": None, "syn_events":None,
#             "save_res" : True, 
#             "tacs_params":None
#             }
#             # "iclamp":{"amp":0.2, "dur":80, "del":10}
#             # }
#     cell = BBPcell(**Params)
#     # cell.run()

#     ## SHAM run
#     ts = time.time()
#     dt = 60e3 ## each 30s = 30e3 ms
#     time2save = np.arange(dt, cell.tstop+1, dt)
#     cell.init(-65)
#     # cell.get_steady_state()
#     fn_dir = f"{cell.resdir}/States/"
#     for t in time2save:
#         neuron.h.continuerun(t)
#         cell.save_state(f"{fn_dir}/SHAM_state_{t//1e3}s_cutaxon.pkl", savesyn =True)
#         print(f"Done simulating {t//1e3} s and saving current voltages...")

#     cell.save_results(f"{cell.resdir}/control_dur360s_cutaxon")
    # print(f'Took {round(time.time() - ts)} s to simulate')

    # cell.tacs_params["tacs_amp"] = 2.0
    # cell.set_tACS()
    # cell.init()
    # cell.get_steady_state()
    # # cell.restore_state("test_state.pkl")
    # cell.run()


