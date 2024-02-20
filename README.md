# tACS Single Cell SSCx
This repository is the code developped to perform the simulations associated with Gaugain G, Al Harrach M, Yochum M, Wendling F, Bikson M, Modolo J, et al. Frequency-dependent phase entrainment of cortical cell types during tACS: Converging modeling evidence. https://doi.org/10.1101/2023.12.15.571874 [[1]].
This code was developp to ease the use of single cell models reconstructed and published by the Blue Brain Project from rats somato-sensory cortex (SSCx), which are available in [NMCPortal](https://bbp.epfl.ch/nmc-portal/welcome.html). To use the code, dowload first the cell to model from this link of the archive containing all the available cells.

## BBPcell class
The class to pilot single cell models is coded in `BBPCell.py`, and functionnalities can be added as needed. Notably, you can instanciate a cell with or without the set of synapses reconstructed by the BBP, add a curent clamp, add extracellular stimulation (which was the aim of the paper). Also, steady state can be computed to initiate the state of the cell to rest (without fluctuation in the membrane voltage in all teh compartments, which also requires the steady state of all ionic mechanisms).
Here is a common syntax to instanciate an object of this class:
`
tacs_params = { "tacs_amp": 1, "tacs_freq": 10, 
                "tacs_phi":0, "tacs_dur": 120e3, "tacs_del": 1e3}
Params = {"temperature": 36, "tstop":121e3,"cell_name": ''  , 
        "use_BBPsynapses":True,  "fullAxon":False, "save_res":False,
        "tacs_params": tacs_params}

cell = BBPcell(**Params)
cell.set_tACS()
`
To run the simulation you can simply add:
`
cell.initialize(-65)
cell.run()
`
Alternativelly, if one want to initiate the simulation at the steady state (see 1_polarization_length as instance):
`
cell.initialize()
cell.get_steady_state()
cell.run()
`

### AberraCell class
Alternativelly, this class allows to use the model built in [[2]]. This class was less developped since the study chose to use the initial model from which were adapted myelinated axonal mechanisms. It also requires to download the appropriate models into the `aberra_cells` directory (you can found them [here](https://github.com/Aman-A/TMSsim_Aberra2019/tree/master/nrn/cells) ). The instanciation of the model uses the script originally developped in [[2]] by A. Aberra (in the `aberra_scripts` folder).


### List of used cells and associated synaptic parameters
The list of considered cells in the study is stored as Python `List` to ease instanciation (long and more complex name) and run simulation for particullar cell type. The synaptic weight scaling parameters are also sotred in the `alpha10Hz` list to generate the close-to-10Hz activity. 

### Main: running multiple simulations
Since different stimulation frequencies and amplitudes were used in this work, **MPI** and its Python API **mpi4pi** were used to speeded up the process and launch simulations in parallel (one cell instanciated and simulated by core). The `main.py` script was used for this purpose. This script calls `run_tacs.py` with required parameters (*cell_name*, *tacs_amp*, or *tacs_freq* and so on). 

### BBPresults and post-processing
The script to read and process the saved results from simulation is available in `BBPresults.py`. I chose to convert all the data in dataframe for further use.

### Figures and Table
The script used the figures related to the anatomically detailed models in [[1]] are attached. The Table S2 of the supplementary can be constructed and exported as a `.csv` file with the script `table1.py`.



## References

Gaugain G, Al Harrach M, Yochum M, Wendling F, Bikson M, Modolo J, et al. Frequency-dependent phase entrainment of cortical cell types during tACS: Converging modeling evidence. https://doi.org/10.1101/2023.12.15.571874 
Aberra AS, Peterchev AV, Grill WM. Biophysically realistic neuron models for simulation of cortical stimulation. J Neural Eng 2018;15:066023. https://doi.org/10.1088/1741-2552/aadbb1.


