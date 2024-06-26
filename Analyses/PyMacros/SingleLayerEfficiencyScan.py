'''
Samuel Grant 2024
Scan PE thresholds get the module and single layer efficiency for the test module

This should be in a standalone directory along with the write script. 

'''

import sys
import numpy as np
import h5py

import Utils as ut

def MakeListsPerLayerDict(): 

    dict_ = {}
    for layer in range(4):
        dict_[f"Layer{layer}"] = []

    return dict_


# ------------------------------------------------
#                      Run
# ------------------------------------------------

def Run(particle, reproc, module):

    tag = "001205_00000000" # loop through these 

    finName = f"../Txt/{reproc}/PEsPerLayer/{tag}/PEsPerLayer_{particle}.h5"

    data_ = {}

    # Open the HDF5 file for reading
    with h5py.File(finName, "r") as fin:
        # Iterate over the groups in the HDF5 file
        for sector in fin.keys():
            # Initialize dictionary for the sector
            data_[sector] = {}
            # Get the dataset within the group
            dataset_ = fin[sector]["PEsPerLayer"][:]
            # Store the dataset in the dictionary
            data_[sector]["PEsPerLayer"] = dataset_

    # Filter module 
    data_ = data_["Sector1"]

    # Vertically slice the array by layer
    data_ = [data_['PEsPerLayer'][:, i] for i in range(4)]

    labels_ = [] 
    [labels_.append(f"Layer {layer}") for layer in range(4)]

    # Now plot 
    # ut.Plot1DOverlay(hists_ = data_ , nbins=125, xmin=0, xmax=125, title=f"Middle module, {particle}", xlabel="PEs per layer", ylabel="Hits", label_ = labels_, fout=f"../Images/{reproc}/ThresholdScan/h1_overlay_PEsPerLayer_module{module}_{particle}.png") 

    # Scan threshold can get the single layer efficiency 
    thresholds_ = np.arange(10, 126, 2)

    # Sanity histograms
    # scanData_ = MakeListsPerLayerDict()

    # Efficiency 
    eff_ = MakeListsPerLayerDict()

    # Efficiency error
    effErr_ = MakeListsPerLayerDict()

    for layer, data in enumerate(data_):
        total = len(data)

        effList = []
        effErrList = []
        
        for threshold in thresholds_:
            data = data[data >= threshold]
            above = len(data)
            below = total - above # inefficiency! 
            eff = (below / total) # * 100
            # This doesn't really makes since above is a subset of total... 
            effErr = eff * np.sqrt( (np.sqrt(above)/above)**2 + (np.sqrt(total)/total)**2)

            effList.append(eff)
            effErrList.append(effErr)


        eff_[f"Layer{layer}"] = effList 
        effErr_[f"Layer{layer}"] = effErrList 

    # Plot
    # labels_ = []
    # [labels_.append(f"{threshold} PEs") for threshold in thresholds_]
    # ut.Plot1DOverlay(hists_ = scanData_["Layer0"], nbins=125, xmin=0, xmax=125, title=f"Middle module, layer 0, {particle}", xlabel="PEs per layer", ylabel="Hits", label_ = labels_, fout=f"../Images/{reproc}/ThresholdScan/h1_overlay_PEsPerLayerThresholds_module{module}_{particle}_layer0.png") 

    # ut.PlotGraphErrors(x=thresholds_, y=)

    # Format data for graphs again 

    graphs_ = {

        "Layer 0" :  [thresholds_, [], eff_["Layer0"], effErr_["Layer0"] ]
        ,"Layer 1" :  [thresholds_, [], eff_["Layer1"], effErr_["Layer2"] ]
        ,"Layer 2" :  [thresholds_, [], eff_["Layer2"], effErr_["Layer2"] ]
        ,"Layer 3" :  [thresholds_, [], eff_["Layer3"], effErr_["Layer3"] ]

    }


    ut.PlotGraphOverlay2(graphs_=graphs_, title=f"Middle module, {particle}", xlabel="PEs per layer threshold", ylabel="Single layer inefficiency", fout=f"../Images/{reproc}/ThresholdScan/gr_overlay_single_layer_ineff_per_threshold_module{module}_{particle}.png") 
    ut.PlotGraphOverlay2(graphs_=graphs_, title=f"Middle module, {particle}", xlabel="PEs per layer threshold", ylabel="Single layer inefficiency", fout=f"../Images/{reproc}/ThresholdScan/gr_log_overlay_single_layer_ineff_per_threshold_module{module}_{particle}.png", log=True) 



    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    # Take command-line arguments
    # finName = sys.argv[1] if len(sys.argv) > 1 else "../Txt/reprocessed/PEsPerLayer/001205_00000000/PEsPerLayer_all.h5" 
    particle = sys.argv[1] if len(sys.argv) > 1 else "all"
    reproc = sys.argv[2] if len(sys.argv) > 2 else "reprocessed" 
    module = sys.argv[3] if len(sys.argv) > 3 else 1

    print("\n--->Running with inputs:\n")
    print("\tparticle:", particle)
    print("\treproc:", reproc)
    print("\tmodule:", module)

    Run(particle=particle, reproc=reproc, module=module) 

    return

if __name__ == "__main__":
    main()
