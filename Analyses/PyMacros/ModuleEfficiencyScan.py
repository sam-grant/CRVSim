'''
Samuel Grant 2024
Plot the results of the PE threshold scan 
'''

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import Utils as ut
import CoincidenceConditions as cc

# May not be strictly accurate.
def GetRatioUncertainty(dataAll):
  n_failures = np.sum(dataAll[" Failures"])
  n_total = np.sum(dataAll[" Total"])
  ratio = n_failures / n_total
  se_ratio = np.sqrt(ratio * (1 - ratio) / n_total)
  return se_ratio

def PlotGraphOverlay(graphs_, title=None, xlabel=None, ylabel=None, ymin=1, ymax=-1, labels_=[], fout="scatter.png", effLine=False, log=False, includeBlack=False, NDPI=300):
    
    # Create figure and axes
    fig, ax = plt.subplots()

    # Iterate over each pair of xy lists
    for i, (label, data_) in enumerate(graphs_.items()):

        x = data_[0]
        xerr = data_[1]
        y = data_[2]
        yerr = data_[3]

         # Plot scatter with error bars
        if len(xerr)==0: xerr = [0] * len(x) # Sometimes we only use yerr
        if len(yerr)==0: yerr = [0] * len(y) # Sometimes we only use yerr

        colour = ut.colours[i+1]
        if includeBlack: colour = ut.colours[i]

        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=colour, markersize=4, ecolor=colour, capsize=2, elinewidth=1, linestyle='-',label=label)

    if log: 
        ax.set_yscale("log")

    # Set title, xlabel, and ylabel
    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel(xlabel, fontsize=13, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=13, labelpad=10)

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=13)  
    ax.tick_params(axis='y', labelsize=13)  

    if (ymax > ymin):
        ax.set_ylim(ymin, ymax)

    # Add a line at 99.99% efficiency
    if effLine:
        ax.text(100, 1.2e-4, "99.99% efficiency", color="gray") #, transform=ax.transAxes, ha='right')
        ax.axhline(y=1e-4, color='gray', linestyle='--')

    ax.legend(loc="best", frameon=False, fontsize=13) # , markerscale=5)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(): 

    dataset = "MDC2020ae"
    # PEs_ = list(range(10, 151, 2)) 
    # PEs_ = np.arange(10, 132.5, 2.5)
    PEs_ = np.arange(10.0, 135.0, 5.0)
    layers_ = [2, 3] 

    # Create new dictionaries with independent keys and values
    # ineff_ = {key: value.copy() for key, value in data_.items()}
    # eff_ = {key: value.copy() for key, value in data_.items()}

    # Tyler's data
    finNameWB = f"../Txt/Wideband/Tyler2024.csv"
    dataWB = pd.read_csv(finNameWB)
    wb_thresholds = dataWB["Thresholds"]
    wb_ineff = dataWB["Observed 3/4 inefficiency"]
    wb_ineff_err = dataWB["Error bars, observed 3/4 inefficiency"]

    # shared_thresholds = []
    # shared_wb_ineff = []
    # shared_wb_ineff_err = []

    # print(wb_thresholds)

    for layer in layers_: 

        ineff_ = { 
            # What's the deal with this structure? 
            # x, xerr, y, yerr
            "All" : [ PEs_, [], [], [] ] 
            ,"Muons" : [ PEs_, [], [], [] ] 
            ,"Non-muons" : [ PEs_, [], [], [] ]     
        }

        nfailures_ = { 
            "All" : [ PEs_, [], [], [] ]
            ,"Muons" : [ PEs_, [], [], [] ] 
            ,"Non-muons" : [ PEs_, [], [], [] ]     
        }

        for PE in PEs_:
                
            finNameAll = f"../Txt/{dataset}/concatenated/results_all_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameMuons = f"../Txt/{dataset}/concatenated/results_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameNonMuons = f"../Txt/{dataset}/concatenated/results_non_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            
            dataAll = pd.read_csv(finNameAll) 
            dataMuons = pd.read_csv(finNameMuons) 
            dataNonMuons = pd.read_csv(finNameNonMuons) 

            ineff_["All"][2].append(np.sum(dataAll[" Failures"])/np.sum(dataAll[" Total"]))
            ineff_["All"][3].append(GetRatioUncertainty(dataAll))

            ineff_["Muons"][2].append(np.sum(dataMuons[" Failures"])/np.sum(dataMuons[" Total"]))
            ineff_["Muons"][3].append(GetRatioUncertainty(dataMuons))

            ineff_["Non-muons"][2].append(np.sum(dataNonMuons[" Failures"])/np.sum(dataNonMuons[" Total"]))
            ineff_["Non-muons"][3].append(GetRatioUncertainty(dataNonMuons))
    
            nfailures_["All"][2].append(np.sum(dataAll[" Failures"]))
            nfailures_["Muons"][2].append(np.sum(dataMuons[" Failures"]))
            nfailures_["Non-muons"][2].append(np.sum(dataNonMuons[" Failures"]))
        
        print("PEs:", PEs_)
        print("Failures (muons):", nfailures_["Muons"][2])
        
        PlotGraphOverlay(nfailures_, xlabel="Threshold [PEs/layer]", ylabel=f"Failures ({layer}/4 layers)", fout=f"../Images/{dataset}/ThresholdScan/gr_log_overlay_nfailures_{layer}layers.png", log=True)
        PlotGraphOverlay(nfailures_, xlabel="Threshold [PEs/layer]", ylabel=f"Failures ({layer}/4 layers)", fout=f"../Images/{dataset}/ThresholdScan/gr_overlay_nfailures_{layer}layers.png", log=False)
        
        PlotGraphOverlay(ineff_, xlabel="Threshold [PEs/layer]", ymin=5e-7, ymax=1e0, ylabel=f"{layer}/4 inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_overlay_ineff_{layer}layers.png", log=False)
        PlotGraphOverlay(ineff_, xlabel="Threshold [PEs/layer]", ymin=5e-7, ymax=1e0, ylabel=f"{layer}/4 inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_log_overlay_ineff_{layer}layers.png", log=True)
        
        # Comparison with data
        if layer != 3: # Tyler only has 3/4 layers
            continue 

        data_ = { 
            # What's the deal with this structure? 
            # x, xerr, y, yerr
            "Data" : [ PEs_, np.array([]), wb_ineff, wb_ineff_err]   
            ,"Sim (all)" : [ PEs_, np.array([]), ineff_["All"][2], ineff_["All"][3] ]
            ,"Sim (muons)" : [ PEs_, np.array([]), ineff_["Muons"][2], ineff_["Muons"][3] ] 
            ,"Sim (non-muons)" : [ PEs_, np.array([]), ineff_["Non-muons"][2], ineff_["Non-muons"][3] ] 
              
        }

        # print(PEs_)
        # print(ineff_["Muons"][2])
        
        PlotGraphOverlay(data_, xlabel="Threshold [PEs/layer]", ymin=5e-7, ymax=1e0, ylabel=f"{layer}/4 inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_data_sim_overlay_ineff_{layer}layers.png", log=False, effLine=True, includeBlack=True)
        PlotGraphOverlay(data_, xlabel="Threshold [PEs/layer]", ymin=5e-7, ymax=1e0, ylabel=f"{layer}/4 inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_log_data_sim_overlay_ineff_{layer}layers.png", log=True, effLine=True, includeBlack=True)
        
    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    Run() 

    return

if __name__ == "__main__":
    main()
