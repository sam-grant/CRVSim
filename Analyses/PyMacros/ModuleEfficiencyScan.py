'''
Samuel Grant 2024
Plot the results of the PE threshold scan 
'''

import sys
import numpy as np
import pandas as pd

import Utils as ut
import CoincidenceConditions as cc

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(): 

    PEs_ = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]
    layers_ = [2]

    # finNames_ = []

    ineff_ = { 
        "All" : [ PEs_, [], [], [] ]
        ,"Muons" : [ PEs_, [], [], [] ] 
        ,"Non-muons" : [ PEs_, [], [], [] ]     
    }

    nfailures_ = { 
        "All" : [ PEs_, [], [], [] ]
        ,"Muons" : [ PEs_, [], [], [] ] 
        ,"Non-muons" : [ PEs_, [], [], [] ]     
    }

    # Create new dictionaries with independent keys and values
    # ineff_ = {key: value.copy() for key, value in data_.items()}
    # eff_ = {key: value.copy() for key, value in data_.items()}
    for layer in layers_: 

        for PE in PEs_:

            finNameAll = f"../Txt/reprocessed/concatenated/results_all_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameMuons = f"../Txt/reprocessed/concatenated/results_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameNonMuons = f"../Txt/reprocessed/concatenated/results_non_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"


            # eff_["All"][1].append(pd.read_csv(finNameAll)[" Efficiency [%]"].mean())
            # eff_["Muons"][1].append(pd.read_csv(finNameMuons)[" Efficiency [%]"].mean())
            # eff_["Non-muons"][1].append(pd.read_csv(finNameNonMuons)[" Efficiency [%]"].mean())
            
            dataAll = pd.read_csv(finNameAll) # [" Inefficiency [%]"]
            dataMuons = pd.read_csv(finNameMuons) # [" Inefficiency [%]"]
            dataNonMuons = pd.read_csv(finNameNonMuons) # [" Inefficiency [%]"]

            ineff_["All"][2].append(np.mean(dataAll[" Inefficiency [%]"]))
            ineff_["Muons"][2].append(np.mean(dataMuons[" Inefficiency [%]"]))
            ineff_["Non-muons"][2].append(np.mean(dataNonMuons[" Inefficiency [%]"]))

            nfailures_["All"][2].append(np.sum(dataAll[" Failures"]))
            nfailures_["Muons"][2].append(np.sum(dataMuons[" Failures"]))
            nfailures_["Non-muons"][2].append(np.sum(dataNonMuons[" Failures"]))


            # print(dataAll)
            # data_["muons"].append(pd.read_csv(finNameMuons))
            # data_["non_muons"].append(pd.read_csv(finNameNonMuons))

            # data_["all"].append(pd.read_csv(finNameAll))
            # data_["muons"].append(pd.read_csv(finNameMuons))
            # data_["non_muons"].append(pd.read_csv(finNameNonMuons))
    # 
        # print(eff_)
        # print(ineff_)
        # ineff_ = {}
        ut.PlotGraph(nfailures_["Non-muons"][0], nfailures_["Non-muons"][2], title="Non-muons", xlabel="PEs per layer threshold", ylabel="Failures", fout=f"../Images/reprocessed/ThresholdScan/gr_overlay_nfailures_non_muons_{layer}layers.png")
        ut.PlotGraphOverlay2(nfailures_, xlabel="PEs per layer threshold", ylabel="Failures", fout=f"../Images/reprocessed/ThresholdScan/gr_log_overlay_nfailures_{layer}layers.png", log=True)
        ut.PlotGraphOverlay2(ineff_, xlabel="PEs per layer threshold", ylabel="Inefficiency [%]", fout=f"../Images/reprocessed/ThresholdScan/gr_log_overlay_ineff_{layer}layers.png", log=True)
        ut.PlotGraphOverlay2(nfailures_, xlabel="PEs per layer threshold", ylabel="Failures", fout=f"../Images/reprocessed/ThresholdScan/gr_overlay_nfailures_{layer}layers.png", log=False)
        ut.PlotGraphOverlay2(ineff_, xlabel="PEs per layer threshold", ylabel="Inefficiency [%]", fout=f"../Images/reprocessed/ThresholdScan/gr_overlay_ineff_{layer}layers.png", log=False)



    # ../Txt/reprocessed/concatenated/failures_concise_non_muons_20PEs3Layers_one_coincidence_per_trigger_sector.csv


    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    Run() 

    return

if __name__ == "__main__":
    main()
