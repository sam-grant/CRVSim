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
#                       Run
# ------------------------------------------------

def Run(): 

    dataset = "MDC2020ae"
    PEs_ = list(range(10, 151, 2)) 
    layers_ = [3] # 2, 3] 

    # Create new dictionaries with independent keys and values
    # ineff_ = {key: value.copy() for key, value in data_.items()}
    # eff_ = {key: value.copy() for key, value in data_.items()}

    # Tyler's data
    finNameWB = f"../Txt/Wideband/Tyler2024.csv"
    dataWB = pd.read_csv(finNameWB)
    wb_thresholds = dataWB["Thresholds"]
    wb_ineff = dataWB["Observed 3/4 inefficiency"]
    wb_ineff_err = dataWB["Error bars, observed 3/4 inefficiency"]

    shared_thresholds = []
    shared_wb_ineff = []
    shared_wb_ineff_err = []

    print(wb_thresholds)

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

        # count = 0
        count = 0

        for PE in PEs_:

            if PE not in np.array([wb_thresholds]): 
                continue
            else:
                print(PE)
                shared_wb_ineff.append(wb_ineff[count])
                shared_wb_ineff_err.append(wb_ineff_err[count])
                shared_thresholds.append(PE)
                count+=1

            # print(PE)
            # print(wb_thresholds[count])
            # print()
            # count += 1
            # continue

            # if (PE != wb_thresholds[count]): 
            #     # count+=1
            #     continue
            # else:
            #     print(PE)
            #     print(wb_thresholds[count])
            #     print()
            #     count+=1

            # continue
                
            finNameAll = f"../Txt/{dataset}/concatenated/results_all_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameMuons = f"../Txt/{dataset}/concatenated/results_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            finNameNonMuons = f"../Txt/{dataset}/concatenated/results_non_muons_{PE}PEs{layer}Layers_one_coincidence_per_trigger_sector.csv"
            
            dataAll = pd.read_csv(finNameAll) 
            dataMuons = pd.read_csv(finNameMuons) 
            dataNonMuons = pd.read_csv(finNameNonMuons) 

            ineff_["All"][2].append(np.sum(dataAll[" Failures"])/np.sum(dataAll[" Total"]))
            ineff_["Muons"][2].append(np.sum(dataMuons[" Failures"])/np.sum(dataMuons[" Total"])) 
            ineff_["Non-muons"][2].append(np.sum(dataNonMuons[" Failures"])/np.sum(dataNonMuons[" Total"])) 

            nfailures_["All"][2].append(np.sum(dataAll[" Failures"]))
            nfailures_["Muons"][2].append(np.sum(dataMuons[" Failures"]))
            nfailures_["Non-muons"][2].append(np.sum(dataNonMuons[" Failures"]))

            # count += 1

        # ut.PlotGraphOverlay2(nfailures_, title=f"{layer}/4 layers", xlabel="PEs per layer threshold", ylabel="Failures", fout=f"../Images/{dataset}/ThresholdScan/gr_log_overlay_nfailures_{layer}layers.png", log=True)
        # ut.PlotGraphOverlay2(ineff_, title=f"{layer}/4 layers", xlabel="PEs per layer threshold", ylabel="Inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_log_overlay_ineff_{layer}layers.png", log=True)
        # ut.PlotGraphOverlay2(nfailures_, title=f"{layer}/4 layers", xlabel="PEs per layer threshold", ylabel="Failures", fout=f"../Images/{dataset}/ThresholdScan/gr_overlay_nfailures_{layer}layers.png", log=False)
        # ut.PlotGraphOverlay2(ineff_, title=f"{layer}/4 layers", xlabel="PEs per layer threshold", ylabel="Inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_overlay_ineff_{layer}layers.png", log=False)
        
        # if (layer==3):
        print(ineff_["All"][2])
        print(wb_ineff)
        data_ = { 
            # What's the deal with this structure? 
            # x, xerr, y, yerr
            "Sim" : [ shared_thresholds, [], ineff_["All"][2], [] ] 
            ,"Data" : [ shared_thresholds, [], shared_wb_ineff, shared_wb_ineff_err]     
        }
        ut.PlotGraphOverlay2(data_, title=f"3/4 layers", xlabel="PEs per layer threshold", ylabel="Inefficiency", fout=f"../Images/{dataset}/ThresholdScan/gr_data_sim_overlay_ineff_3layers.png", log=True)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    Run() 

    return

if __name__ == "__main__":
    main()
