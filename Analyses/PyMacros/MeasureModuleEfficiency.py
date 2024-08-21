'''
Samuel Grant 2024

Analyse CRV module efficiency from TrkAna in KPP geometry. 

KPP: three tiers, top and bottom tiers are trigger tiers, middle is the measurement tier.
* Two CRV-L-end modules parallel to the DS axis (sector 3).
* Four CRV-T modules perpendicular to the DS axis (sector 1).
* Two CRV-DS modules perpendicular to the DS axis (sector 2).

Procedure: 

Use top and bottom modules to trigger on cosmics, look for the failure in the measurement module. 
 
First pass: trigger on events which have ONE coincidience in the top and bottom tier (2/4 layers hit and 10 PEs per layer), store failures.
Second pass: analyses failures using time window and spatial grouping between trigger modules. 
Third pass: analyses failure using the tracker. 

In each case: 
Count for successful and unsuccessful triggers. 
Calculate the inefficiency, write failules and effiency results to file. 

Notes: 

In future, we can expand to >1 coin in the triggers, pass directly these directly to the grouping algorithm.

Other questions: 
What is the purity of the tracker? 
What is the angular distribtuion that the tracker can handle with straight tracks? 
What is the efficiency of the tracker? 

We are trying to learn how well we can use the tracker to measure the module efficiency in situ. 


'''
# External libraries
import sys
import uproot 
import numpy as np
import awkward as ak
import pandas as pd

# Internal libraries
import Utils as ut
import PrintUtils as pr
import CoincidenceConditions as cc

from Mu2eEAF import ReadData as rd 
from Mu2eEAF import Parallelise as pa

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

# def SanityPlots(data_, recon, filterLevel, foutTag, quiet):
    
#     #TODO: add tracker sanity plots.
#     if not quiet: print("\n---> Making sanity plots")

#     ###################################################################
#     # CRV reco
    
#     sectors_ = ak.flatten(data_["crv"]["crvcoincs.sectorType"]) 
#     t_ = ak.flatten(data_["crv"]["crvcoincs.time"]) 
#     x_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fX"]) * 1e-3 # mm -> m
#     y_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fY"]) * 1e-3 # mm -> m
#     z_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fZ"]) * 1e-3 # mm -> m
#     PEs_ = ak.flatten(data_["crv"]["crvcoincs.PEs"]) 
#     PEsPerLayer_ = ak.flatten(data_["crv"]["crvcoincs.PEsPerLayer[4]"])
#     nHits_ = ak.flatten(data_["crv"]["crvcoincs.nHits"]) 
#     nLayers_ = ak.flatten(data_["crv"]["crvcoincs.nLayers"]) 
#     slopes_ = ak.flatten(data_["crv"]["crvcoincs.angle"])

#     # Coincidence
#     ut.PlotGraphOverlay(graphs_=[(x_[sectors_ == 3], y_[sectors_ == 3]), (x_[sectors_ == 1], y_[sectors_ == 1]), (x_[sectors_ == 2], y_[sectors_ == 2])]
#                         , labels_=["Top", "Middle", "Bottom"], xlabel="x-position [m]", ylabel="y-position [m]", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/gr_XY_{foutTag}.png")
#     ut.PlotGraphOverlay(graphs_=[(z_[sectors_ == 3], y_[sectors_ == 3]), (z_[sectors_ == 1], y_[sectors_ == 1]), (z_[sectors_ == 2], y_[sectors_ == 2])]
#                         , labels_=["Top", "Middle", "Bottom"], xlabel="z-position [m]", ylabel="y-position [m]", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/gr_ZY_{foutTag}.png")
#     ut.Plot1DOverlay(hists_=[t_[sectors_ == 3], t_[sectors_ == 1], t_[sectors_ == 2]], nbins=1000, xmin = np.min(t_), xmax = np.max(t_)
#                      , xlabel="Average hit time [ns]", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_times_{foutTag}.png") 
#     ut.Plot1DOverlay(hists_=[nHits_[sectors_ == 3], nHits_[sectors_ == 1], nHits_[sectors_ == 2]], nbins=41, xmin = -0.5, xmax = 40.5
#                      , xlabel="Number of hits", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_nHits_{foutTag}.png") 
#     ut.Plot1DOverlay(hists_=[nLayers_[sectors_ == 3], nLayers_[sectors_ == 1], nLayers_[sectors_ == 2]], nbins=5, xmin = -0.5, xmax = 4.5
#                      , xlabel="Number of layers hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_nLayers_{foutTag}.png") 
#     ut.Plot1DOverlay(hists_=[slopes_[sectors_ == 3], slopes_[sectors_ == 1], slopes_[sectors_ == 2]], nbins=1000, xmin = -2, xmax = 2
#                      , xlabel="Slope", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_slopes_{foutTag}.png") 

#     # Calculate the sum of each length-4 array in PEsPerLayer_
#     PEsPerLayerSum_ = ak.sum(PEsPerLayer_, axis=-1)
#     PEsDiff_ = PEs_ - PEsPerLayerSum_

#     ut.Plot1DOverlay(hists_=[PEs_[sectors_ == 3], PEs_[sectors_ == 1], PEs_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_PEs_{foutTag}.png")
#     ut.Plot1DOverlay(hists_=[PEsPerLayerSum_[sectors_ == 3], PEsPerLayerSum_[sectors_ == 1], PEsPerLayerSum_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs per layer", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_PEsPerLayer_{foutTag}.png")

#     ut.Plot1DOverlay(hists_=[PEsPerLayer_[:, 0], PEsPerLayer_[:, 1], PEsPerLayer_[:, 2], PEsPerLayer_[:, 3]], nbins=500, xmin = 0, xmax = 500, title="PEs per layer", ylabel="Coincidences", label_=["Layer0", "Layer 1", "Layer 2", "Layer 3"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}//h1_PEsPerLayerAllSectors_{foutTag}.png")
#     ut.Plot1DOverlay(hists_=[PEsDiff_[sectors_ == 3], PEsDiff_[sectors_ == 1], PEsDiff_[sectors_ == 2]], nbins=40, xmin = -0.0002, xmax = 0.0002, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_PEsDiff_{foutTag}.png", logY=True)
#     ut.Plot1DOverlay(hists_=[PEsDiff_], nbins=200, xmin = np.min(PEsDiff_)-1e-4, xmax = np.max(PEsDiff_)+1e-4, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["All modules"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_PEsDiffAll_{foutTag}.png", logY=True, includeBlack=True)

#     # MC 
#     # valid_ = ak.flatten(data_["crv"]["crvcoincsmc.valid"])
#     pdgid_ = ak.flatten(data_["crv"]["crvcoincsmc.pdgId"])
#     primaryE_ = ak.flatten(data_["crv"]["crvcoincsmc.primaryE"]) 

#     ut.Plot1D(primaryE_, nbins=100, xmin = 0, xmax = 1e5, xlabel="Primary energy [unit?]", ylabel="Coincidences", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_primaryE_{foutTag}.png") 

#     # ut.Plot1D(hists_=[valid_[sectors_ == 3], valid_[sectors_ == 1], valid_[sectors_ == 2]], nbins=2, xmin = 0, xmax = 2, xlabel="Valid MC event", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_valid_{foutTag}.png")
    
#     label_dict = {
#         2212: 'proton',
#         211: 'pi+',
#         -211: 'pi-',
#         -13: 'mu+',
#         13: 'mu-',
#         -11: 'e+',
#         11: 'e-',
#         "other": "other"
#         # Add more particle entries as needed
#     }

#     ut.BarChart(data_=pdgid_, label_dict=ut.particleDict, ylabel="Coincidences [%]", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/bar_pdgid_{foutTag}.png", percentage=True)
#     ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/bar_overlay_pdgid_{foutTag}.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

#     # title=filters_[filterLevel] 


#     # Tracker stuff
#     title=filters_[filterLevel] 
#     # sector1Condition = data_["crv"]["crvcoincs.sectorType"] == 1

#     # # CRV-DS with 826.0 length 2,570
#     # crvDS_len = 2570
#     # crvDS_wid = 826
#     ut.Plot2D(x=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)
#                      , y=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
#                      , nbinsX=50, xmin=-4000, xmax=4000, nbinsY=50, ymin=-4000, ymax=4000
#                      # , xbox_= [-crvDS_len/2, crvDS_len/2], ybox_ =  [-crvDS_wid/2, crvDS_wid/2]
#                      , title=title, xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]"
#                      , fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h2_ZX_{foutTag}.png")


#     # Requires singles and track cuts
#     if (filterLevel == 2): 
        
#         # Direct comparison between tracker and CRV1 variables 
    
#         sector1Condition = data_["crv"]["crvcoincs.sectorType"] == 1
#         # oneSector1Condition = ak.count(data_["crv"]["crvcoincs.sectorType"], axis=1) == 1 
    
#         # ut.Plot1D(ak.flatten(data_["trkfit"]["klfit"]["time"], axis=None) - ak.flatten(data_["crv"]["crvcoincs.time"][sector1Condition], axis=None)
#         #           , nbins=250, xmin=-30, xmax=30 #, norm=300, mu=0, sigma=3, fitMin=-30, fitMax=30
#         #           , xlabel="$T_{KKInter} - T_{CRV}$ [ns]", ylabel="Counts", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_deltaT_{foutTag}.png"
#         #           , stats=True) 
    
#         ut.Plot1DWithGaussFit(ak.flatten(data_["trkfit"]["klfit"]["time"], axis=None) - ak.flatten(data_["crv"]["crvcoincs.time"][sector1Condition], axis=None)
#                   , nbins=250, xmin=-30, xmax=30, norm=300, mu=0, sigma=3, fitMin=-30, fitMax=30
#                   , title=title, xlabel="$T_{KKInter} - T_{CRV}$ [ns]", ylabel="Counts", fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h1_fit_deltaT_{foutTag}.png"
#                   , stats=True) 
    
    
#         ut.Plot2D(x=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)
#                   , y=ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fZ"][sector1Condition], axis=None)
#                   , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
#                   , title=title, xlabel="KKInter Z [mm]", ylabel="CRV Z [mm]"
#                   , fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h2_ZZ_singles_track_cuts.png")
    
        
    
#         ut.Plot2D(x=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
#                   , y=ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fX"][sector1Condition], axis=None)
#                   , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
#                   , title=title, xlabel="KKInter X [mm]", ylabel="CRV X [mm]"
#                   , fout=f"../Images/{recon}/Sanity/{filters_[filterLevel]}/h2_XX_singles_track_cuts.png")
        
#     if not quiet: print("Done!")

#     return
    
# ------------------------------------------------
# Coincidence finding in measurement sector
# ------------------------------------------------

# Impose stricter conditions than the default, if desired
def FindCoincidences(data_, coincidenceConditions, quiet): 
    
    # print("\n---> Marking coincidences")
    if not quiet: print("\n---> Marking coincidences in measurement sector.")

    # Get conditions from file
    coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]

    # Measurement sector condition 
    # print("* Measurement sector condition")
    sectors_ = data_["crv"]["crvcoincs.sectorType"]
    sectorCondition = sectors_ == 1
    
    # PE threshold per layer
    # print("* PE threshold condition")
    PEs_ = data_["crv"]["crvcoincs.PEsPerLayer[4]"]
    
    # PE condition per layer
    PECondition = PEs_ >= coincidenceConditions_["PEs"]

    # Number of layers hit 
    # (based on number of coincidences with PEs per layer above threshold, not nLayers which is baked into the reconstruction in mcs)
    # print("* Layers condition")
    layerCondition = ak.count(PEs_[PECondition], axis=2) >= coincidenceConditions_["nLayers"] # This can be adjusted on the level of nts 

    # Hit slope: horizontal / vertical direction 
    # print("* Angle condition")
    angleCondition = (abs(data_["crv"]["crvcoincs.angle"]) <= coincidenceConditions_["maxSlope"]) & (abs(data_["crv"]["crvcoincs.angle"]) >= coincidenceConditions_["minSlope"])

    # Combine conditions to mark per-module coincidences
    coincidenceMask = sectorCondition & layerCondition & angleCondition # PECondition is implicit in the layer condition!

    # Add a new field 'is_coincidence' to mark coincidences
    data_["CRVT_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    data_["sector_condition"] = sectorCondition 
    data_["PE_condition"] = PECondition 
    data_["layer_condition"] = layerCondition 
    data_["angle_condition"] = angleCondition 

    if not quiet: print("Done!")

    return 

# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

# "all", "muons", "non_muons"
def FilterParticles(data_, particle, quiet):
    if not quiet: print(f"\n---> Filtering particles, keeping {particle}")
    
    muonCondition = ak.any((data_["crv"]["crvcoincsmc.pdgId"] == 13) | (data_["crv"]["crvcoincsmc.pdgId"] == -13), axis=1)
    if particle == "all":
        return data_
    elif particle == "muons": 
        return data_[muonCondition] 
    elif particle == "non_muons":
        return data_[~muonCondition] 
    else:
        raise ValueError(f"Particle string {particle} not valid!")

#  Events with ONE coincidence in sectors 2 & 3 
#  AND no more than one coincidence in sector 1 (with default coin conditions) 
def FilterSingles(data_, fail, quiet):
    
    if not quiet: print(f"\n---> Filtering singles") 

    sector1Condition = data_["crv"]["crvcoincs.sectorType"] == 1
    sector2Condition = data_["crv"]["crvcoincs.sectorType"] == 2
    sector3Condition = data_["crv"]["crvcoincs.sectorType"] == 3

    oneOrZeroCoincInMeasurementSector = ak.count(data_["crv"]["crvcoincs.sectorType"][sector1Condition], axis=1) < 2
    oneCoincInSector2Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector2Condition], axis=1) == 1
    oneCoincInSector3Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector3Condition], axis=1) == 1
    
    data_["oneOrZeroCoincInMeasurementSector"] = oneOrZeroCoincInMeasurementSector 
    data_["oneCoinInTriggerSectors"] = (oneCoincInSector2Condition & oneCoincInSector3Condition)
    data_["pass_singles"] = (oneOrZeroCoincInMeasurementSector & oneCoincInSector2Condition & oneCoincInSector3Condition)

    if not quiet: print("Done!")
    
    # Cut on event level
    if not fail: 
        return data_[data_["pass_singles"]]
    else: 
        return data_[~data_["pass_singles"]]
        
# Tracker cuts
def ApplyTrackerCuts(arrays_, fail=False, quiet=False):
    
    if not quiet: print(f"\n---> Applying tracker cuts") 

    # Mark cuts on the track and track fit level
    arrays_["trkfit_KLCRV1"] = ( 
        (arrays_["trkfit"]["klfit"]["sid"] == 200) 
        & (arrays_["trkfit"]["klfit"]["sindex"] == 1) )

    arrays_["trk_bestFit"] = ( 
        (arrays_["trk"]["kl.ndof"] >= 10)
        & (arrays_["trk"]["kl.fitcon"] > 0.1)
        & ((arrays_["trk"]["kl.nactive"]/arrays_["trk"]["kl.nhits"]) > 0.99)
        & (arrays_["trk"]["kl.nplanes"] >= 4)
        & ((arrays_["trk"]["kl.nnullambig"]/arrays_["trk"]["kl.nhits"]) < 0.2) )
    
    arrays_["trkfit_bestFit"] = ( 
        (arrays_["trkfit"]["klkl"]["z0err"] < 1) 
        & (arrays_["trkfit"]["klkl"]["d0err"] < 1) 
        & (arrays_["trkfit"]["klkl"]["thetaerr"] < 0.004)
        & (arrays_["trkfit"]["klkl"]["phi0err"] < 0.001) )

    if not fail: 
        # Create masks
        arrays_["trkfit"] = arrays_["trkfit"][(arrays_["trkfit_bestFit"] & arrays_["trkfit_KLCRV1"])]
        arrays_["trk"] = arrays_["trk"][arrays_["trk_bestFit"]]
    else: 
        # Create masks
        arrays_["trkfit"] = arrays_["trkfit"][~(arrays_["trkfit_bestFit"] & arrays_["trkfit_KLCRV1"])]
        arrays_["trk"] = arrays_["trk"][~arrays_["trk_bestFit"]]

    # Check for a track in the event after cuts.
    trkCut = ak.any(arrays_["trk"]["kl.status"], axis=1, keepdims=False) > 0 
    # Check for a track fit in the event after cuts
    trkFitCut = (
        (ak.count(arrays_["trkfit"]["klfit"]["sid"], axis=-1, keepdims=False) > 0) 
        & (ak.count(arrays_["trkfit"]["klkl"]["z0err"], axis=-1, keepdims=False) > 0) )
    
    # Reset to event level
    trkFitCut = ak.any(trkFitCut, axis=-1, keepdims=False) == True 

    # Both do the same thing, but mark them pass/fail for bookkeeping.
    if not fail: 
        arrays_["pass_track_cuts"] = (trkCut & trkFitCut)
        return arrays_[arrays_["pass_track_cuts"]]
    else: 
        arrays_["fail_track_cuts"] = (trkCut & trkFitCut)
        return arrays_[arrays_["fail_track_cuts"]]

# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_, quiet): 
    if not quiet: print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_["crv"]["crvcoincs.timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   

# Basic trigger condition
def Trigger(data_, fail=False, quiet=False): 

    if not quiet: print(f"\n---> Triggering")
        
    # Enforce trigger condititon 
    triggerCondition = (
        ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
        ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1)
    )

    data_["pass_trigger"] = triggerCondition
    
    if not quiet: print("Done!")

    if not fail: 
        return data_[data_["pass_trigger"]]
    else:
        return data_[~data_["pass_trigger"]]
        

# Needs to come after the initial trigger 
def SuccessfulTriggers(data_, success, quiet):

    successStr = ""
    if success: successStr += "successful"
    else: successStr += "unsuccessful"

    if not quiet: print(f"\n---> Getting {successStr} triggers")

    # Ensure at least one passing coincidence the measurement sector
    # Coincidence conditions are set in FindCoincidences()
    successCondition = ak.any(data_["CRVT_coincidence"], axis=1) 
    
    if not quiet: print("Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

# ------------------------------------------------
#                     Output 
# ------------------------------------------------ 

# For debugging single files only!
def WriteSuccessInfo(successes_, recon, finTag, foutTag, coincidenceConditions, quiet):

    # Define the output file path
    foutNameConcise = f"../Txt/{recon}/successes_concise/{finTag}/successes_concise_{foutTag}.csv" 
    foutNameVerbose = f"../Txt/{recon}/successes_verbose/{finTag}/successes_verbose_{foutTag}.csv" 
    
    if not quiet: print(f"\n---> Writing failure info to:\n{foutNameConcise}\n{foutNameVerbose}", flush=True) 
        
   # Concise form
    with open(foutNameConcise, "w") as fout:
        # Write the header
        fout.write("evtinfo.run,evtinfo.subrun,evtinfo.event\n")
        # Write the events
        for event in successes_:
            fout.write(
                f"{event['evt']['evtinfo.run']}, {event['evt']['evtinfo.subrun']}, {event['evt']['evtinfo.event']}\n"
            )
            
    # Verbose form
    if True: 
        with open(foutNameVerbose, "w") as fout:
            # Write the events
            for event in successes_:
                fout.write(
                    PrintEvent(event)+"\n" 
                )

    return

    
def WriteFailureInfo(failures_dict_, recon, finTag, foutTag, coincidenceConditions, quiet):

    for trkTag, failures_ in failures_dict_.items():
        
        # Define the output file path
        foutNameConcise = f"../Txt/{recon}/failures_concise/{finTag}/failures_concise_{foutTag}_{trkTag}.csv" 
        foutNameVerbose = f"../Txt/{recon}/failures_verbose/{finTag}/failures_verbose_{foutTag}_{trkTag}.csv" 
        
        if not quiet: print(f"\n---> Writing failure info to:\n{foutNameConcise}\n{foutNameVerbose}", flush=True) 
    
        # Concise form
        with open(foutNameConcise, "w") as fout:
            # Write the header
            # TODO: changed this to have no spaces!
            fout.write("evtinfo.run, evtinfo.subrun, evtinfo.event\n")
            # Write the events
            for event in failures_:
                fout.write(
                    f"{event['evt']['evtinfo.run']},{event['evt']['evtinfo.subrun']},{event['evt']['evtinfo.event']}\n"
                )
                
        # Verbose form

        # Setup mask string
        masks_ = ["pass_singles", "pass_trigger", "CRVT_coincidence"]
        if trkTag == "track_cuts":
            masks_ += ["trk_bestFit", "trkfit_bestFit", "trkfit_KLCRV1", "pass_track_cuts"]
        
        with open(foutNameVerbose, "w") as fout:
            # Write the events
            for event in failures_:
                maskStr = ''.join([f"{mask}: {event[f'{mask}']}\n" for mask in masks_])
                fout.write(
                    pr.PrintEvent(event, maskStr)+"\n" 
                )

    return

def WriteResults(data_, successes_, failures_, failures_track_cuts_, recon, finTag, foutTag, quiet):

    foutNameNoTrk = f"../Txt/{recon}/results/{finTag}/results_{foutTag}_no_track_cuts.csv"
    foutNameTrk = f"../Txt/{recon}/results/{finTag}/results_{foutTag}_track_cuts.csv"
    if not quiet: print(f"\n---> Writing results to:\n{foutNameNoTrk}\n{foutNameTrk}")

    # Calculate efficiency 
    # Handle edge cases where the remaining subset has length zero.
    tot = len(data_)
    inefficiency = len(failures_) / tot * 100 if tot else np.nan
    inefficiency_track_cuts = len(failures_track_cuts_) / tot * 100 if tot else np.nan

    outputStr = f"""
    ****************************************************
    Input tag: {finTag}
    Output tag: {foutTag}

    Number of successes: {len(successes_)}
    
    No track cuts:
    Failures: {len(failures_)}
    Inefficiency: {len(failures_)}/{tot} = {inefficiency:.2f}%
    
    With track cuts:
    Failures: {len(failures_track_cuts_)}
    Inefficiency: {len(failures_track_cuts_)}/{tot} = {inefficiency_track_cuts:.2f}%
    ****************************************************
    """

    with open(foutNameNoTrk, "w") as foutNoTrk:
        foutNoTrk.write("Total,Successes,Failures,Inefficiency [%]\n") # no spaces!
        foutNoTrk.write(f"{tot}, {len(successes_)}, {len(failures_)}, {inefficiency}\n")

    with open(foutNameTrk, "w") as foutTrk:
        foutTrk.write("Total,Successes,Failures,Inefficiency [%]\n") # no spaces!
        foutTrk.write(f"{tot}, {len(successes_)}, {len(failures_track_cuts_)}, {inefficiency_track_cuts}\n")

    if not quiet: print(outputStr)

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(file, recon, particle, PE, layer, finTag, quiet): 
    
    # Output strings
    coincidenceConditions = f"{PE}PEs{layer}Layers"
    foutTag = particle + "_" + coincidenceConditions 

    # Get data as a set of awkward arrays
    data_ = ut.GetData(file) # finName)

    # Apply start time cut
    data_ = CutOnStartTime(data_, quiet)

    # Filter particle species
    data_ = FilterParticles(data_, particle, quiet)

    # Enforce basic trigger
    data_ = Trigger(data_, fail=False, quiet=quiet)

    # Singles cut
    data_ = FilterSingles(data_, fail=False, quiet=quiet)

    # Mark coincidences in the measurement sector for the remaining subset
    FindCoincidences(data_, coincidenceConditions, quiet)

    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True, quiet=quiet)
    failures_ = SuccessfulTriggers(data_, success=False, quiet=quiet)

    # Apply track cuts to failures
    failures_track_cuts_ =  ApplyTrackerCuts(failures_, fail=False, quiet=quiet)

    # Write failures to file
    WriteFailureInfo( {"no_track_cuts" : failures_, "track_cuts" : failures_track_cuts_}, recon, finTag, foutTag, coincidenceConditions, quiet) 
                       
    # Write results to file
    WriteResults(data_, successes_, failures_, failures_track_cuts_, recon, finTag, foutTag, quiet) 

    return

from concurrent.futures import ThreadPoolExecutor, as_completed

# ------------------------------------------------
#              Multithread 
# ------------------------------------------------
            
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product

# Multithread on one file with many configs
# def Multithread2(processFunction2, fileName, particles_, layers_, PEs_):
    
#     print("\n---> Starting multithreading...")
    
#     totalConfigurations = len(particles_) * len(layers_) * len(PEs_)
#     completedConfigurations = 0
    
#     with ThreadPoolExecutor() as executor:
        
#         # Generate all combinations of (particle, layer, PE)
#         configurations = product(particles_, layers_, PEs_)
        
#         # Submit tasks to the executor
#         futures = {
#             executor.submit(processFunction2, fileName, particle, layer, PE): (particle, layer, PE)
#             for particle, layer, PE in configurations
#         }

#         # Process results as they complete
#         for future in as_completed(futures):
#             (particle, layer, PE) = futures[future]
#             try:
#                 future.result()  # Retrieve the result or raise an exception if occurred
#                 completedConfigurations += 1
#                 percentComplete = (completedConfigurations / totalConfigurations) * 100
#                 print(f'---> Configuration ({particle}, {layer}, {PE}) processed successfully! ({percentComplete:.1f}% complete)')
#             except Exception as exc:
#                 print(f'---> Configuration ({particle}, {layer}, {PE}) generated an exception!\n{exc}')
                
#     print("\nMultithreading completed!")
#     return

def Multithread(processFunction, fileList_, max_workers=381):
    
    print("\n---> Starting multithreading...\n")

    completedFiles = 0
    totalFiles = len(fileList_)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        # Prepare a list of futures and map them to file names
        futures_ = {executor.submit(processFunction, fileName): (fileName) for fileName in fileList_}
        
        # Process results as they complete
        for future in as_completed(futures_):
            fileName = futures_[future]  # Get the file name associated with this future
            future.result()  # Retrieve the result 
            completedFiles += 1
            percentComplete = (completedFiles / (totalFiles)) * 100
            print(f'\n---> {fileName} processed successfully! ({percentComplete:.1f}% complete)')
            # Exceptions are handled in the processFunction. 
                
    print("\n---> Multithreading completed!")
    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def TestMain():
    
    fileName="/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00038/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000006.root"
    finTag = fileName.split('.')[-2] 

    with uproot.open(fileName) as file:
        Run(file, recon="MDC2020ae", particle="all", PE="100", layer="3", finTag=finTag, quiet=False)

    return
    
    
def main():

    # TestMain()
    # return

    #########################################################

    defname = "nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.root"
    recon = "MDC2020ae"
    
    particles_ = ["all", "muons", "non_muons"]
    layers_ = [3, 2]
    PEs_ = np.arange(10, 135, 5)
    quiet = True

    # def processFunction(fileName):
    #     # Always open the file in the processFunction 
    #     file = rd.ReadFile(fileName, quiet)
    #     finTag = fileName.split('.')[-2] 
    #     # Scan PE thresholds
    #     for PE in PEs_: 
    #         # Scan particles
    #         for particle in particles_: 
    #             # Scan layers
    #             for layer in layers_: 
    #                 outputStr = (
    #                     "\n---> Running with:\n"
    #                     f"fileName: {fileName}\n"
    #                     f"recon: {recon}\n"
    #                     f"particle: {particle}\n"
    #                     f"PEs: {PE}\n"
    #                     f"layers: {layer}/4\n"
    #                     f"finTag: {finTag}\n"
    #                     f"quiet: {quiet}\n"
    #                 )
    #                 if not quiet: print(outputStr) 
    #                 try:
    #                     Run(file, recon, particle, PE, layer, finTag, quiet)
    #                 except Exception as exc:
    #                     print(f'---> Exception!\n{row}\n{exc}')
    #                 # Uncomment to test
    #                 # return
    #     return

    # fileList_ = rd.GetFileList(defname) #[:2]
    
    # print(f"---> Got {len(fileList_)} files.")

    # Multithread(processFunction, fileList_)

    # Resubmission of failures.
    def processFunction(fileName):

        # Get failed jobs
        failedJobsFile = "../Txt/MDC2020ae/FailedJobs/failures.csv"
        df_failedJobs_ = pd.read_csv(failedJobsFile)
        
        # Always open the file in the processFunction 
        file = rd.ReadFile(fileName, quiet)
        finTag = fileName.split('.')[-2] 
    
        df_failedJobs_ = df_failedJobs_[df_failedJobs_["Tag"] == finTag]
    
        # helper to get filter level from dict
        def get_key(d, val):
          keys = [k for k, v in d.items() if v == val]
          return keys[0] if keys else None
            
        # Submit failed configs
        for index, row in df_failedJobs_.iterrows():
    
            # print(row["PEs"], row["PEs"], row["Layer"], row["Particle"], row["Filter"]) 
            particle = row["Particle"]
            PE = row["PEs"]
            layer = row["Layer"]
    
            outputStr = (
                "\n---> Running with:\n"
                f"fileName: {fileName}\n"
                f"recon: {recon}\n"
                f"particle: {particle}\n"
                f"PEs: {PE}\n"
                f"layers: {layer}/4\n"
                f"finTag: {finTag}\n"
                f"quiet: {quiet}\n"
            )
    
            if not quiet: print(outputStr) 
    
            try:
                Run(file, recon, particle, PE, layer, finTag, quiet)
            except Exception as exc:
                print(f'---> Exception!\n{row}\n{exc}')
    
            # Testing
            # return
            
        return

    # fileName="nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root"
    # processFunction(fileName)

    # return

    fileList_ = rd.GetFileList(defname) #[:2]
    
    print(f"---> Got {len(fileList_)} files.")

    Multithread(processFunction, fileList_)
    
    # TODO: wrap this in a proper function.
    # # Get failed jobs
    
    # failedJobsFile = "../Txt/MDC2020ae/FailedJobs/failures.csv"
    # df_failedJobs_ = pd.read_csv(failedJobsFile)

    # # Set of tags
    # tags_ = list(set(df_failedJobs_["Tag"]))
    
    # # Extract the tag from the file name
    # def ExtractTag(fileName):
    #     parts = fileName.split('.')
    #     if len(parts) > 1:
    #         return parts[-2]
    #     return None

    # # Filter and sort files based on tags
    # fileList_ = sorted(
    #     [file for file in fileList_ if ExtractTag(file) in tags_]
    #     , key=lambda file: tags_.index(ExtractTag(file))
    # )

    

    # Multithread5(processFunction3, fileList_) 
    # Multithread5(processFunction4, fileList_)
    # Multithread6(processFunction4, fileList_)
    
    # Multithread3(processFunction2, fileList_[:2], particles_, layers_, PEs_) 
    
    # Multithread4(processFunction3, ["nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root"]) 
    # Multithread3(processFunction2, fileList_[:20], particles_, layers_, PEs_) 
    # Multithread2(processFunction2, fileName, particles_, layers_, PEs_) 

    return

if __name__ == "__main__":
    main()
