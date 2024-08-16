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
import CoincidenceConditions as cc

from Mu2eEAF import ReadData as rd 
from Mu2eEAF import Parallelise as pa
    
# ------------------------------------------------
#                     Input 
# ------------------------------------------------ 

#TODO: remove this and use the version in Utils.py 

evtBranchNames_ = [ 
                # Event info
                "evtinfo.run"
                , "evtinfo.subrun"
                , "evtinfo.event"
]

crvBranchNames_ = [
                # Coincidences 
                "crvcoincs.sectorType" 
                , "crvcoincs.pos.fCoordinates.fX" 
                , "crvcoincs.pos.fCoordinates.fY" 
                , "crvcoincs.pos.fCoordinates.fZ"
                , "crvcoincs.time" 
                , "crvcoincs.timeStart"
                , "crvcoincs.PEs" 
                , "crvcoincs.nHits" 
                , "crvcoincs.nLayers" 
                , "crvcoincs.PEsPerLayer[4]"
                , "crvcoincs.angle" 
                # Coincidences (truth)
                , "crvcoincsmc.valid" 
                , "crvcoincsmc.pdgId" 
]

trkBranchNames_ = [
                # Tracks
                "kl.status"
                , "kl.nactive"
                , "kl.nhits"
                , "kl.nplanes"
                , "kl.nnullambig"
                , "kl.ndof"
                , "kl.fitcon"
]

trkFitBranchNames_ = [
                # Track fits (vector of vector like objects)
                "klfit"
                , "klkl"
] 

allBranchNames_ = { "evt" : evtBranchNames_
                   ,"crv" : crvBranchNames_
                   ,"trk" : trkBranchNames_ 
                   ,"trkfit" : trkFitBranchNames_
                  }

# headers_ = ["evt", "crv", "trk", "trkfit"]
# allBranchNames_ = evtBranchNames_ + crvBranchNames_ + trkBranchNames_ + trkFitBranchNames_

def GetData(file, quiet): # finName):

    # 1. Coincidence masks should ONLY be applied to coincidences.
    # 2. Global track masks should be applied to tracks and track fits.
    # 3. Local track masks should be applied to track fits.
    # ... We need to split the data array up into these parts.

    # Get data
    # finName = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00089/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000015.root"
    # treeName = "TrkAnaExt/trkana"
    data_dict_ = {}

    # Open tree
    # with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
    with file["TrkAnaExt/trkana"] as tree:
        # Seperate event info, coincidnces, tracks, and track fits. 
        # This way we can apply masks independently. 
        for field, branch in allBranchNames_.items():
            data_dict_[field] = tree.arrays(branch)
            
        # evtData_ = tree.arrays(evtBranchNames_) 
        # crvData_ = tree.arrays(crvBranchNames_) 
        # trkData_ = tree.arrays(trkBranchNames_)
        # trkFitData_ = tree.arrays(trkFitBranchNames_)

    # Zip them together 
    return ak.zip(data_dict_) 
    # return ak.zip({"evt" : evtData_, "crv" : crvData_, "trk" : trkData_, "trkfit" : trkFitData_}) 

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(data_, recon, foutTag, quiet):
    
    #TODO: add tracker sanity plots.
    if not quiet: print("\n---> Making sanity plots")

    # Reco
    sectors_ = ak.flatten(data_["crv"]["crvcoincs.sectorType"]) 
    t_ = ak.flatten(data_["crv"]["crvcoincs.time"]) 
    x_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fX"]) * 1e-3 # mm -> m
    y_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fY"]) * 1e-3 # mm -> m
    z_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fZ"]) * 1e-3 # mm -> m
    PEs_ = ak.flatten(data_["crv"]["crvcoincs.PEs"]) 
    PEsPerLayer_ = ak.flatten(data_["crv"]["crvcoincs.PEsPerLayer[4]"])
    nHits_ = ak.flatten(data_["crv"]["crvcoincs.nHits"]) 
    nLayers_ = ak.flatten(data_["crv"]["crvcoincs.nLayers"]) 
    slopes_ = ak.flatten(data_["crv"]["crvcoincs.angle"])

    ut.PlotGraphOverlay(graphs_=[(x_[sectors_ == 3], y_[sectors_ == 3]), (x_[sectors_ == 1], y_[sectors_ == 1]), (x_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="x-position [m]", ylabel="y-position [m]", fout=f"../Images/{recon}/Sanity/gr_XY_{foutTag}.png")
    ut.PlotGraphOverlay(graphs_=[(z_[sectors_ == 3], y_[sectors_ == 3]), (z_[sectors_ == 1], y_[sectors_ == 1]), (z_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="z-position [m]", ylabel="y-position [m]", fout=f"../Images/{recon}/Sanity/gr_ZY_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[t_[sectors_ == 3], t_[sectors_ == 1], t_[sectors_ == 2]], nbins=1000, xmin = np.min(t_), xmax = np.max(t_), xlabel="Average hit time [ns]", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_times_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nHits_[sectors_ == 3], nHits_[sectors_ == 1], nHits_[sectors_ == 2]], nbins=41, xmin = -0.5, xmax = 40.5, xlabel="Number of hits", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_nHits_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nLayers_[sectors_ == 3], nLayers_[sectors_ == 1], nLayers_[sectors_ == 2]], nbins=5, xmin = -0.5, xmax = 4.5, xlabel="Number of layers hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_nLayers_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[slopes_[sectors_ == 3], slopes_[sectors_ == 1], slopes_[sectors_ == 2]], nbins=1000, xmin = -2, xmax = 2, xlabel="Slope", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_slopes_{foutTag}.png") 

    # Calculate the sum of each length-4 array in PEsPerLayer_
    PEsPerLayerSum_ = ak.sum(PEsPerLayer_, axis=-1)
    PEsDiff_ = PEs_ - PEsPerLayerSum_

    ut.Plot1DOverlay(hists_=[PEs_[sectors_ == 3], PEs_[sectors_ == 1], PEs_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_PEs_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsPerLayerSum_[sectors_ == 3], PEsPerLayerSum_[sectors_ == 1], PEsPerLayerSum_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs per layer", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_PEsPerLayer_{foutTag}.png")

    ut.Plot1DOverlay(hists_=[PEsPerLayer_[:, 0], PEsPerLayer_[:, 1], PEsPerLayer_[:, 2], PEsPerLayer_[:, 3]], nbins=500, xmin = 0, xmax = 500, title="PEs per layer", ylabel="Coincidences", label_=["Layer0", "Layer 1", "Layer 2", "Layer 3"], fout=f"../Images/{recon}/Sanity//h1_PEsPerLayerAllSectors_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsDiff_[sectors_ == 3], PEsDiff_[sectors_ == 1], PEsDiff_[sectors_ == 2]], nbins=40, xmin = -0.0002, xmax = 0.0002, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_PEsDiff_{foutTag}.png", logY=True)
    ut.Plot1DOverlay(hists_=[PEsDiff_], nbins=200, xmin = np.min(PEsDiff_)-1e-4, xmax = np.max(PEsDiff_)+1e-4, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["All modules"], fout=f"../Images/{recon}/Sanity/h1_PEsDiffAll_{foutTag}.png", logY=True, includeBlack=True)

    # MC 
    valid_ = ak.flatten(data_["crv"]["crvcoincsmc.valid"])
    pdgid_ = ak.flatten(data_["crv"]["crvcoincsmc.pdgId"])

    ut.Plot1DOverlay(hists_=[valid_[sectors_ == 3], valid_[sectors_ == 1], valid_[sectors_ == 2]], nbins=2, xmin = 0, xmax = 2, xlabel="Valid MC event", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{recon}/Sanity/h1_valid_{foutTag}.png") 
    
    label_dict = {
        2212: 'proton',
        211: 'pi+',
        -211: 'pi-',
        -13: 'mu+',
        13: 'mu-',
        -11: 'e+',
        11: 'e-',
        "other": "other"
        # Add more particle entries as needed
    }

    ut.BarChart(data_=pdgid_, label_dict=ut.particle_dict, ylabel="Coincidences [%]", fout=f"../Images/{recon}/Sanity/bar_pdgid_{foutTag}.png", percentage=True)
    ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout=f"../Images/{recon}/Sanity/bar_overlay_pdgid_{foutTag}.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

    if not quiet: print("Done!")

    return

# Printout helper functions
def PrintEvent(event, showMasks=False):
    
    eventStr = (
        f"-------------------------------------------------------------------------------------\n"
        f"***** evt *****\n"
        f"evtinfo.run: {event['evt']['evtinfo.run']}\n" 
        f"evtinfo.subrun: {event['evt']['evtinfo.subrun']}\n" 
        f"evtinfo.eventid: {event['evt']['evtinfo.event']}\n"
        f"***** crv *****\n"
        f"crvcoincs.sectorType: {event['crv']['crvcoincs.sectorType']}\n"
        f"crvcoincs.nLayers {event['crv']['crvcoincs.nLayers']}\n"
        f"crvcoincs.angle: {event['crv']['crvcoincs.angle']}\n"
        f"crvcoincs.pos.fCoordinates: ({event['crv']['crvcoincs.pos.fCoordinates.fX']}, {event['crv']['crvcoincs.pos.fCoordinates.fY']}, {event['crv']['crvcoincs.pos.fCoordinates.fZ']})\n"
        f"crvcoincs.timeStart: {event['crv']['crvcoincs.timeStart']}\n"
        f"crvcoincs.time: {event['crv']['crvcoincs.time']}\n"
        f"crvcoincs.PEs: {event['crv']['crvcoincs.PEs']}\n"
        f"crvcoincs.PEsPerLayer[4]: {event['crv']['crvcoincs.PEsPerLayer[4]']}\n"
        f"crvcoincs.nHits: {event['crv']['crvcoincs.nHits']}\n"
        f"crvcoincsmc.pdgId: {event['crv']['crvcoincsmc.pdgId']}\n"
        f"crvcoincsmc.valid: {event['crv']['crvcoincsmc.valid']}\n"
        f"***** trk *****\n"
        f"kl.status: {event['trk']['kl.status']}\n"
        f"kl.nactive: {event['trk']['kl.nactive']}\n"
        f"kl.nhits: {event['trk']['kl.nhits']}\n"
        f"kl.nplanes: {event['trk']['kl.nplanes']}\n"
        f"kl.nnullambig: {event['trk']['kl.nnullambig']}\n"
        f"kl.ndof: {event['trk']['kl.ndof']}\n"
        f"kl.kl.fitcon: {event['trk']['kl.fitcon']}\n"
        f"***** trkfit *****\n"
        f"klfit: {event['trkfit']['klfit']}\n"
        f"klfit.sid: {event['trkfit']['klfit']['sid']}\n"
        f"klfit.sindex: {event['trkfit']['klfit']['sindex']}\n"
        # f"Example variables from klfit...\n"
        # f"klfit.sid: {event['trkfit']['klfit']['sid']}\n"
        # f"klfit.sindex: {event['trkfit']['klfit']['sindex']}\n"
        # f"klfit: {event['trkfit']['klfit']}\n"
        f"klkl: {event['trkfit']['klkl']}\n"
        f"klkl.z0err: {event['trkfit']['klkl']['z0err']}\n"
        f"klkl.d0err: {event['trkfit']['klkl']['d0err']}\n"
        f"klkl.thetaerr: {event['trkfit']['klkl']['thetaerr']}\n"
        f"klkl.phi0err: {event['trkfit']['klkl']['phi0err']}\n"
        f"-------------------------------------------------------------------------------------\n"
    )

    if showMasks: 
        eventStr = eventStr[:-86]
        eventStr += f"***** masks *****\n" 
        eventStr += f"sector_condition: {event['sector_condition']}\n" 
        eventStr += f"PE_condition: {event['PE_condition']}\n" 
        eventStr += f"layer_condition: {event['layer_condition']}\n" 
        eventStr += f"angle_condition: {event['angle_condition']}\n" 
        eventStr += f"CRVT_coincidence: {event['CRVT_coincidence']}\n" 
        eventStr += f"trigger: {event['trigger']}\n" 
        eventStr += f"oneOrZeroCoincInMeasurementSector: {event['oneOrZeroCoincInMeasurementSector']}\n"
        eventStr += f"oneCoinInTriggerSectors: {event['oneCoinInTriggerSectors']}\n"
        eventStr += f"filter_zero: {event['filter_zero']}\n"
        # eventStr += f"trk_bestFit: {event['trk_bestFit']}\n"
        # eventStr += f"trkfit_bestFit: {event['trkfit_bestFit']}\n"
        # eventStr += f"trkfit_KLCRV1: {event['trkfit_KLCRV1']}\n"
        # eventStr += f"goodTrk: {event['goodTrk']}\n"
        # eventStr += f"goodTrkFit: {event['goodTrkFit']}\n"

        eventStr += f"-------------------------------------------------------------------------------------\n"
    
    return eventStr

def PrintNEvents(data_, nEvents=10, showMasks=False):
     # Iterate event-by-event
    for i, event in enumerate(data_, start=1):
        print(PrintEvent(event, showMasks))
        if i >= nEvents: 
            return

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

    return # data_ # ["crv"][data_["is_coincidence"]]

# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

# Not needed.
# def RemoveEmptyEvents(data_):
#     print(f"\n---> Removing empty events")
#     emptyEventCondition = ak.num(data_["crv"]["crvcoincs.nHits"]) == 0
#     return data_[~emptyEventCondition]

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

    # Use ak.mask here, since the output from ak.count will not match the dimensions of sectorType array. 
    # That is, a element like [2, 3] gets masked with "True" for the whole array
    # data_["crv"] = ak.mask(data_["crv"], (oneHitInSector2Condition & oneHitInSector3Condition))
    
    data_["oneOrZeroCoincInMeasurementSector"] = oneOrZeroCoincInMeasurementSector 
    data_["oneCoinInTriggerSectors"] = (oneCoincInSector2Condition & oneCoincInSector3Condition)
    data_["filter_zero"] = (oneOrZeroCoincInMeasurementSector & oneCoincInSector2Condition & oneCoincInSector3Condition)
    # data_["filter_zero"] = (oneCoincInSector2Condition & oneCoincInSector3Condition)
    
    if not quiet: print("Done!")
    
    # Cut on event level
    if not fail: 
        return data_[data_["filter_zero"]]
    else: 
        return data_[~data_["filter_zero"]]
        

# Events that pass the tracker cuts 
def ApplyTrackerCuts(data_, fail, quiet):
    
    if not quiet: print(f"\n---> Applying tracker cuts") 
    
    data_["trkfit_KLCRV1"] = ( 
        (data_["trkfit"]["klfit"]["sid"] == 200) 
        & (data_["trkfit"]["klfit"]["sindex"] == 1) )

    data_["trk_bestFit"] = ( 
        (data_["trk"]["kl.ndof"] >= 10)
        & (data_["trk"]["kl.fitcon"] > 0.1)
        & ((data_["trk"]["kl.nactive"]/data_["trk"]["kl.nhits"]) > 0.99)
        & (data_["trk"]["kl.nplanes"] >= 4)
        & ((data_["trk"]["kl.nnullambig"]/data_["trk"]["kl.nhits"]) < 0.2) )
    
    data_["trkfit_bestFit"] = ( 
        (data_["trkfit"]["klkl"]["z0err"] < 1) 
        & (data_["trkfit"]["klkl"]["d0err"] < 1) 
        & (data_["trkfit"]["klkl"]["thetaerr"] < 0.004)
        & (data_["trkfit"]["klkl"]["phi0err"] < 0.001) )

    # Apply cuts on the track and track fit level
    if not fail: 
        data_["trkfit"] = data_["trkfit"][(data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"])]
        data_["trk"] = data_["trk"][data_["trk_bestFit"]]
    else: 
        data_["trkfit"] = data_["trkfit"][~(data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"])]
        data_["trk"] = data_["trk"][~data_["trk_bestFit"]]

    # Mark events which now have empty tracks or track fits after cuts
    data_["goodTrk"] = ak.any(data_["trk"]["kl.status"], axis=1, keepdims=False) > 0 
    # Could be handled better? 
    # data_["goodTrkFit"] = (
    #     (ak.any(data_["trkfit"]["klfit"]["sid"], axis=-1, keepdims=False) > 0) 
    #     & (ak.any(data_["trkfit"]["klkl"]["z0err"], axis=-1, keepdims=False) > 0) )
    data_["goodTrkFit"] = (
        (ak.count(data_["trkfit"]["klfit"]["sid"], axis=-1, keepdims=False) > 0) 
        & (ak.count(data_["trkfit"]["klkl"]["z0err"], axis=-1, keepdims=False) > 0) )
    # Flatten (needed for masking)
    data_["goodTrkFit"] = ak.any(data_["goodTrkFit"], axis=-1, keepdims=False) == True 
        # & (ak.any(data_["trkfit"]["klkl"]["z0err"], axis=-1, keepdims=False) > 0) )
    
    # (ak.any(data_["trkfit_KLCRV1"], axis=-1, keepdims=False) == True))

    # # Apply cut on event level
    # data_["goodTrk"] = ak.any(data_["trk"]["kl.status"], axis=1, keepdims=False) > 0 
    # # Check for a good track fit after cuts
    # data_["goodTrkFit"] = (ak.all(data_["trkfit"]["klfit"]["sid"], axis=1, keepdims=False) > 0) & (ak.count(data_["trkfit"]["klkl"], axis=1, keepdims=False) > 0) 
    # Return events passing track fits
    return data_[(data_["goodTrk"] & data_["goodTrkFit"])]
    
    # if not fail:  
    #     return data_[data_["goodTrk"]]
    # else:
    #     return data_[~data_["goodTrk"]]

# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_, quiet): 
    if not quiet: print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_["crv"]["crvcoincs.timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   

# Basic trigger condition
def Trigger(data_, fail, quiet): 

    if not quiet: print(f"\n---> Triggering")

    # Enforce basic trigger condition
    # triggerCondition = (
    #     ak.any(data_["is_coincidence"] & (data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
    #     ak.any(data_["is_coincidence"] & (data_["crv"]["crvcoincs.sectorType"] == 3), axis=1)
    # )

    # Leave the trigger sectors at 2 layers and 10 PEs. 
    triggerCondition = (
        ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
        ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1)
    )

    data_["trigger"] = triggerCondition
    
    if not quiet: print("Done!")

    # return data_[~triggerCondition]

    if not fail: 
        return data_[triggerCondition]
    else:
        return data_[~triggerCondition]
        

# Needs to come after the initial trigger 
def SuccessfulTriggers(data_, success, quiet):

    successStr = ""
    if success: successStr += "successful"
    else: successStr += "unsuccessful"

    if not quiet: print(f"\n---> Getting {successStr} triggers")

    # Ensure at least one passing coincidence the measurement sector
    # Coincidence conditions are set in FindCoincidences()
    successCondition = ak.any(data_["CRVT_coincidence"], axis=1) 
    
    # ONE coincidence in the measurement sector 
    # successCondition = ak.count(data_["CRVT_coincidence"], axis=1) == 1
    
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
        fout.write("evtinfo.run, evtinfo.subrun, evtinfo.event\n")
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

    
def WriteFailureInfo(failures_, recon, finTag, foutTag, coincidenceConditions, quiet):

    # Define the output file path
    foutNameConcise = f"../Txt/{recon}/failures_concise/{finTag}/failures_concise_{foutTag}.csv" 
    foutNameVerbose = f"../Txt/{recon}/failures_verbose/{finTag}/failures_verbose_{foutTag}.csv" 
    
    if not quiet: print(f"\n---> Writing failure info to:\n{foutNameConcise}\n{foutNameVerbose}", flush=True) 

# allBranchNames_ = evtBranchNames_ + crvBranchNames_ + trkBranchNames_ + trkFitBranchNames_

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
    if True: 
        with open(foutNameVerbose, "w") as fout:
            # Write the events
            for event in failures_:
                fout.write(
                    PrintEvent(event)+"\n" 
                )

    return

def WriteResults(data_, successes_, failures_, recon, finTag, foutTag, quiet):

    foutName = f"../Txt/{recon}/results/{finTag}/results_{foutTag}.csv"
    if not quiet: print(f"\n---> Writing results to {foutName}")
    tot = len(data_)
    efficiency = len(successes_) / tot * 100
    inefficiency = len(failures_) / tot * 100

    outputStr = f"""
    ****************************************************
    Input tag: {finTag}
    Output tag: {foutTag}
    Number of failures: {len(failures_)}
    Number of successes: {len(successes_)}
    Efficiency: {len(successes_)}/{tot} = {efficiency:.2f}%
    Inefficiency: {len(failures_)}/{tot} = {inefficiency:.2f}%
    ****************************************************
    """

    with open(foutName, "w") as fout:
        fout.write("Total,Successes,Failures,Efficiency [%],Inefficiency [%]\n") # no spaces!
        fout.write(f"{tot}, {len(successes_)}, {len(failures_)}, {efficiency}, {inefficiency}\n")
        # fout.write(outputStr)

    if True: print(outputStr)

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

# Best practice?
filters_ = { 
    0 : "singles"
    , 1 : "singles_track_cuts"
}

# def Run(file, recon, particle, coincidenceConditions, finTag, foutTag, sanityPlots, quiet): 
def Run(file, recon, particle, PE, layer, filterLevel, finTag, sanityPlots, quiet): 

    # A bit ugly?
    coincidenceConditions = f"{PE}PEs{layer}Layers"
    foutTag = particle + "_" + coincidenceConditions + "_" + filters_[filterLevel]
    
    # Placeholder
    fail = False

    # Get data as a set of awkward arrays
    data_ = GetData(file, quiet) # finName)
    
    # Sanity plots 
    if (sanityPlots): SanityPlots(data_, recon, foutTag, quiet)

    # Apply start time cut
    data_ = CutOnStartTime(data_, quiet)

    data_ = FilterParticles(data_, particle, quiet)

    # Basic trigger
    data_ = Trigger(data_, fail, quiet)

    # Pass 1 
    # This limits the analysis to a subset of the data

    # Is this the best way to handle this??? 
    # if coincidenceFilterLevel == "pass0":
    if filterLevel == 0:
        data_ = FilterSingles(data_, fail, quiet)
    elif filterLevel == 1: 
        data_ = FilterSingles(data_, fail, quiet)
        data_ = ApplyTrackerCuts(data_, fail, quiet)
        
    
    # elif coincidenceFilterLevel == "pass1":
    #     print()
    # elif coincidenceFilterLevel == "pass2":
    #     print()
    # elif coincdienceFilterLevel == "pass3":
    #     data_ = PassThree(data_, fail=False)

    # This needs to come after all the filtering!
    # Well it needs to come after the trigger I guess. 
    # Find coincidences in measurement sector
    # This could actually be written into SuccessfulTriggers no? 
    FindCoincidences(data_, coincidenceConditions, quiet)

    # PrintNEvents(data_, 10, True)
    # thisEventCondition = ((data_["evt"]["evtinfo.run"] == 1205) & (data_["evt"]["evtinfo.subrun"] == 1497) & (data_["evt"]["evtinfo.event"] == 19255))
    # data_["evt"] = data_["evt"][thisEventCondition]
    
    # PrintEvent(data_["evt"], showMasks=True)
    
    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True, quiet=quiet)
    failures_ = SuccessfulTriggers(data_, success=False, quiet=quiet)

    # Write failures to file
    # WriteSuccessInfo(successes_, recon, finTag, foutTag, coincidenceConditions, quiet) 
    WriteFailureInfo(failures_, recon, finTag, foutTag, coincidenceConditions, quiet) 

    # Write results to file
    WriteResults(data_, successes_, failures_, recon, finTag, foutTag, quiet) 

    return

from concurrent.futures import ThreadPoolExecutor, as_completed

# ------------------------------------------------
#              Multithread 
# ------------------------------------------------
import subprocess

# Read a single file
def ReadFile2(fileName, quiet=False): 
    try:
        # Setup commands
        commands = "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh; muse setup ops;"
        commands += f"echo {fileName} | mdh print-url -s root -"
        if not quiet:
            print(f"---> Reading file:\n{fileName}")
        # Execute commands 
        fileName = subprocess.check_output(commands, shell=True, universal_newlines=True)
        if not quiet:
            print(f"\n---> Created xroot url:\n{fileName}")
            print("---> Opening file with uproot...") 
        # Open the file 
        with uproot.open(fileName) as file:
            if not quiet: 
                print("Done!")
            return file 
    except OSError as e:
        # Setup alternative commands 
        print(f"\n----> Exception timeout while opening file with xroot, retrying locally: {fileName}")                    
        commands = "source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh; muse setup ops;"
        commands += f"echo {fileName} | mdh copy-file -s tape -l local -" 
        # Execute commands
        subprocess.check_output(commands, shell=True, universal_newlines=True)
        # Return the opened file 
        with uproot.open(fileName) as file:
            return file
            
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product

# Multithread on one file with many configs
def Multithread2(processFunction2, fileName, particles_, layers_, PEs_):
    
    print("\n---> Starting multithreading...")
    
    totalConfigurations = len(particles_) * len(layers_) * len(PEs_)
    completedConfigurations = 0
    
    with ThreadPoolExecutor() as executor:
        
        # Generate all combinations of (particle, layer, PE)
        configurations = product(particles_, layers_, PEs_)
        
        # Submit tasks to the executor
        futures = {
            executor.submit(processFunction2, fileName, particle, layer, PE): (particle, layer, PE)
            for particle, layer, PE in configurations
        }

        # Process results as they complete
        for future in as_completed(futures):
            (particle, layer, PE) = futures[future]
            try:
                future.result()  # Retrieve the result or raise an exception if occurred
                completedConfigurations += 1
                percentComplete = (completedConfigurations / totalConfigurations) * 100
                print(f'---> Configuration ({particle}, {layer}, {PE}) processed successfully! ({percentComplete:.1f}% complete)')
            except Exception as exc:
                print(f'---> Configuration ({particle}, {layer}, {PE}) generated an exception!\n{exc}')
                
    print("\nMultithreading completed!")
    return

# Multithread on many files with many configs
# This seems to use huge amount of memory! Like 35 GB per file... 
def Multithread3(processFunction2, fileList_, particles_, layers_, PEs_, max_workers=254):
    
    print("\n---> Starting multithreading...\n")
    
    totalFiles = len(fileList_)
    totalConfigurations = len(particles_) * len(layers_) * len(PEs_)
    completedTasks = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        # Prepare the list of futures
        futures = []
        
        # Submit tasks for each file
        for fileName in fileList_:
            # Generate all combinations of (particle, layer, PE) for each file
            configurations = product(particles_, layers_, PEs_)
            for particle, layer, PE in configurations:
                futures.append(
                    executor.submit(processFunction2, fileName, particle, layer, PE)
                )
        
        # Process results as they complete
        for future in as_completed(futures):
            try:
                future.result()  # Retrieve the result or raise an exception if occurred
                completedTasks += 1
                percentComplete = (completedTasks / (totalFiles * totalConfigurations)) * 100
                print(f'---> Task {completedTasks} processed successfully! ({percentComplete:.1f}% complete)')
            except Exception as exc:
                print(f'---> Task generated an exception!\n{exc}')
                
    print("\n---> Multithreading completed!")
    return


def Multithread4(processFunction3, fileList_, max_workers=381):
    
    print("\n---> Starting multithreading...\n")

    completedTasks = 0
    totalFiles = len(fileList_)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        # Prepare the list of futures
        futures_ = [executor.submit(processFunction3, fileName) for fileName in fileList_]
        
        # Process results as they complete
        for future in as_completed(futures_):
            try:
                future.result()  # Retrieve the result or raise an exception if occurred
                completedTasks += 1
                percentComplete = (completedTasks / (totalFiles)) * 100
                print(f'\n---> Task {completedTasks} processed successfully! ({percentComplete:.1f}% complete)')
            except Exception as exc:
                print(f'\n---> Task generated an exception!\n{exc}')
                
    print("\n---> Multithreading completed!")
    return
    
# This is good. TODO: add this to Mu2eEAF. 
def Multithread5(processFunction, fileList_, max_workers=381):
    
    print("\n---> Starting multithreading...\n")

    completedFiles = 0
    totalFiles = len(fileList_)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        # Prepare a list of futures and map them to file names
        futures_ = {executor.submit(processFunction, fileName): (fileName) for fileName in fileList_}
        
        # Process results as they complete
        for future in as_completed(futures_):
            fileName = futures_[future]  # Get the file name associated with this future
            # print(f'\n---> Processing {fileName}...') 
            try:
                future.result()  # Retrieve the result or raise an exception if occurred
                completedFiles += 1
                percentComplete = (completedFiles / (totalFiles)) * 100
                print(f'\n---> {fileName} processed successfully! ({percentComplete:.1f}% complete)')
            except Exception as exc:
                print(f'\n---> {fileName} generated an exception!\n{exc}')
                
    print("\n---> Multithreading completed!")
    return

def Multithread6(processFunction, fileList_, max_workers=381):
    
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
    
def main():

    # return
    defname = "nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.root"
    recon = "MDC2020ae"
    
    filterLevels_ = [0, 1]
    particles_ = ["all", "muons", "non_muons"]
    layers_ = [3, 2]
    # PEs_ = np.arange(10, 135, 5)
    PEs_ = np.arange(100, 135, 5)

    # filterLevel_ = [0, 1] # [0, 1]
    # particles_ = ["all"] # , "muons", "non_muons"]
    # layers_ = [2] # , 3]
    # PEs_ = [10] # np.arange(10, 135, 5)

    sanityPlots = False
    quiet = True # True

    def processFunction2(fileName, particle, layer, PE):
        # Always open the file in the processFunction 
        file = rd.ReadFile(fileName, quiet)
        finTag = fileName.split('.')[-2] 
        outputStr = (
                        "\n---> Running with:\n"
                        f"fileName: {fileName}\n"
                        f"recon: {recon}\n"
                        f"particle: {particle}\n"
                        f"PEs: {PE}\n"
                        f"layers: {layer}/4\n"
                        f"filterLevel: {filterLevel}\n"
                        f"finTag: {finTag}\n"
                        f"sanityPlots: {sanityPlots}\n"
                        f"quiet: {quiet}\n"
                    )
        if not quiet: print(outputStr, flush=True)
        Run(file, recon, particle, PE, layer, filterLevel, finTag, sanityPlots, quiet)
        return

    def processFunction3(fileName):
        # Always open the file in the processFunction 
        file = rd.ReadFile(fileName, quiet)
        finTag = fileName.split('.')[-2] 
        
        # Build one datapoint at a time, so PEs need to go on the bottom
        # Scan PE thresholds
        for PE in PEs_: # Data point level
            # Scan particles
            for particle in particles_: # Series level
                # Scan filters
                for filterLevel in filterLevels_: # Plot level
                    # Scan layers
                    for layer in layers_: # Plot level
                    
                        outputStr = (
                            "\n---> Running with:\n"
                            f"fileName: {fileName}\n"
                            f"recon: {recon}\n"
                            f"particle: {particle}\n"
                            f"PEs: {PE}\n"
                            f"layers: {layer}/4\n"
                            f"filterLevel: {filterLevel}\n"
                            f"finTag: {finTag}\n"
                            f"sanityPlots: {sanityPlots}\n"
                            f"quiet: {quiet}\n"
                        )
                        if not quiet: print(outputStr) 
                            
                        Run(file, recon, particle, PE, layer, filterLevel, finTag, sanityPlots, quiet)
        return

    # Resubmission of failures.
    def processFunction4(fileName):

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
            filterLevel = get_key(filters_, row["Filter"])

            outputStr = (
                    "\n---> Running with:\n"
                    f"fileName: {fileName}\n"
                    f"recon: {recon}\n"
                    f"particle: {particle}\n"
                    f"PEs: {PE}\n"
                    f"layers: {layer}/4\n"
                    f"filterLevel: {filterLevel}\n"
                    f"finTag: {finTag}\n"
                    f"sanityPlots: {sanityPlots}\n"
                    f"quiet: {quiet}\n"
                )

            if not quiet: print(outputStr) 

            try:
                Run(file, recon, particle, PE, layer, filterLevel, finTag, sanityPlots, quiet)
            except Exception as exc:
                print(f'---> Exception!\n{row}\n{exc}')
            
        return

        
    fileList_ = rd.GetFileList(defname) 


    # Get failed jobs
    failedJobsFile = "../Txt/MDC2020ae/FailedJobs/failures.csv"
    df_failedJobs_ = pd.read_csv(failedJobsFile)

    # Set of tags
    tags_ = list(set(df_failedJobs_["Tag"]))
    
    # Extract the tag from the file name
    def ExtractTag(fileName):
        parts = fileName.split('.')
        if len(parts) > 1:
            return parts[-2]
        return None

    # Filter and sort files based on tags
    fileList_ = sorted(
        [file for file in fileList_ if ExtractTag(file) in tags_]
        , key=lambda file: tags_.index(ExtractTag(file))
    )

    print(f"---> Got {len(fileList_)} files.")
    # Multithread5(processFunction3, fileList_) 
    Multithread5(processFunction4, fileList_)
    Multithread6(processFunction4, fileList_)

    
    # Multithread3(processFunction2, fileList_[:2], particles_, layers_, PEs_) 
    
    # Multithread4(processFunction3, ["nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root"]) 
    # Multithread3(processFunction2, fileList_[:20], particles_, layers_, PEs_) 
    # Multithread2(processFunction2, fileName, particles_, layers_, PEs_) 

    return

if __name__ == "__main__":
    main()
