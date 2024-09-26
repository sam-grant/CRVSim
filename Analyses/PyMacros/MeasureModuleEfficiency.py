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
from itertools import product
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from statsmodels.stats.proportion import proportion_confint

# Internal libraries
import Utils as ut
import PrintUtils as pr
import CoincidenceConditions as cc

from mu2etools import read_data as rd 
from mu2etools import parallelise as pa
    
# ------------------------------------------------
# Coincidence finding in measurement sector
# ------------------------------------------------

# Impose stricter conditions than the default, if desired
# This works but needs a be a refactor, it's a bit awkard.
def MarkCoincidences(data_, coincidenceConditions, quiet): 
    
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

    # if not quiet: print("Done!")

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

def FilterSingles(data_, triggerMode, fail, quiet):
    
    if not quiet: print(f"\n---> Filtering singles") 

    sector1Condition = data_["crv"]["crvcoincs.sectorType"] == 1
    sector2Condition = data_["crv"]["crvcoincs.sectorType"] == 2
    sector3Condition = data_["crv"]["crvcoincs.sectorType"] == 3

    oneOrZeroCoincInMeasurementSector = ak.count(data_["crv"]["crvcoincs.sectorType"][sector1Condition], axis=1) <= 1 # Allows for failures
    oneCoincInSector2Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector2Condition], axis=1) == 1
    oneCoincInSector3Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector3Condition], axis=1) == 1
    
    # data_["oneOrZeroCoincInMeasurementSector"] = oneOrZeroCoincInMeasurementSector 
    # data_["oneCoinInTriggerSectors"] = (oneCoincInSector2Condition & oneCoincInSector3Condition)

    if triggerMode in ["crv_trigger", "trk_crv_trigger", "crv_2layers_trigger", "crv_3layers_trigger", "trk_crv_2layers_trigger", "trk_crv_3layers_trigger"]:
        data_["pass_singles"] = (oneOrZeroCoincInMeasurementSector & oneCoincInSector2Condition & oneCoincInSector3Condition)
    elif triggerMode in ["crv2_trigger", "trk_crv2_trigger", "trk_crv2_2layers_trigger", "trk_crv2_3layers_trigger"]:
        data_["pass_singles"] = (oneOrZeroCoincInMeasurementSector & oneCoincInSector2Condition)
    elif triggerMode in ["crv3_trigger", "trk_crv3_trigger"]:
        data_["pass_singles"] = (oneOrZeroCoincInMeasurementSector & oneCoincInSector3Condition)        
    elif triggerMode in ["trk_trigger"]:
        data_["pass_singles"] = (oneOrZeroCoincInMeasurementSector)  
    else:
        raise ValueError(f"---> triggerMode '{triggerMode}' not valid...")
            
    # if not quiet: print("Done!")
    
    # Cut on event level
    if not fail: 
        return data_[data_["pass_singles"]]
    else: 
        return data_[~data_["pass_singles"]]
    
def ApplyTrackerCuts(data_, triggerMode="default", fail=False, quiet=False): 

    if triggerMode in ["crv_trigger", "crv2_trigger", "crv3_trigger", "crv_2layers_trigger", "crv_3layers_trigger"]:
        raise ValueError(f"Improper trigger mode {triggerMode} called in ApplyTrackerCuts()...")
    
    if not quiet: print(f"\n---> Marking tracker cuts") 

    # Mark cuts on the track and track fit level
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

    ###### Fiducial area ######

    # Dimensions from crv_parameters.xlsx
    # CRV-T module overall length 6100 mm
    # CRV-T module overall width 951.0 mm 
    # CRV-DS length = 2570 mm
    # CRV-DS width = 826 mm 
    # CRV-L-end length = 3388 mm
    # CRV-L-end width = 951 mm 
    
    # Total layer offset 42.00*(4-1) = 127 mm
    # Four modules side-by-side w_tot = w + (N_mod-1)*(w_mod-off) 
        
    # CRV-T
    # width ---> 951 + (4-1)*(951-127) = 3388
    # z: 3388/2 (offset by -500 mm). Need to subtract the layer offset as well. 3388/2 - 127. 
    # x: 6100/2 
    # (z, x)
    # min_box_coords = (-3388/2-500, -6100/2)
    # max_box_coords = (3388/2-500, 6100/2)
    
    # CRV-DS
    # width ---> 826 + (2-1)*(826-127) = 1525
    # It is rotated. length -> width. 
    # z-coords is length for DS 
    # z: 2570/2 (offset by -500 mm) 
    # x: 1525/2 
    # (z, x)
    # min_box_coords = (-2570/2-500, -1525/2)
    # max_box_coords = (2570/2-500, +1525/2)

    # CRV-L-end
    # width ---> 951 + (2-1)*(951-127) = 1775
    # z-coords is length for DS 
    # z: 1775/2 (offset by -500 mm) 
    # x: 3388/2  
    # (z, x)
    # min_box_coords = (-(1775/2)-500, -(3388/2))
    # max_box_coords = (+(1775/2)-500, +(3388/2))

    # Trigger modules 
    # In z:
    # CRV-L-end width < CRV-DS length
    
    # In x:
    # CRV-L-end length > CRV-DS width
    
    # So the box is defined by:
    # z: CRV-DS length (2570/2) 
    # x: CRV-L-end length (3388/2) 
    # (z, x)
    # min_box_coords = (-(2570/2)-500, -(3388/2))s
    # max_box_coords = (+(2570/2)-500, +(3388/2))

    # Tested in TrackCuts.ipynb and CompCRV.ipynb
    
    # Measurement module
    data_["trkfit_CRV1Fiducial"] = ( 
        (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"]) < 6100/2) 
        & (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < (3388/2 - 127*2)) ) # Make it 2*layer offset to be  sure. 

    # Trigger modules
    data_["trkfit_CRV23Fiducial"] = ( 
        (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"]) < 3388/2)
        & (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 2570/2) ) 

    # CRV-DS 
    data_["trkfit_CRV2Fiducial"] = ( 
        (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"]) < 1525/2)
        & (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 2570/2) ) 

    # CRV-L-end
    data_["trkfit_CRV3Fiducial"] = ( 
        (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"]) < 3388/2)
        & (abs(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 1775/2) ) 

    ################################################
    
    # Track condition 
    data_["pass_trk"] = data_["trk_bestFit"]
    
    # Track fit (segments) condition (default has no area cut)
    if triggerMode == "default":
        data_["pass_trkfit"] = (data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"])
    elif triggerMode in ["trk_crv_trigger", "trk_crv_2layers_trigger", "trk_crv_3layers_trigger"]:
        data_["pass_trkfit"] = (data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"] & data_["trkfit_CRV23Fiducial"])
    elif triggerMode in ["trk_crv2_trigger", "trk_crv2_2layers_trigger", "trk_crv2_3layers_trigger"]:
        data_["pass_trkfit"] = (data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"] & data_["trkfit_CRV2Fiducial"])
    elif triggerMode == "trk_crv3_trigger":
        data_["pass_trkfit"] = (data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"] & data_["trkfit_CRV3Fiducial"]) 
    elif triggerMode == "trk_trigger":
        data_["pass_trkfit"] = (data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"] & data_["trkfit_CRV1Fiducial"]) 
    else:
        raise ValueError(f"---> triggerMode '{triggerMode}' not valid...")
    
    if not quiet: 
        # print("Done!") 
        print(f"\n---> Applying tracker cuts") 

    # Mask the data on trk/trkfit level
    if not fail: 
        data_["trk"] = data_["trk"][data_["pass_trk"]] 
        data_["trkfit"] = data_["trkfit"][data_["pass_trkfit"]]
    else: 
        data_["trk"] = data_["trk"][~data_["pass_trk"]] 
        data_["trkfit"] = data_["trkfit"][~data_["pass_trkfit"]]

    # Now clean up empty events after cuts
    goodTrk = ak.any(data_["trk"]["kl.status"], axis=1, keepdims=False) == True
    goodTrkFit = (
        (ak.any(data_["trkfit"]["klfit"]["sid"], axis=-1, keepdims=False) == True) 
        & (ak.any(data_["trkfit"]["klkl"]["z0err"], axis=-1, keepdims=False) == True) 
    )
    # Reset to event level
    goodTrkFit = ak.any(goodTrkFit, axis=-1, keepdims=False) == True 

    # Mark passing events
    # We need to include this flag in the trigger condition 
    data_["pass_track_cuts"] = (goodTrk & goodTrkFit)

    # Apply the event level cut
    data_ = data_[ data_["pass_track_cuts"] ] 

    # Sanity check: this flag should always be True
    if ak.any(data_["pass_track_cuts"] == False, axis=0) == True:
        raise ValueError("pass_track_cuts is False!") 

    # if not quiet: print("Done!") 
        
    return data_ 


# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_, quiet): 
    if not quiet: print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_["crv"]["crvcoincs.timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   

# Trigger condition
def Trigger(data_, triggerMode, fail, quiet): 

    if not quiet: print(f"\n---> Triggering")
        
    # Trigger condititon with CRV modules
    triggerCondition = False 

    # Modes: 
    # 1. CRV-DS and CRV-L trigger, crv_trigger 
    # X 2. CRV-DS trigger, crv2_trigger 
    # X 3. CRV-L-end trigger, crv3_trigger 
    # 4. Tracker trigger, trk_trigger
    # 5. CRV and tracker trigger, trk_crv_trigger
    # 6. CRV-DS and tracker trigger, trk_crv2_trigger
    # 7. CRV-L-end and tracker trigger, trk_crv3_trigger 
    # 8. CRV-DS (3 layers active) and tracker trigger, trk_crv2_3layers_trigger
    # 9. CRV-DS (2 layers active) and tracker trigger, trk_crv2_2layers_trigger
    # 10. CRV-DS and CRV-L trigger (2 layers active), crv_2layers_trigger
    # 11. CRV-DS and CRV-L trigger (3 layers active), crv_3layers_trigger 
    # 12. CRV and tracker trigger (2 layers active), trk_crv_2layers_trigger
    # 13. CRV and tracker trigger (3 layers active), trk_crv_3layers_trigger 
    
    if triggerMode == "crv_trigger": 
        triggerCondition = (
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1)
        )
    elif triggerMode == "crv2_trigger": 
        triggerCondition = ( 
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) 
        )
    elif triggerMode == "crv3_trigger": 
        triggerCondition = ( 
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1) 
        )
    elif triggerMode == "trk_trigger": 
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet)             
        triggerCondition = ( 
            data_["pass_track_cuts"]
            # This is always True after we apply the track cuts. 
            # Somehow it's nice to include it in the trigger condition anyway, for readability. 
            # We need to apply fit-level masks to the data in order to get the event level flag, so it's bit odd. 
            # Track cuts are tricky. 
        )
    elif triggerMode == "trk_crv_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        triggerCondition = (
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1) &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "trk_crv2_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        triggerCondition = (
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "trk_crv3_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        triggerCondition = (
            ak.any((data_["crv"]["crvcoincs.sectorType"] == 3), axis=1) &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "trk_crv2_2layers_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        # Look for PEsPerLayer >= 10 in the BOTTOM two layers (indices 2 and 3) 
        # layerCondition = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -2:] >= 10
        layerCondition = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, :2] >= 10 # bottom two layers
        layerCondition = ak.flatten((ak.all(layerCondition, axis=-1, keepdims=False)==True), axis=None)
        triggerCondition = (
            layerCondition &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "trk_crv2_3layers_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        # Look for PEsPerLayer >= 10 in 2/3 of the top three layers 
        layerCondition = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -3:] >= 10  
        layerCondition = ak.flatten(ak.sum(layerCondition == True, axis=-1, keepdims=False), axis=None) >= 2 
        triggerCondition = (
            layerCondition &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "crv_2layers_trigger":  
        # layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -2:] >= 10 # top two layers
        # layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, :2] >= 10 # bottom two layers
        layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, :2] >= 10 # bottom two layers
        layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, -2:] >= 10 # top two layers
        layerConditionUpper = ak.flatten((ak.all(layerConditionUpper, axis=-1, keepdims=False)==True), axis=None)
        layerConditionLower = ak.flatten((ak.all(layerConditionLower, axis=-1, keepdims=False)==True), axis=None)
        layerCondition = layerConditionUpper & layerConditionLower
        triggerCondition = (
            layerCondition 
        )
    elif triggerMode == "crv_3layers_trigger":
        # Look for PEsPerLayer >= 10 in the 2/3 of layers of the top and bottom modules 
        layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -3:] >= 10 #  three layers
        layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, :3] >= 10 # bottom three layers        
        layerConditionUpper = ak.flatten(ak.sum(layerConditionUpper == True, axis=-1, keepdims=False), axis=None) >= 2 
        layerConditionLower = ak.flatten(ak.sum(layerConditionLower == True, axis=-1, keepdims=False), axis=None) >= 2 
        layerCondition = layerConditionUpper & layerConditionLower
        triggerCondition = (
            layerCondition 
        )
    elif triggerMode == "trk_crv_2layers_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        # Look for PEsPerLayer >= 10 in the top two layers (indices 2 and 3) 
        # layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -2:] >= 10 # top two layers
        # layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, :2] >= 10 # bottom two layers
        layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, :2] >= 10 # bottom two layers
        layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, -2:] >= 10 # top two layers
        layerConditionUpper = ak.flatten((ak.all(layerConditionUpper, axis=-1, keepdims=False)==True), axis=None)
        layerConditionLower = ak.flatten((ak.all(layerConditionLower, axis=-1, keepdims=False)==True), axis=None)
        layerCondition = layerConditionUpper & layerConditionLower
        triggerCondition = (
            layerCondition &
            data_["pass_track_cuts"]
        )
    elif triggerMode == "trk_crv_3layers_trigger":
        ApplyTrackerCuts(data_, triggerMode=triggerMode, fail=fail, quiet=quiet) 
        # Look for PEsPerLayer >= 10 in the 2/3 of layers of the top and bottom modules
        layerConditionUpper = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 2][:, :, -3:] >= 10 # top three layers
        layerConditionLower = data_["crv"]["crvcoincs.PEsPerLayer[4]"][data_["crv"]["crvcoincs.sectorType"] == 3][:, :, :3] >= 10 # bottom three layers
        layerConditionUpper = ak.flatten(ak.sum(layerConditionUpper == True, axis=-1, keepdims=False), axis=None) >= 2 
        layerConditionLower = ak.flatten(ak.sum(layerConditionLower == True, axis=-1, keepdims=False), axis=None) >= 2 
        layerCondition = layerConditionUpper & layerConditionLower
        triggerCondition = (
            layerCondition &
            data_["pass_track_cuts"]
        )
    else:
        raise ValueError(f"---> triggerMode '{triggerMode}' not valid...")

    data_["pass_trigger"] = triggerCondition
    
    # if not quiet: print("Done!")

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
    
    # if not quiet: print("Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

# ------------------------------------------------
#                     Output 
# ------------------------------------------------ 

def WriteFailureInfo(failures_dict_, recon, finTag, foutTag, coincidenceConditions, quiet):

    for trkTag, failures_ in failures_dict_.items():
        
        # Define the output file path
        foutNameConcise = f"../Txt/{recon}/failures_concise/{finTag}/failures_concise_{foutTag}_{trkTag}.csv" 
        foutNameVerbose = f"../Txt/{recon}/failures_verbose/{finTag}/failures_verbose_{foutTag}_{trkTag}.txt" 
        
        if not quiet: print(f"\n---> Writing failure info to:\n{foutNameConcise}\n{foutNameVerbose}", flush=True) 
    
        # Concise form
        with open(foutNameConcise, "w") as fout:
            # Write the header
            fout.write("evtinfo.run,evtinfo.subrun,evtinfo.event\n")
            # Write the events
            for event in failures_:
                fout.write(
                    f"{event['evt']['evtinfo.run']},{event['evt']['evtinfo.subrun']},{event['evt']['evtinfo.event']}\n"
                )
                
        # Verbose form

        # Setup mask string
        masks_ = ["pass_singles", "pass_trigger", "CRVT_coincidence"]
        if trkTag == "track_cuts":
            masks_ += ["trk_bestFit", "trkfit_bestFit", "trkfit_KLCRV1", "trkfit_CRV1Fiducial", "pass_track_cuts"]
        
        with open(foutNameVerbose, "w") as fout:
            # Write the events
            for event in failures_:
                maskStr = ''.join([f"{mask}: {event[f'{mask}']}\n" for mask in masks_])
                # print(maskStr)
                fout.write(
                    pr.PrintEvent(event, maskStr)+"\n" 
                )

    return

# Write inefficiency results
# Successes and failures are dicts mapping cut labels to a tuple containing a success/failure pair
def WriteResults(results_, recon, finTag, foutTag, quiet):

    # Start output string 
    resultStr = f"""
    ****************************************************
    Input tag: {finTag}
    Output tag: {foutTag}
    """

    foutStr = f"\n----> Written results to:"

    # Loop through tuple
    for cut, result in results_.items():
    
        foutName = f"../Txt/{recon}/results/{finTag}/results_{foutTag}_{cut}.csv"

        foutStr += f"\n{foutName}"

        # if not quiet: print(f"\n---> Writing results to:\n{foutName}")

        # Calculate efficiency 
        successes_ = result[0]
        failures_ = result[1]
        
        # len should give the same value as ak.count(failures_, axis=0)
        # Check this though...
        if ak.num(failures_, axis=0) != len(failures_):
            raise ValueError("---> Missing values in failure array!")

        if ak.num(successes_, axis=0) != len(successes_):
            raise ValueError("---> Missing values in successes array!")
            
        nfailures = len(failures_)
        nsuccesses = len(successes_)
        ntotal = nfailures + nsuccesses
        
        inefficiency = (nfailures / ntotal) * 100 if ntotal > 0 else np.nan

        # # Error from Wilson interval
        # ineffErr = -1 # np.nan
        # if ntotal > 0:
        #     lower, upper = proportion_confint(nfailures, ntotal, method="wilson")
        #     point = nfailures/ntotal
        #     ineffErr = abs((upper - point) / 2) * 100
        # else:
        #     ineffErr = np.nan 
            
        # Append output string 
        resultStr += f"""\n
        Cut: {cut}
        Successes: {nsuccesses}
        Failures: {nfailures}
        Total: {ntotal}
        Inefficiency: {nfailures}/{ntotal} = {inefficiency:.3f}% 
        """
        #+/-{ineffErr:.3f}%
        
        # Write to file
        with open(foutName, "w") as fout:
            fout.write("Total,Successes,Failures\n") # no spaces!
            fout.write(f"{ntotal},{nsuccesses},{nfailures}\n")

    # Finish output string 
    resultStr += f"""
    ****************************************************
    """

    if not quiet: print(foutStr)
    if not quiet: print(resultStr)

    return

# ------------------------------------------------
#   Measure module eff with CRV and tracker
# ------------------------------------------------

def Run(file, recon, particle, PE, layer, finTag, triggerMode, quiet): 
    
    # Output strings
    coincidenceConditions = f"{PE}PEs{layer}Layers"
    foutTag = particle + "_" + coincidenceConditions 

    # Get data as a set of awkward arrays
    data_ = ut.GetData(file) # finName)

    # Apply start time cut
    data_ = CutOnStartTime(data_, quiet)

    # Filter particle species
    data_ = FilterParticles(data_, particle, quiet)

    # Singles cut
    data_ = FilterSingles(data_, triggerMode=triggerMode, fail=False, quiet=quiet)

    # Enforce trigger
    data_ = Trigger(data_, triggerMode=triggerMode, fail=False, quiet=quiet)

    # Sanity check
    if False:

        # trk_crv_trigger
        # min_box_coords = (-(2570/2)-500, -(3388/2))
        # max_box_coords = (+(2570/2)-500, +(3388/2))

        # trk_crv2_trigger
        min_box_coords = (-2570/2-500, -1525/2)
        max_box_coords = (2570/2-500, +1525/2)

        # # # trk_crv3_trigger
        # min_box_coords = (-(1775/2)-500, -(3388/2))
        # max_box_coords = (+(1775/2)-500, +(3388/2))

        # # trk_trigger
        # min_box_coords = (-3388/2-500, -6100/2)
        # max_box_coords = (3388/2-500, 6100/2)

        # Plot the track distribution in zx as a sanity check 
        ut.Plot2D(x=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)#[:1000]
                  , y=ak.flatten(data_["trkfit"]["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)#[:1000]
                  , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
                  , min_box_coords=min_box_coords, max_box_coords=max_box_coords, box_colour="w"
                  , title="Track cuts", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]"
                  , fout=f"../Images/{recon}/Sanity/h2_trkZX.png")
                  # , log=True)

        ut.Plot2D(x=ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fZ"][data_["crv"]["crvcoincs.sectorType"] == 1], axis=None)#[:1000]
                  , y=ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fX"][data_["crv"]["crvcoincs.sectorType"] == 1], axis=None)#[:1000]
                  , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
                  , min_box_coords=min_box_coords, max_box_coords=max_box_coords, box_colour="w"
                  , title="Track cuts", xlabel="CRV Z [mm]", ylabel="CRV X [mm]"
                  , fout=f"../Images/{recon}/Sanity/h2_crvZX.png")
                  # , log=True)


    # Mark coincidences in the measurement sector for the remaining subset
    MarkCoincidences(data_, coincidenceConditions, quiet)

    # Successful and unsuccessful triggers 
    successes_ = SuccessfulTriggers(data_, success=True, quiet=quiet)
    failures_ = SuccessfulTriggers(data_, success=False, quiet=quiet)

    WriteFailureInfo({triggerMode : failures_}, recon, finTag, foutTag, coincidenceConditions, quiet) 
    WriteResults({triggerMode : (successes_, failures_)}, recon, finTag, foutTag, quiet) 

    return
    
# ------------------------------------------------
#                      Main
# ------------------------------------------------

def TestMain():
    
    fileName="/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00023/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root" 
    
    finTag = fileName.split('.')[-2] 

    with uproot.open(fileName) as file:
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv3_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv2_trigger", quiet=False) # tested
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv3_trigger", quiet=False)vtrk_crv2_2layers_trigger
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_2layers_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_3layers_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv_2layers_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv_3layers_trigger", quiet=False)

        Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_trigger", quiet=False)
        
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="crv_2layers_trigger", quiet=False)

        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv_2layers_trigger", quiet=False)

        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_trigger", quiet=False)
        # Run(file, recon="MDC2020ae", particle="all", PE="10", layer="3", finTag=finTag, triggerMode="trk_crv2_2layers_trigger", quiet=False)


        

    return

def main():

    # TestMain()
    # return

    #########################################################
    # Try multiprocessing, called from an external script
    
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Measure Module Efficiency")

    # Define positional or optional arguments
    parser.add_argument('fileName', type=str, help="Dile name")
    parser.add_argument('recon', type=str, help="Dataset tag (e.g. 'MDC2020ae')")
    parser.add_argument('particles_', type=str, help="List of particles (e.g. ['all', 'muons', 'non_muons']")
    parser.add_argument('PEs_', type=str, help="List of PE thresholds (e.g. [10, 15, ...]'")
    parser.add_argument('layers_', type=str, help="List of layers (e.g. [2, 3])")
    parser.add_argument('triggerModes_', type=str, help="Trigger mode, e.g. ['crv_trigger', 'crv_trk_trigger']")
    parser.add_argument('quiet', type=str, help="Quiet mode flag")
    parser.add_argument('resub', type=str, help="Resubmission mode flag")

    # Parse arguments
    args = parser.parse_args()

    # Convert string arguments back to lists and booleans
    fileName = args.fileName
    recon = args.recon
    particles_ = args.particles_.split(',')
    PEs_ = [int(PE) for PE in args.PEs_.split(',')]
    layers_ = [int(layer) for layer in args.layers_.split(',')]
    triggerModes_ = [str(triggerMode) for triggerMode in args.triggerModes_.split(',')]
    quiet = args.quiet.lower() == 'true'
    resub = args.resub.lower() == 'true'
    
    # Open the file
    # vomsCert not working
    # file = rd.read_file(fileName, quiet)

    with uproot.open(fileName) as file:
        # with uproot.open(fileName) as file:
        finTag = fileName.split('.')[-2]
        if not resub:
            # Scan PE thresholds
            for PE in PEs_: 
                # Scan particles
                for particle in particles_: 
                    # Scan layers
                    for layer in layers_: 
                        # Scan trigger modes
                        for triggerMode in triggerModes_:
                            outputStr = (
                                "\n---> Running with:\n"
                                f"fileName: {fileName}\n"
                                f"recon: {recon}\n"
                                f"particle: {particle}\n"
                                f"PEs: {PE}\n"
                                f"layers: {layer}/4\n"
                                f"finTag: {finTag}\n"
                                f"triggerMode: {triggerMode}\n"
                                f"quiet: {quiet}\n"
                            )
                            if not quiet: print(outputStr) 
                            try:
                                Run(file, recon, particle, PE, layer, finTag, triggerMode, quiet)
                            except Exception as exc:
                                print(f'---> Exception!\n{exc}')
        else: # resubmit jobs
            for PE, particle, layer, triggerMode in zip(PEs_, particles_, layers_, triggerModes_):
                outputStr = (
                    "\n---> Running with:\n"
                    f"fileName: {fileName}\n"
                    f"recon: {recon}\n"
                    f"particle: {particle}\n"
                    f"PEs: {PE}\n"
                    f"layers: {layer}/4\n"
                    f"finTag: {finTag}\n"
                    f"triggerMode: {triggerMode}\n"
                    f"quiet: {quiet}\n"
                )
                if not quiet: print(outputStr) 
                try:
                    Run(file, recon, particle, PE, layer, finTag, triggerMode, quiet)
                except Exception as exc:
                    print(f'---> Exception!\n{exc}')
        return

if __name__ == "__main__":
    main()
    
   #  #########################################################

   #  defname = "nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.root"
   #  recon = "MDC2020ae"
   #  particles_ = ["all"] #, "muons", "non_muons"]
   #  layers_ = [3] #, 2]
   #  PEs_ = np.arange(15, 135, 5)
   #  trkOnly = True
   #  quiet = False
    
   #  def processFunctionA(fileName):
   #      # Always open the file in the processFunction 
   #      file = rd.read_file(fileName, quiet)
   #      # with uproot.open(fileName) as file:
   #      finTag = fileName.split('.')[-2] 
   #      # Scan PE thresholds
   #      for PE in PEs_: 
   #          # Scan particles
   #          for particle in particles_: 
   #              # Scan layers
   #              for layer in layers_: 
   #                  outputStr = (
   #                      "\n---> Running with:\n"
   #                      f"fileName: {fileName}\n"
   #                      f"recon: {recon}\n"
   #                      f"particle: {particle}\n"
   #                      f"PEs: {PE}\n"
   #                      f"layers: {layer}/4\n"
   #                      f"finTag: {finTag}\n"
   #                      f"quiet: {quiet}\n"
   #                  )
   #                  if not quiet: print(outputStr) 
   #                  try:
   #                      Run(file, recon, particle, PE, layer, finTag, trkOnly, quiet)
   #                  except Exception as exc:
   #                      print(f'---> Exception!\n{exc}')
   #                      # Uncomment to test
   #                      # return
    
   #  fileList_ = rd.get_file_list(defname) #[:2]
   #  # fileList_ = ut.ReadFileList("../Txt/FileLists/MDC2020aeOnExpData.txt")
    
   #  print(f"---> Got {len(fileList_)} files.")

   #  # Multithread(processFunctionA, [fileList_[0]])
   #  Multithread(processFunctionA, fileList_)
   #  # Multiprocess(processFunctionA, fileList_)
    
   #  return

   # # Resubmission of failures.
   #  def processFunctionB(fileName):

   #      # Get failed jobs
   #      failedJobsFile = "../Txt/MDC2020ae/FailedJobs/failures.csv"
   #      df_failedJobs_ = pd.read_csv(failedJobsFile)
        
   #      # Always open the file in the processFunction 
   #      # Aren't you wasting resources by opening the file first? 
   #      file = rd.read_file(fileName, quiet)
   #      finTag = fileName.split('.')[-2] 
    
   #      df_failedJobs_ = df_failedJobs_[df_failedJobs_["Tag"] == finTag]
    
   #      # helper to get filter level from dict
   #      # def get_key(d, val):
   #      #   keys = [k for k, v in d.items() if v == val]
   #      #   return keys[0] if keys else None
            
   #      # Submit failed configs
   #      for index, row in df_failedJobs_.iterrows():
    
   #          # print(row["PEs"], row["PEs"], row["Layer"], row["Particle"], row["Filter"]) 
   #          particle = row["Particle"]
   #          PE = row["PEs"]
   #          layer = row["Layer"]
    
   #          outputStr = (
   #              "\n---> Running with:\n"
   #              f"fileName: {fileName}\n"
   #              f"recon: {recon}\n"
   #              f"particle: {particle}\n"
   #              f"PEs: {PE}\n"
   #              f"layers: {layer}/4\n"
   #              f"finTag: {finTag}\n"
   #              f"quiet: {quiet}\n"
   #          )
    
   #          if not quiet: print(outputStr) 
    
   #          try:
   #              Run(file, recon, particle, PE, layer, finTag, quiet)
   #          except Exception as exc:
   #              print(f'---> Exception!\n{row}\n{exc}')
    
   #          # Testing
   #          # return
            
   #      return

   #  fileList_ = rd.get_file_list(defname) #[:2]
    
   #  print(f"---> Got {len(fileList_)} files.")

   #  Multithread(processFunctionB, fileList_)

   #  #########################################################

   #  return
    




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
        
#     # if not quiet: print("Done!")

#     return
