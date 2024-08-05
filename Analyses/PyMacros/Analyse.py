'''
Samuel Grant 2024.
Analyse CRV efficiency from TrkAna in KPP geometry.

KPP: three tiers, top and bottom tiers are trigger tiers, middle is the measurement tier.
* Two CRV-L-end modules parallel to the DS axis (sector 3).
* Four CRV-T modules perpendicular to the DS axis (sector 1).
* Two CRV-DS modules perpendicular to the DS axis (sector 2).

Procedure: 
 
1. Find coincidences in all sectors for each unique eventID/runid/subrunid. TrkAna V4 crvcoincs are already concidences, but we can tighten the conditions if desired.
2. Filter the dataset by particle and the number of cosmics per trigger sector -- simplify things and remove the need for clustering.
3. Trigger on cosmics in the top and bottom tier.
4. Check for successful and unsuccessful triggers. 
5. Calculate the inefficiency and write the failures to file. 

'''

import sys
import numpy as np
import awkward as ak

import Utils as ut
import CoincidenceConditions as cc

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(data_, reproc, foutTag):
    
    print("\n---> Making sanity plots")

    # Reco
    sectors_ = ak.flatten(data_[ut.coincsBranchName+".sectorType"]) 
    t_ = ak.flatten(data_[ut.coincsBranchName+".time"]) 
    x_ = ak.flatten(data_[ut.coincsBranchName+".pos.fCoordinates.fX"]) * 1e-3 # mm -> m
    y_ = ak.flatten(data_[ut.coincsBranchName+".pos.fCoordinates.fY"]) * 1e-3 # mm -> m
    z_ = ak.flatten(data_[ut.coincsBranchName+".pos.fCoordinates.fZ"]) * 1e-3 # mm -> m
    PEs_ = ak.flatten(data_[ut.coincsBranchName+".PEs"]) 
    PEsPerLayer_ = ak.flatten(data_[ut.coincsBranchName+".PEsPerLayer[4]"])
    nHits_ = ak.flatten(data_[ut.coincsBranchName+".nHits"]) 
    nLayers_ = ak.flatten(data_[ut.coincsBranchName+".nLayers"]) 
    slopes_ = ak.flatten(data_[ut.coincsBranchName+".angle"])

    ut.PlotGraphOverlay(graphs_=[(x_[sectors_ == 3], y_[sectors_ == 3]), (x_[sectors_ == 1], y_[sectors_ == 1]), (x_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="x-position [m]", ylabel="y-position [m]", fout=f"../Images/{reproc}/Sanity/gr_XY_{foutTag}.png")
    ut.PlotGraphOverlay(graphs_=[(z_[sectors_ == 3], y_[sectors_ == 3]), (z_[sectors_ == 1], y_[sectors_ == 1]), (z_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="z-position [m]", ylabel="y-position [m]", fout=f"../Images/{reproc}/Sanity/gr_ZY_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[t_[sectors_ == 3], t_[sectors_ == 1], t_[sectors_ == 2]], nbins=1000, xmin = np.min(t_), xmax = np.max(t_), xlabel="Average hit time [ns]", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_times_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nHits_[sectors_ == 3], nHits_[sectors_ == 1], nHits_[sectors_ == 2]], nbins=41, xmin = -0.5, xmax = 40.5, xlabel="Number of hits", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_nHits_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nLayers_[sectors_ == 3], nLayers_[sectors_ == 1], nLayers_[sectors_ == 2]], nbins=5, xmin = -0.5, xmax = 4.5, xlabel="Number of layers hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_nLayers_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[slopes_[sectors_ == 3], slopes_[sectors_ == 1], slopes_[sectors_ == 2]], nbins=1000, xmin = -2, xmax = 2, xlabel="Slope", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_slopes_{foutTag}.png") 

    # Calculate the sum of each length-4 array in PEsPerLayer_
    PEsPerLayerSum_ = ak.sum(PEsPerLayer_, axis=-1)
    PEsDiff_ = PEs_ - PEsPerLayerSum_

    ut.Plot1DOverlay(hists_=[PEs_[sectors_ == 3], PEs_[sectors_ == 1], PEs_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEs_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsPerLayerSum_[sectors_ == 3], PEsPerLayerSum_[sectors_ == 1], PEsPerLayerSum_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs per layer", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEsPerLayer_{foutTag}.png")

    ut.Plot1DOverlay(hists_=[PEsPerLayer_[:, 0], PEsPerLayer_[:, 1], PEsPerLayer_[:, 2], PEsPerLayer_[:, 3]], nbins=500, xmin = 0, xmax = 500, title="PEs per layer", ylabel="Coincidences", label_=["Layer0", "Layer 1", "Layer 2", "Layer 3"], fout=f"../Images/{reproc}/Sanity//h1_PEsPerLayerAllSectors_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsDiff_[sectors_ == 3], PEsDiff_[sectors_ == 1], PEsDiff_[sectors_ == 2]], nbins=40, xmin = -0.0002, xmax = 0.0002, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEsDiff_{foutTag}.png", logY=True)
    ut.Plot1DOverlay(hists_=[PEsDiff_], nbins=200, xmin = np.min(PEsDiff_)-1e-4, xmax = np.max(PEsDiff_)+1e-4, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["All modules"], fout=f"../Images/{reproc}/Sanity/h1_PEsDiffAll_{foutTag}.png", logY=True, includeBlack=True)

    # MC 
    valid_ = ak.flatten(data_[ut.coincsBranchName+"mc.valid"])
    pdgid_ = ak.flatten(data_[ut.coincsBranchName+"mc.pdgId"])

    ut.Plot1DOverlay(hists_=[valid_[sectors_ == 3], valid_[sectors_ == 1], valid_[sectors_ == 2]], nbins=2, xmin = 0, xmax = 2, xlabel="Valid MC event", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_valid_{foutTag}.png") 
    
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

    ut.BarChart(data_=pdgid_, label_dict=ut.particle_dict, ylabel="Coincidences [%]", fout=f"../Images/{reproc}/Sanity/bar_pdgid_{foutTag}.png", percentage=True)
    ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout=f"../Images/{reproc}/Sanity/bar_overlay_pdgid_{foutTag}.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

    print("...Done!")

    return

def PrintEvent(event, coincidenceConditions):

    coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]
    
    eventStr = (
        f"evtinfo.run: {event['evtinfo.run']}\n" 
        f"evtinfo.subrun: {event['evtinfo.subrun']}\n" 
        f"evtinfo.eventid: {event['evtinfo.event']}\n"
        f"coincidenceConditions_['PEthreshold']: {coincidenceConditions_['PEthreshold']}\n"
        f"coincidenceConditions_['nLayers']: {coincidenceConditions_['nLayers']}\n"
        f"coincidenceConditions_['minSlope']: {coincidenceConditions_['minSlope']}\n"
        f"coincidenceConditions_['maxSlope']: {coincidenceConditions_['maxSlope']}\n"
        f"is_coincidence: {event['is_coincidence']}\n"
        f"PE_condition: {event['PE_condition']}\n"
        f"layer_condition: {event['layer_condition']}\n"
        f"angle_condition: {event['angle_condition']}\n"
        f"{ut.coincsBranchName}.nLayers {event[f'{ut.coincsBranchName}.nLayers']}\n"
        f"{ut.coincsBranchName}.angle: {event[f'{ut.coincsBranchName}.angle']}\n"
        f"{ut.coincsBranchName}.sectorType: {event[f'{ut.coincsBranchName}.sectorType']}\n"
        f"{ut.coincsBranchName}.pos.fCoordinates: ({event[f'{ut.coincsBranchName}.pos.fCoordinates.fX']}, {event[f'{ut.coincsBranchName}.pos.fCoordinates.fY']}, {event[f'{ut.coincsBranchName}.pos.fCoordinates.fZ']})\n"
        f"{ut.coincsBranchName}.timeStart: {event[f'{ut.coincsBranchName}.timeStart']}\n"
        f"{ut.coincsBranchName}.timeEnd: {event[f'{ut.coincsBranchName}.timeEnd']}\n"
        f"{ut.coincsBranchName}.time: {event[f'{ut.coincsBranchName}.time']}\n"
        f"{ut.coincsBranchName}.PEs: {event[f'{ut.coincsBranchName}.PEs']}\n"
        f"{ut.coincsBranchName}.PEsPerLayer[4]: {event[f'{ut.coincsBranchName}.PEsPerLayer[4]']}\n"
        f"{ut.coincsBranchName}.nHits: {event[f'{ut.coincsBranchName}.nHits']}\n"
        f"{ut.coincsBranchName}mc.valid: {event[f'{ut.coincsBranchName}mc.valid']}\n"
        f"{ut.coincsBranchName}mc.pdgId: {event[f'{ut.coincsBranchName}mc.pdgId']}\n"
    )
    

    return eventStr

# This seems to crash the mu2e node! 
# def PrintEvent(event, coincidenceConditions):

#     coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]
    
#     eventStr = (
#         f"evtinfo.runid: {event['evtinfo.runid']}\n" 
#         f"evtinfo.subrunid: {event['evtinfo.subrunid']}\n" 
#         f"evtinfo.eventid: {event['evtinfo.eventid']}\n"
#         f"coincidenceConditions_['PEthreshold']: {coincidenceConditions_['PEthreshold']}\n"
#         f"coincidenceConditions_['nLayers']: {coincidenceConditions_['nLayers']}\n"
#         f"coincidenceConditions_['minSlope']: {coincidenceConditions_['minSlope']}\n"
#         f"coincidenceConditions_['maxSlope']: {coincidenceConditions_['maxSlope']}\n"
#         f"is_coincidence: {event['is_coincidence']}\n"
#         f"PE_condition: {ak.to_list(event['PE_condition'])}\n"
#         f"layer_condition: {ak.to_list(event['layer_condition'])}\n"
#         f"angle_condition: {ak.to_list(event['angle_condition'])}\n"
#         f"{ut.coincsBranchName}.nLayers {ak.to_list(event[f'{ut.coincsBranchName}.nLayers'])}\n"
#         f"{ut.coincsBranchName}.angle: {ak.to_list(event[f'{ut.coincsBranchName}.angle'])}\n"
#         f"{ut.coincsBranchName}.sectorType: {ak.to_list(event[f'{ut.coincsBranchName}.sectorType'])}\n"
#         f"{ut.coincsBranchName}.pos.fCoordinates: ({ak.to_list(event[f'{ut.coincsBranchName}.pos.fCoordinates.fX'])}, {ak.to_list(event[f'{ut.coincsBranchName}.pos.fCoordinates.fY'])}, {ak.to_list(event[f'{ut.coincsBranchName}.pos.fCoordinates.fZ'])})\n"
#         f"{ut.coincsBranchName}.timeStart: {ak.to_list(event[f'{ut.coincsBranchName}.timeStart'])}\n"
#         f"{ut.coincsBranchName}.timeEnd: {ak.to_list(event[f'{ut.coincsBranchName}.timeEnd'])}\n"
#         f"{ut.coincsBranchName}.time: {ak.to_list(event[f'{ut.coincsBranchName}.time'])}\n"
#         f"{ut.coincsBranchName}.PEs: {ak.to_list(event[f'{ut.coincsBranchName}.PEs'])}\n"
#         f"{ut.coincsBranchName}.PEsPerLayer[4]: {ak.to_list(event[f'{ut.coincsBranchName}.PEsPerLayer[4]'])}\n"
#         f"{ut.coincsBranchName}.nHits: {ak.to_list(event[f'{ut.coincsBranchName}.nHits'])}\n"
#         f"{ut.coincsBranchName}mc.valid: {ak.to_list(event[f'{ut.coincsBranchName}mc.valid'])}\n"
#         f"{ut.coincsBranchName}mc.pdgId: {ak.to_list(event[f'{ut.coincsBranchName}mc.pdgId'])}\n"
#     )
    
#     return eventStr

def PrintRawEvent(event):
    
    eventStr = (
        f"evtinfo.runid: {event['evtinfo.run']}\n"
        f"evtinfo.subrunid: {event['evtinfo.subrun']}\n"
        f"evtinfo.event: {event['evtinfo.event']}\n"
        f"{ut.coincsBranchName}.nLayers {event[f'{ut.coincsBranchName}.nLayers']}\n"
        f"{ut.coincsBranchName}.angle: {event[f'{ut.coincsBranchName}.angle']}\n"
        f"{ut.coincsBranchName}.sectorType: {event[f'{ut.coincsBranchName}.sectorType']}\n"
        f"{ut.coincsBranchName}.pos.fCoordinates: ({event[f'{ut.coincsBranchName}.pos.fCoordinates.fX']}, {event[f'{ut.coincsBranchName}.pos.fCoordinates.fY']}, {event[f'{ut.coincsBranchName}.pos.fCoordinates.fZ']})\n"
        f"{ut.coincsBranchName}.timeStart: {event[f'{ut.coincsBranchName}.timeStart']}\n"
        f"{ut.coincsBranchName}.timeEnd: {event[f'{ut.coincsBranchName}.timeEnd']}\n"
        f"{ut.coincsBranchName}.time: {event[f'{ut.coincsBranchName}.time']}\n"
        f"{ut.coincsBranchName}.PEs: {event[f'{ut.coincsBranchName}.PEs']}\n"
        f"{ut.coincsBranchName}.nHits: {event[f'{ut.coincsBranchName}.nHits']}\n"
        f"{ut.coincsBranchName}mc.valid: {event[f'{ut.coincsBranchName}mc.valid']}\n"
        f"{ut.coincsBranchName}mc.pdgId: {event[f'{ut.coincsBranchName}mc.pdgId']}\n"
    )

    return

def PrintNEvents(data_, nEvents=10, coincidenceConditions="default", raw=False):

     # Iterate event-by-event
    for i, event in enumerate(data_):
        
        if raw: print(PrintRawEvent(event))
        else: print(PrintEvent(event, coincidenceConditions))

        if i >= nEvents: 
            return

# ------------------------------------------------
#               Coincidence finding 
# ------------------------------------------------

# OLD! 
# def FindCoincidences(data_, coincidenceConditions="default"): 
    
#     print("\n---> Marking coincidences")

#     # Get conditions from file
#     coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]

#     # PE threshold per layer
#     print("* PE threshold condition")
#     PE_ = data_[ut.coincsBranchName+".PEs"]
    
#     # PE condition per layer
#     PECondition = PE_ >= coincidenceConditions_["PEthreshold"]

#     # Number of layers hit 
#     # (based on number of coincidences with PEs per layer above threshold, not nLayers which is baked into the reconstruction in mcs)
#     print("* Layers condition")
#     nLayers_ = data_[ut.coincsBranchName+".nLayers"]
#     layerCondition = nLayers_ >= coincidenceConditions_["nLayers"] # This is baked into the reconstruction
#     # layerCondition = ak.count(PE_[PECondition], axis=2) > coincidenceConditions_["nLayers"] # This can be adjusted on the level of nts 

#     # Hit slope: horizontal / vertical direction 
#     print("* Angle condition")
#     angleCondition = (abs(data_[ut.coincsBranchName+".angle"]) <= coincidenceConditions_["maxSlope"]) & (abs(data_[ut.coincsBranchName+".angle"]) >= coincidenceConditions_["minSlope"])

#     # Combine conditions to mark per-module coincidences
#     coincidenceMask = layerCondition & angleCondition & PECondition 

#     # Add a new field 'is_coincidence' to mark coincidences
#     data_["is_coincidence"] = coincidenceMask

#     # Mark the individual coniditions for debugging 
#     data_["PE_condition"] = PECondition 
#     data_["layer_condition"] = layerCondition 
#     data_["angle_condition"] = angleCondition 

#     # print(data_)

#     print("...Done!")

#     return data_

# Just used to impose stricter conditions than the default, if desired
def FindCoincidences(data_, coincidenceConditions="default"): 
    
    print("\n---> Marking coincidences")

    # Get conditions from file
    coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]

    # PE threshold per layer
    print("* PE threshold condition")
    PE_ = data_[ut.coincsBranchName+".PEsPerLayer[4]"]
    
    # PE condition per layer
    PECondition = PE_ >= coincidenceConditions_["PEthreshold"]

    # Number of layers hit 
    # (based on number of coincidences with PEs per layer above threshold, not nLayers which is baked into the reconstruction in mcs)
    print("* Layers condition")
    nLayers_ = data_[ut.coincsBranchName+".nLayers"]
    # layerCondition = nLayers_ >= coincidenceConditions_["nLayers"] # This is baked into the reconstruction

    # You could also sort the array and check if the 3rd or 2nd element is above threshold 
    # Might be more efficient

    layerCondition = ak.count(PE_[PECondition], axis=2) >= coincidenceConditions_["nLayers"] # This can be adjusted on the level of nts 

    # Hit slope: horizontal / vertical direction 
    print("* Angle condition")
    angleCondition = (abs(data_[ut.coincsBranchName+".angle"]) <= coincidenceConditions_["maxSlope"]) & (abs(data_[ut.coincsBranchName+".angle"]) >= coincidenceConditions_["minSlope"])

    # Combine conditions to mark per-module coincidences
    coincidenceMask = layerCondition & angleCondition # PECondition is implicit in the layer condition!

    # Add a new field 'is_coincidence' to mark coincidences
    data_["is_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    data_["PE_condition"] = PECondition 
    data_["layer_condition"] = layerCondition 
    data_["angle_condition"] = angleCondition 

    # print(data_)

    print("...Done!")

    return data_

# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

def RemoveEmptyEvents(data_):

    print(f"\n---> Removing empty events")

    emptyEventCondition = ak.num(data_[ut.coincsBranchName+".nHits"]) == 0
    return data_[~emptyEventCondition]

# Currently just "all", "muons", "non_muons"
# Can generalise if needed
def FilterParticles(data_, particle):

    print(f"\n---> Filtering particles, keeping {particle}")

    # I think this really should be in trigger sectors only FIXME
    muonCondition = ak.any((data_[ut.coincsBranchName+"mc.pdgId"] == 13) | (data_[ut.coincsBranchName+"mc.pdgId"] == -13), axis=1)
    # muonCondition = ak.any(abs(data_[ut.coincsBranchName+"mc.pdgId"]), axis=1) == 13

    if particle == "all":
        return data_
    elif particle == "muons": 
        return data_[muonCondition] 
    elif particle == "non_muons":
        return data_[~muonCondition] 
    else:
        raise ValueError(f"Particle string {particle} not valid!")

# Filter coincidences 
def FilterCoincidences(data_, coincidenceFilter):
    
    print(f"\n---> Filtering event with condition {coincidenceFilter}")

    if coincidenceFilter == "no_filter": 
        return data_
    else:
        oneHitInSector1Condition = ak.sum(data_[ut.coincsBranchName+".sectorType"] == 1, axis=1) == 1
        oneHitInSector2Condition = ak.sum(data_[ut.coincsBranchName+".sectorType"] == 2, axis=1) == 1
        oneHitInSector3Condition = ak.sum(data_[ut.coincsBranchName+".sectorType"] == 3, axis=1) == 1
        if coincidenceFilter == "one_coincidence_per_sector":
            return data_[oneHitInSector1Condition & oneHitInSector2Condition & oneHitInSector3Condition]
        elif coincidenceFilter == "one_coincidence_per_trigger_sector":
            return data_[oneHitInSector2Condition & oneHitInSector3Condition]
        else: 
            raise ValueError(f"!!! Invalid filter condition {coincidenceFilter} !!!") 

    print("...Done!")

# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_): 
    print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_[ut.coincsBranchName+".timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   
 
def Trigger(data_):

    print(f"\n---> Triggering (at least one passing coincidence in each trigger sector)")

    # Enforce trigger condition
    triggerCondition = (
        ak.any(data_["is_coincidence"] & (data_[ut.coincsBranchName+".sectorType"] == 2), axis=1) &
        ak.any(data_["is_coincidence"] & (data_[ut.coincsBranchName+".sectorType"] == 3), axis=1)
    )

    print("...Done!")

    # Mask data
    return data_[triggerCondition]

# Needs to come after the initial trigger 
def SuccessfulTriggers(data_, success):

    successStr = ""
    if success: successStr += "successful"
    else: successStr += "unsuccessful"

    print(f"\n---> Getting {successStr} triggers")

    # Ensure at least one passing coincidence in each trigger sector
    successCondition = (
        ak.any(data_["is_coincidence"] & (data_[ut.coincsBranchName+".sectorType"] == 1), axis=1) 
    )

    print("...Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

# ------------------------------------------------
#                     Output 
# ------------------------------------------------ 

def WriteFailuresToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions):

    # Define the output file path
    foutName = f"../Txt/{reproc}/failures_ntuple/{doutTag}/failures_ntuple_{foutTag}.csv" 

    print(f"\n---> Writing failures to:\n{foutName}") 

    with open(foutName, "w") as fout:
        # Write the header
        header = "\t".join(ut.branchNamesTrkAna_) + "\n"
        fout.write(header)

        for event in failures_:
            data = "\t".join(str(event[name]) for name in ut.branchNamesTrkAna_) + "\n"
            fout.write(data)

    return

def WriteFailureInfoToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions, verbose):

    # Define the output file path
    foutNameConcise = f"../Txt/{reproc}/failures_concise/{doutTag}/failures_concise_{foutTag}.csv" 
    foutNameVerbose = f"../Txt/{reproc}/failures_verbose/{doutTag}/failures_verbose_{foutTag}.csv" 

    print(f"\n---> Writing failure info to:\n{foutNameConcise}") #\n{foutNameVerbose}")
    if verbose: print(foutNameVerbose)

    # Concise form
    with open(foutNameConcise, "w") as fout:
        # Write the header
        fout.write("evtinfo.run, evtinfo.subrun, evtinfo.event\n")
        # Write the events
        for event in failures_:
            fout.write(
                f"{event['evtinfo.run']}, {event['evtinfo.subrun']}, {event['evtinfo.event']}\n"
            )

    # Verbose form
    if verbose: 
        with open(foutNameVerbose, "w") as fout:
            # Write the header
            # fout.write(f"evtinfo.runid,evtinfo.subrunid,evtinfo.eventid,is_coincidence,PE_condition,layer_condition,angle_condition,{ut.coincsBranchName}.nLayers,{ut.coincsBranchName}.angle,{ut.coincsBranchName}.sectorType,{ut.coincsBranchName}.pos.fCoordinates.fX,{ut.coincsBranchName}.pos.fCoordinates.fY,{ut.coincsBranchName}.pos.fCoordinates.fZ,{ut.coincsBranchName}.timeStart,{ut.coincsBranchName}.timeEnd,{ut.coincsBranchName}.time,{ut.coincsBranchName}.PEs,{ut.coincsBranchName}mc.valid,{ut.coincsBranchName}mc.pdgId\n")
            # Write the events
            for event in failures_:
                fout.write(
                    PrintEvent(event, coincidenceConditions)+"\n" #f"{event['evtinfo.runid']},{event['evtinfo.subrunid']},{event['evtinfo.eventid']},{event['is_coincidence']},{event['PE_condition']},{event['layer_condition']},{event['angle_condition']},{event[f'{ut.coincsBranchName}.nLayers']},{event[f'{ut.coincsBranchName}.angle']},{event[f'{ut.coincsBranchName}.sectorType']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fX']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fY']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fZ']},{event[f'{ut.coincsBranchName}.timeStart']},{event[f'{ut.coincsBranchName}.timeEnd']},{event[f'{ut.coincsBranchName}.time']},{event[f'{ut.coincsBranchName}.PEs']},{event[f'{ut.coincsBranchName}.nHits']},{event[f'{ut.coincsBranchName}mc.valid']},{event[f'{ut.coincsBranchName}mc.pdgId']}\n"
                )

    return

def WriteResultsToFile(data_, successes_, failures_, doutTag, foutTag, reproc):

    foutName = f"../Txt/{reproc}/results/{doutTag}/results_{foutTag}.csv"
    print(f"\n---> Writing results to {foutName}")
    tot = len(data_)
    efficiency = len(successes_) / tot * 100
    inefficiency = len(failures_) / tot * 100

    outputStr = f"""
    ****************************************************
    Number of failures: {len(failures_)}
    Number of successes: {len(successes_)}
    Efficiency: {len(successes_)}/{tot} = {efficiency:.2f}%
    Inefficiency: {len(failures_)}/{tot} = {inefficiency:.2f}%
    ****************************************************
    """

    with open(foutName, "w") as fout:
        fout.write("Total, Successes, Failures, Efficiency [%], Inefficiency [%]\n")
        fout.write(f"{tot}, {len(successes_)}, {len(failures_)}, {efficiency}, {inefficiency}\n")
        # fout.write(outputStr)

    print(outputStr)

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, particle, coincidenceConditions, reproc, coincidenceFilter, sanityPlots, verbose):

    # Get data as a set of awkward arrays
    data_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Remove empty events
    data_ = RemoveEmptyEvents(data_)

    # Output dir/file tag 
    doutTag = finName.split('.')[-2] 
    foutTag = particle + "_" + coincidenceConditions + "_" + coincidenceFilter # + "_" +finName.split('.')[-2] 

    # Sanity plots 
    if (sanityPlots): SanityPlots(data_, reproc, foutTag)

    # Filter particles
    data_ = FilterParticles(data_, particle)

    # Apply start time cut
    data_ = CutOnStartTime(data_)

    # Find coincidences
    data_ = FindCoincidences(data_, coincidenceConditions)

    # Useful debugging tool
    # PrintNEvents(data_, 100, coincidenceConditions)

    # Filter dataset 
    data_ = FilterCoincidences(data_, coincidenceFilter)

    # Trigger
    data_ = Trigger(data_)

    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True)
    failures_ = SuccessfulTriggers(data_, success=False)

    # PrintNEvents(successes_, 10, coincidenceConditions)

    #return 

    # Write failures to file
    WriteFailuresToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions) # write ntuple to table
    WriteFailureInfoToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions, verbose) #  

    # Write results to file
    WriteResultsToFile(data_, successes_, failures_, doutTag, foutTag, reproc) 

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    # Take command-line arguments
    finName = sys.argv[1] if len(sys.argv) > 1 else "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatTriggered/MDC2020ae_best_v1_3/root/c4/15/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000231.root" # "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/40/73/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # /pnfs/mu2e/scratch/users/sgrant/workflow/CosmicCRYExtractedTrk.MDC2020z2_best_v1_1/outstage/67605881/00/00000/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    particle = sys.argv[2] if len(sys.argv) > 2 else "all"
    coincidenceConditions = sys.argv[3] if len(sys.argv) > 3 else "10PEs2Layers" # "ana1" # "default"
    reproc = sys.argv[4] if len(sys.argv) > 4 else "MDC2020ae" # "original"
    coincidenceFilter = sys.argv[5] if len(sys.argv) > 5 else "one_coincidence_per_trigger_sector"
    sanityPlots = bool(sys.argv[6]) if len(sys.argv) > 6 else False
    verbose = bool(sys.argv[7]) if len(sys.argv) > 7 else True # False

    print("\n---> Running with inputs:\n")
    print("\tfinName:", finName)
    print("\tparticle:", particle)
    print("\tcoincidenceConditions:", particle)
    print("\tcoincidenceFilter:", coincidenceFilter)
    print("\tsanityPlots:", sanityPlots)
    print("\tverbose:", verbose, "\n")

    Run(finName=finName, particle=particle, coincidenceConditions=coincidenceConditions, reproc=reproc, coincidenceFilter=coincidenceFilter, sanityPlots=sanityPlots, verbose=verbose) 

    return

if __name__ == "__main__":
    main()
