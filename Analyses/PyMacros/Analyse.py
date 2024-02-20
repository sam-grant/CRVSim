'''
Samuel Grant 2024
Analyse CRV efficiency from TrkAna in KPP geometry 

KPP: three tiers, top and bottom tiers are trigger tiers, middle is the measurement tier
* Two CRV-D modules with scintillator bars parallel to the DS axis (sector 3) 
* Four CRV-T modules with scintillator bars perpendicular to the DS axis (sector 1) 
* Two CRV side end cap modules with scintillator bars perpendicular to the DS axis (sector 2)

Procedure: 
 
1. Find coincidences in all sectors for each unique eventID/runid/subrunid. TrkAna V4 crvhits are already concidences, but we can tighten the conditions if desired.
2. Filter the dataset by particle and the number of cosmics per trigger sector -- simplify things and remove the need for clustering.
3. Trigger on cosmics in the top and bottom tier.
4. Check for successful and unsuccessful triggers. 
5. Calculate the inefficiency and write the failures to file. 

'''

import sys
import numpy as np
import awkward as ak

# from collections import Counter
# from itertools import combinations

import Utils as ut
import CoincidenceConditions as cc

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(data_):
    
    print("\n---> Making sanity plots")

    # Reco
    sectors_ = ak.flatten(data_["crvhit.sectorType"]) 
    t_ = ak.flatten(data_["crvhit.time"]) 
    x_ = ak.flatten(data_["crvhit.pos.fCoordinates.fX"]) * 1e-3 # mm -> m
    y_ = ak.flatten(data_["crvhit.pos.fCoordinates.fY"]) * 1e-3 # mm -> m
    z_ = ak.flatten(data_["crvhit.pos.fCoordinates.fZ"]) * 1e-3 # mm -> m
    PEs_ = ak.flatten(data_["crvhit.PEs"]) 
    nHits_ = ak.flatten(data_["crvhit.nHits"]) 
    nLayers_ = ak.flatten(data_["crvhit.nLayers"]) 
    slopes_ = ak.flatten(data_["crvhit.angle"])
    # PEsPerHit_ = PEs_ / nHits_ 

    ut.PlotGraphOverlay(graphs_=[(x_[sectors_ == 3], y_[sectors_ == 3]), (x_[sectors_ == 1], y_[sectors_ == 1]), (x_[sectors_ == 2], y_[sectors_ == 2]) ], labels_=["Top", "Middle", "Bottom"], xlabel="x-position [m]", ylabel="y-position [m]", fout="../Images/Sanity/gr_XY.png")
    ut.PlotGraphOverlay(graphs_=[(z_[sectors_ == 3], y_[sectors_ == 3]), (z_[sectors_ == 1], y_[sectors_ == 1]), (z_[sectors_ == 2], y_[sectors_ == 2]) ], labels_=["Top", "Middle", "Bottom"], xlabel="z-position [m]", ylabel="y-position [m]", fout="../Images/Sanity/gr_ZY.png")
    ut.Plot1DOverlay(hists_=[t_[sectors_ == 3], t_[sectors_ == 1], t_[sectors_ == 2]], nbins=1000, xmin = np.min(t_), xmax = np.max(t_), xlabel="Average hit time [ns]", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_times.png") 
    ut.Plot1DOverlay(hists_=[PEs_[sectors_ == 3], PEs_[sectors_ == 1], PEs_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_PEs.png") 
    # ut.Plot1DOverlay(hists_=[PEsPerHit_[sectors_ == 3], PEsPerHit_[sectors_ == 1], PEsPerHit_[sectors_ == 2]], nbins=100, xmin = 0, xmax = 100, xlabel="Average PEs per hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_PEs_per_hit.png") 
    ut.Plot1DOverlay(hists_=[nHits_[sectors_ == 3], nHits_[sectors_ == 1], nHits_[sectors_ == 2]], nbins=41, xmin = -0.5, xmax = 40.5, xlabel="Number of hits", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_nHits.png") 
    ut.Plot1DOverlay(hists_=[nLayers_[sectors_ == 3], nLayers_[sectors_ == 1], nLayers_[sectors_ == 2]], nbins=5, xmin = -0.5, xmax = 4.5, xlabel="Number of layers hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_nLayers.png") 
    ut.Plot1DOverlay(hists_=[slopes_[sectors_ == 3], slopes_[sectors_ == 1], slopes_[sectors_ == 2]], nbins=1000, xmin = -2, xmax = 2, xlabel="Slope", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_slopes.png") 

    # MC 
    valid_ = ak.flatten(data_["crvhitmc.valid"])
    pdgid_ = ak.flatten(data_["crvhitmc.pdgId"])

    ut.Plot1DOverlay(hists_=[valid_[sectors_ == 3], valid_[sectors_ == 1], valid_[sectors_ == 2]], nbins=2, xmin = 0, xmax = 2, xlabel="Valid MC event", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout="../Images/Sanity/h1_valid.png") 
    
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

    ut.BarChart(data_=pdgid_, label_dict=ut.particle_dict, ylabel="Coincidences [%]", fout="../Images/Sanity/bar_pdgid.png", percentage=True)
    ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout="../Images/Sanity/bar_overlay_pdgid.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

    print("...Done!")

    return

def PrintEvent(event):
    
    print(
        f"evtinfo.runid: {event['evtinfo.runid']}\n"
        f"evtinfo.subrunid: {event['evtinfo.subrunid']}\n"
        f"evtinfo.eventid: {event['evtinfo.eventid']}\n"
        f"is_coincidence: {event['is_coincidence']}\n"
        f"PE_condition: {event['PE_condition']}\n"
        f"layer_condition: {event['layer_condition']}\n"
        f"angle_condition: {event['angle_condition']}\n"
        f"crvhit.nLayers {event['crvhit.nLayers']}\n"
        f"crvhit.angle: {event['crvhit.angle']}\n"
        f"crvhit.sectorType: {event['crvhit.sectorType']}\n"
        f"crvhit.pos.fCoordinates: ({event['crvhit.pos.fCoordinates.fX']}, {event['crvhit.pos.fCoordinates.fY']}, {event['crvhit.pos.fCoordinates.fZ']})\n"
        f"crvhit.timeStart: {event['crvhit.timeStart']}\n"
        f"crvhit.timeEnd: {event['crvhit.timeEnd']}\n"
        f"crvhit.time: {event['crvhit.time']}\n"
        f"crvhit.PEs: {event['crvhit.PEs']}\n"
        f"crvhit.nHits: {event['crvhit.nHits']}\n"
        f"crvhitmc.valid: {event['crvhitmc.valid']}\n"
        f"crvhitmc.pdgId: {event['crvhitmc.pdgId']}\n"
    )

    return

def PrintRawEvent(event):
    
    print(
        f"evtinfo.runid: {event['evtinfo.runid']}\n"
        f"evtinfo.subrunid: {event['evtinfo.subrunid']}\n"
        f"evtinfo.eventid: {event['evtinfo.eventid']}\n"
        f"crvhit.nLayers {event['crvhit.nLayers']}\n"
        f"crvhit.angle: {event['crvhit.angle']}\n"
        f"crvhit.sectorType: {event['crvhit.sectorType']}\n"
        f"crvhit.pos.fCoordinates: ({event['crvhit.pos.fCoordinates.fX']}, {event['crvhit.pos.fCoordinates.fY']}, {event['crvhit.pos.fCoordinates.fZ']})\n"
        f"crvhit.timeStart: {event['crvhit.timeStart']}\n"
        f"crvhit.timeEnd: {event['crvhit.timeEnd']}\n"
        f"crvhit.time: {event['crvhit.time']}\n"
        f"crvhit.PEs: {event['crvhit.PEs']}\n"
        f"crvhit.nHits: {event['crvhit.nHits']}\n"
        f"crvhitmc.valid: {event['crvhitmc.valid']}\n"
        f"crvhitmc.pdgId: {event['crvhitmc.pdgId']}\n"
    )

    return

def PrintNEvents(data_, nEvents=10, raw=False):

     # Iterate event-by-event
    for i, event in enumerate(data_):
        
        if raw: PrintRawEvent(event)
        else: PrintEvent(event)

        if i >= nEvents: 
            return

# ------------------------------------------------
#               Coincidence finding 
# ------------------------------------------------

# Just used to impose stricter conditions than the default, if desired
def FindCoincidences(data_, coincidenceConditions="default"): # debug is just a placeholder here
    
    print("\n---> Marking coincidences")

    # Get conditions from file
    coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]

    # PE threshold
    print("* PE threshold condition")
    PE_ = data_["crvhit.PEs"]
    PECondition = PE_ >= coincidenceConditions_["PEthreshold"]

    # Layers hit 
    print("* Layers condition")
    nLayers_ = data_["crvhit.nLayers"]
    layerCondition = nLayers_ >= coincidenceConditions_["nLayers"]

    # Hit slope: horizontal / vertical direction 
    print("* Angle condition")
    angleCondition = (abs(data_["crvhit.angle"]) <= coincidenceConditions_["maxSlope"]) & (abs(data_["crvhit.angle"]) >= coincidenceConditions_["minSlope"])

    # We do not have access to time difference in TrkAna... 
    # Hit must have (timeEnd - timeStart) <= 15 ns
    # print("* Time difference condition")
    # data_["crvhit.timeDiff"] = data_["crvhit.timeEnd"] - data_["crvhit.timeStart"]
    # timeCondition = abs(data_["crvhit.timeDiff"]) <= 20 # 20 # (15.0 * 1e6) # ns -> ms

    # Combine conditions to mark per-module coincidences
    coincidenceMask = PECondition & layerCondition & angleCondition # & timeCondition 

    # Add a new field 'is_coincidence' to mark coincidences
    data_["is_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    data_["PE_condition"] = PECondition 
    data_["layer_condition"] = layerCondition 
    data_["angle_condition"] = angleCondition 
    # data_["time_condition"] = timeCondition 

    print("...Done!")

    return data_

# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

def RemoveEmptyEvents(data_):

    print(f"\n---> Removing empty events")

    emptyEventCondition = ak.num(data_["crvhit.nHits"]) == 0
    return data_[~emptyEventCondition]

# Currently just "all", "muons", "non_muons"
# Can generalise if needed
def FilterParticles(data_, particle):

    print(f"\n---> Filtering particles, keeping {particle}")

    # I think this really should be in trigger sectors only FIXME
    muonCondition = ak.any((data_["crvhitmc.pdgId"] == 13) | (data_["crvhitmc.pdgId"] == -13), axis=1)

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
        oneHitInSector1Condition = ak.sum(data_["crvhit.sectorType"] == 1, axis=1) == 1
        oneHitInSector2Condition = ak.sum(data_["crvhit.sectorType"] == 2, axis=1) == 1
        oneHitInSector3Condition = ak.sum(data_["crvhit.sectorType"] == 3, axis=1) == 1
        if coincidenceFilter == "one_coincidence_per_sector":
            return data_[oneHitInSector1Condition & oneHitInSector2Condition & oneHitInSector3Condition]
        elif coincidenceFilter == "one_coincidence_per_trigger_sector":
            return data_[oneHitInSector2Condition & oneHitInSector3Condition]
        else: 
            raise ValueError(f"!!! Invalid filter condition {coincidenceFilter} !!!") 

    print("...Done!")

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   
 
def Trigger(data_):

    print(f"\n---> Triggering (at least one passing coincidence in each trigger sector)")

    # Enforce trigger condition
    triggerCondition = (
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 2), axis=1) &
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 3), axis=1)
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
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 1), axis=1) 
    )

    print("...Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

# ------------------------------------------------
#                     Output 
# ------------------------------------------------ 

def WriteFailuresToFile(failures_, foutTag):

    # Define the output file path
    foutNameConcise = "../Txt/failures_concise_" + foutTag + ".csv" 
    foutNameVerbose = "../Txt/failures_verbose_" + foutTag + ".csv" 

    print(f"\n---> Writing failures to:\n{foutNameConcise}\n{foutNameVerbose}")

    # Concise form
    with open(foutNameConcise, "w") as fout:
        # Write the header
        fout.write("evtinfo.runid, evtinfo.subrunid, evtinfo.eventid\n")
        # Write the events
        for event in failures_:
            fout.write(
                f"{event['evtinfo.runid']}, {event['evtinfo.subrunid']}, {event['evtinfo.eventid']}\n"
            )

    # Verbose form
    with open(foutNameVerbose, "w") as fout:
        # Write the header
        fout.write("evtinfo.runid,evtinfo.subrunid,evtinfo.eventid,is_coincidence,PE_condition,layer_condition,angle_condition,crvhit.nLayers,crvhit.angle,crvhit.sectorType,crvhit.pos.fCoordinates.fX,crvhit.pos.fCoordinates.fY,crvhit.pos.fCoordinates.fZ,crvhit.timeStart,crvhit.timeEnd,crvhit.time,crvhit.PEs,crvhitmc.valid,crvhitmc.pdgId\n")
        # Write the events
        for event in failures_:
            fout.write(
                f"{event['evtinfo.runid']},{event['evtinfo.subrunid']},{event['evtinfo.eventid']},{event['is_coincidence']},{event['PE_condition']},{event['layer_condition']},{event['angle_condition']},{event['crvhit.nLayers']},{event['crvhit.angle']},{event['crvhit.sectorType']},{event['crvhit.pos.fCoordinates.fX']},{event['crvhit.pos.fCoordinates.fY']},{event['crvhit.pos.fCoordinates.fZ']},{event['crvhit.timeStart']},{event['crvhit.timeEnd']},{event['crvhit.time']},{event['crvhit.PEs']},{event['crvhit.nHits']},{event['crvhitmc.valid']},{event['crvhitmc.pdgId']}\n"
            )

    return

def WriteResultsToFile(data_, successes_, failures_, foutTag):

    foutName = f"../Txt/results_{foutTag}.csv"
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

def Run(finName, particle, coincidenceConditions, coincidenceFilter, sanityPlots):

    # Get data as a set of awkward arrays
    data_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Remove empty events
    data_ = RemoveEmptyEvents(data_)

    # Sanity plots 
    if (sanityPlots): SanityPlots(data_)

    # Filter particles
    data_ = FilterParticles(data_, particle)

    # Find coincidences
    data_ = FindCoincidences(data_, coincidenceConditions)

    # Filter dataset 
    data_ = FilterCoincidences(data_, coincidenceFilter)

    # Trigger
    data_ = Trigger(data_)

    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True)
    failures_ = SuccessfulTriggers(data_, success=False)

    # Output file tag 
    foutTag = particle + "_" + coincidenceConditions + "_" + coincidenceFilter + "_" +finName.split('.')[-2] 

    # Write failures to file
    WriteFailuresToFile(failures_, foutTag)

    # Write results to file
    WriteResultsToFile(data_, successes_, failures_, foutTag) 

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    # Take command-line arguments
    finName = sys.argv[1] if len(sys.argv) > 1 else "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    particle = sys.argv[2] if len(sys.argv) > 2 else "all"
    coincidenceConditions = sys.argv[3] if len(sys.argv) > 3 else "default"
    coincidenceFilter = sys.argv[4] if len(sys.argv) > 4 else "one_coincidence_per_trigger_sector"
    sanityPlots = bool(sys.argv[5]) if len(sys.argv) > 5 else False

    print("\n--->Running with inputs:\n")
    print("\tfinName:", finName)
    print("\tparticle:", particle)
    print("\tcoincidenceConditions:", particle)
    print("\tcoincidenceFilter:", coincidenceFilter)
    print("\tsanityPlots:", sanityPlots, "\n")

    Run(finName=finName, particle=particle, coincidenceConditions=coincidenceConditions, coincidenceFilter=coincidenceFilter, sanityPlots=sanityPlots) 

    return

if __name__ == "__main__":
    main()
