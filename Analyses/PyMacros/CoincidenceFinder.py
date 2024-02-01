"""
Samuel Grant
Jan 2024

* Take input from simulation CRV in KPP configuration (three sectors, sectors 3 1 2 top to bottom, where 3 and 2 are the trigger sectors)
* Read TrkAna ROOT trees to Python Awkward arrays 
* ID coincidences in all sectors, based on certain conditions 
* Filter and simplify the dataset, based on requirements
* Calculate the inefficiency and store failed event IDs for event display 

Pass 1: 

* Filter events that have more than one coincidence (passing or failing) in the trigger sectors 
* Calculate the inefficiency as (no successful coincidence in sector 1) / (# successful trigger)
# Write out the events that fulfill this condition
# Will need to adjust to run over multiple files

"""

import sys
import numpy as np
import awkward as ak
import h5py
import os
from collections import Counter
from itertools import combinations

import Utils as ut

# ------------------------------------------------
#                Helper functions 
# ------------------------------------------------

# Stolen from Yuri
def GetPosition3D(data_, branch="crvhit"):
    
    pos_ = ak.zip({"px": ak.flatten(data_['%s.pos.fCoordinates.fX'%branch]), 
                    "py": ak.flatten(data_['%s.pos.fCoordinates.fY'%branch]), 
                    "pz": ak.flatten(data_['%s.pos.fCoordinates.fZ'%branch]),}, with_name="Position3D")

    return pos_

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(data_):
    
    print("\n---> Making sanity plots")

    sectors = ak.flatten(data_["crvhit.sectorType"])
    times = ak.flatten(data_["crvhit.time"])
    startTimes = ak.flatten(data_["crvhit.timeStart"])
    endTimes = ak.flatten(data_["crvhit.timeEnd"])

    ut.Plot1DOverlayOriginal([times[sectors == 3], times[sectors == 1], times[sectors == 2]], nbins=1000, xmin = np.min(times), xmax = np.max(times), xlabel="Average hit time [ns]", ylabel="Hits", labels = ["Sector 3", "Sector 1", "Sector 2"], fout="../Images/Sanity/h1_overlay_time.png")
    ut.Plot1DOverlayOriginal([startTimes[sectors == 3], startTimes[sectors == 1], startTimes[sectors == 2]], nbins=1000, xmin = np.min(startTimes), xmax = np.max(startTimes), xlabel="Hit start time [ns]", ylabel="Hits", labels = ["Sector 3", "Sector 1", "Sector 2"], fout="../Images/Sanity/h1_overlay_startTime.png")
    ut.Plot1DOverlayOriginal([endTimes[sectors == 3], endTimes[sectors == 1], endTimes[sectors == 2]], nbins=1000, xmin = np.min(endTimes), xmax = np.max(endTimes), xlabel="Hit end time [ns]", ylabel="Hits", labels = ["Sector 3", "Sector 1", "Sector 2"], fout="../Images/Sanity/h1_overlay_endTime.png")

    pos_ = GetPosition3D(data_)

    ut.PlotGraph(pos_["px"], pos_["py"], xlabel="x-position [mm]", ylabel="y-position [mm]", fout="../Images/Sanity/gr_XY.png")
    ut.PlotGraph(pos_["pz"], pos_["py"], xlabel="z-position [mm]", ylabel="y-position [mm]", fout="../Images/Sanity/gr_ZY.png")

    print("...Done!")

    return

def PrintEvent(event):
    
    print(
        f"\n"
        f"evtinfo.eventid: {event['evtinfo.eventid']}\n"
        f"is_coincidence: {event['is_coincidence']}\n"
        f"PE_condition: {event['PE_condition']}\n"
        f"layer_condition: {event['layer_condition']}\n"
        f"angle_condition: {event['angle_condition']}\n"
        f"time_condition: {event['time_condition']}\n"
        f"crvhit.nLayers {event['crvhit.nLayers']}\n"
        f"crvhit.angle: {event['crvhit.angle']}\n"
        f"crvhit.sectorType: {event['crvhit.sectorType']}\n"
        f"crvhit.pos.fCoordinates: ({event['crvhit.pos.fCoordinates.fX']}, {event['crvhit.pos.fCoordinates.fY']}, {event['crvhit.pos.fCoordinates.fZ']})\n"
        f"crvhit.timeStart: {event['crvhit.timeStart']}\n"
        f"crvhit.timeEnd: {event['crvhit.timeEnd']}\n"
        f"crvhit.time: {event['crvhit.time']}\n"
        f"crvhit.PEs: {event['crvhit.PEs']}\n"
        f"crvhit.nHits: {event['crvhit.nHits']}\n"
        # f"cluster_ID: {event['cluster_ID']}\n"
    )

    return

def PrintNEvents(data_, nEvents=10):

     # Iterate event-by-event
    for i, event in enumerate(data_):
        
        # print(f"\n ----> Event {i}:")
        PrintEvent(event)

        if i >= nEvents: 
            return

def PrintSpecificEvent(data_, eventID=0): 

    return
# ------------------------------------------------
#               Coincidence finding 
# ------------------------------------------------

# Get coincindences 
# These will need some tuning! 
def MarkCoincidences(data_, debug=False): # debug is just a placeholder here
    
    print("\n---> Marking coincidences")

    # PE threshold of 10 for now
    print("* PE threshold condition")
    PE_ = data_["crvhit.PEs"]
    PECondition = PE_ >= 10

    # Hit must have >=3 layers hits per module
    print("* Layers condition")
    nLayers_ = data_["crvhit.nLayers"]
    layerCondition = nLayers_ >= 3 

    # Hit must not have an angle greater than tan(2)
    print("* Angle condition")
    # "Angle" is just the slope
    angleCondition = abs(data_["crvhit.angle"]) <= 2.0

    # Hit must have (timeEnd - timeStart) <= 15 ms
    print("* Time difference condition")
    data_["crvhit.timeDiff"] = data_["crvhit.timeEnd"] - data_["crvhit.timeStart"]
    timeCondition = data_["crvhit.timeDiff"] <= (15.0 * 1e6) # ns -> ms

    # Combine conditions to mark per-module coincidences
    coincidenceMask = PECondition & layerCondition & angleCondition & timeCondition 

    # Add a new field 'is_coincidence' to mark coincidences
    data_["is_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    data_["PE_condition"] = PECondition 
    data_["layer_condition"] = layerCondition 
    data_["angle_condition"] = angleCondition 
    data_["time_condition"] = timeCondition

    print("...Done!")

    return data_

# ------------------------------------------------
#                     Filter 
# ------------------------------------------------   

# Filter coincidences 
def Filter(data_, filterCondition, debug):
    
    print(f"\n---> Filtering event with condition {filterCondition}")

    if filterCondition == "no_filter": 
        return data_
    else:
        oneHitInSector1 = ak.sum(data_["crvhit.sectorType"] == 1, axis=1) == 1
        oneHitInSector2 = ak.sum(data_["crvhit.sectorType"] == 2, axis=1) == 1
        oneHitInSector3 = ak.sum(data_["crvhit.sectorType"] == 3, axis=1) == 1
        if filterCondition == "one_coincidence_per_sector":
            return data_[oneHitInSector1 & oneHitInSector2 & oneHitInSector3]
        elif filterCondition == "one_coincidence_per_trigger_sector":
            return data_[oneHitInSector2 & oneHitInSector3]
        else: 
            raise ValueError(f"!!! Invalid filter condition {filterCondition} !!!") 

    print("...Done!")

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   
 
def Trigger(data_, debug):

    print(f"\n---> Triggering (at least one passing coincidence in each trigger sector)")
    # Pass 1: at least one = one 

    # We should ID distinct triggers within each event as well...

    # Ensure at least one passing coincidence in each trigger sector
    triggerCondition = (
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 2), axis=1) &
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 3), axis=1)
    )

    print("...Done!")

    return data_[triggerCondition]

# Needs to come after the initial trigger 
def SuccessfulTriggers(data_, success, debug):

    successStr = ""
    if success: successStr += "successful"
    else: successStr += "unsuccessful"

    print(f"\n---> Getting {successStr} triggers")

    # Pass 1: is there a true coincidence in sector 1?

    # We should ID distinct triggers within each event as well...

    # Ensure at least one passing coincidence in each trigger sector
    successCondition = (
        ak.any(data_["is_coincidence"] & (data_["crvhit.sectorType"] == 1), axis=1) 
    )

    print("...Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

def WriteToAwkd(data_, foutName="out.awkd"):

    print(f"\n---> Writing {foutName}")

    with open(foutName, "wb") as fout:
        ak.to_parquet(fout, data_)

    print("...Done!")

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, filterCondition, sanityPlots, debug):

    # Get data as a set of awkward arrays
    data_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # PrintSpecificEvent(data_, eventID=8830)

    # return

    # Sanity plots 
    if (sanityPlots): SanityPlots(data_)

    # Mark coincidences
    data_ = MarkCoincidences(data_)

    # Filter dataset 
    data_ = Filter(data_, filterCondition, debug)

    # Trigger
    data_ = Trigger(data_, debug)

    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True, debug=debug)
    failures_ = SuccessfulTriggers(data_, success=False, debug=debug)

    if debug: PrintNEvents(failures_)
    PrintNEvents(failures_, len(failures_))

    print("\n****************************************************")

    print("File:", finName)

    print("Number of failures:", len(failures_))
    print("Number of successes:", len(successes_))

    tot = len(data_) # after filtering

    print(f"Efficiency: {len(successes_)}/{tot} = {len(successes_)/tot*100:.2f}%")

    inefficiency = (len(failures_) / tot) * 100
    print(f"Inefficiency: {len(failures_)}/{tot} = {inefficiency:.2f}%")

    print("****************************************************\n")

    # Write to counts to file

    # Strip path name and extension from input file string 
    # foutNameSuccess = "../awkd/"+filterCondition+"/success_"+os.path.splitext(os.path.basename(finName))[0]+".awkd"  
    # foutNameFailures = "../awkd/"+filterCondition+"/failures_"+os.path.splitext(os.path.basename(finName))[0]+".awkd"  
    # # Write
    # WriteToAwkd(successes_, foutNameSuccess)
    # WriteToAwkd(failures_, foutNameFailures)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
  
    # # Take input file name command-line argument
    # if len(sys.argv) != 3:
    #     # print("Input and outname file names required as arguments")
    #     print("Error: filter condition and input file name required as argument, example:")
    #     print("python CoincidenceFinder.py one_coincidence_per_trigger_sector /pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka")
    #     sys.exit(1)
    
    # Take input file name command-line argument
    if len(sys.argv) != 2:
        # print("Input and outname file names required as arguments")
        print("Input file name required as argument, example:")
        print("python CoincidenceFinder.py /pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka")
        sys.exit(1)

    # filterCondition = sys.argv[1]
    finName = sys.argv[1] 
    # finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka" # sys.argv[1] # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
 
    # filterCondition = "no_filter"
    # filterCondition = "one_coincidence_per_sector"
    filterCondition = "one_coincidence_per_trigger_sector"
    sanityPlots = False
    debug = False

    # I haven't thought of a sensible way to store the output from each file, so just do the whole dataset in one step
    # h5 doesn't work well with awkward arrays
    # awkd just doesn't seem to work at all
    # Maybe as the analysis matures I will figure out exactly what I need to store, and can be more selective 
    # Or could just fill histograms 
    #  /pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/dc/2e/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00002077.tka
    finName_ = [
        "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/2a/d4/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00001415.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/e4/1b/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00001697.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/dc/2e/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00002077.tka"
        # ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/8f/d1/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00004890.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/0e/04/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00005361.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/6a/e1/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00005673.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/58/a1/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00006238.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/3b/37/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00009634.tka"
        ,"/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/2e/9c/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00009752.tka"
        ]

    # for finName in finName_:
    #     Run(finName=finName, filterCondition=filterCondition, sanityPlots=sanityPlots, debug=debug) 

    Run(finName=finName, filterCondition=filterCondition, sanityPlots=sanityPlots, debug=debug) 

    return

if __name__ == "__main__":
    main()
