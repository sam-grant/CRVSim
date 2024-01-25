# Samuel Grant 2024
# Count CRV KPP coincidences
# Output lists of event IDs for true/false coincidences per sector
# Filter conditions: no_filter, one_coincidence_per_sector, one_coincidence_per_trigger_sector

import sys
import numpy as np
import awkward as ak
import h5py
import os
from collections import Counter

import Utils as ut

# ------------------------------------------------
#                Helper functions 
# ------------------------------------------------

# Stolen from Yuri
def GetPosition3D(arr_, branch="crvhit"):
    
    pos_ = ak.zip({"px": ak.flatten(arr_['%s.pos.fCoordinates.fX'%branch]), 
                    "py": ak.flatten(arr_['%s.pos.fCoordinates.fY'%branch]), 
                    "pz": ak.flatten(arr_['%s.pos.fCoordinates.fZ'%branch]),}, with_name="Position3D")

    return pos_

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(arr_):
    
    print("\n---> Making sanity plots")

    pos_ = GetPosition3D(arr_)

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
        f"crvhit.nHit: {event['crvhit.nHits']}\n"
    )

    return

# ------------------------------------------------
#      Coincidence finding and filtering
# ------------------------------------------------

# Get per module coincindences 
# These will need some tuning 
def MarkCoincidences(arr_, debug=False): # debug is just a placeholder here
    
    print("\n---> Marking coincidences")

    # PE threshold of 10 for now
    print("* PE threshold condition")
    PE_ = arr_["crvhit.PEs"]
    PECondition = PE_ >= 10

    # Hit must have >=3 layers hits per module
    print("* Layers condition")
    nLayers_ = arr_["crvhit.nLayers"]
    layerCondition = nLayers_ >= 3 

    # Hit must not have an angle greater than tan(2)
    print("* Angle condition")
    # angleCondition = abs(np.arctan(arr_["crvhit.angle"])) <= 2.0
    # "Angle" is just the slope
    angleCondition = abs(arr_["crvhit.angle"]) <= 2.0

    # Hit must have (timeEnd - timeStart) <= 15 ms
    print("* Time difference condition")
    arr_["crvhit.timeDiff"] = arr_["crvhit.timeEnd"] - arr_["crvhit.timeStart"]
    timeCondition = arr_["crvhit.timeDiff"] <= (15.0 * 1e6) # ns -> ms

    # Combine conditions to mark per-module coincidences
    coincidenceMask = PECondition & layerCondition & angleCondition & timeCondition 

    # Add a new field 'is_coincidence' to mark coincidences
    arr_["is_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    arr_["PE_condition"] = PECondition 
    arr_["layer_condition"] = layerCondition 
    arr_["angle_condition"] = angleCondition 
    arr_["time_condition"] = timeCondition

    print("...Done!")

    return arr_

# Filter coincidences 
def Filter(event, filterCondition, sectors_, debug):
    
    if debug: print(f"\n---> Filtering event with condition {filterCondition}")

    if filterCondition == "no_filter": 
        return False
    elif filterCondition == "one_coincidence_per_sector":
        # More than one instance of any sector in the event
        if len(sectors_) != len(set(sectors_)): 
            return True
    elif filterCondition == "one_coincidence_per_trigger_sector":
        # More than one instance of sectors 2 or 3 in the event
        if np.count_nonzero(sectors_ == 2) > 1 or np.count_nonzero(sectors_ == 3) > 1: 
            return True
    else: 
        print("!!! Error: invalid filterCondition. !!!")
        return True        

    return False

def CountCoincidences(arr_, filterCondition, debug):

    print("\n---> Counting coincidences.")
    
    # Setup counting dicts, element is a list of event IDS
    trueCoincidences_ = { "S1" : [], "S2" : [],  "S3" : [] }
    falseCoincidences_ = { "S1" : [], "S2" : [],  "S3" : [] }

    # Start iterating through coincidences
    

    # For tracking progress inside the loop
    totEvents = len(arr_)

    # Iterate event-by-event
    for i, event in enumerate(arr_):

        if debug: 
            print("\n****************************")
            print(f"---> Next event i: {i}.")

        # Coincidences
        coincidences_ = event["is_coincidence"]
        # Event IDs 
        eventID = event["evtinfo.eventid"]
        # Sectors
        sectors_ = event["crvhit.sectorType"]

        # Skip empty events (ones which miss the CRV)
        # nCoin = len(event["is_coincidence"])
        if len(coincidences_) < 1: 
            if debug: print("---> Skipping empty event.")
            continue

        if debug: PrintEvent(event)
        
        # Filtering events
        # 0. "no_filter"
        # 1. "one_coincidence_per_sector"
        # 2. "one_coincidence_per_trigger_sector" 
        if Filter(event, filterCondition=filterCondition, sectors_=sectors_, debug=debug):  
            if debug: print("---> Sector list is", sectors_, "... filtering!")
            continue

        # Start counting 
        for j, coincidence in enumerate(coincidences_):
            sectorKey = "S" + str(sectors_[j])

            if coincidence == True:
                if debug: print(f"\n---> True coincidence at index {j}")
                trueCoincidences_[sectorKey].append(eventID)
            elif coincidence == False:
                if debug: print(f"\n---> False coincidence at index {j}")
                falseCoincidences_[sectorKey].append(eventID)
            else: 
                # This should never happen
                print("!!! Error, coincidence not found !!!")
                break

        progress = (i + 1) / totEvents * 100
        print(f"Progress: {progress:.2f}%", end='\r', flush=True)

        # Testing
        # if progress > 5.0: 
        #     break

    print("\n...Done!")

    print("\nNumber of True Coincidences:")
    for key, value in trueCoincidences_.items():
        print(f"{key}: {len(value)}")

    print("\nNumber of False Coincidences:")
    for key, value in falseCoincidences_.items():
        print(f"{key}: {len(value)}")

    return trueCoincidences_, falseCoincidences_

def WriteCoincidences(trueCoincidences_, falseCoincidences_, foutName, debug):

    print(f"\n---> Writing to h5.")

    # Write to HDF5 file
    with h5py.File(foutName, "w") as file:
        # Write data 
        for key, value in trueCoincidences_.items():
            file.create_dataset(f'trueCoincidences/{key}', data=np.array(value))
        # Write data from dict2
        for key, value in falseCoincidences_.items():
            file.create_dataset(f'falseCoincidences/{key}', data=np.array(value))
    print("...Done")    

    return

def CheckWrite(trueCoincidences_, falseCoincidences_, foutName, debug):

    print("\n---> Checking write!")

    # Read the h5 file
    readDict1 = {}
    readDict2 = {}

    with h5py.File(foutName, 'r') as file:
        # Read data into dict1
        for key in file['trueCoincidences']:
            readDict1[key] = list(file[f'trueCoincidences/{key}'])
        # Read data into dict2
        for key in file['falseCoincidences']:
            readDict2[key] = list(file[f'falseCoincidences/{key}'])

    if readDict1 != trueCoincidences_ or readDict2 != falseCoincidences_:
        print("\n!!! Error writing coincidences to file !!!")
        return
    if readDict2 != falseCoincidences_:
        print("\n!!! Error writing falseCoincidences_ to file !!!")

    print("...Done!")

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, filterCondition="no_filter", sanityPlots=False, debug=False):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Sanity plots 
    if (sanityPlots): SanityPlots(arr_)

    # Mark coincidences
    arr_ = MarkCoincidences(arr_)

    # Count coincidences 
    trueCoincidences_, falseCoincidences_ = CountCoincidences(arr_, filterCondition, debug)

    # Write to counts to h5

    # Strip path name and extension from input file  string (input file name without path and extension)
    foutName = "../h5/TrueAndFalseCoincidences/"+os.path.splitext(os.path.basename(finName))[0]+".h5"
    # Write
    WriteCoincidences(trueCoincidences_, falseCoincidences_, foutName, debug)
    # Check that the write worked as expected
    CheckWrite(trueCoincidences_, falseCoincidences_, foutName, debug)

    print(f"\n---> Output file is:\n{foutName}")

    # Analysis is done in CoincidenceAna.py

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
  
    # Take input file name command-line argument
    if len(sys.argv) != 3:
        # print("Input and outname file names required as arguments")
        print("Error: filter condition and input file name required as argument, example:")
        print("python CoincidenceFinder.py one_coincidence_per_trigger_sector /pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka")
        sys.exit(1)
    
    filterCondition = sys.argv[1]
    finName = sys.argv[2] 

    # finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka" # sys.argv[1] # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    # filterCondition = "no_filter"
    # filterCondition = "one_coincidence_per_sector"
    # filterCondition = "one_coincidence_per_trigger_sector"

    Run(finName=finName, filterCondition=filterCondition, debug=False) 

    return

if __name__ == "__main__":
    main()
