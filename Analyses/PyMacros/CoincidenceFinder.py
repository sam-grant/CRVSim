# Samuel Grant 2024
# Find CRV KPP coincidences
# Trigger on coincidences in the top and bottom sectors
# Count the number of coincidences in the middle sector / number of triggers
# Ignore events with more than one coincience in the same trigger sector 

import sys
import numpy as np
import awkward as ak
import h5py
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

    print("Done!")

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
def MarkCoincidences(arr_):

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

    print("Done!")

    return arr_

# Filter coincidences 
def Filter(event, filterCondition, verbose=False):

    if verbose: print(f"\n---> Filter event with condition {filterCondition}")

    if filterCondition=="no_filter": 
        return False

    # return False

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, filterCondition="no_filter", sanityPlots=False, verbose=False):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Sanity plots 
    if (sanityPlots): SanityPlots(arr_)

    # Mark coincidences
    arr_ = MarkCoincidences(arr_)

    # Setup counting dicts, element is a list of event IDS
    trueCoincidences_ = { "S1" : [], "S2" : [],  "S3" : [] }
    falseCoincidences_ = { "S1" : [], "S2" : [],  "S3" : [] }

    # Start iterating through coincidences
    print("\n---> Iterating through coincidences.\n")

    # For tracking progress inside the loop
    totEvents = len(arr_)

    # Iterate event-by-event
    for i, event in enumerate(arr_):

        if verbose: 
            print("\n****************************\n")
            print(f"---> Next event i: {i}.")
            PrintEvent(event)

        # Get coincidence list
        coincidences_ = event["is_coincidence"]

        # Skip empty events (ones which miss the CRV)
        nCoin = len(event["is_coincidence"])
        if nCoin < 1: 
            if verbose: print("---> Skipping empty event.")
            continue

        # Filtering events
        # 0. "no_filter"
        # 1. "one_coincidence_per_sector"
        # 2. "one_coincidence_per_trigger_sector" 
        # if Filter(event, filterCondition=filterCondition, verbose=verbose):
        #     continue

        # Event IDs and sectors and 
        eventID = event["evtinfo.eventid"]
        sectors_ = event["crvhit.sectorType"]

        # Start counting 
        for j, coincidence in enumerate(coincidences_):
            sectorKey = "S" + str(sectors_[j])

            if coincidence == True:
                if verbose: print(f"\n---> True coincidence at index {j}")
                trueCoincidences_[sectorKey].append(eventID)
            elif coincidence == False:
                if verbose: print(f"\n---> False coincidence at index {j}")
                falseCoincidences_[sectorKey].append(eventID)
            else: 
                # This should never happen
                print("!!! Error, no coincidence found !!!")
                break

        progress = (i + 1) / totEvents * 100
        print(f"Progress: {progress:.2f}%", end='\r', flush=True)

    print("Done!")

    if verbose:
        print("\ntrueCoincidences_:\n", trueCoincidences_)
        print("falseCoincidences_:\n", falseCoincidences_)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
    
    # Take input file name command-line argument
    # if len(sys.argv) != 2:
    #     # print("Input and outname file names required as arguments")
    #     print("Input file name required as argument")
    #     sys.exit(1)
    
    # finName = sys.argv[1] # 
    finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka" # sys.argv[1] # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"

    Run(finName=finName, verbose=False) 

    return

if __name__ == "__main__":
    main()



# How often do you get a coincidence in sectors 2 & 3 when you also have on in sector 1
# We need to handle multiple hits

def FilterAndCountCoincidences(arr_):

    print("\n---> Filtering and counting coincidences")

    totEvents = len(arr_)

    # Dictionary to count coincidences event-by-event

    # Count coincidences event-by-event, append the event IDs to lists 
    singles_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }
    doubles_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }
    triples_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }   
    
    # Count false coincidences 
    falses_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }  

    # Iterate event-by-event
    for i, entry in enumerate(arr_):

        # Check if any of the arrays in the entry are empty
        # This just slows everything down!
        # if any(not entry[field].tolist() for field in entry.fields):
        #     continue

        # Number of coincidences, event ID and sectors
        nCoin = len(entry["is_coincidence"])
        eventID = entry['evtinfo.eventid']
        sectors_ = entry["crvhit.sectorType"]

        # Check for false coincidences
        if any(entry["is_coincidence"] == False):
            print("\n!!! False coincidence found !!!")
            PrintEntry(entry)
            # Store them 
            for sector in sectors_:  
                falses_["sector_"+str(sector)].append(eventID)
            # If there are any, we can cycle back later


        # Handle multiple coincidences per trigger sector 
        sectorCounts = Counter(sectors_)
        # print(f"\nEvent {eventID}")
        for sector, n in sectorCounts.items():
            if n>1: 
                print(f"Sector {sector} has {n} counts")
                PrintEntry(entry)

        # break
        
        continue


        # You've been getting away with it because most of them are true.


        # Count single coincidences
        if (nCoin == 1):
            for sector in sector_: # this should always be the top sector 
                singles_["sector_"+str(sector)].append(eventID)
        # Count double coincidences 
        elif (len(entry["is_coincidence"]) == 2):
            for sector in sector_:
                doubles_["sector_"+str(sector)].append(eventID)
        # Count triple coincidences
        elif (len(entry["is_coincidence"]) == 3):
            for sector in sector_:
                triples_["sector_"+str(sector)].append(eventID)
        else: 
            # Is this the correct place to be doing this?

            beyond_.append(eventID)
            print(
                f"*** WARNING: CountCoincidences() ***\n"
            )

            PrintEntry(entry)

            # print(
            #     f"*** WARNING: CountCoincidences() ***\n"
            #     f"Event ID: {eventID}\n"
            #     f"Number of coincidences: {len(entry['is_coincidence'])}\n"
            #     f"Coincidences: {entry['is_coincidence']}\n"
            #     f"Entry: {entry}"
            # )

        progress = (i + 1) / totEvents * 100
        print(f"Progress: {progress:.2f}%", end='\r', flush=True)

        # if (i > 500): break

    # print(
    #     f"\nSingles: {singles_}\n"
    #     f"Doubles: {doubles_}\n"
    #     f"Triples: {triples_}\n"
    # )

    # Make bar charts 
    # TODO: handle this better!
    singlesCounts_ = { "sector_1" : len(singles_["sector_1"]), "sector_2" : len(singles_["sector_2"]),  "sector_3" : len(singles_["sector_3"]) }
    doublesCounts_ = { "sector_1" : len(doubles_["sector_1"]), "sector_2" : len(doubles_["sector_2"]),  "sector_3" : len(doubles_["sector_3"]) }
    triplesCounts_ = { "sector_1" : len(triples_["sector_1"]), "sector_2" : len(triples_["sector_2"]),  "sector_3" : len(triples_["sector_3"]) }

    ut.BarChart2(data_dict=singlesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_singles.png")
    ut.BarChart2(data_dict=doublesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_doubles.png")
    ut.BarChart2(data_dict=triplesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_triples.png")

    print("Done!")

    return
