# Samuel Grant 2023
# Test coincidence finding
import sys
import numpy as np
import awkward as ak

import Utils as ut

# ------------------------------------------------
#                Analysis functions 
# ------------------------------------------------

# Get per module coincindences 
def MarkCoincidences(arr_):

    print("---> Getting coincidences")

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
    angleCondition = abs(np.arctan(arr_["crvhit.angle"])) <= 2.0

    # Hit must have (timeEnd - timeStart) <= 15 ms
    print("* Time difference condition")
    arr_["crvhit.timeDiff"] = arr_["crvhit.timeEnd"] - arr_["crvhit.timeStart"]
    timeCondition = arr_["crvhit.timeDiff"] <= (15.0 * 1e6) # ns -> ms

    # Should also be a PE threshold 

    # They should also reside in different layers 

    # Not sure how nHits plays into this

    # Combine conditions to mark per-module coincidences
    coincidenceMask = PECondition & layerCondition & angleCondition & timeCondition 

    # Add a new field 'is_coincidence' to mark coincidences
    arr_["is_coincidence"] = coincidenceMask

    arr_["PE_condition"] = PECondition 
    arr_["layer_condition"] = layerCondition 
    arr_["angle_condition"] = angleCondition 
    arr_["time_condition"] = timeCondition

    print("Done.")

    return arr_ 

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, foutName):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Remove empty arrays (some tracks have no CRV hits)
    # This just slows everything down massively
    # arr_ = RemoveEmptyArrays(arr_)

    # Mark coincidences
    arr_ = MarkCoincidences(arr_)

    print("Writing to", foutName)

    # Wrap this in a function
    # Entries
    with open(foutName, mode='w') as fout:

        # Write headers
        fout.write("evtinfo.eventid, is_coincidence, PE_condition, layer_condition, angle_condition, time_condition, crvhit.nLayers, crvhit.angle, crvhit.sectorType, crvhit.pos.fCoordinates.fX, crvhit.pos.fCoordinates.fY, crvhit.pos.fCoordinates.fZ, crvhit.timeStart, crvhit.timeEnd, crvhit.time, crvhit.PEs, crvhit.nHit\n")

        # Write entries
        for entry in arr_:

            # Check if any of the arrays in the entry are empty
            if any(not entry[field].tolist() for field in entry.fields):
                continue

            # Format the data as a CSV line
            line = ",".join([str(entry[field]) for field in entry.fields]) + "\n"
            fout.write(line)

    return


    # Make printout

    # Headers
    # print("evtinfo.eventid, is_coincidence, PE_condition, layer_condition, angle_condition, time_condition, crvhit.nLayers, crvhit.angle, crvhit.sectorType, crvhit.pos.fCoordinates, crvhit.timeStart, crvhit.timeEnd, crvhit.time, crvhit.PEs, crvhit.nHit")

    # Entries
    for i, entry in enumerate(arr_):

        # Check if any of the arrays in the entry are empty
        if any(not entry[field].tolist() for field in entry.fields):
            continue

        # Format the data as a CSV line
        # line = f"{entry['evtinfo.eventid']}, {entry['is_coincidence']}, {entry['PE_condition']}, {entry['layer_condition']}, {entry['angle_condition']}, {entry['time_condition']}, {entry['crvhit.nLayers']}, {entry['crvhit.angle']}, {entry['crvhit.sectorType']}, {entry['crvhit.pos.fCoordinates.fX']}, {entry['crvhit.pos.fCoordinates.fY']}, {entry['crvhit.pos.fCoordinates.fZ']}, {entry['crvhit.timeStart']}, {entry['crvhit.timeEnd']}, {entry['crvhit.time']}, {entry['crvhit.PEs']}, {entry['crvhit.nHits']}"
    
        print(
        # f"\n"
        f"evtinfo.eventid: {entry['evtinfo.eventid']}\n"
        f"is_coincidence: {entry['is_coincidence']}\n"
        f"PE_condition: {entry['PE_condition']}\n"
        f"layer_condition: {entry['layer_condition']}\n"
        f"angle_condition: {entry['angle_condition']}\n"
        f"time_condition: {entry['time_condition']}\n"
        f"crvhit.nLayers {entry['crvhit.nLayers']}\n"
        f"crvhit.angle: {entry['crvhit.angle']}\n"
        f"crvhit.sectorType: {entry['crvhit.sectorType']}\n"
        f"crvhit.pos.fCoordinates: ({entry['crvhit.pos.fCoordinates.fX']}, {entry['crvhit.pos.fCoordinates.fY']}, {entry['crvhit.pos.fCoordinates.fZ']})\n"
        f"crvhit.timeStart: {entry['crvhit.timeStart']}\n"
        f"crvhit.timeEnd: {entry['crvhit.timeEnd']}\n"
        f"crvhit.time: {entry['crvhit.time']}\n"
        f"crvhit.PEs: {entry['crvhit.PEs']}\n"
        f"crvhit.nHit: {entry['crvhit.nHits']}\n"
        )

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
    
    # # Take output file name command-line argument
    # if len(sys.argv) != 3:
    #     print("Input and outname file names required as arguments")
    #     sys.exit(1)
    
    finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka" # sys.argv[1] # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    foutName = "test.txt" # sys.argv[2] 
    # treeName = "TrkAnaExt/trkana"

    Run(finName, foutName) 

    return

if __name__ == "__main__":
    main()
