# Samuel Grant 2023
# Test coincidence finding

import numpy as np
import awkward as ak

import Utils as ut

# ------------------------------------------------
#                Analysis functions 
# ------------------------------------------------

# Get per module coincindences 
def MarkCoincidences(arr_):

    print("---> Getting coincidences")
    
    # Hit must have >=3 layers hits per module
    print("* Layers condition")
    nLayers_ = ak.num(arr_["crvhit.nLayers"])  
    layerCondition = nLayers_ >= 3 

    # Hit must not have an angle greater than 2 rad
    print("* Angle condition")
    angleCondition = np.arctan(arr_["crvhit.angle"]) <= 2.0

    # Hit must have (timeEnd - timeStart) <= 15 ms
    print("* Time difference condition")
    arr_["crvhit.timeDiff"] = arr_["crvhit.timeEnd"] - arr_["crvhit.timeStart"]
    timeCondition = arr_["crvhit.timeDiff"] <= (15.0 * 1e6) # ns -> ms

    # Should also be a PE threshold 

    # They should also reside in different layers 

    # Not sure how nHits plays into this

    # Combine conditions to mark per-module coincidences
    coincidenceMask = layerCondition & angleCondition & timeCondition

    # Add a new field 'is_coincidence' to mark coincidences
    arr_["is_coincidence"] = coincidenceMask

    print("Done.")

    return arr_ 

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, treeName):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Remove empty arrays (some tracks have no CRV hits)
    # This just slows everything down massively
    # arr_ = RemoveEmptyArrays(arr_)

    # Mark coincidences
    arr_ = MarkCoincidences(arr_)

    for i, entry in enumerate(arr_):

        # Check if any of the arrays in the entry is empty
        if any(not entry[field].tolist() for field in entry.fields):
            continue
    
        print(
        f"\n"
        f"evtinfo.eventid: {entry['evtinfo.eventid']}\n"
        f"is_coincidence: {entry['is_coincidence']}\n"
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

        if (i > 100): break

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
    
    finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    treeName = "TrkAnaExt/trkana"

    Run(finName, treeName) 

    return

if __name__ == "__main__":
    main()
