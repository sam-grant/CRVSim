# Samuel Grant 2023
# Test coincidence finding
import sys
import numpy as np
import awkward as ak
import h5py

import Utils as ut

# ------------------------------------------------
#                Helper functions 
# ------------------------------------------------

def GetPosition3D(arr_, branch="crvhit"):

    pos_ = ak.zip({"px": ak.flatten(arr_['%s.pos.fCoordinates.fX'%branch]), 
                    "py": ak.flatten(arr_['%s.pos.fCoordinates.fY'%branch]), 
                    "pz": ak.flatten(arr_['%s.pos.fCoordinates.fZ'%branch]),}, with_name="Position3D")

    return pos_

# ------------------------------------------------
#                Analysis functions 
# ------------------------------------------------

def SanityPlots(arr_):

    print("---> Making sanity plots")

    pos_ = GetPosition3D(arr_)

    ut.PlotGraph(pos_["px"], pos_["py"], xlabel="x-position [mm]", ylabel="y-position [mm]", fout="Images/Sanity/gr_XY.png")
    ut.PlotGraph(pos_["pz"], pos_["py"], xlabel="z-position [mm]", ylabel="y-position [mm]", fout="Images/Sanity/gr_ZY.png")

    print("Done!")

    return

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

# Print coincidences in human readable format
def PrintCoincidence(arr_, n):

    for i, entry in enumerate(arr_):

        # Check if any of the arrays in the entry are empty
        if any(not entry[field].tolist() for field in entry.fields):
            continue

        print(
        f"\n"
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

        if(i == n):
            break

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Sanity plots 
    # SanityPlots(arr_)

    # return

    # Mark coincidences
    arr_ = MarkCoincidences(arr_)

    # Printout
    PrintCoincidence(arr_, 100)

    # Write to file
    # WriteCoincidencesToHDF5(arr_, foutName)

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
    # foutName = "test.h5" # sys.argv[2] 
    # treeName = "TrkAnaExt/trkana"

    Run(finName) 

    return

if __name__ == "__main__":
    main()
