import numpy as np
import pandas as pd
import awkward as ak
import random

import Utils as ut

# ------------------------------------------------
#              Sanity check functions
# ------------------------------------------------
def FilterEmptyArrays(arr_):

    return arr_[len(arr_["crvhit.nLayers"] > 0)] 


def PrintAwkwardInfo(arr_, indent=""):

    print(f"{indent}---> Awkward arrays info")

    if len(arr_) > 0:

        print(f"{indent}Length: {len(arr_)}")

        # Access the keys of the first dictionary to get branch names
        branchNames_ = ak.fields(arr_[0])

        # Print the branch names and types
        for branchName in branchNames_:
            print(f"\nBranch name: {branchName}")
            print(f"Array: {arr_[branchName]}")
            # print(f"Array length: {len(arr_[branchName])}")
            PrintType(ak.type(arr_[0][branchName]))

    else:
        print(f"{indent}No events in the Awkward Array.")

    print("\n")

    return

def PrintType(type_):
    if isinstance(type_, ak.types.RecordType):
        # print(f"{indent}RecordType:")
        for field in type_.fields:
            PrintType(type_.field(field))
    elif isinstance(type_, ak.types.ArrayType):
        # print(f"Sub-branches: {len(type_.fields)})")
        PrintType(type_.content)
    else:
        print(f"Type: {type_}")

def SampleCoincidences(arr_, num_samples=10):

    print("\n---> Sampling coincidences")

    # Get n random indices
    random_indices = random.sample(range(len(arr_)), num_samples)

    # Print corresponding entries
    for i in random_indices:
        entry = arr_[i]
        timeDiff = entry['crvhit.timeEnd'] - entry['crvhit.timeStart']
        # time_diff = [end - start for start, end in zip(entry['crvhit.timeStart'], entry['crvhit.timeEnd'])]
        print(f"\nEntry {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {timeDiff}")


    # Mark coincidences and filter entries
    # arr_ = MarkCoincidences(arr_)

    # # Separate coincidences and non-coincidences
    # coincidences = arr_[arr_["is_coincidence"]]
    # non_coincidences = arr_[~arr_["is_coincidence"]]

    # # Convert Awkward Array to regular Python list
    # coincidences = ak.to_list(coincidences)
    # non_coincidences = ak.to_list(non_coincidences)

    # # Sample 10 coincidences and 10 non-coincidences
    # sampled_coincidences = random.sample(coincidences, min(num_samples, len(coincidences)))
    # sampled_non_coincidences = random.sample(non_coincidences, min(num_samples, len(non_coincidences)))

    # # Print sampled coincidences
    # print("Sampled Coincidences:")
    # for i, entry in enumerate(sampled_coincidences):
    #     print(f"\nSampled Coincidence {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")

    # # Print sampled non-coincidences
    # # print("\nSampled Non-Coincidences:")
    # # for i, entry in enumerate(sampled_non_coincidences):
    # #     print(f"\nSampled Non-Coincidence {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")


    return  

    # # Sample 10 coincidences and 10 non-coincidences
    # sampled_coincidences = random.sample(coincidences, min(num_samples, len(coincidences)))
    # sampled_non_coincidences = random.sample(non_coincidences, min(num_samples, len(non_coincidences)))

    # # Print sampled coincidences
    # print("Sampled Coincidences:")
    # for i, entry in enumerate(sampled_coincidences):
    #     print(f"\nSampled Coincidence {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")

    # # Print sampled non-coincidences
    # print("\nSampled Non-Coincidences:")
    # for i, entry in enumerate(sampled_non_coincidences):
    #     print(f"\nSampled Non-Coincidence {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")

# ------------------------------------------------
#                Analysis functions 
# ------------------------------------------------

def GetPosition3D(arr_, branch="crvhit"):

    pos_ = ak.zip({"px": ak.flatten(arr_['%s.pos.fCoordinates.fX'%branch]), 
                    "py": ak.flatten(arr_['%s.pos.fCoordinates.fY'%branch]), 
                    "pz": ak.flatten(arr_['%s.pos.fCoordinates.fZ'%branch]),}, with_name="Position3D")

    return pos_

# Get per module coincindences 
def MarkCoincidences(arr_):

    print("---> Getting coincidences")
    
    # Hit must have >=3 layers hits per module
    print("* Layers condition")
    nLayers_ = ak.num(arr_["crvhit.nLayers"])  
    layerCondition = nLayers_ >= 3 

    # Hit must not have an angle greater than 2 rad
    print("* Angle condition")
    angleCondition = arr_["crvhit.angle"] <= 2.0

    # Hit must have (timeEnd - timeStart) <= 15 ms
    print("* Time difference condition")
    arr_["crvhit.timeDiff"] = arr_["crvhit.timeEnd"] - arr_["crvhit.timeStart"]
    timeCondition = arr_["crvhit.timeDiff"]*1e-6 <= 15.0 

    # Combine conditions to mark per-module coincidences
    coincidenceMask = layerCondition & angleCondition & timeCondition

    # Add a new field 'is_coincidence' to mark coincidences
    arr_["is_coincidence"] = coincidenceMask

    # Drop empty arrays 
    arr_ = arr_[arr_["is_coincidence"]]

    return arr_ # [coincidenceMask]

# ------------------------------------------------
#                       Run
# ------------------------------------------------
def Run(finName, treeName):

    # Get data as a set of awkward arrays
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNamesTrkAna_)

    # Print array info (sanity check)
    PrintAwkwardInfo(arr_)

    arr_ = FilterEmptyArrays(arr_)

    # Get coincidences 
    arr_ = MarkCoincidences(arr_[0])

    for i, entry in enumerate(arr_):
        
        print(f"\nEntry {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")


        # is_coincidence = entry['is_coincidence']

        # if (all(is_coincidence)): 
        # # if (entry['is_coincidence'][0]):
        #     print(f"\nEntry {i+1}\n* Coincidence: {entry['is_coincidence']}\n* Number of Layers: {entry['crvhit.nLayers']}\n* Angle: {entry['crvhit.angle']}\n* Time diff: {entry['crvhit.timeEnd'] - entry['crvhit.timeStart']}")


    # Print out a sample of coincidences 
    # SampleCoincidences(arr_[arr_["is_coincidence"]], 10)
    # SampleCoincidences(arr_[~arr_["is_coincidence"]], 10)


    # print(ak.flatten(arr_["is_coincidence"]))
    # print(ak.flatten(arr_[arr_["is_coincidence"]]["crvhit.time"]))
    # Three of more hits where evtinfo.eventid is the same, but crvhit.sectorType is different hit

    # print(arr_["crvhit.time"])
    # print(ak.flatten(arr_["crvhit.time"]))

    return
    

    # Get time, position, momentum 
    # eventID_ = arr_["evtinfo.eventid"]
    # time_ = ak.flatten(arr_["crvhit.time"])
    # pos_ = GetPosition3D(arr_)

    # print(len(eventID_), len(time_))

    # print(eventID_[:10], time_[:10])

    # return

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
