'''
Samuel Grant 2024

Write PEs per layer information to file before analysing the single layer inefficiency. 

Do this so I can parallelise, although one file seems to be fine to me. 

'''

import sys
import numpy as np
import uproot
import awkward as ak

import Utils as ut

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, particle, dataset): 

    # Load data.
    print("---> Loading data...")
    data_ = ak.Array([])
    with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
        data_ = tree.arrays(["crvcoincs.sectorType", "crvcoincs.PEsPerLayer[4]", "crvcoincsmc.pdgId"]) 
    # print("Done.")

    # Remove empty events.
    print("---> Removing empty events...")
    emptyEventCondition = ak.num(data_["crvcoincs.sectorType"], axis=1) == 0
    data_ = data_[~emptyEventCondition]
    # print("Done.")

    # Filter particles 
    print(f"---> Filtering particles ({particle})...")
    muonCondition = ak.any((data_[ut.coincsBranchName+"mc.pdgId"] == 13) | (data_[ut.coincsBranchName+"mc.pdgId"] == -13), axis=1)
    if (particle=="muons"):
        data_ = data_[muonCondition]
    elif (particle=="non_muons"):
        data_ = data_[~muonCondition]

    # print("PDG ID:", data_["crvcoincsmc.pdgId"])

    # Output dir/file tag 
    doutTag = finName.split('.')[-2] 
    foutTag = particle 

    # Flatten arrays
    print("---> Flattening arrays.")
    sectors_ = ak.flatten(data_["crvcoincs.sectorType"]) 
    PEsPerLayer_ = ak.flatten(data_["crvcoincs.PEsPerLayer[4]"])

    results = { 
        "Sector1" : PEsPerLayer_[sectors_ == 1] 
        ,"Sector2" : PEsPerLayer_[sectors_ == 2] 
        ,"Sector3" : PEsPerLayer_[sectors_ == 3] 
    }

    # Write data
    # TODO wrap this in a function 

    import h5py

    foutName = f"../Txt/{dataset}/PEsPerLayer/{doutTag}/PEsPerLayer_{foutTag}.h5"
    print(f"---> Writing to {foutName}.")
    # Create a new HDF5 file
    with h5py.File(foutName, "w") as fout:
        # Iterate over the sectors in the results dictionary
        for sector, data in results.items():
            # Create a group for each sector
            sector_group = fout.create_group(sector)
            # Store the data for each sector
            fout.create_dataset(sector + "/PEsPerLayer", data=data)

    print(f"---> Written {foutName}\n")

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
    finName="/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatTriggered/MDC2020ae_best_v1_3/root/c4/15/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000231.root" 
    dataset="MDC2020ae"

    # Run(finName=finName, particle="all", dataset=dataset)
    Run(finName=finName, particle="muons", dataset=dataset)
    Run(finName=finName, particle="non_muons", dataset=dataset)
    
    return

    # Take command-line arguments finName="/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00033/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000080.root"
    finName = sys.argv[1] if len(sys.argv) > 1 else "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatTriggered/MDC2020ae_best_v1_3/root/c4/15/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000231.root" # "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/40/73/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" 
    particle = sys.argv[2] if len(sys.argv) > 2 else "all"
    dataset = sys.argv[3] if len(sys.argv) > 3 else "MDC2020ae" 

    print("\n--->Running with inputs:\n")
    print("\tfinName:", finName)
    print("\tparticle:", particle)
    print("\tdataset:", dataset)

    Run(finName=finName, particle=particle, dataset=dataset) 

    return

if __name__ == "__main__":
    main()
