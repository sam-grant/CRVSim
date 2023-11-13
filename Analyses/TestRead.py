import uproot
import pandas as pd

branchNames_ = [
    # ---> evtinfo
    "evtinfo.eventid"
    , "evtinfo.runid"
    , "evtinfo.subrunid"
    # , "evtinfo.nprotons"
    # , "evtinfo.pbtime"
    # , "evtinfo.pbterr"
    # ---> crvhit
    , "crvhit.sectorType"
    , "crvhit.pos.fCoordinates.fX"
    , "crvhit.pos.fCoordinates.fY"
    , "crvhit.pos.fCoordinates.fZ"
    , "crvhit.timeStart"
    , "crvhit.timeEnd"
    , "crvhit.time"
    , "crvhit.PEs"
    , "crvhit.nHits"
    , "crvhit.nLayers"
    , "crvhit.angle"
    , "@size"
]

def TTreeToDataFrame(finName, treeName, branchNames):
    
    print("---> Reading", treeName, "in", finName)

    # Open the file
    fin = uproot.open(finName)

    # Get the tree
    tree = fin[treeName]

    if len(tree) == 0:
        return

    # Create an empty dictionary to store the selected columns as NumPy arrays
    branchData = {}

    # Iterate over the specified branch names
    for branchName in branchNames:
        # Check if the branch name exists in the TTree
        if branchName in tree:
            # Load values into an array
            branchData[branchName] = tree[branchName].array(library="np")

    # Create the DataFrame directly from the dictionary of column data
    df = pd.DataFrame(branchData)

    # Close the ROOT file
    fin.file.close()

    # Print the DataFrame
    print(df)

    # Return the DataFrame
    return df

# Call the function with the correct branch names list
df = TTreeToDataFrame(
    "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka",
    "TrkAnaExt/trkana",
    branchNames_
)

# Make a map (is that what it's called) of coincidences 
# Plot these guys in space and time 
# find events with coincidences in top and bottom CRV sectors.
# Then estimate the CRV efficiency of the middle module from that.
# Can we get to 99.99% efficiency from that?

print(df)