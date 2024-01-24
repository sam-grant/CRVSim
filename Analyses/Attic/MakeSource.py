import numpy as np
import pandas as pd
import h5py

import Utils as ut

# Function to flatten arrays in the specified column
def Flatten(arr):
    return [item for sublist in arr for item in sublist]

# Get data, flatten and convert DataFrame column to numpy array
def GetPar(df, parName, flatten=False):
    data = df[parName] 
    if flatten: data = Flatten(data)
    data = np.array(data)    
    return data

def Run(finName, foutName):

    # Get data
    df = ut.TTreeToDataFrame(finName, "TrkAnaExt/trkana", ut.branchNames_)

    # Load and format parameters
    eventID_ = GetPar(df, "evtinfo.eventid") 
    runID_ = GetPar(df, "evtinfo.runid") # just one (per file?)
    subRunID_ = GetPar(df, "evtinfo.subrunid") # seems somewhat arbitrary
    time_ = GetPar(df, "crvhit.time", flatten=True)
    timeStart_ = GetPar(df, "crvhit.timeStart", flatten=True)
    timeEnd_ = GetPar(df, "crvhit.timeEnd", flatten=True)
    x_ = GetPar(df, "crvhit.pos.fCoordinates.fX", flatten=True)
    y_ = GetPar(df, "crvhit.pos.fCoordinates.fY", flatten=True)
    z_ = GetPar(df, "crvhit.pos.fCoordinates.fZ", flatten=True)
    PEs_ = GetPar(df, "crvhit.PEs", flatten=True)
    nLayers_ = GetPar(df, "crvhit.nLayers", flatten=True)
    nHits_ = GetPar(df, "crvhit.nHits", flatten=True)
    angle_ = GetPar(df, "crvhit.angle", flatten=True)

    # Write to h5 
    with h5py.File(foutName, 'w') as h5file:
        h5file.create_dataset('eventID', data=eventID_)
        h5file.create_dataset('runID', data=runID_)
        h5file.create_dataset('subRunID', data=subRunID_)
        h5file.create_dataset('time', data=time_)
        h5file.create_dataset('timeStart', data=timeStart_)
        h5file.create_dataset('timeEnd', data=timeEnd_)
        h5file.create_dataset('x', data=x_)
        h5file.create_dataset('y', data=y_)
        h5file.create_dataset('z', data=z_)
        h5file.create_dataset('PEs', data=PEs_)
        h5file.create_dataset('nLayers', data=nLayers_)
        h5file.create_dataset('nHits', data=nHits_)
        h5file.create_dataset('angle', data=angle_)

    print("---> Written ", foutName)

    return

def main():

    Run("/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka", "../h5/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka.h5")

    return

if __name__ == "__main__":
    main()