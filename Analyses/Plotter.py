import numpy as np
import pandas as pd

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

def Run(finName):

    SanityPlots = True
    CRVHitPlots = False

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

    # Sanity plots (plot the parameters being used)
    if SanityPlots: 

        # ut.Plot1D(data=eventID_, nBins=int(np.max(eventID_)/5e3), xmin=0, xmax=np.max(eventID_), xlabel="Event ID", ylabel="Events", fout="Images/Sanity/h1_eventID.png") 
        # # ut.Plot1D(data=runID_, nBins=int(np.max(runID_)), xmin=0, xmax=np.max(runID_), xlabel="Run ID", ylabel="Runs", fout="Images/h1_runID.png") 
        # # ut.Plot1D(data=subRunID_, nBins=int(np.max(subRunID_)), xmin=0, xmax=np.max(subRunID_), xlabel="Subrun ID", ylabel="Subruns", fout="Images/h1_subrunID.png") 
        # ut.Plot1D(data=time_, nBins=int(np.max(time_)/100), xmin=0, xmax=np.max(time_), xlabel="Time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_time.png") 
        # ut.Plot1D(data=timeStart_, nBins=int(np.max(timeStart_)/100), xmin=0, xmax=np.max(timeStart_), xlabel="Start time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_startTime.png") 
        # ut.Plot1D(data=timeEnd_, nBins=int(np.max(timeEnd_)/100), xmin=0, xmax=np.max(timeEnd_), xlabel="End time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_endTime.png") 
        # ut.Plot1D(data=x_, nBins=int((np.max(x_)+abs(np.min(x_)))/100), xmin=np.min(x_), xmax=np.max(x_), xlabel="x [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_x.png") 
        # ut.Plot1D(data=y_, nBins=int((np.max(y_)-np.min(y_))/100), xmin=np.min(y_), xmax=np.max(y_), xlabel="y [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_y.png") 
        # ut.Plot1D(data=z_, nBins=int((np.max(z_)-np.min(z_))/100), xmin=np.min(z_), xmax=np.max(z_), xlabel="z [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_z.png") 
        # ut.Plot1D(data=PEs_, nBins=int(np.max(PEs_)/100), xmin=0, xmax=np.max(PEs_), xlabel="PE", ylabel="CRV hits", fout="Images/Sanity/h1_PEs.png") # I think this is PE integral
        ut.Plot1D(data=nLayers_, nBins=int(np.max(nLayers_)+1), xmin=-0.5, xmax=np.max(nLayers_)+0.5, xlabel="Layers hit", ylabel="CRV hits", fout="Images/Sanity/h1_nLayers.png") 
        # ut.Plot1D(data=nHits_, nBins=int(np.max(nHits_)), xmin=0, xmax=np.max(nHits_), xlabel="Hits", ylabel="CRV hits", fout="Images/Sanity/h1_hits.png")
        # ut.Plot1D(data=angle_*(180/(2*np.pi)), nBins=int(360), xmin=-180, xmax=180, xlabel="Angle [deg]", ylabel="CRV hits", fout="Images/Sanity/h1_angle.png") 

    # CRV hits
    if CRVHitPlots:

        ut.PlotGraph(x=z_, y=y_, xlabel="z [mm]", ylabel="y [mm]", fout="Images/CRVHits/gr_ZY.png") 
        ut.PlotGraph(x=z_, y=x_, xlabel="z [mm]", ylabel="x [mm]", fout="Images/CRVHits/gr_ZX.png")
        ut.PlotGraph(x=x_, y=y_, xlabel="x [mm]", ylabel="y [mm]", fout="Images/CRVHits/gr_XY.png")

    return

    # Get true coincidences

    # tree->Scan("crvhit.pos.fCoordinates.fY:crvhit.pos.fCoordinates.fZ:evtinfo.eventid:crvhit.sectorType");
    # ...
    # *       33 *        0 * 4774.9931 * 12846.164 *      9297 *         1 *
    # *       33 *        1 * 4622.1709 * 13562.849 *      9297 *         2 *
    # *       33 *        2 * 4910.0312 * 12823.271 *      9297 *         3 *

    # Conditions: 
    # 1) evtinfo.eventid must be the same in each hit 
    # 2) 
    
    true_coincidences = df[df.groupby('evtinfo.uniqueid')['evtinfo.uniqueid'].transform('count') >= 3]


    print(true_coincidences)


    # Count the occurrences of each unique combination
    df['count'] = df.groupby(['evtinfo.eventid', 'evtinfo.runid', 'evtinfo.subrunid'])['evtinfo.eventid'].transform('count')

    # Filter rows with a count of 3 or more
    coincidences = df[df['count'] >= 3].copy()

    # Drop the 'count' column if not needed
    coincidences.drop(columns='count', inplace=True)

    print(coincidences)

    return


    # Define true coincidences, group rows of >=3 that have the same eventid, runid, and subrunid
    # df_true_coin = df.groupby(['evtinfo.eventid', 'evtinfo.runid', 'evtinfo.subrunid'])
    # df_true_coin = df_true_coin.filter(lambda group: len(group) >= 3)

    # print(df_true_coin)

    df_true_coin = df.merge(df, on=["evtinfo.eventid", "evtinfo.runid", "evtinfo.subrunid"], suffixes=("", "_duplicate"), how="inner")

    print(df_true_coin)

    return
    # Group by the specified columns and count the occurrences
    grouped = df.groupby(['evtinfo.eventid', 'evtinfo.runid', 'evtinfo.subrunid'])
    print(grouped)
    count = grouped.size().reset_index(name='count')

    # Filter the rows with a count of 3 or more
    filtered_df = df.merge(count[count['count'] >= 3], on=['evtinfo.eventid', 'evtinfo.runid', 'evtinfo.subrunid'])

    print(filtered_df)

    return

def main():

    Run("/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka")

    return

if __name__ == "__main__":
    main()