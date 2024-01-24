import numpy as np
import pandas as pd
import awkward as ak

import Utils as ut

# Sanity check
def PrintAwkwardInfo(arr_):

    print("\n*** Awkward arrays info ***")

    if len(arr_) > 0:

        print(f"Number of arrays: {len(arr_)}")

        # Access the keys of the first dictionary to get branch names
        branchNames_ = ak.fields(arr_[0])

        # Print the branch names
        for branchName in branchNames_: 
            print(f"Branch name: {branchName}")
            print(f"Type: {ak.type(arr_[0][branchName])}")

    else:
        print("No events in the Awkward Array.")

    print("*** ***\n")

    return  

def GetPosition3D(arr_, branch="crvhit"):

    pos_ = ak.zip({"px": ak.flatten(arr_['%s.pos.fCoordinates.fX'%branch]), 
                    "py": ak.flatten(arr_['%s.pos.fCoordinates.fY'%branch]), 
                    "pz": ak.flatten(arr_['%s.pos.fCoordinates.fZ'%branch]),}, with_name="Position3D")

    return pos_

def Run(finName, treeName):

    # Get data  
    arr_ = ut.TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", ut.branchNames_)

    

    # PrintAwkwardInfo(arr_)

    print(arr_["crvhit.time"])
    print(ak.flatten(arr_["crvhit.time"]))

    return
    

    # Get time, position, momentum 
    eventID_ = arr_["evtinfo.eventid"]
    time_ = ak.flatten(arr_["crvhit.time"])
    pos_ = GetPosition3D(arr_)

    print(len(eventID_), len(time_))

    print(eventID_[:10], time_[:10])

    return

    # pos_ = GetPosition3D(arr_) 

    # eventID_ = GetPar(df, "evtinfo.eventid") 
    # runID_ = GetPar(df, "evtinfo.runid") # just one (per file?)
    # subRunID_ = GetPar(df, "evtinfo.subrunid") # seems somewhat arbitrary
    # time_ = GetPar(df, "crvhit.time", flatten=True)
    # timeStart_ = GetPar(df, "crvhit.timeStart", flatten=True)
    # timeEnd_ = GetPar(df, "crvhit.timeEnd", flatten=True)
    # x_ = GetPar(df, "crvhit.pos.fCoordinates.fX", flatten=True)
    # y_ = GetPar(df, "crvhit.pos.fCoordinates.fY", flatten=True)
    # z_ = GetPar(df, "crvhit.pos.fCoordinates.fZ", flatten=True)
    # PEs_ = GetPar(df, "crvhit.PEs", flatten=True)
    # nLayers_ = GetPar(df, "crvhit.nLayers", flatten=True)
    # nHits_ = GetPar(df, "crvhit.nHits", flatten=True)
    # angle_ = GetPar(df, "crvhit.angle", flatten=True)







    # SanityPlots = True
    # CRVHitPlots = False


    # # Read data from the HDF5 file
    # df = pd.read_hdf(finName)

    # return

    # # Load and format parameters
    # eventID_ = GetPar(df, "evtinfo.eventid") 
    # runID_ = GetPar(df, "evtinfo.runid")
    # subRunID_ = GetPar(df, "evtinfo.subrunid")
    # time_ = GetPar(df, "crvhit.time", flatten=True)
    # timeStart_ = GetPar(df, "crvhit.timeStart", flatten=True)
    # timeEnd_ = GetPar(df, "crvhit.timeEnd", flatten=True)
    # x_ = GetPar(df, "crvhit.pos.fCoordinates.fX", flatten=True)
    # y_ = GetPar(df, "crvhit.pos.fCoordinates.fY", flatten=True)
    # z_ = GetPar(df, "crvhit.pos.fCoordinates.fZ", flatten=True)
    # PEs_ = GetPar(df, "crvhit.PEs", flatten=True)
    # nLayers_ = GetPar(df, "crvhit.nLayers", flatten=True)
    # nHits_ = GetPar(df, "crvhit.nHits", flatten=True)
    # angle_ = GetPar(df, "crvhit.angle", flatten=True)


    # # Sanity plots (plot the parameters being used)
    # if SanityPlots: 

    #     # ut.Plot1D(data=eventID_, nBins=int(np.max(eventID_)/5e3), xmin=0, xmax=np.max(eventID_), xlabel="Event ID", ylabel="Events", fout="Images/Sanity/h1_eventID.png") 
    #     # # ut.Plot1D(data=runID_, nBins=int(np.max(runID_)), xmin=0, xmax=np.max(runID_), xlabel="Run ID", ylabel="Runs", fout="Images/h1_runID.png") 
    #     # # ut.Plot1D(data=subRunID_, nBins=int(np.max(subRunID_)), xmin=0, xmax=np.max(subRunID_), xlabel="Subrun ID", ylabel="Subruns", fout="Images/h1_subrunID.png") 
    #     # ut.Plot1D(data=time_, nBins=int(np.max(time_)/100), xmin=0, xmax=np.max(time_), xlabel="Time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_time.png") 
    #     # ut.Plot1D(data=timeStart_, nBins=int(np.max(timeStart_)/100), xmin=0, xmax=np.max(timeStart_), xlabel="Start time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_startTime.png") 
    #     # ut.Plot1D(data=timeEnd_, nBins=int(np.max(timeEnd_)/100), xmin=0, xmax=np.max(timeEnd_), xlabel="End time [ns]", ylabel="CRV hits", fout="Images/Sanity/h1_endTime.png") 
    #     # ut.Plot1D(data=x_, nBins=int((np.max(x_)+abs(np.min(x_)))/100), xmin=np.min(x_), xmax=np.max(x_), xlabel="x [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_x.png") 
    #     # ut.Plot1D(data=y_, nBins=int((np.max(y_)-np.min(y_))/100), xmin=np.min(y_), xmax=np.max(y_), xlabel="y [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_y.png") 
    #     # ut.Plot1D(data=z_, nBins=int((np.max(z_)-np.min(z_))/100), xmin=np.min(z_), xmax=np.max(z_), xlabel="z [mm]", ylabel="CRV hits", fout="Images/Sanity/h1_z.png") 
    #     # ut.Plot1D(data=PEs_, nBins=int(np.max(PEs_)/100), xmin=0, xmax=np.max(PEs_), xlabel="PE", ylabel="CRV hits", fout="Images/Sanity/h1_PEs.png") # I think this is PE integral
    #     ut.Plot1D(data=nLayers_, nBins=int(np.max(nLayers_)+1), xmin=-0.5, xmax=np.max(nLayers_)+0.5, xlabel="Layers hit", ylabel="CRV hits", fout="Images/Sanity/h1_nLayers.png") 
    #     # ut.Plot1D(data=nHits_, nBins=int(np.max(nHits_)), xmin=0, xmax=np.max(nHits_), xlabel="Hits", ylabel="CRV hits", fout="Images/Sanity/h1_hits.png")
    #     # ut.Plot1D(data=angle_*(180/(2*np.pi)), nBins=int(360), xmin=-180, xmax=180, xlabel="Angle [deg]", ylabel="CRV hits", fout="Images/Sanity/h1_angle.png") 

    # # CRV hits
    # if CRVHitPlots:

    #     ut.PlotGraph(x=z_, y=y_, xlabel="z [mm]", ylabel="y [mm]", fout="Images/CRVHits/gr_ZY.png") 
    #     ut.PlotGraph(x=z_, y=x_, xlabel="z [mm]", ylabel="x [mm]", fout="Images/CRVHits/gr_ZX.png")
    #     ut.PlotGraph(x=x_, y=y_, xlabel="x [mm]", ylabel="y [mm]", fout="Images/CRVHits/gr_XY.png")

    # return

    # return

def main():
    # Specify the name of the HDF5 file to read from
    # h5_file_path = "output.h5"
    
    # Run the analysis
    Run("/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka", "TrkAnaExt/trkana") # "../h5/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka.h5")

    return

if __name__ == "__main__":
    main()
