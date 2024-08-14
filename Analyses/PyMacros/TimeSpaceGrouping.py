'''
Samuel Grant 2024

Find times/space grouping cuts for coincidences in trigger sector. 

'''

# External libraries
import sys
import uproot 
import numpy as np
import awkward as ak

# Internal libraries
import Utils as ut
import CoincidenceConditions as cc

from Mu2eEAF import ReadData as rd 
from Mu2eEAF import Parallelise as pa
    
# ------------------------------------------------
#                     Input 
# ------------------------------------------------ 
evtBranchNames_ = [ 
                # Event info
                "evtinfo.run"
                , "evtinfo.subrun"
                , "evtinfo.event"
]

crvBranchNames_ = [
                # Coincidences 
                "crvcoincs.sectorType" 
                , "crvcoincs.pos.fCoordinates.fX" 
                , "crvcoincs.pos.fCoordinates.fY" 
                , "crvcoincs.pos.fCoordinates.fZ"
                , "crvcoincs.time" 
                , "crvcoincs.timeStart"
                , "crvcoincs.PEs" 
                , "crvcoincs.nHits" 
                , "crvcoincs.nLayers" 
                , "crvcoincs.PEsPerLayer[4]"
                , "crvcoincs.angle" 
                # Coincidences (truth)
                , "crvcoincsmc.valid" 
                , "crvcoincsmc.pdgId" 
]

# allBranchNames_ = { "evt" : evtBranchNames_
#                    ,"crv" : crvBranchNames_
#                   }

allBranchNames_ = { "evt" : evtBranchNames_
                   ,"crv" : crvBranchNames_
                  }

def GetData(file, quiet): # finName):

    data_dict_ = {}

    # Open tree
    # with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
    with file["TrkAnaExt/trkana"] as tree:
        # Seperate event info, coincidnces, tracks, and track fits. 
        # This way we can apply masks independently. 
        for field, branch in allBranchNames_.items():
            data_dict_[field] = tree.arrays(branch)

    # Zip them together 
    return ak.zip(data_dict_) 

# ------------------------------------------------
#                Plotting functions 
# ------------------------------------------------

# Fill histogram objects
# delta T 
# delta X
# delta Y
# delta Z 

# Fill histogram object.
# This function is only possible with arrays with one coincidence in sectors 2 & 3 in every event.
# Filtering beforehand is essential.
def FillHistogram(data_, recon, foutTag):

    # Filter muons
    data_muons_ = ak.copy(data_)
    data_muons_ = FilterParticles(data_muons_, "muons", quiet=False)

    # Deep copies
    sector2_ = ak.copy(data_)
    sector3_ = ak.copy(data_)
    sector2_muons_ = ak.copy(data_muons_)
    sector3_muons_ = ak.copy(data_muons_)
    
    
    # Sector 2 & 3 arrays
    sector2_["crv"] = sector2_["crv"][sector2_["crv"]["crvcoincs.sectorType"] == 2]
    sector3_["crv"] = sector3_["crv"][sector3_["crv"]["crvcoincs.sectorType"] == 3]
    sector2_muons_["crv"] = sector2_muons_["crv"][sector2_muons_["crv"]["crvcoincs.sectorType"] == 2]
    sector3_muons_["crv"] = sector3_muons_["crv"][sector3_muons_["crv"]["crvcoincs.sectorType"] == 3]

    print("\nExample sector 2 event:")
    PrintNEvents(sector2_, 1)
    print("Example sector 3 event:")
    PrintNEvents(sector3_, 1)

    # # Filter muons
    # sector2_muons_ = FilterParticles(sector2_, "muons", quiet=False)
    # sector3_muons_ = FilterParticles(sector3_, "muons", quiet=False)

    # Time/space difference 
    # 3 is top, 2 is bottom
    dT_ = ak.flatten(sector2_["crv"]["crvcoincs.time"]) - ak.flatten(sector3_["crv"]["crvcoincs.time"])
    dX_ = ak.flatten(sector2_["crv"]["crvcoincs.pos.fCoordinates.fX"]) - ak.flatten(sector3_["crv"]["crvcoincs.pos.fCoordinates.fX"])
    dY_ = ak.flatten(sector2_["crv"]["crvcoincs.pos.fCoordinates.fY"]) - ak.flatten(sector3_["crv"]["crvcoincs.pos.fCoordinates.fY"])
    dY_ = ak.flatten(sector2_["crv"]["crvcoincs.pos.fCoordinates.fY"]) - ak.flatten(sector3_["crv"]["crvcoincs.pos.fCoordinates.fY"])
    dPE_ = ak.flatten(sector2_["crv"]["crvcoincs.PEs"]) - ak.flatten(sector3_["crv"]["crvcoincs.PEs"]) 
    dSlope_ = ak.flatten(sector2_["crv"]["crvcoincs.angle"]) - ak.flatten(sector3_["crv"]["crvcoincs.angle"]) 

    dT_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.time"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.time"])
    dX_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.pos.fCoordinates.fX"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.pos.fCoordinates.fX"])
    dY_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.pos.fCoordinates.fY"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.pos.fCoordinates.fY"])
    dY_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.pos.fCoordinates.fY"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.pos.fCoordinates.fY"])
    dPE_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.PEs"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.PEs"]) 
    dSlope_muons_ = ak.flatten(sector2_muons_["crv"]["crvcoincs.angle"]) - ak.flatten(sector3_muons_["crv"]["crvcoincs.angle"]) 
    
    # No resolution in z on the bottom layer. 
    # dZ_ = ak.flatten(sector2_["crv"]["crvcoincs.pos.fCoordinates.fZ"]) - ak.flatten(sector3_["crv"]["crvcoincs.pos.fCoordinates.fZ"])

    print("dT_:", dT_)
    print("dX_:", dX_)
    print("dY_:", dY_)
    print("dY_:", dY_)
    print("dPE_:", dPE_)
    print("dSlope_:", dSlope_)
    print()

    # Plot(dT_, dX_, dY_, dPE_, dSlope_, recon, foutTag)

    # ut.Plot1D(dT_, nbins=100, xmin=-15, xmax=25, xlabel=r"$\Delta t$ [ns]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dT_{foutTag}.png", underOver=True)
    # ut.Plot1D(dX_, nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta X$ [mm]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dX_{foutTag}.png", underOver=True)
    # ut.Plot1D(dY_, nbins=100, xmin=-400, xmax=-200, xlabel=r"$\Delta Y$ [mm]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dY_{foutTag}.png", underOver=True)
    # ut.Plot1D(dPE_, nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta$ PEs", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dPE_{foutTag}.png", underOver=True)
    # ut.Plot1D(dSlope_, nbins=100, xmin=-2.5, xmax=2.5, xlabel=r"$\Delta$ slope", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dSlope_{foutTag}.png", underOver=True)

    # def Plot1DOverlay(hists_, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", label_=None, legPos="best", NDPI=300, includeBlack=False, logY=False, legFontSize=12):
    
    label_ = ["All", "Muons"]
    
    ut.Plot1DOverlay([dT_, dT_muons_], nbins=100, xmin=-15, xmax=25, xlabel=r"$\Delta t$ [ns]", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dT_{foutTag}.png")
    ut.Plot1DOverlay([dX_, dX_muons_], nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta X$ [mm]", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dX_{foutTag}.png")
    ut.Plot1DOverlay([dY_, dY_muons_], nbins=100, xmin=-400, xmax=-200, xlabel=r"$\Delta Y$ [mm]", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dY_{foutTag}.png")
    ut.Plot1DOverlay([dPE_, dPE_muons_], nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta$ PEs", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dPE_{foutTag}.png")
    ut.Plot1DOverlay([dSlope_, dSlope_muons_], nbins=100, xmin=-2.5, xmax=2.5, xlabel=r"$\Delta$ slope", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dSlope_{foutTag}.png")

    # Hail mary?
    # Normalise the arrays by there mean, multiply them together? 
    dT_ = dT_ / np.mean(dT_)
    dX_ = dX_ / np.mean(dX_)
    dY_ = dY_ / np.mean(dY_)
    dPE_ = dPE_ / np.mean(dPE_)
    dSlope_ = dSlope_ / np.mean(dSlope_) 

    dTot_ = dT_ * dX_ * dY_ * dPE_ * dSlope_
    dTot_ = dTot_ / np.mean(dTot_)

    dT_muons_ = dT_muons_ / np.mean(dT_muons_)
    dX_muons_ = dX_muons_ / np.mean(dX_muons_)
    dY_muons_ = dY_muons_ / np.mean(dY_muons_)
    dPE_muons_ = dPE_muons_ / np.mean(dPE_muons_)
    dSlope_muons_ = dSlope_muons_ / np.mean(dSlope_muons_) 

    dTot_muons_ = dT_muons_ * dX_muons_ * dY_muons_ * dPE_muons_ * dSlope_muons_
    dTot_muons_ = dTot_muons_ / np.mean(dTot_muons_)

    ut.Plot1DOverlay([dTot_, dTot_muons_], nbins=100, xmin=-100, xmax=100, xlabel=r"$\Delta$ total [normalised]", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dTot_{foutTag}.png")

    ut.Plot1DOverlay([dTot_, dTot_muons_], nbins=100, xmin=-100, xmax=100, xlabel=r"$\Delta$ total [normalised]", ylabel="Coincidences", label_=label_, fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_overlay_dTot_log_{foutTag}.png", logY=True)

    



    return
    
# def Plot(dT_, dX_, dY_, dPE_, dSlope_, recon, foutTag):

#     # def Plot1D(data, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", legPos="best", stats=True, underOver=False, errors=False, NDPI=300):
    
#     ut.Plot1D(dT_, nbins=100, xmin=-15, xmax=25, xlabel=r"$\Delta T$ [ns]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dT_{foutTag}.png", underOver=True)
#     ut.Plot1D(dX_, nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta X$ [mm]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dX_{foutTag}.png", underOver=True)
#     ut.Plot1D(dY_, nbins=100, xmin=-400, xmax=-200, xlabel=r"$\Delta Y$ [mm]", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dY_{foutTag}.png", underOver=True)
#     ut.Plot1D(dPE_, nbins=100, xmin=-1000, xmax=1000, xlabel=r"$\Delta$ PEs", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dPE_{foutTag}.png", underOver=True)
#     ut.Plot1D(dSlope_, nbins=100, xmin=-2.5, xmax=2.5, xlabel=r"$\Delta$ slope", ylabel="Coincidences", fout=f"../Images/{recon}/TimeAndSpaceGrouping/h1_dSlope_{foutTag}.png", underOver=True)

#     return

# Printout helper functions
def PrintEvent(event):
    
    eventStr = (
        f"-------------------------------------------------------------------------------------\n"
        f"***** evt *****\n"
        f"evtinfo.run: {event['evt']['evtinfo.run']}\n" 
        f"evtinfo.subrun: {event['evt']['evtinfo.subrun']}\n" 
        f"evtinfo.eventid: {event['evt']['evtinfo.event']}\n"
        f"***** crv *****\n"
        f"crvcoincs.sectorType: {event['crv']['crvcoincs.sectorType']}\n"
        f"crvcoincs.nLayers {event['crv']['crvcoincs.nLayers']}\n"
        f"crvcoincs.angle: {event['crv']['crvcoincs.angle']}\n"
        f"crvcoincs.pos.fCoordinates: ({event['crv']['crvcoincs.pos.fCoordinates.fX']}, {event['crv']['crvcoincs.pos.fCoordinates.fY']}, {event['crv']['crvcoincs.pos.fCoordinates.fZ']})\n"
        f"crvcoincs.timeStart: {event['crv']['crvcoincs.timeStart']}\n"
        f"crvcoincs.time: {event['crv']['crvcoincs.time']}\n"
        f"crvcoincs.PEs: {event['crv']['crvcoincs.PEs']}\n"
        f"crvcoincs.PEsPerLayer[4]: {event['crv']['crvcoincs.PEsPerLayer[4]']}\n"
        f"crvcoincs.nHits: {event['crv']['crvcoincs.nHits']}\n"
        f"crvcoincsmc.pdgId: {event['crv']['crvcoincsmc.pdgId']}\n"
        f"crvcoincsmc.valid: {event['crv']['crvcoincsmc.valid']}\n"
        f"-------------------------------------------------------------------------------------\n"
    )

    return eventStr

# Print N events in verbose format. 
def PrintNEvents(data_, nEvents=10, showMasks=False):
     # Iterate event-by-event
    for i, event in enumerate(data_, start=1):
        print(PrintEvent(event))
        if i >= nEvents: 
            return


# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

def FilterParticles(data_, particle, quiet):
    if not quiet: print(f"\n---> Filtering particles, keeping {particle}")
    
    muonCondition = ak.any((data_["crv"]["crvcoincsmc.pdgId"] == 13) | (data_["crv"]["crvcoincsmc.pdgId"] == -13), axis=1)
    if particle == "all":
        return data_
    elif particle == "muons": 
        return data_[muonCondition] 
    elif particle == "non_muons":
        return data_[~muonCondition] 
    else:
        raise ValueError(f"Particle string {particle} not valid!")
        
#  Events with ONE coincidence in sectors 2 & 3
def PassZero(data_, fail, quiet):
    if not quiet: print(f"\n---> Running pass zero filter") 

    sector2Condition = data_["crv"]["crvcoincs.sectorType"] == 2
    sector3Condition = data_["crv"]["crvcoincs.sectorType"] == 3

    oneCoincInSector2Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector2Condition], axis=1) == 1
    oneCoincInSector3Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector3Condition], axis=1) == 1
    
    data_["oneCoinInTriggerSectors"] = (oneCoincInSector2Condition & oneCoincInSector3Condition)

    if not quiet: print("Done!")
    
    # Cut on event level
    if not fail: 
        return data_[data_["oneCoinInTriggerSectors"]]
    else: 
        return data_[~data_["oneCoinInTriggerSectors"]]


# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_, quiet): 
    if not quiet: print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_["crv"]["crvcoincs.timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(file, recon, doutTag, foutTag, quiet): 

    # Placeholder
    fail = False

    # Get data as a set of awkward arrays
    data_ = GetData(file, quiet) # finName)

    # Apply start time cut
    data_ = CutOnStartTime(data_, quiet)

    # Filter particle species
    # data_all_ = FilterParticles(data_, "all", quiet)
    # data_muons_ = FilterParticles(data_, "muons", quiet)

    # Run pass zero
    data_ = PassZero(data_, fail, quiet)
    
    # Make plots
    # Plot(data_, recon, foutTag, quiet)
    FillHistogram(data_, recon, foutTag)
    # PrintNEvents(data_)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    # Testing
    TestMain()
    return
    
    ##################################
    # Input parameters
    ##################################

    defname = "nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.root"
    recon = "MDC2020ae"
    # coincidenceConditions = "10PEs2Layers" # should be in a loop
    coincidenceFilterLevel = "pass0" 
    sanityPlots = False
    quiet = True
    particles_ = ["all", "muon", "non_muon"]
    layers_ = [2, 3]
    PEs_ = np.arange(10.0, 135.0, 5) # Same steps as Tyler
    
    # Not sure about this 
    coincidenceFilters = {
        "pass0" : "one_coincidence_per_trigger_sector"
        , "pass1" : "coincidence_grouping"
    }

    # Get file list
    fileList = rd.GetFileList(defname) 

    trials = len(fileList) + len(particles_) + len(layers_) + len(PEs_)
    
    ##################################
    # Setup
    ##################################
    
    # Wrap Run() in a process function
    def processFunction(fileName):
        file = rd.ReadFile(fileName, quiet)
        doutTag = fileName.split('.')[-2] 
        count = 0 
        # Scan particles
        for particle in particles_:
            # Scan layers
            for layer in layers_:
                # Scan PE thresholds
                for PE in PEs_: 
                    coincidenceConditions = f"{PE}PEs{layer}Layers"
                    foutTag = particle + "_" + coincidenceConditions + "_" + coincidenceFilters[coincidenceFilterLevel]
                    percent_done = (count / trials) * 100
                    
                    outputStr = (
                        "\n---> Running with:\n"
                        f"fileName: {fileName}\n"
                        f"recon: {recon}\n"
                        f"particle: {particle}\n"
                        f"layers: {layer}/4\n"
                        f"PEs: {PE}\n"
                        f"coincidenceConditions: {coincidenceConditions}\n"
                        f"coincidenceFilterLevel: {coincidenceFilterLevel}\n"
                        f"doutTag: {doutTag}\n"
                        f"foutTag: {foutTag}\n"
                        f"sanityPlots: {sanityPlots}\n"
                        f"Percent done: {percent_done:.2f}%\n" # Doesn't really work with multithreading. Worth fixing?
                    )
                    
                    print(outputStr, flush=True)
                
                    # Run script
                    Run(file, recon, particle, coincidenceConditions, doutTag, foutTag, sanityPlots, quiet) 
      
                    count += 1
                    # print() 
                    # Uncomment for testing 
                    # return   
        return 

    ##################################
    # Submit 
    ##################################    

    # Testing 
    # processFunction(fileList[0])
    # processFunction(fileList[1])

    # Submit jobs
    pa.Multithread(fileList, processFunction) 

    print("---> Analysis complete!", flush=True)
    
    return

# TestParallelise()

def TestMain():

    fileName = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00023/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root"
    file = uproot.open(fileName)
    coincidenceConditions = "10PEs2Layers"
    recon = "MDC2020ae"
    # particle = "muons"
    # particles_ = ["all", "muon", "non_muon"]
    # foutTag = particle + "_" + coincidenceConditions + "_one_coincidence_per_trigger_sector"
    foutTag = coincidenceConditions + "_one_coincidence_per_trigger_sector"
    doutTag = fileName.split('.')[-2]
    quiet = False
    
    outputStr = (
        "\n---> Running with:\n"
        f"fileName: {fileName}\n"
        f"recon: {recon}\n"
        # f"particle: {particle}\n"
        f"doutTag: {doutTag}\n"
        f"foutTag: {foutTag}\n"
        f"quiet: {quiet}\n"
    )
    print(outputStr)

    Run(file, recon, doutTag, foutTag, quiet) 
    # [Run(file, recon, particle, doutTag, foutTag, quiet) for particle in particles_]

    return

if __name__ == "__main__":
    main()
