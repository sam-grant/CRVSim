# Samuel Grant 2024
# Replicate Dave Brown's CompCRV.C in Python 

import sys
import awkward as ak
import uproot

import Utils as ut


def EventStr(event, showCutMasks):
    
    eventStr = (
        f"evtinfo.runid: {event[f'evtinfo.'][f'evtinfo.runid']}\n"
        f"evtinfo.subrunid: {event[f'evtinfo.'][f'evtinfo.subrunid']}\n"
        f"evtinfo.eventid: {event[f'evtinfo.'][f'evtinfo.eventid']}\n"
        # f"crvsummary.nHitCounters: {event[f'crvsummary.'][f'crvsummary.nHitCounters']}\n"
        f"crvcoincs.nLayers {event[f'crvcoincs'][f'crvcoincs.nLayers']}\n"
        f"crvcoincs.angle: {event[f'crvcoincs'][f'crvcoincs.angle']}\n"
        f"crvcoincs.sectorType: {event[f'crvcoincs'][f'crvcoincs.sectorType']}\n"
        f"crvcoincs.pos.fCoordinates: ({event[f'crvcoincs'][f'crvcoincs.pos.fCoordinates.fX']}, {event[f'crvcoincs'][f'crvcoincs.pos.fCoordinates.fY']}, {event[f'crvcoincs'][f'crvcoincs.pos.fCoordinates.fZ']})\n"
        f"crvcoincs.timeStart: {event[f'crvcoincs'][f'crvcoincs.timeStart']}\n"
        f"crvcoincs.timeEnd: {event[f'crvcoincs'][f'crvcoincs.timeEnd']}\n"
        f"crvcoincs.time: {event[f'crvcoincs'][f'crvcoincs.time']}\n"
        f"crvcoincs.PEs: {event[f'crvcoincs'][f'crvcoincs.PEs']}\n"
        f"crvcoincs.nHits: {event[f'crvcoincs'][f'crvcoincs.nHits']}\n"
        # f"crvcoincsmc.valid: {event[f'crvcoincsmc'][f'crvcoincsmc.valid']}\n"
        # f"crvcoincsmc.pdgId: {event[f'crvcoincsmc'][f'crvcoincsmc.pdgId']}\n"
        f"kl.status: {event[f'kl'][f'kl.status']}\n"
        f"kl.nactive: {event[f'kl'][f'kl.nactive']}\n"
        f"kl.nhits: {event[f'kl'][f'kl.nhits']}\n"
        f"kl.nplanes: {event[f'kl'][f'kl.nplanes']}\n"
        f"kl.nnullambig: {event[f'kl'][f'kl.nnullambig']}\n"
        f"kl.ndof: {event[f'kl'][f'kl.ndof']}\n"
        f"kl.fitcon: {event[f'kl'][f'kl.fitcon']}\n"
        f"klfit.sid: {event[f'klfit']['sid']}\n"
        f"klfit.sindex: {event[f'klfit']['sindex']}\n"
        f"klfit.pos.X() {event[f'klfit']['pos']['fCoordinates']['fX']}\n"
        f"klfit.pos.Y() {event[f'klfit']['pos']['fCoordinates']['fY']}\n"
        f"klfit.pos.Z() {event[f'klfit']['pos']['fCoordinates']['fZ']}\n"
        f"klkl.z0err: {event[f'klkl']['z0err']}\n"
        f"klkl.d0err: {event[f'klkl']['d0err']}\n"
        f"klkl.thetaerr: {event[f'klkl']['thetaerr']}\n"
        f"klkl.phi0err: {event[f'klkl']['phi0err']}\n"
    )


    if showCutMasks:

        eventStr += (  
            f"goodTrk: {event[f'goodTrk']}\n"
            f"goodCRV: {event[f'goodCRV']}\n"
            f"noCRV: {event[f'noCRV']}\n"
            f"CRV1: {event[f'CRV1']}\n"
            f"KLCRV1: {event[f'KLCRV1']}\n"
            f"kl.bestFit: {event[f'kl.bestFit']}\n"
            f"klkl.bestFit: {event[f'klkl.bestFit']}\n"
            f"L1Fiducial: {event[f'L1Fiducial']}\n"
        )

    return eventStr

def PrintNEvents(data_, nEvents=5, showCutMasks=False): 
    # Iterate event-by-event
    for i, event in enumerate(data_):
        if event is None: continue
        # event[f'evtinfo.'][f'evtinfo.runid']
        print(EventStr(event, showCutMasks)) 
        if i >= nEvents - 1: return


def MarkCuts(arrays): 

    arrays["goodTrk"] = ak.any(arrays["kl"]["kl.status"], axis=1, keepdims=False) > 0
                        
    arrays["goodCRV"] = ak.any(arrays["crvcoincs"]["crvcoincs.nHits"], axis=1, keepdims=False) > 0
    arrays["noCRV"] = ~arrays["goodCRV"]

    # CRV1: are there CRV hits in sector 1? 
    # arrays["CRV1"] = ak.any(arrays["crvcoincs"]["crvcoincs.sectorType"], axis=-1, keepdims=True) == 1
    arrays["CRV1"] = arrays["crvcoincs"]["crvcoincs.sectorType"] == 1 

    # KLCRV1: are there tracks intersecting the CRV?
    arrays["KLCRV1"] = ( (arrays["klfit"]["sid"] == 200) 
                        & (arrays["klfit"]["sindex"] == 1) )

    # bestFit: best possible track fit
    # needs to be seperated by kl and klkl 
    arrays["kl.bestFit"] = ( (arrays["kl"]["kl.ndof"] >= 10)
                            & (arrays["kl"]["kl.fitcon"] > 0.1)
                            & ((arrays["kl"]["kl.nactive"]/arrays["kl"]["kl.nhits"]) > 0.99)
                            & (arrays["kl"]["kl.nplanes"] >= 4)
                            & ((arrays["kl"]["kl.nnullambig"]/arrays["kl"]["kl.nhits"]) < 0.2) )

    arrays["klkl.bestFit"] = ( (arrays["klkl"]["z0err"] < 1) 
                            & (arrays["klkl"]["d0err"] < 1) 
                            & (arrays["klkl"]["thetaerr"] < 0.004)
                            & (arrays["klkl"]["phi0err"] < 0.001) )

    # Hits within L1 fidicuial area
    arrays["L1Fiducial"] = ( (abs(arrays["klfit"]["pos"]["fCoordinates"]["fX"]) < 2500)
                          & (abs(arrays["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 1500) ) 

    return 

def ApplyCuts(arraysMain, cuts): 

    # Make a deep copy of main array(no shared memory)
    arrays = ak.copy(arraysMain)

    # Does the order matter? I don't think so? 
    if "goodTrk" in cuts:
        mask = arrays["goodTrk"] 
        arrays = arrays[mask] # Event level
    if "goodCRV" in cuts:
        mask = arrays["goodCRV"] 
        arrays = arrays[mask] # Event level
    if "noCRV" in cuts:
        mask = arrays["noCRV"] 
        arrays = arrays[mask]  # Event level
    if "CRV1" in cuts: 
        mask = arrays["CRV1"] #  Global level
        arrays["crvcoincs"] = arrays["crvcoincs"][mask] # This works. 
    if "KLCRV1" in cuts:
        mask = arrays["KLCRV1"] 
        # arrays["klfit"] = ak.mask(arrays["klfit"], mask) # This actually does work. 
        arrays["klfit"] = arrays["klfit"][mask] # Local track level
    if "L1Fiducial" in cuts:
        mask = arrays["L1Fiducial"]
        arrays["klfit"] = arrays["klfit"][mask] # Local track level
    if "bestFit" in cuts: 
        mask1 = arrays["kl.bestFit"] # Be cognizant! This will wipe out the entries in kl if the fit is no good, even if it's a "goodTrk".
        mask2 = arrays["klkl.bestFit"]
        arrays["kl"] = arrays["kl"][mask1] # Global track level
        arrays["klkl"] = arrays["klkl"][mask2] # Local track level
        
    # if "kl.bestFit" in cuts: 
    #     mask = arrays["kl.bestFit"]
    #     arrays["kl"] = arrays["kl"][mask] # Local track level
    # if "klkl.bestFit" in cuts: # Doesn't seem to work?
    #     mask = arrays["klkl.bestFit"]
    #     arrays["klkl"] = arrays["klkl"][mask] # Local track level

    return arrays

def Plot1(arrays):

    # TH2D* tcrvpg = new TH2D("tcrvpg","KKInter TCRV Layer 1 Position, Has CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    # TH2D* tcrvpb = new TH2D("tcrvpb","KKInter TCRV Layer 1 Position,  No CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);

    # ta->Project("tcrvpg","klfit.pos.X():klfit.pos.Z()",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);
    # ta->Project("tcrvpb","klfit.pos.X():klfit.pos.Z()",goodtrk+KLCRV1+bestfit+noCRV);

    tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "CRV1", "KLCRV1", "bestFit"])
    tcrvpb = ApplyCuts(arrays, ["noCRV", "goodTrk", "CRV1", "KLCRV1", "bestFit"])

    # Plot
    ut.Plot2D(x=ak.flatten(tcrvpg["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None), y=ak.flatten(tcrvpg["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
            , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
            , title="Has coincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_hasCoin.png")

    ut.Plot2D(x=ak.flatten(tcrvpb["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None), y=ak.flatten(tcrvpb["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
            , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
            , title="No coincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_noCoin.png")

    # tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "CRV1", "KLCRV1", "bestFit"])

    #######################################

    # This works! 
    # tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "KLCRV1", "CRV1"])

    # This doesn't. 
    # tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "KLCRV1", "CRV1", "bestFit"])
    # We get output like this 

    # evtinfo.runid: 1205
    # evtinfo.subrunid: 0
    # evtinfo.eventid: 31673
    # crvcoincs.nLayers [3]
    # crvcoincs.angle: [0.327]
    # crvcoincs.sectorType: [1]
    # crvcoincs.pos.fCoordinates: ([1.82e+03], [4.76e+03], [12])
    # crvcoincs.timeStart: [9.08e+04]
    # crvcoincs.timeEnd: [9.1e+04]
    # crvcoincs.time: [9.08e+04]
    # crvcoincs.PEs: [162]
    # crvcoincs.nHits: [13]
    # kl.status: []
    # kl.nactive: []
    # kl.nhits: []
    # kl.nplanes: []
    # kl.nnullambig: []
    # kl.ndof: []
    # kl.fitcon: []
    # klfit.sid: [[200]]
    # klfit.sindex: [[1]]
    # klfit.pos.X() [[-1.56e+03]]
    # klfit.pos.Y() [[4.78e+03]]
    # klfit.pos.Z() [[-588]]
    # klkl.z0err: [[]]
    # goodTrk: True
    # goodCRV: True
    # noCRV: False
    # CRV1: [True, False]
    # KLCRV1: [[False, False, False, True, False]]
    # kl.bestFit: [False]
    # klkl.bestFit: [[False, False, False, False, False]]
    # L1Fiducial: [[True, True, True, True, True]]

    # Which doesn't make sense, since 
    # PrintNEvents(arrays[arrays["evtinfo."]["evtinfo.eventid"] == 31673], nEvents=1, showCutMasks=False)
    # evtinfo.runid: 1205
    # evtinfo.subrunid: 0
    # evtinfo.eventid: 31673
    # crvcoincs.nLayers [3, 4]
    # crvcoincs.angle: [0.327, 0.326]
    # crvcoincs.sectorType: [1, 3]
    # crvcoincs.pos.fCoordinates: ([1.82e+03, 1.41e+03], [4.76e+03, 4.92e+03], [12, 17])
    # crvcoincs.timeStart: [9.08e+04, 9.08e+04]
    # crvcoincs.timeEnd: [9.1e+04, 9.09e+04]
    # crvcoincs.time: [9.08e+04, 9.08e+04]
    # crvcoincs.PEs: [162, 368]
    # crvcoincs.nHits: [13, 16]
    # kl.status: [1]
    # kl.nactive: [11]
    # kl.nhits: [17]
    # kl.nplanes: [1]
    # kl.nnullambig: [9]
    # kl.ndof: [23]
    # kl.fitcon: [0.768]
    # klfit.sid: [[4, 4, 200, 200, 200]]
    # klfit.sindex: [[0, 0, 0, 1, 2]]
    # klfit.pos.X() [[390, 835, -1.49e+03, -1.56e+03, -1.63e+03]]
    # klfit.pos.Y() [[755, -162, 4.62e+03, 4.78e+03, 4.92e+03]]
    # klfit.pos.Z() [[-418, -379, -582, -588, -595]]
    # klkl.z0err: [[1.29, 1.29, 1.31, 1.31, 1.31]]

    #######################################

    # # What about this
    # tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "KLCRV1", "CRV1", "kl.bestFit"])
    # PrintNEvents(tcrvpg, nEvents=25, showCutMasks=True)

    # # I don't think so? 
    # # Maybe because kl is actually a global parameter?  
    # evtinfo.runid: 1205
    # evtinfo.subrunid: 0
    # evtinfo.eventid: 891
    # crvcoincs.nLayers [4]
    # crvcoincs.angle: [0.354]
    # crvcoincs.sectorType: [1]
    # crvcoincs.pos.fCoordinates: ([1.61e+03], [4.77e+03], [380])
    # crvcoincs.timeStart: [2.88e+04]
    # crvcoincs.timeEnd: [2.89e+04]
    # crvcoincs.time: [2.88e+04]
    # crvcoincs.PEs: [252]
    # crvcoincs.nHits: [16]
    # kl.status: []
    # kl.nactive: []
    # kl.nhits: []
    # kl.nplanes: []
    # kl.nnullambig: []
    # kl.ndof: []
    # kl.fitcon: []
    # klfit.sid: [[200]]
    # klfit.sindex: [[1]]
    # klfit.pos.X() [[1.41e+03]]
    # klfit.pos.Y() [[4.78e+03]]
    # klfit.pos.Z() [[362]]
    # klkl.z0err: [[0.175, 0.166, 0.166, 0.175, 0.166, 0.166, 0.166]]
    # goodTrk: True
    # goodCRV: True
    # noCRV: False
    # CRV1: [True]
    # KLCRV1: [[False, False, False, False, False, True, False]]
    # kl.bestFit: [False]
    # klkl.bestFit: [[False, False, False, False, False, False, False]]
    # L1Fiducial: [[True, False, True, True, True, True, True]]

    # # What about this
    # tcrvpg = ApplyCuts(arrays, ["goodCRV", "goodTrk", "KLCRV1", "CRV1", "klkl.bestFit"])
    # PrintNEvents(tcrvpg, nEvents=25, showCutMasks=True)

    # evtinfo.runid: 1205
    # evtinfo.subrunid: 0
    # evtinfo.eventid: 31673
    # crvcoincs.nLayers [3]
    # crvcoincs.angle: [0.327]
    # crvcoincs.sectorType: [1]
    # crvcoincs.pos.fCoordinates: ([1.82e+03], [4.76e+03], [12])
    # crvcoincs.timeStart: [9.08e+04]
    # crvcoincs.timeEnd: [9.1e+04]
    # crvcoincs.time: [9.08e+04]
    # crvcoincs.PEs: [162]
    # crvcoincs.nHits: [13]
    # kl.status: [1]
    # kl.nactive: [11]
    # kl.nhits: [17]
    # kl.nplanes: [1]
    # kl.nnullambig: [9]
    # kl.ndof: [23]
    # kl.fitcon: [0.768]
    # klfit.sid: [[200]]
    # klfit.sindex: [[1]]
    # klfit.pos.X() [[-1.56e+03]]
    # klfit.pos.Y() [[4.78e+03]]
    # klfit.pos.Z() [[-588]]
    # klkl.z0err: [[]]
    # klkl.d0err: [[]]
    # klkl.thetaerr: [[]]
    # klkl.phi0err: [[]]
    # goodTrk: True
    # goodCRV: True
    # noCRV: False
    # CRV1: [True, False]
    # KLCRV1: [[False, False, False, True, False]]
    # kl.bestFit: [False]
    # klkl.bestFit: [[False, False, False, False, False]]
    # L1Fiducial: [[True, True, True, True, True]]

    # I think I understand now. 
    # The cuts are working fine. 
    # It's just that the kl.bestFit wipes out the entries in kl, leaving the event intact. 
    # So, if you were to flatten everything you would get some odd results if you didn't know what you were doing.  
    
    #######################################

    return

def Run(finName):


    # TrkAna tutorial:
    # https://github.com/Mu2e/TrkAna/blob/main/tutorial

    # We have event level, global coincidence level, global track level, and local track level
    # Trying to mask the entire array based on, for example, a local track level mask, is tricky. 
    # One approach is to split the array into the different levels/trees and mask the appropirate tree
    arrays = ak.Array([])
    with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
        arrays = tree.arrays(["evtinfo.", "crvcoincs", "kl", "klfit", "klkl"]) 

    # # Convert awkward array to nested list
    # nested_list = ak.to_list(arrays[:10])

    # print(nested_list)

    # # return
    # # Define a function to convert nested list to a string
    # def nested_list_to_string(nested_list):
    #     lines = []
    #     for sublist in nested_list:
    #         line = " ".join(map(str, sublist))
    #         lines.append(line)
    #     return "\n".join(lines)

    # # Convert nested list to a string
    # output_string = nested_list_to_string(nested_list)

    # # Write the string to a text file
    # with open("output.txt", "w") as file:
    #     file.write(output_string)

    # return

    # print(arrays.fields)
    # PrintNEvents(arrays)

    # print(arrays["crvcoincs"])
    # print(arrays["klfit"])

    # return

    # Mark cuts 
    MarkCuts(arrays)

    # PrintNEvents(arrays, showCutMasks=True)
   
    # First plot 
    Plot1(arrays)

    return



def main():

    # Take command-line arguments
    finName = sys.argv[1] if len(sys.argv) > 1 else "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/40/73/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # /pnfs/mu2e/scratch/users/sgrant/workflow/CosmicCRYExtractedTrk.MDC2020z2_best_v1_1/outstage/67605881/00/00000/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"

    print("\n---> Running with inputs:\n")
    print("\tfinName:", finName)

    Run(finName=finName) 

    return

if __name__ == "__main__":
    main()