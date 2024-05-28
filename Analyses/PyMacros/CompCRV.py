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
        f"crvcoincsmc.pdgId: {event[f'crvcoincsmc'][f'crvcoincsmc.pdgId']}\n"
        f"kl.status: {event[f'kl'][f'kl.status']}\n"
        f"kl.nactive: {event[f'kl'][f'kl.nactive']}\n"
        f"kl.nhits: {event[f'kl'][f'kl.nhits']}\n"
        f"kl.nplanes: {event[f'kl'][f'kl.nplanes']}\n"
        f"kl.nnullambig: {event[f'kl'][f'kl.nnullambig']}\n"
        f"kl.ndof: {event[f'kl'][f'kl.ndof']}\n"
        f"kl.fitcon: {event[f'kl'][f'kl.fitcon']}\n"
        f"klfit.sid: {event[f'klfit']['sid']}\n"
        f"klfit.sindex: {event[f'klfit']['sindex']}\n"
        f"klfit.time: {event[f'klfit']['time']}\n"
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
            f"singleCRV1: {event[f'singleCRV1']}\n"
            f"muonsInCRV: {event[f'muonsInCRV']}\n"
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


def MarkCuts(arrays, cuts=[]): 

    if len(cuts) == 0: 
        cuts = ["goodCRV", "noCRV", "goodTrk", "CRV1", "singleCRV1", "muonsInCRV", "KLCRV1", "bestFit", "L1Fiducial"]
    
    if "goodCRV" or "noCRV" in cuts: # Hits in the CRV
        arrays["goodCRV"] = ak.any(arrays["crvcoincs"]["crvcoincs.nHits"], axis=1, keepdims=False) > 0
        # arrays["goodCRV"] = arrays["crvcoincs"]["crvcoincs.nHits"] > 0
        arrays["noCRV"] = ~arrays["goodCRV"]

    if "goodTrk" in cuts: # Tracks in the tracker.
        # arrays["goodTrk"] = arrays["kl"]["kl.status"] > 0
        arrays["goodTrk"] = ak.any(arrays["kl"]["kl.status"], axis=1, keepdims=False) > 0

    if "CRV1" in cuts: # Coincidences in CRV sector 1
        # arrays["CRV1"] = ak.any(arrays["crvcoincs"]["crvcoincs.sectorType"], axis=-1, keepdims=True) == 1
        arrays["CRV1"] = arrays["crvcoincs"]["crvcoincs.sectorType"] == 1 

    if "singleCRV1" in cuts: # Single coincidences in CRV sector 1
        # arrays["CRV1"] = ak.any(arrays["crvcoincs"]["crvcoincs.sectorType"], axis=1, keepdims=True) == 1
        arrays["singleCRV1"] = ak.sum(arrays["crvcoincs"]["crvcoincs.sectorType"] == 1, axis=1, keepdims=False) == 1

    # if "singleCRV1" in cuts: # Single coincidences in CRV sector 1
    #     # arrays["CRV1"] = ak.any(arrays["crvcoincs"]["crvcoincs.sectorType"], axis=1, keepdims=True) == 1
    #     arrays["singleCRV1"] = ak.sum(arrays["crvcoincs"]["crvcoincs.sectorType"] == 1, axis=1, keepdims=False) == 1

    if "muonsInCRV" in cuts: # Coincidences in CRV sector 1
        # arrays["CRV1"] = ak.any(arrays["crvcoincs"]["crvcoincs.sectorType"], axis=1, keepdims=True) == 1
        arrays["muonsInCRV"] = ( (arrays["crvcoincsmc"]["crvcoincsmc.pdgId"] == 13) 
                                    | (arrays["crvcoincsmc"]["crvcoincsmc.pdgId"] == -13) )

    if "KLCRV1" in cuts: # Tracks intersecting the CRV
        arrays["KLCRV1"] = ( (arrays["klfit"]["sid"] == 200) 
                            & (arrays["klfit"]["sindex"] == 1) )

        # For some reason this just doesn't work properly 
        # I think klfit resets the array axis indexing, so it's axis=1 not 2 
        # arrays["KLCRV1"] = ( (ak.any(arrays["klfit"]["sid"], axis=1, keepdims=True) == 200) 
        #                     & (ak.any(arrays["klfit"]["sindex"], axis=1, keepdims=True) == 1) )

    if "bestFit" in cuts: # Good quality track fit 
        # needs to be seperated by kl and klkl in this scheme
        arrays["kl.bestFit"] = ( (arrays["kl"]["kl.ndof"] >= 10)
                                & (arrays["kl"]["kl.fitcon"] > 0.1)
                                & ((arrays["kl"]["kl.nactive"]/arrays["kl"]["kl.nhits"]) > 0.99)
                                & (arrays["kl"]["kl.nplanes"] >= 4)
                                & ((arrays["kl"]["kl.nnullambig"]/arrays["kl"]["kl.nhits"]) < 0.2) )
        arrays["klkl.bestFit"] = ( (arrays["klkl"]["z0err"] < 1) 
                                & (arrays["klkl"]["d0err"] < 1) 
                                & (arrays["klkl"]["thetaerr"] < 0.004)
                                & (arrays["klkl"]["phi0err"] < 0.001) )

    if "L1Fiducial" in cuts: # Hits within L1 fidicuial area
        arrays["L1Fiducial"] = ( (abs(arrays["klfit"]["pos"]["fCoordinates"]["fX"]) < 2500)
                            & (abs(arrays["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 1500) ) 

    return 

def ApplyCuts(arraysMain, cuts): 

    # Make a deep copy of main array(no shared memory)
    arrays = ak.copy(arraysMain)

    if "goodTrk" in cuts:
        mask = arrays["goodTrk"] 
        # mask = ak.any(arrays["goodTrk"], axis=-1, keepdims=True)
        arrays = arrays[mask] # Event level
        # arrays = ak.mask(arrays, mask)
        # arrays["kl"] = arrays["kl"][mask] # Event level
    if "goodCRV" in cuts:
        mask = arrays["goodCRV"] 
        # mask = ak.any(arrays["goodCRV"], axis=-1, keepdims=True)
        arrays = arrays[mask] # Event level
        # arrays = ak.mask(arrays, mask)
        # arrays["crvcoincs"] = arrays["crvcoincs"][mask] # Event level
    if "noCRV" in cuts:
        mask = arrays["noCRV"] 
        # mask = ak.any(arrays["noCRV"], axis=-1, keepdims=True)
        arrays = arrays[mask]  # Event level
        # arrays = ak.mask(arrays, mask)
    if "CRV1" in cuts: 
        mask = arrays["CRV1"] #  Global level
        arrays["crvcoincs"] = arrays["crvcoincs"][mask]
        # Fancier masks lead to some really weird behaviour
        # arrays["crvcoincs"] = ak.mask(arrays["crvcoincs"], ak.any(mask, keepdims=True))
        # arrays["crvcoincs"] = ak.mask(arrays["crvcoincs"], mask) # , keepdims=True)
        # mask = ak.any(arrays["CRV1"], axis=1, keepdims=True)
    if "singleCRV1" in cuts:
        mask = arrays["singleCRV1"] #  Global level
        arrays = arrays[mask] 
    if "KLCRV1" in cuts:
        mask = arrays["KLCRV1"] 
        arrays["klfit"] = arrays["klfit"][mask] # Local track level
    if "L1Fiducial" in cuts:
        mask = arrays["L1Fiducial"]
        # arrays["klfit"] = arrays["klfit"][mask] # Local track level
        arrays["klfit"] = ak.mask(arrays["klfit"], mask)
    if "bestFit" in cuts: 
        mask1 = arrays["kl.bestFit"] 
        mask2 = arrays["klkl.bestFit"]
        arrays["kl"] = arrays["kl"][mask1] # Global track level
        arrays["klfit"] = arrays["klfit"][mask1] # Local track level, is this implicit in Dave's code? 
        arrays["klkl"] = arrays["klkl"][mask2] # Local track level

    # Cuts to klfit will wipe out the entries in kl if the fit is no good, even if it's a "goodTrk".
    # Need to reapply the goodTrk markers 
    if ("bestFit" or "goodTrk") in cuts:
        # This is a possbile fix: reapply the goodTrk markers and cut on goodTrk 
        MarkCuts(arrays, ["goodTrk"])
        arrays = ApplyCuts(arrays, ["goodTrk"])

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

    return

def Plot2(arrays):

    # TH1D* tcrvts = new TH1D("tcrvts","KKInter Time - CRV Time, Layer 1 #Delta T;T_{KKInter} - T_{CRV} (ns)",250,-30,30);
    # ta->Project("tcrvts","klfit.time-crvcoincs.time",goodtrk+KLCRV1+CRV1+bestfit+goodCRV);

    tcrvts = ApplyCuts(arrays, ["goodCRV", "CRV1", "goodTrk", "KLCRV1", "bestFit", "singleCRV1"]) # , "singleKL"]) # , "singleTrack"]) 


    # y_klfit = ak.flatten(tcrvy["klfit"]["pos"]["fCoordinates"]["fY"], axis=None)
    # y_crvcoincs = ak.flatten(tcrvy["crvcoincs"]["crvcoincs.pos.fCoordinates.fY"], axis=None)

    t_klfit = ak.flatten(tcrvts["klfit"]["time"], axis=None)
    t_crvcoincs = ak.flatten(tcrvts["crvcoincs"]["crvcoincs.time"], axis=None)

    # nY_klfit = len(ak.flatten(tcrvy["klfit"]["pos"]["fCoordinates"]["fY"], axis=None))
    # nY_crvcoincs = len(ak.flatten(tcrvy["crvcoincs"]["crvcoincs.pos.fCoordinates.fY"], axis=None))


    if (len(t_klfit) != len(t_crvcoincs)): 

        print("\nTimes:")
        print("klfit = ", len(t_klfit))
        print("crvcoincs = ", len(t_crvcoincs)) 

        # Print out to find the discrepency 

        # Iterate event-by-event
        count = 0
        for i, event in enumerate(tcrvts):
            if event is None: continue
            # if len(ak.flatten(event["crvcoincs"]["crvcoincs.pos.fCoordinates.fY"], axis=None)) == len(ak.flatten(event["klfit"]["pos"]["fCoordinates"]["fY"],axis=None)): continue
            if len(ak.flatten(event["crvcoincs"]["crvcoincs.time"], axis=None)) == len(ak.flatten(event["klfit"]["time"],axis=None)): continue
            print(EventStr(event, showCutMasks=True))
            count += 1 

    # OK. So we get an event which has an empty list in klfit? 
    # I guess this is this goodTrk thing, which needs to be reapplied... 
    
    # 1
    # 60851
    # [[]]
    # [4.76e+03]
    # True
    # [[False, False, False, False, False]]
    # [True]

    # ut.Plot1D(data=ak.flatten(tcrvy["klfit"]["pos"]["fCoordinates"]["fY"], axis=None) - ak.flatten(tcrvy["crvcoincs"]["crvcoincs.pos.fCoordinates.fY"], axis=None)
    #         , nbins=250, xmin=-30, xmax=30
    #         , title="$\Delta T$", xlabel="$T_{KKInter} - T_{CRV}$ [ns]", ylabel="Counts", fout="../Images/CompCRV/h1_deltaT.png")
        
    ut.Plot1D(data=ak.flatten(tcrvts["klfit"]["time"], axis=None) - ak.flatten(tcrvts["crvcoincs"]["crvcoincs.time"], axis=None)
            , nbins=250, xmin=-30, xmax=30
            , title="$\Delta T$", xlabel="$T_{KKInter} - T_{CRV}$ [ns]", ylabel="Counts", fout="../Images/CompCRV/h1_deltaT.png")
        
    
    # With 1000 events there is one discrenpancy 
    # One klfit less that coincidences.   
    # 231
    # 232

    # 1205
    # 0
    # 216980
    # [[4.78e+03]]
    # [4.77e+03, 4.76e+03]

    # So in this case we have multiple coincidences 
    # Can we look at the event display? Done.


        # event[f'evtinfo.'][f'evtinfo.runid']
        # if i >= 10 - 1: return

    # PrintNEvents(tcrvy["crvcoincs"][multipleCoincidence], showCutMasks=True)
    # print(len(ak.flatten(tcrvy)

    # 47941
    # 48047

    # 8637
    # 8628

    # OK. With single CRV 1 hits, we reduce to 

    # 47795
    # 47833

    # I think the next thing is electron coincidences, can I find an example? 

    # With an muon only cut, we have

    # 47795
    # 47833
    # Length difference of 38

    # So that's not it. It's also not multiple hits? 


    # How to do an apples to apples comparison of klfit and crvcoincs? 


    # tcrvy["klfit"]["pos"]["fCoordinates"]["fY"] 



    # # Plot
    # ut.Plot2D(x=ak.flatten(tcrvpg["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None), y=ak.flatten(tcrvpg["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    #         , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
    #         , title="Has coincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_hasCoin.png")

    # ut.Plot2D(x=ak.flatten(tcrvpb["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None), y=ak.flatten(tcrvpb["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    #         , nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000
    #         , title="No coincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_noCoin.png")

def Run(finName):

    # TrkAna tutorial:
    # https://github.com/Mu2e/TrkAna/blob/main/tutorial

    # We have event level, global coincidence level, global track level, and local track level
    # Trying to mask the entire array based on, for example, a local track level mask, is tricky. 
    # One approach is to split the array into the different levels/trees and mask the appropirate tree
    arrays = ak.Array([])
    with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
        arrays = tree.arrays(["evtinfo.", "crvcoincs", "crvcoincsmc", "kl", "klfit", "klkl"]) 

    # arrays = arrays[:2500]

    # Mark cuts 
    MarkCuts(arrays)

    # arrays = ApplyCuts(arrays, ["goodTrk", "goodCRV", "CRV1", "KLCRV1"] ) # , "CRV1"]) # , "KLCRV1", "bestFit"])

    # PrintNEvents(arrays, nEvents=20, showCutMasks=True)

    # return

    # return

    # PrintNEvents(arrays, showCutMasks=True)
    
    # Plot 
    # Plot1(arrays) # KKInter TCRV Layer 1 Position, with/without coincidence
    Plot2(arrays) # KKInter Time - CRV Time Layer 1 

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