# Samuel Grant 2024
# Replicate Dave Brown's CompCRV.C in Python 

import sys
import awkward as ak
import uproot
import math

import Utils as ut

coincsBranchName = "crvcoincs"

branchNames_ = [ 

    # ---> evtinfo
    "evtinfo.runid" # run ID 
    , "evtinfo.subrunid" # sub-run ID 
    , "evtinfo.eventid" # event ID 

    # ---> crvhit (reco)
    , "crvcoincs.sectorType" # CRV sector hit
    , "crvcoincs.pos.fCoordinates.fX" # Reconstructed position of the cluster in X 
    , "crvcoincs.pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
    , "crvcoincs.pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
    , "crvcoincs.timeStart" # Earliest time recorded at either end of all bars in the hit
    , "crvcoincs.timeEnd" # Latest time recorded at either end of all bars in the hit
    , "crvcoincs.time" # average reconstructed hit time of the cluster.
    , "crvcoincs.PEs" # total number of photoelectrons in this cluser
    , "crvcoincs.nHits" # Number of individual bar hits combined in this hit
    , "crvcoincs.nLayers" # Number of CRV layers that are part of this cluster
    , "crvcoincs.PEsPerLayer[4]"
    , "crvcoincs.angle" # slope (in the plane perpendicular to the bar axis of a sector) of the track assumed to be responsible for the cluster (=change in the "layer direction" / change in the "thickness direction")

    # ---> crvhitmc (truth)
    , "crvcoincsmc.valid" # Records if there is a valid MC match to this CRV reco hit
    , "crvcoincsmc.pdgId" # PDG ID of the track mostly likely responsible for this cluster

    # ----> kl (tracks)
    , "kl.status"

]

branchNamesTrkAna_ = [

    # ---> evtinfo
    "evtinfo.runid" # run ID 
    ,"evtinfo.subrunid" # sub-run ID 
    ,"evtinfo.eventid" # event ID 

    # ---> crvsummary
    ,"crvsummary.nHitCounters"
    
    # ---> crvhit (reco)
    , "crvcoincs.sectorType" # CRV sector hit
    , "crvcoincs.pos.fCoordinates.fX" # Reconstructed position of the cluster in X 
    , "crvcoincs.pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
    , "crvcoincs.pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
    , "crvcoincs.timeStart" # Earliest time recorded at either end of all bars in the hit
    , "crvcoincs.timeEnd" # Latest time recorded at either end of all bars in the hit
    , "crvcoincs.time" # average reconstructed hit time of the cluster.
    , "crvcoincs.PEs" # total number of photoelectrons in this cluser
    , "crvcoincs.nHits" # Number of individual bar hits combined in this hit
    , "crvcoincs.nLayers" # Number of CRV layers that are part of this cluster
    , "crvcoincs.PEsPerLayer[4]"
    , "crvcoincs.angle" # slope (in the plane perpendicular to the bar axis of a sector) of the track assumed to be responsible for the cluster (=change in the "layer direction" / change in the "thickness direction")

    # ---> crvhitmc (truth)
    , "crvcoincsmc.valid" # Records if there is a valid MC match to this CRV reco hit
    , "crvcoincsmc.pdgId" # PDG ID of the track mostly likely responsible for this cluster
    # , "crvcoincsmc.pdgId"
    # , "crvhitmc.primaryPdgId" # PDG ID of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.primaryE" # energy of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.primary.pos.fCoordinates.fX" # start position of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parentPdgId" # PDG ID of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parentE" # start energy of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parent.fCoordinates.fX" # X start position of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parent.fCoordinates.fY" # Y
    # , "crvhitmc.parent.fCoordinates.fZ" # Z
    # , "crvhitmc.gparentPdgId" # grandparent info..
    # , "crvhitmc.gparentE" # "
    # , "crvhitmc.gparent" # "
    # , "crvhitmc.depositedEnergy" # total deposited energy of the cluster based on the CrvSteps

    # ----> kl (tracks)
    , "kl.status"
    , "kl.nactive"
    , "kl.nhits"
    , "kl.nplanes"
    , "kl.nnullambig"
    , "kl.ndof"
    , "kl.fitcon"
    # # arrays of structs, handled slightly differently 
    , "klfit" 
    , "klkl"
    # , "kl.fitcon"
    # , "klfit.sid"
    # , "klfit.sindex"
    # , "klfit.pos.X"
    # , "klfit.pos.Y"
    # , "klfit.pos.Z"
    # , "klfit.time"
    # , "klkl.z0err"
    # , "klkl.d0err"
    # , "klkl.thetaerr"
    # , "klkl.phi0err"
]


def EventStr(event):
    
    eventStr = (
        f"evtinfo.runid: {event['evtinfo.runid']}\n"
        f"evtinfo.subrunid: {event['evtinfo.subrunid']}\n"
        f"evtinfo.eventid: {event['evtinfo.eventid']}\n"
        f"crvsummary.nHitCounters: {event['crvsummary.nHitCounters']}\n"
        f"crvcoincs.nLayers {event[f'crvcoincs.nLayers']}\n"
        f"crvcoincs.angle: {event[f'crvcoincs.angle']}\n"
        f"crvcoincs.sectorType: {event[f'crvcoincs.sectorType']}\n"
        f"crvcoincs.pos.fCoordinates: ({event[f'crvcoincs.pos.fCoordinates.fX']}, {event[f'crvcoincs.pos.fCoordinates.fY']}, {event[f'crvcoincs.pos.fCoordinates.fZ']})\n"
        f"crvcoincs.timeStart: {event[f'crvcoincs.timeStart']}\n"
        f"crvcoincs.timeEnd: {event[f'crvcoincs.timeEnd']}\n"
        f"crvcoincs.time: {event[f'crvcoincs.time']}\n"
        f"crvcoincs.PEs: {event[f'crvcoincs.PEs']}\n"
        f"crvcoincs.nHits: {event[f'crvcoincs.nHits']}\n"
        f"crvcoincsmc.valid: {event[f'crvcoincsmc.valid']}\n"
        f"crvcoincsmc.pdgId: {event[f'crvcoincsmc.pdgId']}\n"
        f"kl.status: {event[f'kl.status']}\n"
        f"kl.nactive: {event[f'kl.nactive']}\n"
        f"kl.nhits: {event[f'kl.nhits']}\n"
        f"kl.nplanes: {event[f'kl.nplanes']}\n"
        f"kl.nnullambig: {event[f'kl.nnullambig']}\n"
        f"kl.ndof: {event[f'kl.ndof']}\n"
        f"kl.fitcon: {event[f'kl.fitcon']}\n"
        # f"klfit: {event[f'klfit']}\n"
        # f"klfit.sid: {event[f'klfit.sid']}\n"
        # f"klfit.sindex: {event[f'klfit.sindex']}\n"
        f"klfit.sid: {event[f'klfit']['sid']}\n"
        f"klfit.sindex: {event[f'klfit']['sindex']}\n"
        f"klfit.pos.X() {event[f'klfit']['pos']['fCoordinates']['fX']}\n"
        f"klfit.pos.Z() {event[f'klfit']['pos']['fCoordinates']['fZ']}\n"
        f"klkl.z0err: {event['klkl']['z0err']}\n"
        f"goodTrk: {event[f'goodTrk']}\n"
        f"CRV1: {event[f'CRV1']}\n"
        f"KLCRV1: {event[f'KLCRV1']}\n"
        f"goodCRV: {event[f'goodCRV']}\n"
        f"noCRV: {event[f'noCRV']}\n"
        f"bestFit1: {event[f'bestFit1']}\n"
        f"bestFit2: {event[f'bestFit2']}\n"
        f"L1Fiducial: {event[f'L1Fiducial']}\n"
    )

    # # bestFit needs to be seperated by kl and klkl 
    # data_["bestFit1"] = ( (data_["kl.ndof"] >= 10)
    #                     & (data_["kl.fitcon"] > 0.1)
    #                     & ((data_["kl.nactive"]/data_["kl.nhits"]) > 0.99)
    #                     & (data_["kl.nplanes"] >= 4)
    #                     & ((data_["kl.nnullambig"]/data_["kl.nhits"]) < 0.2) )

    # data_["bestFit2"] =  ( (data_["klkl"]["z0err"] < 1) 
    #                     & (data_["klkl"]["d0err"] < 1) 
    #                     & (data_["klkl"]["thetaerr"] < 0.004)
    #                     & (data_["klkl"]["phi0err"] < 0.001) )
    return eventStr

def PrintNEvents(data_, nEvents=10): 
    # Iterate event-by-event
    for i, event in enumerate(data_):
        if event is None: continue
        print(EventStr(event)) 
        if i >= nEvents - 1: return

def TTreeToAwkwardArray(finName, treeName, branchNames):

    print("\n---> Converting TTree to awkward array")

    # Open the ROOT file and access the TTree
    file = uproot.open(finName)
    tree = file[treeName]

    # Read the data into Awkward Arrays, if you exclude branchNames it reads all of them 
    arrays = tree.arrays(branchNames, library="ak")

    # Open the ROOT file and access the TTree
    with uproot.open(finName) as file:
        tree = file[treeName]
        # Read the data into Awkward Arrays
        arrays = tree.arrays(branchNames, library="ak")
        
    print("Done!")

    return arrays

# def MarkTrackerCuts(data_): 

#     data_["goodCRV"] = data_["crvsummary.nHitCounters"] > 0 # ak.num(data_["crvcoincs.nHits"]) > 0

#     # data_["noCRV"] = data_["crvsummary.nHitCounters"] == 0 
    
#     data_["goodTrk"] = data_["kl.status"] >= 0

#     data_["CRV1"] = (ak.num(data_["crvcoincs.sectorType"]) > 0) & (ak.any(data_["crvcoincs.sectorType"] == 1)) # data_["crvcoincs.sectorType"] == 1 if ak.num(data_["crvcoincs.sectorType"], axis=1) > 0 else False

#     # data_["KLCRV1"] =   ( (data_["klfit"]["sid"] == 200) 
#     #                     & (data_["klfit"]["sindex"] == 1) )

#     data_["KLCRV1"] =   ( (data_["klfit.sid"] == 200) 
#                         & (data_["klfit.sindex"] == 1) )

#     data_["bestFit1"] =  ( (data_["klkl"]["z0err"] < 1) 
#                         & (data_["klkl"]["d0err"] < 1) 
#                         & (data_["klkl"]["thetaerr"] < 0.004)
#                         & (data_["klkl"]["phi0err"] < 0.001) )

#     data_["bestFit2"] = ( (data_["kl.ndof"] >= 10)
#                         & (data_["kl.fitcon"] > 0.1)
#                         & ((data_["kl.nactive"]/data_["kl.nhits"]) > 0.99)
#                         & (data_["kl.nplanes"] >= 4)
#                         & ((data_["kl.nnullambig"]/data_["kl.nhits"]) < 0.2) )



#     data_["L1Fiducial"] = ( (abs(data_["klfit"]["pos"]["fCoordinates"]["fX"]) < 2500)
#                           & (abs(data_["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 1500) ) 

#     return 

def MarkTrackerCuts(data_): 

    # Tracks?
    data_["goodTrk"] = data_["kl.status"] >= 0

    # CRV hits?
    data_["goodCRV"] = ak.num(data_["crvcoincs.nHits"], axis=1) > 0
    data_["noCRV"] = ~data_["goodCRV"] # "goodCRV"]ak.num(data_["crvcoincs.nHits"], axis=1) > 0

    # CRV hits in sector 1? 
    data_["CRV1"] = data_["crvcoincs.sectorType"] == 1

    # Tracks intersecting the CRV?
    data_["KLCRV1"] =   ( (data_["klfit"]["sid"] == 200) 
                        & (data_["klfit"]["sindex"] == 1) )

    # Best fit? 
    # bestFit needs to be seperated by kl and klkl 
    data_["bestFit1"] = ( (data_["kl.ndof"] >= 10)
                        & (data_["kl.fitcon"] > 0.1)
                        & ((data_["kl.nactive"]/data_["kl.nhits"]) > 0.99)
                        & (data_["kl.nplanes"] >= 4)
                        & ((data_["kl.nnullambig"]/data_["kl.nhits"]) < 0.2) )

    data_["bestFit2"] =  ( (data_["klkl"]["z0err"] < 1) 
                        & (data_["klkl"]["d0err"] < 1) 
                        & (data_["klkl"]["thetaerr"] < 0.004)
                        & (data_["klkl"]["phi0err"] < 0.001) )

    data_["L1Fiducial"] = ( (abs(data_["klfit"]["pos"]["fCoordinates"]["fX"]) < 2500)
                          & (abs(data_["klfit"]["pos"]["fCoordinates"]["fZ"] + 500) < 1500) ) 

    return 

# Maybe a reasonable solution? Still ugly though. 
# Apply cuts sequentially according to if statements
def ApplyTrackerCuts(data_, cuts_): 

    # Check if MarkTrackerCutsTest has been applied here

    # Make a deep copy (no shared memory)
    dataCopy_ = ak.copy(data_)

    if "goodCRV" in cuts_:
        mask = dataCopy_["goodCRV"] # Outer array mask
        # dataCopy_ = ak.mask(dataCopy_, mask) 
        dataCopy_ = dataCopy_[mask] # Filter everything 
    if "noCRV" in cuts_:
        mask = dataCopy_["noCRV"] # Outer array mask
        # dataCopy_ = ak.mask(dataCopy_, mask) 
        dataCopy_ = dataCopy_[mask] # Filter everything 
    if "goodTrk" in cuts_: 
        mask = ak.any(dataCopy_["goodTrk"], axis=1) # Inner array mask
        dataCopy_ = ak.mask(dataCopy_, mask) # Mark non-matching elements as None
    if "CRV1" in cuts_: 
        mask = ak.any(dataCopy_["CRV1"], axis=1) # Inner array mask
        dataCopy_ = ak.mask(dataCopy_, mask)
    if "KLCRV1" in cuts_:
        mask = dataCopy_["KLCRV1"] # Outer array mask for klfit only
        dataCopy_["klfit"] = ak.mask(dataCopy_["klfit"], mask) 
    if "bestFit" in cuts_:
        mask = ak.any(dataCopy_["bestFit1"], axis=1) # Inner array mask 
        dataCopy_ = ak.mask(dataCopy_, mask)
        mask = dataCopy_["bestFit2"] # Outer array mask for klkl
        # dataCopy_ = ak.mask(dataCopy_["klkl"], mask)
        dataCopy_["klkl"] = dataCopy_["klkl"][mask]
        # dataCopy_ = data_[ak.any(dataCopy_["bestFit1"], axis=1)]
        # dataCopy_["klkl"] = ak.mask(dataCopy_["klkl"], dataCopy_["bestFit2"]) 
    if "L1Fiducial" in cuts_:
        mask = ak.any(dataCopy_["L1Fiducial"], axis=1) # Inner array mask 
        dataCopy_ = ak.mask(dataCopy_, mask)
    
    return dataCopy_


def Run(finName):
    # Load data
    data_ = TTreeToAwkwardArray(finName, "TrkAnaExt/trkana", branchNamesTrkAna_)

    # Mark the tracker cuts 
    MarkTrackerCuts(data_)

    # Apply cuts 
    # goodtrk+KLCRV1+CRV1+bestfit+goodCRV
    tcrvpg_ = ApplyTrackerCuts(data_, ["goodCRV", "CRV1", "goodTrk", "bestFit", "KLCRV1"]) 
    # goodtrk+KLCRV1+bestfit+noCRV
    tcrvpb_ = ApplyTrackerCuts(data_, ["noCRV", "goodTrk", "bestFit", "KLCRV1"])

    # Flatten
    # See https://awkward-array.org/doc/main/user-guide/how-to-restructure-flatten.html
    x_hasCoin = ak.flatten(tcrvpg_["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    z_hasCoin = ak.flatten(tcrvpg_["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)

    x_noCoin = ak.flatten(tcrvpb_["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    z_noCoin = ak.flatten(tcrvpb_["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)

    # Plot
    ut.Plot2D(x=z_hasCoin, y=x_hasCoin, nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000, xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_hasCoin.png")
    ut.Plot2D(x=z_noCoin, y=x_noCoin, nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000, xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_noCoin.png")

    return

    # # Cut 1 
    # goodTrk = ak.any(data_["goodTrk"], axis=1) 
    # CRV1 = data_["CRV1"] 
    # bestFit2 = ak.any(data_["bestFit2"], axis=1)
    # # Apply regular masks 
    # tcrvpg_ = data_[goodTrk & CRV1 & bestFit2]

    # # Now klkl cuts, this is unbelievably awkward!
    # tcrvpg_["klkl"] = ak.mask(tcrvpg_["klkl"], tcrvpg_["bestFit1"]) 
    # # And klfit cuts
    # tcrvpg_["klfit"] = ak.mask(tcrvpg_["klfit"], tcrvpg_["KLCRV1"])  

    # # PrintNEvents(tcrvpg_, 50)

    # # See https://awkward-array.org/doc/main/user-guide/how-to-restructure-flatten.html
    # x = ak.flatten(tcrvpg_["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    # z = ak.flatten(tcrvpg_["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)

    # ut.Plot2D(z, x, nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000, title="KKInter TCRV Layer 1 Position, Has CRVCoincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_hasCoin.png")

    # # TH2D* tcrvpb = new TH2D("tcrvpb","KKInter TCRV Layer 1 Position,  No CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    # # ta->Project("tcrvpb","klfit.pos.X():klfit.pos.Z()",goodTrk+KLCRV1+bestFit+noCRV);
    # noCRV = data_["noCRV"]
    # tcrvpb_ = data_[goodTrk & noCRV & bestFit2]
    # tcrvpb_["klkl"] = ak.mask(tcrvpb_["klkl"], tcrvpb_["bestFit1"]) 
    # tcrvpb_["klfit"] = ak.mask(tcrvpb_["klfit"], tcrvpb_["KLCRV1"]) 

    # x = ak.flatten(tcrvpb_["klfit"]["pos"]["fCoordinates"]["fX"], axis=None)
    # z = ak.flatten(tcrvpb_["klfit"]["pos"]["fCoordinates"]["fZ"], axis=None)

    # ut.Plot2D(z, x, nbinsX=100, xmin=-8000, xmax=8000, nbinsY=100, ymin=-8000, ymax=8000, title="KKInter TCRV Layer 1 Position, No CRVCoincidence", xlabel="KKInter Z [mm]", ylabel="KKInter X [mm]", fout="../Images/CompCRV/h2_zx_noCoin.png")

    # return

    # print(data_)

    # return

    # data_["klfit.sindex"] = data_["klfit"]["sindex"]
    # data_["klfit.sid"] = data_["klfit"]["sid"]

    # data_ = data_[[column for column in data_.fields if column != "klfit"]]
    # data_ = data_[[column for column in data_.fields if column != "kl"]]
    # PrintNEvents(data_, 1)

    # This is extremely awkward!
    # klfit and klkl masks must be apply directly to the fields 
    # data_["klfit"] = ak.mask(data_["klfit"], data_["KLCRV1"]) 

    # mask1 = ak.mask(data_["klfit"], data_["KLCRV1"])#  & data_["CRV1"] 
    # # mask2 = data_["CRV1"] # ak.any(data_["KLCRV1"], axis=2)

    # PrintNEvents(data_[mask1], 1)

    # Apply the mask to specific fields
    # mask1 = ak.mask(array["kl"], data_["KLCRV1"])

    # # Combine field1 with the masked field2 using zip
    # data_ = ak.zip({"field1": array["field1"], "field2": masked_field2})

    # return

    # 
    # mask = data_["KLCRV1"]
    # # masked_data = data_[mask]

    # flat_mask = ak.flatten(mask)
    # masked_data = ak.flatten(data_)[flat_mask]

    # evtinfo.runid: 1205
    # evtinfo.subrunid: 0
    # evtinfo.eventid: 2582
    # crvcoincs.nLayers [4]
    # crvcoincs.angle: [-0.39]
    # crvcoincs.sectorType: [1]
    # crvcoincs.pos.fCoordinates: ([2.74e+03], [4.76e+03], [-1.4e+03])
    # crvcoincstimeStart: [6.01e+04]
    # crvcoincs.timeEnd: [6.02e+04]
    # crvcoincs.time: [6.01e+04]
    # crvcoincs.PEs: [372]
    # crvcoincs.nHits: [16]
    # crvcoincsmc.valid: [True]
    # crvcoincsmc.pdgId: [11]
    # kl.status: [1]
    # kl.nactive: [24]
    # kl.nhits: [24]
    # kl.nplanes: [6]
    # kl.nnullambig: [4]
    # kl.ndof: [49]
    # kl.fitcon: [0.996]
    # klfit.sid: [[0, 2, 4, 4, 200, 200, 200]]
    # klfit.sindex: [[0, 0, 0, 0, 0, 1, 2]]
    # KLCRV1: [[False, False, False, False, False, True, False]]

    # I think that klfit and klkl cuts can't be treated on an event level. They only apply to parameters on their level. 
    # We want the 200, 1 fit within this event. If it contains at least one of these then we keep the event. 
    # Isn't it just ak.any(axis=2) then? I'm confused. 

    # MarkTrackerCuts(data_)

    # PrintNEvents(data_, 1)

    # mask1 = (
    #     ak.any(data_["klfit"]["sid"] == 200, axis=1)
    #     & ak.any(data_["klfit"]["sindex"] == 1, axis=1)
    # )

    # print(data_) 

    # print(data_.fields)

    # ta->Project("tcrvpg","klfit.pos.X():klfit.pos.Z()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("tcrvpb","klfit.pos.X():klfit.pos.Z()",goodTrk+KLCRV1+bestFit+noCRV);
    # ta->Project("tcrvy","klfit.pos.Y()-crvcoincs.pos.Y()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("tcrvts","klfit.time-crvcoincs.time",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("dtvz","klfit.time-crvcoincs.time:klfit.pos.Z()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("dtvx","klfit.time-crvcoincs.time:klfit.pos.X()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("xvx","crvcoincs.pos.X():klfit.pos.X()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);
    # ta->Project("zvz","crvcoincs.pos.Z():klfit.pos.Z()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);

    # axis = 0 is the outer array structure, axis = 1 is the inner array structure, axis = 2 is the inner-inner array structure
    # mask1 = ak.any(data_["goodTrk"], axis=1) # & ak.any(data_["KLCRV1"], axis=1) & data_["CRV1"] & ak.any(data_["bestFit"], axis=1) & data_["goodCRV"]
    # mask1 = data_["KLCRV1"] # ak.any(data_["KLCRV1"], axis=2)
    # print("data_['KLCRV1']\n", data_["KLCRV1"])
    # mask1 = ak.any(data_["goodTrk"], axis=1
    # print("mask\n", mask1)
    # fuckyou = data_[mask1]
    # PrintNEvents(data_, 1)
    # print(data_["KLCRV1"])
    # mask1 = data_["KLCRV1"]
    # masked_data = data_[mask1]
    # PrintNEvents(data_, 1)
    
    # return

    # data_["goodTrk"]]
    # print(data_["goodTrk"])
    # data_ = data_[data_["goodTrk"]]

    # ut.Plot2D(data_[mask1]["klfit"]["pos"]["fCoordinates"]["fX"], data_[mask1]["klfit"]["pos"]["fCoordinates"]["fZ"], nbinsX=100,xmin=-8000,xmax=8000,nBinsY=100,ymin=-8000,ymax=8000)
    # TH2D* tcrvpg = new TH2D("tcrvpg","KKInter TCRV Layer 1 Position, Has CRVCoincidence;KKInter Z (mm);KKInter X (mm)",100,-8000,8000,100,-8000,8000);
    # ta->Project("tcrvpg","klfit.pos.X():klfit.pos.Z()",goodTrk+KLCRV1+CRV1+bestFit+goodCRV);

    # print(data_["L1Fiducial"])

    # Test cut

    # return 


def main():

    # Take command-line arguments
    finName = sys.argv[1] if len(sys.argv) > 1 else "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/40/73/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # /pnfs/mu2e/scratch/users/sgrant/workflow/CosmicCRYExtractedTrk.MDC2020z2_best_v1_1/outstage/67605881/00/00000/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"

    print("\n--->Running with inputs:\n")
    print("\tfinName:", finName)

    Run(finName=finName) 

    # print(data_)
    # print(data_["klfit"]['mom']['fCoordinates']["fY"])

    return

if __name__ == "__main__":
    main()
