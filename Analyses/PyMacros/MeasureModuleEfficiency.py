'''
Samuel Grant 2024

Analyse CRV module efficiency from TrkAna in KPP geometry. 

KPP: three tiers, top and bottom tiers are trigger tiers, middle is the measurement tier.
* Two CRV-L-end modules parallel to the DS axis (sector 3).
* Four CRV-T modules perpendicular to the DS axis (sector 1).
* Two CRV-DS modules perpendicular to the DS axis (sector 2).

Procedure: 

Use top and bottom modules to trigger on cosmics, look for the failure in the measurement module. 
 
First pass: trigger on events which have ONE coincidience in the top and bottom tier (2/4 layers hit and 10 PEs per layer), store failures.
Second pass: analyses failures using time window and spatial grouping between trigger modules. 
Third pass: analyses failure using the tracker. 

In each case: 
Count for successful and unsuccessful triggers. 
Calculate the inefficiency, write failules and effiency results to file. 

Notes: 

In future, we can expand to >1 coin in the triggers, pass directly these directly to the grouping algorithm.

Other questions: 
What is the purity of the tracker? 
What is the angular distribtuion that the tracker can handle with straight tracks? 
What is the efficiency of the tracker? 

We are trying to learn how well we can use the tracker to measure the module efficiency in situ. 


'''
# External libraries
import sys
import uproot 
import numpy as np
import awkward as ak

# Internal libraries
import Utils as ut
import CoincidenceConditions as cc

# ------------------------------------------------
#                Debugging functions 
# ------------------------------------------------

def SanityPlots(data_, reproc, foutTag):
    
    #TODO: add tracker sanity plots.
    print("\n---> Making sanity plots")

    # Reco
    sectors_ = ak.flatten(data_["crv"]["crvcoincs.sectorType"]) 
    t_ = ak.flatten(data_["crv"]["crvcoincs.time"]) 
    x_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fX"]) * 1e-3 # mm -> m
    y_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fY"]) * 1e-3 # mm -> m
    z_ = ak.flatten(data_["crv"]["crvcoincs.pos.fCoordinates.fZ"]) * 1e-3 # mm -> m
    PEs_ = ak.flatten(data_["crv"]["crvcoincs.PEs"]) 
    PEsPerLayer_ = ak.flatten(data_["crv"]["crvcoincs.PEsPerLayer[4]"])
    nHits_ = ak.flatten(data_["crv"]["crvcoincs.nHits"]) 
    nLayers_ = ak.flatten(data_["crv"]["crvcoincs.nLayers"]) 
    slopes_ = ak.flatten(data_["crv"]["crvcoincs.angle"])

    ut.PlotGraphOverlay(graphs_=[(x_[sectors_ == 3], y_[sectors_ == 3]), (x_[sectors_ == 1], y_[sectors_ == 1]), (x_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="x-position [m]", ylabel="y-position [m]", fout=f"../Images/{reproc}/Sanity/gr_XY_{foutTag}.png")
    ut.PlotGraphOverlay(graphs_=[(z_[sectors_ == 3], y_[sectors_ == 3]), (z_[sectors_ == 1], y_[sectors_ == 1]), (z_[sectors_ == 2], y_[sectors_ == 2])], labels_=["Top", "Middle", "Bottom"], xlabel="z-position [m]", ylabel="y-position [m]", fout=f"../Images/{reproc}/Sanity/gr_ZY_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[t_[sectors_ == 3], t_[sectors_ == 1], t_[sectors_ == 2]], nbins=1000, xmin = np.min(t_), xmax = np.max(t_), xlabel="Average hit time [ns]", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_times_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nHits_[sectors_ == 3], nHits_[sectors_ == 1], nHits_[sectors_ == 2]], nbins=41, xmin = -0.5, xmax = 40.5, xlabel="Number of hits", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_nHits_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[nLayers_[sectors_ == 3], nLayers_[sectors_ == 1], nLayers_[sectors_ == 2]], nbins=5, xmin = -0.5, xmax = 4.5, xlabel="Number of layers hit", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_nLayers_{foutTag}.png") 
    ut.Plot1DOverlay(hists_=[slopes_[sectors_ == 3], slopes_[sectors_ == 1], slopes_[sectors_ == 2]], nbins=1000, xmin = -2, xmax = 2, xlabel="Slope", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_slopes_{foutTag}.png") 

    # Calculate the sum of each length-4 array in PEsPerLayer_
    PEsPerLayerSum_ = ak.sum(PEsPerLayer_, axis=-1)
    PEsDiff_ = PEs_ - PEsPerLayerSum_

    ut.Plot1DOverlay(hists_=[PEs_[sectors_ == 3], PEs_[sectors_ == 1], PEs_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEs_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsPerLayerSum_[sectors_ == 3], PEsPerLayerSum_[sectors_ == 1], PEsPerLayerSum_[sectors_ == 2]], nbins=1000, xmin = 0, xmax = 1000, title="PEs per layer", xlabel="PEs", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEsPerLayer_{foutTag}.png")

    ut.Plot1DOverlay(hists_=[PEsPerLayer_[:, 0], PEsPerLayer_[:, 1], PEsPerLayer_[:, 2], PEsPerLayer_[:, 3]], nbins=500, xmin = 0, xmax = 500, title="PEs per layer", ylabel="Coincidences", label_=["Layer0", "Layer 1", "Layer 2", "Layer 3"], fout=f"../Images/{reproc}/Sanity//h1_PEsPerLayerAllSectors_{foutTag}.png")
    ut.Plot1DOverlay(hists_=[PEsDiff_[sectors_ == 3], PEsDiff_[sectors_ == 1], PEsDiff_[sectors_ == 2]], nbins=40, xmin = -0.0002, xmax = 0.0002, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_PEsDiff_{foutTag}.png", logY=True)
    ut.Plot1DOverlay(hists_=[PEsDiff_], nbins=200, xmin = np.min(PEsDiff_)-1e-4, xmax = np.max(PEsDiff_)+1e-4, xlabel="PEs $\minus$ PEs per layer sum", ylabel="Coincidences", label_=["All modules"], fout=f"../Images/{reproc}/Sanity/h1_PEsDiffAll_{foutTag}.png", logY=True, includeBlack=True)

    # MC 
    valid_ = ak.flatten(data_["crv"]["crvcoincsmc.valid"])
    pdgid_ = ak.flatten(data_["crv"]["crvcoincsmc.pdgId"])

    ut.Plot1DOverlay(hists_=[valid_[sectors_ == 3], valid_[sectors_ == 1], valid_[sectors_ == 2]], nbins=2, xmin = 0, xmax = 2, xlabel="Valid MC event", ylabel="Coincidences", label_=["Top", "Middle", "Bottom"], fout=f"../Images/{reproc}/Sanity/h1_valid_{foutTag}.png") 
    
    label_dict = {
        2212: 'proton',
        211: 'pi+',
        -211: 'pi-',
        -13: 'mu+',
        13: 'mu-',
        -11: 'e+',
        11: 'e-',
        "other": "other"
        # Add more particle entries as needed
    }

    ut.BarChart(data_=pdgid_, label_dict=ut.particle_dict, ylabel="Coincidences [%]", fout=f"../Images/{reproc}/Sanity/bar_pdgid_{foutTag}.png", percentage=True)
    ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout=f"../Images/{reproc}/Sanity/bar_overlay_pdgid_{foutTag}.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

    print("...Done!")

    return

# Printout helper functions
def PrintEvent(event, showMasks):
    
    eventStr = (
        f"-------------------------------------------------------------------------------------\n"
        f"evt...\n"
        f"evtinfo.run: {event['evt']['evtinfo.run']}\n" 
        f"evtinfo.subrun: {event['evt']['evtinfo.subrun']}\n" 
        f"evtinfo.eventid: {event['evt']['evtinfo.event']}\n"
        f"crv...\n"
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
        f"trk...\n"
        f"kl.status: {event['trk']['kl.status']}\n"
        f"kl.nactive: {event['trk']['kl.nactive']}\n"
        f"kl.nhits: {event['trk']['kl.nhits']}\n"
        f"kl.nplanes: {event['trk']['kl.nplanes']}\n"
        f"kl.nnullambig: {event['trk']['kl.nnullambig']}\n"
        f"kl.ndof: {event['trk']['kl.ndof']}\n"
        f"kl.kl.fitcon: {event['trk']['kl.fitcon']}\n"
        f"trkfit...\n"
        f"klfit: {event['trkfit']['klfit']}\n"
        # f"Example variables from klfit...\n"
        # f"klfit.sid: {event['trkfit']['klfit']['sid']}\n"
        # f"klfit.sindex: {event['trkfit']['klfit']['sindex']}\n"
        # f"klfit: {event['trkfit']['klfit']}\n"
        f"klkl: {event['trkfit']['klkl']}\n"
        f"-------------------------------------------------------------------------------------\n"
    )

    if showMasks: 
        eventStr = eventStr[:-86]
        eventStr += f"masks...\n" 
        eventStr += f"PE_condition: {event['PE_condition']}\n" 
        eventStr += f"layer_condition: {event['layer_condition']}\n" 
        eventStr += f"angle_condition: {event['angle_condition']}\n" 
        eventStr += f"is_coincidence: {event['is_coincidence']}\n" 
        eventStr += f"trigger: {event['trigger']}\n" 
        # eventStr += f"crv_mask: {event['crv_mask']}\n" 
        # eventStr += f"trk_mask: {event['trk_mask']}\n" 
        # eventStr += f"trkfit_mask: {event['trkfit_mask']}\n" 
        eventStr += f"-------------------------------------------------------------------------------------\n"
    
    return eventStr

def PrintNEvents(data_, nEvents=10, showMasks=False):
     # Iterate event-by-event
    for i, event in enumerate(data_, start=1):
        print(PrintEvent(event, showMasks))
        if i >= nEvents: 
            return

# ------------------------------------------------
#               Coincidence finding 
# ------------------------------------------------

# Impose stricter conditions than the default, if desired
def FindCoincidences(data_, coincidenceConditions="10PEs2Layers"): 
    
    print("\n---> Marking coincidences")

    # Get conditions from file
    coincidenceConditions_ = cc.coincidenceConditions_[coincidenceConditions]

    # PE threshold per layer
    print("* PE threshold condition")
    PEs_ = data_["crv"]["crvcoincs.PEsPerLayer[4]"]
    
    # PE condition per layer
    PECondition = PEs_ >= coincidenceConditions_["PEs"]

    # Number of layers hit 
    # (based on number of coincidences with PEs per layer above threshold, not nLayers which is baked into the reconstruction in mcs)
    print("* Layers condition")
    layerCondition = ak.count(PEs_[PECondition], axis=2) >= coincidenceConditions_["nLayers"] # This can be adjusted on the level of nts 

    # Hit slope: horizontal / vertical direction 
    print("* Angle condition")
    angleCondition = (abs(data_["crv"]["crvcoincs.angle"]) <= coincidenceConditions_["maxSlope"]) & (abs(data_["crv"]["crvcoincs.angle"]) >= coincidenceConditions_["minSlope"])

    # Combine conditions to mark per-module coincidences
    coincidenceMask = layerCondition & angleCondition # PECondition is implicit in the layer condition!

    # Add a new field 'is_coincidence' to mark coincidences
    data_["is_coincidence"] = coincidenceMask

    # Mark the individual coniditions for debugging 
    data_["PE_condition"] = PECondition 
    data_["layer_condition"] = layerCondition 
    data_["angle_condition"] = angleCondition 

    print("...Done!")

    return data_ # ["crv"][data_["is_coincidence"]]

# ------------------------------------------------
#                     Filtering 
# ------------------------------------------------ 

# Not needed.
# def RemoveEmptyEvents(data_):
#     print(f"\n---> Removing empty events")
#     emptyEventCondition = ak.num(data_["crv"]["crvcoincs.nHits"]) == 0
#     return data_[~emptyEventCondition]

# Not needed.
# Currently just "all", "muons", "non_muons"
# Can generalise if needed
# def FilterParticles(data_, particle):
#     print(f"\n---> Filtering particles, keeping {particle}")
#     # This is really just for debugging. 
#     muonCondition = ak.any((data_["crv"]["crvcoincsmc.pdgId"] == 13) | (data_["crv"]["crvcoincsmc.pdgId"] == -13), axis=1)
#     if particle == "all":
#         return data_
#     elif particle == "muons": 
#         return data_[muonCondition] 
#     elif particle == "non_muons":
#         return data_[~muonCondition] 
#     else:
#         raise ValueError(f"Particle string {particle} not valid!")

#  Events with one coincidence in sectors 2 & 3
def PassOneFilter(data_):
    
    print(f"\n---> Running pass one filter") 

    sector2Condition = data_["crv"]["crvcoincs.sectorType"] == 2
    sector3Condition = data_["crv"]["crvcoincs.sectorType"] == 3

    oneHitInSector2Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector2Condition], axis=1) == 1
    oneHitInSector3Condition = ak.count(data_["crv"]["crvcoincs.sectorType"][sector3Condition], axis=1) == 1

    # Use ak.mask here, since the output from ak.count will not match the dimensions of sectorType array. 
    # That is, a element like [2, 3] gets masked with "True" for the whole array
    # data_["crv"] = ak.mask(data_["crv"], (oneHitInSector2Condition & oneHitInSector3Condition))

    return data_

# Events that pass the tracker cuts 
def PassThreeFilter(data_):
    
    print(f"\n---> Running pass three filter") 
    
    data_["trkfit_KLCRV1"] = ( 
        (data_["trkfit"]["klfit"]["sid"] == 200) 
        & (data_["trkfit"]["klfit"]["sindex"] == 1) )

    data_["trk_bestFit"] = ( 
        (data_["trk"]["kl"]["kl.ndof"] >= 10)
        & (data_["trk"]["kl"]["kl.fitcon"] > 0.1)
        & ((data_["trk"]["kl"]["kl.nactive"]/data_["trk"]["kl"]["kl.nhits"]) > 0.99)
        & (data_["trk"]["kl"]["kl.nplanes"] >= 4)
        & ((data_["trk"]["kl"]["kl.nnullambig"]/data_["trk"]["kl"]["kl.nhits"]) < 0.2) )
    
    data_["trkfit_bestFit"] = ( 
        (data_["trkfit"]["klkl"]["z0err"] < 1) 
        & (data_["trkfit"]["klkl"]["d0err"] < 1) 
        & (data_["trkfit"]["klkl"]["thetaerr"] < 0.004)
        & (data_["trkfit"]["klkl"]["phi0err"] < 0.001) )
    
    # data_["trkfit"] = data_["trkfit"][data_["trkfit_bestFit"] & data_["trkfit_KLCRV1"]]
    # data_["trk"] = data_["trk"][data_["trk_bestFit"]]
        
    return data_

# Events at the end of the digitisation window get messed up
def CutOnStartTime(data_): 
    print(f"\n---> Cutting on start time")
    startTimeCondition = ak.all(data_["crv"]["crvcoincs.timeStart"] <= 99500, axis=1)
    return data_[startTimeCondition]

# ------------------------------------------------
#                     Trigger 
# ------------------------------------------------   
 
def Trigger(data_, passes_ = []): # "pass1", "pass2", "pass3"]):

    print(f"\n---> Triggering")

    # Enforce basic trigger condition
    triggerCondition = (
        ak.any(data_["is_coincidence"] & (data_["crv"]["crvcoincs.sectorType"] == 2), axis=1) &
        ak.any(data_["is_coincidence"] & (data_["crv"]["crvcoincs.sectorType"] == 3), axis=1)
    )

    data_["trigger"] = triggerCondition
    # data_["crv"] = ak.mask(data_["crv"], triggerCondition) 

    # Enforce pass one
    if "pass1" in passes_:
        data_ = PassOneFilter(data_)
    # if "pass2" in passes_:
    #     data_ = PassTwoFilter(data_)
    # Enforce pass three
    if "pass3" in passes_:
        data_ = PassThreeFilter(data_)

    print("...Done!")

    # Mask data

    # TODO: wouldn't it be better if we also stored the non-triggers?
    return data_[triggerCondition]

# Needs to come after the initial trigger 
def SuccessfulTriggers(data_, success):

    successStr = ""
    if success: successStr += "successful"
    else: successStr += "unsuccessful"

    print(f"\n---> Getting {successStr} triggers")

    # Ensure at least one passing coincidence in each trigger sector
    successCondition = (
        ak.any(data_["is_coincidence"] & (data_["crv"]["crvcoincs.sectorType"] == 1), axis=1) 
    )

    print("...Done!")

    if success: return data_[successCondition] # successful triggers
    else: return data_[~successCondition] # unsuccessful triggers

# ------------------------------------------------
#                     Input 
# ------------------------------------------------ 

def GetData(finName):

    # 1. Coincidence masks should ONLY be applied to coincidences.
    # 2. Global track masks should be applied to tracks and track fits.
    # 3. Local track masks should be applied to track fits.
    # ... We need to split the data array up into these parts.

    # Get data
    # finName = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00089/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000015.root"
    treeName = "TrkAnaExt/trkana"

    evtBranchNames = [ 
                    # Event info
                    "evtinfo.run"
                    , "evtinfo.subrun"
                    , "evtinfo.event"
    ]

    crvBranchNames = [
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

    trkBranchNames = [
                    # Tracks
                    "kl.status"
                    , "kl.nactive"
                    , "kl.nhits"
                    , "kl.nplanes"
                    , "kl.nnullambig"
                    , "kl.ndof"
                    , "kl.fitcon"
    ]

    trkFitBranchNames = [
                    # Track fits (vector of vector like objects)
                    "klfit"
                    , "klkl"
    ] 

    # Open tree
    with uproot.open(finName+":TrkAnaExt/trkana") as tree: 
        # Seperate event info, coincidnces, tracks, and track fits. 
        # This way we can apply masks independently. 
        evtData_ = tree.arrays(evtBranchNames) 
        crvData_ = tree.arrays(crvBranchNames) 
        trkData_ = tree.arrays(trkBranchNames)
        trkFitData_ = tree.arrays(trkFitBranchNames)

    # Zip them together 
    return ak.zip({"evt" : evtData_, "crv" : crvData_, "trk" : trkData_, "trkfit" : trkFitData_}) 

# ------------------------------------------------
#                     Output 
# ------------------------------------------------ 

def WriteFailuresToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions):

    # Define the output file path
    foutName = f"../Txt/{reproc}/failures_ntuple/{doutTag}/failures_ntuple_{foutTag}.csv" 

    print(f"\n---> Writing failures to:\n{foutName}") 

    with open(foutName, "w") as fout:
        # Write the header
        header = "\t".join(ut.branchNamesTrkAna_) + "\n"
        fout.write(header)

        for event in failures_:
            data = "\t".join(str(event[name]) for name in ut.branchNamesTrkAna_) + "\n"
            fout.write(data)

    return

def WriteFailureInfoToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions, verbose):

    # Define the output file path
    foutNameConcise = f"../Txt/{reproc}/failures_concise/{doutTag}/failures_concise_{foutTag}.csv" 
    foutNameVerbose = f"../Txt/{reproc}/failures_verbose/{doutTag}/failures_verbose_{foutTag}.csv" 

    print(f"\n---> Writing failure info to:\n{foutNameConcise}") #\n{foutNameVerbose}")
    if verbose: print(foutNameVerbose)

    # Concise form
    with open(foutNameConcise, "w") as fout:
        # Write the header
        fout.write("evtinfo.run, evtinfo.subrun, evtinfo.event\n")
        # Write the events
        for event in failures_:
            fout.write(
                f"{event['evtinfo.run']}, {event['evtinfo.subrun']}, {event['evtinfo.event']}\n"
            )

    # Verbose form
    if verbose: 
        with open(foutNameVerbose, "w") as fout:
            # Write the header
            # fout.write(f"evtinfo.runid,evtinfo.subrunid,evtinfo.eventid,is_coincidence,PE_condition,layer_condition,angle_condition,{ut.coincsBranchName}.nLayers,{ut.coincsBranchName}.angle,{ut.coincsBranchName}.sectorType,{ut.coincsBranchName}.pos.fCoordinates.fX,{ut.coincsBranchName}.pos.fCoordinates.fY,{ut.coincsBranchName}.pos.fCoordinates.fZ,{ut.coincsBranchName}.timeStart,{ut.coincsBranchName}.timeEnd,{ut.coincsBranchName}.time,{ut.coincsBranchName}.PEs,{ut.coincsBranchName}mc.valid,{ut.coincsBranchName}mc.pdgId\n")
            # Write the events
            for event in failures_:
                fout.write(
                    PrintEvent(event, coincidenceConditions)+"\n" #f"{event['evtinfo.runid']},{event['evtinfo.subrunid']},{event['evtinfo.eventid']},{event['is_coincidence']},{event['PE_condition']},{event['layer_condition']},{event['angle_condition']},{event[f'{ut.coincsBranchName}.nLayers']},{event[f'{ut.coincsBranchName}.angle']},{event[f'{ut.coincsBranchName}.sectorType']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fX']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fY']},{event[f'{ut.coincsBranchName}.pos.fCoordinates.fZ']},{event[f'{ut.coincsBranchName}.timeStart']},{event[f'{ut.coincsBranchName}.timeEnd']},{event[f'{ut.coincsBranchName}.time']},{event[f'{ut.coincsBranchName}.PEs']},{event[f'{ut.coincsBranchName}.nHits']},{event[f'{ut.coincsBranchName}mc.valid']},{event[f'{ut.coincsBranchName}mc.pdgId']}\n"
                )

    return

def WriteResultsToFile(data_, successes_, failures_, doutTag, foutTag, reproc):

    foutName = f"../Txt/{reproc}/results/{doutTag}/results_{foutTag}.csv"
    print(f"\n---> Writing results to {foutName}")
    tot = len(data_)
    efficiency = len(successes_) / tot * 100
    inefficiency = len(failures_) / tot * 100

    outputStr = f"""
    ****************************************************
    Number of failures: {len(failures_)}
    Number of successes: {len(successes_)}
    Efficiency: {len(successes_)}/{tot} = {efficiency:.2f}%
    Inefficiency: {len(failures_)}/{tot} = {inefficiency:.2f}%
    ****************************************************
    """

    with open(foutName, "w") as fout:
        fout.write("Total, Successes, Failures, Efficiency [%], Inefficiency [%]\n")
        fout.write(f"{tot}, {len(successes_)}, {len(failures_)}, {efficiency}, {inefficiency}\n")
        # fout.write(outputStr)

    print(outputStr)

    return

# ------------------------------------------------
#                       Run
# ------------------------------------------------

def Run(finName, dataset, coincidenceConditions, sanityPlots, verbose):

    # Get data as a set of awkward arrays
    data_ = GetData(finName)

    # Output dir/file tag 
    # doutTag = finName.split('.')[-2] 
    # foutTag = particle + "_" + coincidenceConditions + "_" + coincidenceFilter # + "_" +finName.split('.')[-2] 

    # Sanity plots 
    # if (sanityPlots): SanityPlots(data_, dataset, foutTag)

    # Filter particles
    # data_ = FilterParticles(data_, particle)

    # Apply start time cut
    # data_ = CutOnStartTime(data_)

    # Find coincidences
    data_ = FindCoincidences(data_, coincidenceConditions)

    # Useful debugging tool
    # PrintNEvents(data_, 100, coincidenceConditions)

    # Filter dataset 
    # data_ = FilterCoincidences(data_, coincidenceFilter)

    # Trigger
    data_ = Trigger(data_)


    PrintNEvents(data_, 10, True)

    print(data_["trigger"])

    return

    # Successful and unsuccessful triggers
    successes_ = SuccessfulTriggers(data_, success=True)
    failures_ = SuccessfulTriggers(data_, success=False)

    # PrintNEvents(successes_, 10, coincidenceConditions)

    #return 

    # Write failures to file
    WriteFailuresToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions) # write ntuple to table
    WriteFailureInfoToFile(failures_, doutTag, foutTag, reproc, coincidenceConditions, verbose) #  

    # Write results to file
    WriteResultsToFile(data_, successes_, failures_, doutTag, foutTag, reproc) 

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    # Take command-line arguments
    finName = sys.argv[1] if len(sys.argv) > 1 else "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00089/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000015.root" #  "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatTriggered/MDC2020ae_best_v1_3/root/c4/15/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000231.root" # "/pnfs/mu2e/tape/usr-nts/nts/sgrant/CosmicCRYExtractedCatDigiTrk/MDC2020z2_best_v1_1/root/40/73/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # /pnfs/mu2e/scratch/users/sgrant/workflow/CosmicCRYExtractedTrk.MDC2020z2_best_v1_1/outstage/67605881/00/00000/nts.sgrant.CosmicCRYExtractedCatDigiTrk.MDC2020z2_best_v1_1.001205_00000000.root" # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    dataset = sys.argv[2] if len(sys.argv) > 2 else "MDC2020ae" # "original"
    # particle = sys.argv[2] if len(sys.argv) > 2 else "all"
    coincidenceConditions = sys.argv[3] if len(sys.argv) > 3 else "10PEs2Layers" # "ana1" # "default"
    
    # coincidenceFilter = sys.argv[5] if len(sys.argv) > 5 else "one_coincidence_per_trigger_sector"
    sanityPlots = bool(sys.argv[4]) if len(sys.argv) > 4 else False
    verbose = bool(sys.argv[5]) if len(sys.argv) > 5 else False # False

    print("\n---> Running with inputs:\n")
    print("\tfinName:", finName)
    print("\dataset:", dataset)
    print("\tcoincidenceConditions:", coincidenceConditions)
    # print("\tcoincidenceFilter:", coincidenceFilter)
    print("\tsanityPlots:", sanityPlots)
    print("\tverbose:", verbose, "\n")

    Run(finName=finName, dataset=dataset, coincidenceConditions=coincidenceConditions, sanityPlots=sanityPlots, verbose=verbose) 

    return

if __name__ == "__main__":
    main()
