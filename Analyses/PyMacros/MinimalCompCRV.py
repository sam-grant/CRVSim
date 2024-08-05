''' 
Minimal example demonstrating the application of tracker-CRV cuts. 
Samuel Grant 2024
'''

# Internal libraries 
import awkward as ak
import uproot

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
        eventStr += f"crv_mask: {event['crv_mask']}\n" 
        eventStr += f"trk_mask: {event['trk_mask']}\n" 
        eventStr += f"trkfit_mask: {event['trkfit_mask']}\n" 
        eventStr += f"-------------------------------------------------------------------------------------\n"
    
    return eventStr

def PrintNEvents(data_, nEvents=10, showMasks=False):
     # Iterate event-by-event
    for i, event in enumerate(data_):
        i += 1
        # if len(cuts) == 0: print(PrintEvent(event, showMasks))
        # else:
        print(PrintEvent(event, showMasks))

        if i >= nEvents: 
            return

print("\nSome rules:")
print("1) Coincidence masks should ONLY be applied to coincidences.")
print("2) Global track masks should be applied to tracks and track fits.")
print("3) Local track masks should be applied to track fits.")
print("... We need to split the data array up into these parts.")

# Get data
finName = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00089/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000015.root"
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
    evtArrays = tree.arrays(evtBranchNames) 
    crvArrays = tree.arrays(crvBranchNames) 
    trkArrays = tree.arrays(trkBranchNames)
    trkFitArrays = tree.arrays(trkFitBranchNames)

# Zip them together 
arrays = ak.zip({"evt" : evtArrays, "crv" : crvArrays, "trk" : trkArrays, "trkfit" : trkFitArrays}) 

# Printouts
print(f"\n---> Got tree {treeName} from file {finName}.")
print(f"\n---> Tree fields {arrays.fields}")

print("\n---> First few events:")
PrintNEvents(arrays, nEvents=5)

print("---> One event look in klfit (local track fit level):")

eventMask = arrays["evt"]["evtinfo.event"] == 2221

# Iterate through each event and print the subarrays in `klfit`
for event_index, event in enumerate(arrays[eventMask]):
    for subarray_index, subarray in enumerate(event["trkfit"]["klfit"]):
        for item_index, item in enumerate(subarray):
            print(f"{item_index}: {item}")

# Now, some enforce masks on each level.
# Append the mask conditions to the array so that we can see the structure

arrays["crv_mask"] = arrays["crv"]["crvcoincs.sectorType"] == 1 
       
arrays["trk_mask"] = ( (arrays["trk"]["kl.ndof"] == 20) ) # 10) ) 
                        #  & (arrays["trk"]["kl.fitcon"] > 0.1)
                        #  & ((arrays["trk"]["kl.nactive"]/arrays["trk"]["kl.nhits"]) > 0.99)
                        #  & (arrays["trk"]["kl.nplanes"] >= 4)
                        #  & ((arrays["trk"]["kl.nnullambig"]/arrays["trk"]["kl.nhits"]) < 0.2) )

arrays["trkfit_mask"] = ( (arrays["trkfit"]["klfit"]["sid"] == 200) 
                           & (arrays["trkfit"]["klfit"]["sindex"] == 1) )

PrintNEvents(arrays, nEvents=5, showMasks=True)

# Apply masks 
print("---> Applying masks...")
print(f"Number of events before masks {len(ak.flatten(arrays['evt']['evtinfo.event'], axis=None))}")

arrays["crv"] = arrays["crv"][arrays["crv_mask"]]
arrays["trk"] = arrays["trk"][arrays["trk_mask"]]
arrays["trkfit"] = arrays["trkfit"][arrays["trkfit_mask"]]

PrintNEvents(arrays, nEvents=5, showMasks=True)

print("---> Masks applied.")

print("---> All sectors should be 1, all ndof should be 20, all sid should be 200.")
print(ak.flatten(arrays["crv"]["crvcoincs.sectorType"], axis=None))
print(ak.flatten(arrays["trk"]["kl.ndof"], axis=None))
print(ak.flatten(arrays["trkfit"]["klfit"]["sid"], axis=None))