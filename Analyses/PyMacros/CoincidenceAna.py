
import h5py
import os
import numpy as np

import Utils as ut

# ------------------------------------------------
#                      Read
# ------------------------------------------------

def ReadCoincidences(finName, debug):

    # Read the h5 file
    print(f"\n---> Reading in coincidences from {finName}")

    # True/false coincidences dicts
    trueCoincidences_ = {}
    falseCoincidences_ = {}

    with h5py.File(finName, "r") as file:

        # Read data 
        for key in file["trueCoincidences"]:
            trueCoincidences_[key] = list(file[f"trueCoincidences/{key}"])

        # Read data into dict2
        for key in file["falseCoincidences"]:
            falseCoincidences_[key] = list(file[f"falseCoincidences/{key}"])

    print("...Done!")

    return trueCoincidences_, falseCoincidences_

# ------------------------------------------------
#                      Analyse
# ------------------------------------------------

def AnalyseCoincidences(trueCoincidences_, falseCoincidences_, filterCondition, debug):

    print("\n---> Analysing coincidences.\n")

    # Count true and false coincidences in all sectors
    trueCounts_ = { "Sector 1" : len(trueCoincidences_["S1"]), "Sector 2" : len(trueCoincidences_["S2"]),  "Sector 3" : len(trueCoincidences_["S3"]) }
    falseCounts_ = { "Sector 1" : len(falseCoincidences_["S1"]), "Sector 2" : len(falseCoincidences_["S2"]),  "Sector 3" : len(falseCoincidences_["S3"]) }
    # fractionPassing_ = { "Sector 1" : 1-(falseCounts_["Sector 1"]/trueCounts_["Sector 1"]), "Sector 2" : 1-(falseCounts_["Sector 2"]/trueCounts_["Sector 2"]),  "Sector 3" : 1-(falseCounts_["Sector 3"]/trueCounts_["Sector 3"]) }

    print("Passing coincidences:\n", trueCounts_)
    print("Failing coincidences:\n", falseCounts_)

    ut.BarChart(data_dict=trueCounts_, title="Passing coincidences", ylabel="Counts / sector", fout="../Images/TrueAndFalseCoincidences/"+filterCondition+"/bar_true_coin_counts_per_sector.png")
    ut.BarChart(data_dict=falseCounts_, title="Failing coincidences", ylabel="Counts / sector", fout="../Images/TrueAndFalseCoincidences/"+filterCondition+"/bar_false_coin_counts_per_sector.png")

    # ut.BarChart(data_dict=fractionPassing_, ylabel="Fraction passing / sector", fout="../Images/TrueAndFalseCoincidences/bar_frac_pass_coin_counts_per_sector.png")
    # ut.BarChartOverlay(data_dict_=[trueCounts_, falseCounts_], labels_=["Passing coincidences", "Failing coincidences"], ylabel="Counts / sector", fout="../Images/TrueAndFalseCoincidences/bar_true_false_coin_counts_per_sector.png")

    # Fraction of hits in S1 to triggers, event-by-event
    # print(trueCoincidences_)
    # print(np.min(trueCoincidences_["S1"]))

    ut.Plot1DOverlay(trueCoincidences_, nbins=int(np.max(trueCoincidences_["S1"])), xmin=0, xmax=np.max(trueCoincidences_["S1"]), xlabel="Event ID", ylabel="Passing coincidences / event", fout="../Images/TrueAndFalseCoincidences/"+filterCondition+"/h1_passing_coin_vs_eventID_per_sector.png")
    ut.Plot1DOverlay(falseCoincidences_, nbins=int(np.max(trueCoincidences_["S1"])), xmin=0, xmax=np.max(trueCoincidences_["S1"]), xlabel="Event ID", ylabel="Failing coincidences / event", fout="../Images/TrueAndFalseCoincidences/"+filterCondition+"/h1_failing_coin_vs_eventID_per_sector.png")

    # Triggers? 

    # Probably easier to do this during the event loop, right?


    return

# ------------------------------------------------
#                      Run
# ------------------------------------------------

def Run(finName, filterCondition, debug=True):

    # Read in
    trueCoincidences_, falseCoincidences_ = ReadCoincidences(finName, debug)

    # Analyse
    AnalyseCoincidences(trueCoincidences_, falseCoincidences_, filterCondition, debug)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():
  
    # # Take input file name command-line argument
    # if len(sys.argv) != 3:
    #     # print("Input and outname file names required as arguments")
    #     print("Error: filter condition and input file name required as argument, example:")
    #     print("python CoincidenceFinder.py one_coincidence_per_trigger_sector /pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka")
    #     sys.exit(1)
    
    # filterCondition = sys.argv[1]
    # finName = sys.argv[2] 

    # finName = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka" # sys.argv[1] # "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CosmicCRYExtractedTrk/MDC2020z1_best_v1_1_std_v04_01_00/tka/82/e8/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.tka"
    # filterCondition = "no_filter"
    filterCondition = "one_coincidence_per_sector"
    # filterCondition = "one_coincidence_per_trigger_sector"

    finName = "../h5/TrueAndFalseCoincidences/"+filterCondition+"/nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.001205_00000000.h5"

    Run(finName=finName, filterCondition=filterCondition, debug=False) 

    return

if __name__ == "__main__":
    main()



# How often do you get a coincidence in sectors 2 & 3 when you also have on in sector 1
# We need to handle multiple hits

# def FilterAndCountCoincidences(arr_):
#     print("\n---> Filtering and counting coincidences")

#     totEvents = len(arr_)

#     # Dictionary to count coincidences event-by-event

#     # Count coincidences event-by-event, append the event IDs to lists 
#     singles_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }
#     doubles_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }
#     triples_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }   
    
#     # Count false coincidences 
#     falses_ = { "sector_1" : [], "sector_2" : [],  "sector_3" : [], }  

#     # Iterate event-by-event
#     for i, entry in enumerate(arr_):

#         # Check if any of the arrays in the entry are empty
#         # This just slows everything down!
#         # if any(not entry[field].tolist() for field in entry.fields):
#         #     continue

#         # Number of coincidences, event ID and sectors
#         nCoin = len(entry["is_coincidence"])
#         eventID = entry['evtinfo.eventid']
#         sectors_ = entry["crvhit.sectorType"]

#         # Check for false coincidences
#         if any(entry["is_coincidence"] == False):
#             print("\n!!! False coincidence found !!!")
#             PrintEntry(entry)
#             # Store them 
#             for sector in sectors_:  
#                 falses_["sector_"+str(sector)].append(eventID)
#             # If there are any, we can cycle back later


#         # Handle multiple coincidences per trigger sector 
#         sectorCounts = Counter(sectors_)
#         # print(f"\nEvent {eventID}")
#         for sector, n in sectorCounts.items():
#             if n>1: 
#                 print(f"Sector {sector} has {n} counts")
#                 PrintEntry(entry)

#         # break
        
#         continue


#         # You've been getting away with it because most of them are true.


#         # Count single coincidences
#         if (nCoin == 1):
#             for sector in sector_: # this should always be the top sector 
#                 singles_["sector_"+str(sector)].append(eventID)
#         # Count double coincidences 
#         elif (len(entry["is_coincidence"]) == 2):
#             for sector in sector_:
#                 doubles_["sector_"+str(sector)].append(eventID)
#         # Count triple coincidences
#         elif (len(entry["is_coincidence"]) == 3):
#             for sector in sector_:
#                 triples_["sector_"+str(sector)].append(eventID)
#         else: 
#             # Is this the correct place to be doing this?

#             beyond_.append(eventID)
#             print(
#                 f"*** WARNING: CountCoincidences() ***\n"
#             )

#             PrintEntry(entry)

#             # print(
#             #     f"*** WARNING: CountCoincidences() ***\n"
#             #     f"Event ID: {eventID}\n"
#             #     f"Number of coincidences: {len(entry['is_coincidence'])}\n"
#             #     f"Coincidences: {entry['is_coincidence']}\n"
#             #     f"Entry: {entry}"
#             # )

#         progress = (i + 1) / totEvents * 100
#         print(f"Progress: {progress:.2f}%", end='\r', flush=True)

#         # if (i > 500): break

#     # print(
#     #     f"\nSingles: {singles_}\n"
#     #     f"Doubles: {doubles_}\n"
#     #     f"Triples: {triples_}\n"
#     # )

#     # Make bar charts 
#     # TODO: handle this better!
#     singlesCounts_ = { "sector_1" : len(singles_["sector_1"]), "sector_2" : len(singles_["sector_2"]),  "sector_3" : len(singles_["sector_3"]) }
#     doublesCounts_ = { "sector_1" : len(doubles_["sector_1"]), "sector_2" : len(doubles_["sector_2"]),  "sector_3" : len(doubles_["sector_3"]) }
#     triplesCounts_ = { "sector_1" : len(triples_["sector_1"]), "sector_2" : len(triples_["sector_2"]),  "sector_3" : len(triples_["sector_3"]) }

#     ut.BarChart2(data_dict=singlesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_singles.png")
#     ut.BarChart2(data_dict=doublesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_doubles.png")
#     ut.BarChart2(data_dict=triplesCounts_, ylabel="Counts / sector", fout="../Images/Coincidences/bar_triples.png")

#     print("...Done!")

#     return