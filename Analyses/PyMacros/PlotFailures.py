'''
Samuel Grant 2024
Make plot for failures

This script is an absolute mess. 

''' 

import numpy as np
import pandas as pd
import awkward as ak
import ast

import Utils as ut

# ------------------------------------------------
#           Read text based ntuple as awk 
# ------------------------------------------------

def TextToAwkward(finName):

    print("\n---> Converting text file to awkward array")

    # Initialize an empty dictionary to hold the arrays
    array_ = {}

    # Initialize empty lists for each branch
    data_ = [[] for _ in range(len(ut.branchNamesTrkAna_))]

    # Count the total number of lines in the file
    # with open(finName, 'r') as fin:
    #     totLines = sum(1 for line in fin)

    # Open the text file and loop through each line
    with open(finName, 'r') as fin:
        lineCount = 0 
        for i, line in enumerate(fin):
            # Print status
            # lineCount += 1
            # percentageDone = lineCount / totLines * 100
            # print(f"\rProcessing: {percentageDone:.2f}%", end='', flush=True)
            # Skip header
            if i == 0: continue
            # Skip empty rows
            if len(line) <= 1: continue
            # Split the line into values
            values = line.strip().split('\t')
            
            # Append values to corresponding data lists
            for i, branch in enumerate(ut.branchNamesTrkAna_):
                # print(i, branch) 
                # val = values[i].apply(ast.literal_eval)
                # val_ = ast.literal_eval(values[i])
                # if isinstance(val_, list):
                #     val_ = [ast.literal_eval(v) for v in val_]
                # [val for sublist in col.apply(ast.literal_eval) for val in sublist]
                data_[i].append(values[i])

    # Convert data lists to awkward arrays and store them in the arrays dictionary
    for i, branchName in enumerate(ut.branchNamesTrkAna_):
        array_[branchName] = ak.Array(data_[i])

    print("Done!")

    return array_

def RecastArray(array_):

    # for array_ in data_: 
    #     for subarray_ in array_: 
    #         subarray_ = ast.literal_eval(subarray_)
    #         if isinstance(subarray_, list):
    #             for subsubarray_ in subarray_:
    #                 subsubarray_ = ast.literal_eval(subsubarray_)
                # [ast.literal_eval(val) for val in col]

    return ast.literal_eval(array_)
# ------------------------------------------------
#                Helper functions 
# ------------------------------------------------

# Readable
# def FlattenColumn(col):
#     col = col.tolist()
#     flatCol = [] 
#     for sublist in col:
#         sublist = ast.literal_eval(sublist)
#         for val in sublist:
#             flatCol.append(val)
#     return flatCol


# Fancy
def FlattenColumn(col):
    return [val for sublist in col.apply(ast.literal_eval) for val in sublist]
    # return [val for sublist in ast.literal_eval(col) for val in sublist]

# ------------------------------------------------
#                Analysis functions 
# ------------------------------------------------

def MakeStartTimePlot():

    data1_ = pd.read_csv("../Txt/reprocessed/concatenated/noStartTimeCut/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
    data2_ = pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
    
    timeStart1_ = FlattenColumn(data1_["crvcoincs.timeStart"]) 
    timeStart2_ = FlattenColumn(data2_["crvcoincs.timeStart"]) 

    ut.Plot1DOverlay(hists_ = [np.array(timeStart1_), np.array(timeStart2_)], nbins=len(timeStart1_)+1000, xmin=0, xmax=np.max(timeStart1_)+1000, xlabel="Start time [ns]", ylabel="Failures", label_ = ["No cut", "99500 ns cut"], fout="../Images/reprocessed/Failures/h1_timeStart_overlay.png") 

    return

def MakePDGIDPlot(): 

    # Plot particle IDs of failures
    layers_ = [2,3]
    # df_ = pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
    ak_ = [TextToAwkward(f"../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs{layer}Layers_one_coincidence_per_trigger_sector.csv") for layer in layers_]
        # TextToAwkward("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv") 
        # TextToAwkward("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv") 

    PDGIDsPerLayer_ = [] 

    for layer, data_ in zip(layers_, ak_):

        # TODO: wrap this in a function and generalise for the whole array
        # Alternatively, stop using text files for this!!! 
        PDGIDs_ = [ast.literal_eval(event) for event in data_["crvcoincsmc.pdgId"]]

        # Get the unique particle in the event
        PDGIDs_ = [set(sublist) for sublist in PDGIDs_]
        # Get all particles 
        # PDGIDs_ = [sublist for sublist in PDGIDs_]

        # Flatten 
        PDGIDs_ = [particle for sublist in PDGIDs_ for particle in sublist]

        # PDGIDs_ = [array for subarray in RecastArray(ak_["crvcoincsmc.pdgId"])
        ut.BarChart(data_=PDGIDs_, label_dict=ut.latexParticleDict, title=f"{layer}/4 layers", ylabel="Unique particles / failure", fout=f"../Images/reprocessed/Failures/bar_all_failures_uniquePDGid_{layer}layers.png")
        # ut.BarChart(data_dict=PDGIDs_[PDGIDs_ != ], title="Failures coincidences", ylabel="Failures / sector", fout="../Images/reprocessed/Failures/bar_all_PDDIDs.png")

        # ut.BarChartOverlay(data_=[pdgid_[sectors_ == 3], pdgid_[sectors_ == 1], pdgid_[sectors_ == 2]], label_dict=label_dict, ylabel="Coincidences [%]", fout=f"../Images/{reproc}/Sanity/bar_overlay_pdgid_{foutTag}.png", percentage=True, label_= ["Top", "Middle", "Bottom"])

        PDGIDsPerLayer_.append(PDGIDs_)

    ut.BarChartOverlay(data_=[PDGIDsPerLayer_[0], PDGIDsPerLayer_[1]], label_dict=ut.latexParticleDict, ylabel="Unique particle species / failure", fout=f"../Images/reprocessed/Failures/bar_overlay_all_failures_uniquePDGid.png", percentage=False, label_= ["2/3 layers", "3/4 layers"])


    return

def MakeAngularDistributionPlot():

    data_ = pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")

    slope_ = FlattenColumn(data_["crvcoincs.angle"]) 
    slope_ = np.array(slope_)

    # Suppress zeros
    slope_ = slope_[slope_ != 0]

    print(slope_)

    ut.Plot1D(slope_, nbins=1000, xmin=np.min(slope_), xmax=np.max(slope_), fout=f"../Images/reprocessed/Failures/h1_failures_slope.png")

    return

# ------------------------------------------------
#                      Run
# ------------------------------------------------

def Run(): 

    # MakeStartTimePlot()
    # MakePDGIDPlot()
    MakeAngularDistributionPlot()


    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    Run() 

    return

if __name__ == "__main__":
    main()
