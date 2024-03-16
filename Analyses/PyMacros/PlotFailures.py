'''
Samuel Grant 2024
Make plot for failures
''' 

import numpy as np
import pandas as pd
import ast

import Utils as ut

# ------------------------------------------------
#                Helper functions 
# ------------------------------------------------

# def FlattenColumn(col):
#     # Subarrays in DataFrame as interpreted as strings, evaluate them literally 
#     col = col.apply(lambda x: ast.literal_eval(x))
#     # Convert to numpy array
#     col = np.array(col)
#     # Flatten the array
#     col = np.concatenate([np.array(sublist) for sublist in col])
#     # More nonsense
#     col = ast.literal_eval(col) 
#     return col

# def FlattenColumn(col):
#     # Subarrays in DataFrame as interpreted as strings, evaluate them literally 
#     col = col.apply(lambda x: ast.literal_eval(x))
#     # Convert to numpy array
#     col = np.array(col)
#     # Flatten the array
#     col = np.concatenate([np.array(sublist) for sublist in col])
#     return col

# def FlattenColumn(col):
#     # Subarrays in DataFrame as interpreted as strings, evaluate them literally 
#     col = col.apply(lambda x: ast.literal_eval(x))
#     # Convert to numpy array
#     col = np.concatenate([np.array(sublist) for sublist in col])
#     return col

# Columns containing arrays are intepreted as strings 
# Convert them into 1D lists

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

# ------------------------------------------------
#                      Run
# ------------------------------------------------

def Run(): 

    data_ = { 
        "All" : pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
        ,"Muons" : pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_muons_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
        ,"Non-muons" : pd.read_csv("../Txt/reprocessed/concatenated/failures_ntuple_non_muons_10PEs2Layers_one_coincidence_per_trigger_sector.csv", sep="\t")
    }

    timeStart = FlattenColumn(data_["All"]["crvcoincs.timeStart"]) # .apply(lambda x: ast.literal_eval(x))

    ut.Plot1D(np.array(timeStart), nbins=len(timeStart)+1000, xmin=0, xmax=np.max(timeStart)+1000, xlabel="Start time [ns]", ylabel="Failures", fout="../Images/reprocessed/Failures/h1_timeStart.png", stats=True)

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    Run() # finName=finName, particle=particle, coincidenceConditions=coincidenceConditions, reproc=reproc, coincidenceFilter=coincidenceFilter) 

    return

if __name__ == "__main__":
    main()
