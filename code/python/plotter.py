# Somewhat experimental, switching to Python only analysis of trees using uproot and pandas
import uproot
import pandas as pd

def read(finName):

    print("---> Reading TTree...")

    # Read ROOT file
    fin = uproot.open(finName)

    print("---> Got input file ", finName, ", ", (fin))

    # Access a specific TTree within the ROOT file
    tree = fin["TrkAnaExt/trkana"]

    print("---> Got tree ", str(tree))

    # Convert the TTree to a pandas DataFrame
    branches = tree.arrays(library="np")

    # For each branch, we convert the corresponding array data into a pandas Series 
    # Store these Series in a dictionary called series, then concatenate the Series into a DataFrame
    # This is needed to handle columns of varying length
    series = {}
    for branch_name, branch_data in branches.items():
        series[branch_name] = pd.Series(branch_data)

    # Concatenate the Series into a DataFrame
    df = pd.DataFrame(series)

    print("---> Loaded events into DataFrame:\n", df)

    # Close the file
    fin.close()

    return df

def main():

    # Get data 
    finName = "../../data/nts.owner.trkana-reco.version.sequencer.root"
    df = read(finName)


if __name__ == "__main__":
    main()



# import uproot
# import pandas as pd

# # Read a ROOT file
# file = uproot.open("../../data/nts.owner.trkana-reco.version.sequencer.root")

# print("---> Opened file " + str(file))

# # Access a specific TTree within the ROOT file
# tree = file["TrkAnaExt/trkana"]

# print("---> Got tree " + str(tree))

# # Convert the TTree to a pandas DataFrame
# branches = tree.arrays(library="np")


# # For each branch, we convert the corresponding array data into a pandas Series 
# # Store these Series in a dictionary called series, then concatenate the Series into a DataFrame
# # This is needed to handle columns of varying length
# # Otherwise you could just do 
# # arrays = tree.arrays()
# # df = pd.DataFrame(arrays)

# series = {}
# for branch_name, branch_data in branches.items():
#     series[branch_name] = pd.Series(branch_data)

# # Concatenate the Series into a DataFrame
# df = pd.DataFrame(series)

# # 

# # Close the file when you're done
# file.close()

# # import argparse
# # import re, glob, os, sys, re, subprocess
# # import pandas as pd
# # from root_pandas import read_root # see https://github.com/scikit-hep/root_pandas 
# # import h5py # https://github.com/h5py/h5py

# # arg_parser = argparse.ArgumentParser()
# # arg_parser.add_argument("--fin", type=str, default="../../data/nts.owner.trkana-reco.version.sequencer.root") # input file

# # def main():

# # 	run()

# #     # if(args.sim_skim==True):
# #     #     sim_skim()

# #     # #default is to skim many Tree
# #     # if(args.add==False):
# #     #     skim()

# #     # #add HDF5s into one based on label (name)
# #     # if(args.add==True):
# #     #     add()

# #     # #add HDF5s into one based on label (name)
# #     # if(args.add_R1==True):
# #     #     add_R1()

# #    # tree->SetBranchAddress("klmcsim.time", &klmcsim_time);
# #    # tree->SetBranchAddress("klmcsim.pos.fCoordinates.fX", &klmcsim_pos_fCoordinates_fX);
# #    # tree->SetBranchAddress("klmcsim.pos.fCoordinates.fY", &klmcsim_pos_fCoordinates_fY);
# #    # tree->SetBranchAddress("klmcsim.pos.fCoordinates.fZ", &klmcsim_pos_fCoordinates_fZ);

# #    # tree->SetBranchAddress("crvinfomcplane._time", &crvinfomcplane_time); 
# #    # tree->SetBranchAddress("crvinfomcplane._x", &crvinfomcplane_x); 
# #    # tree->SetBranchAddress("crvinfomcplane._y", &crvinfomcplane_y); 
# #    # tree->SetBranchAddress("crvinfomcplane._z", &crvinfomcplane_z); 

# # def run():

# #     # Only open required columns 
# #     for i_file, file in enumerate(sorted(os.listdir(args.trees))):
# #         for i_key, key in enumerate(keys):
# #             print("Opening", key, "data in", args.trees+"/"+file)
# #             data_all = read_root(args.trees+"/"+file, key, columns=["klmcsim.time"]) #, "trackT0", "trackMomentum", "trackMomentumY", "trackPValue"])
# #             # data_all['trackT0']=data_all['trackT0']*1e-3   # ns -> us
# #             total_tv[i_key]+=data_all.shape[0] # add to total from each file 
# #             print("Total of", data_all.shape[0], "entries")
           
# #             #define the time and energy cuts (otherwise the wiggle looks wobbly!)
# #             # time_cut= (data_all['trackT0'] > args.t_cut)  # us, define a time_cut with time > 30 us 
# #             # mom_cut = (data_all['trackMomentum'] > args.p_cut) # MeV 
# #             #Apply cuts in one go! 
# #             # data = data_all[time_cut & mom_cut]
# #             # total_tv_skim[i_key]+=data.shape[0] # add to total from each file 
# #             # print("Total of", data.shape[0], "entries after energy and momentum cuts")

# #             #save the skimmed dataframe in a compressed HDF5 format
# #             # here we are appending to the file over tracks and then vertices
# #             print("Saving compressed data...")
# #             data=data.reset_index() # reset index from 0 
# #             cols_to_keep = ["klmcsim.time"] #, "trackPValue"] # only write for time and station 
# #             # cols_to_keep = ["station", "trackT0", "trackMomentum", "trackMomentumY"] # only write for time and station 
# #             data[cols_to_keep].to_hdf(args.df+"_"+str(i_file)+".h5", key=key, mode='a', complevel=9, complib="zlib", format="fixed")
# #             print("Dataframe saved to disk", args.df, "\n")

# #     # print("Grand total of (M)", total_tv[0]/1e6, "tracks,", total_tv[1]/1e6, "vertices")
# #     # print("After the cuts (M)", total_tv_skim[0]/1e6, "tracks,", total_tv_skim[1]/1e6, "vertices")

# # if __name__ == "__main__":
# #     main()
