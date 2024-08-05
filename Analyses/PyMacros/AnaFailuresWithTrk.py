# Analyse failures with the tracker. 
# Just take the concatenated failure file for testing.  

# External libraries
import pandas as pd
import uproot
import awkward as ak

# Internal libraries
import Utils as ut
import CompCRVUtils as trkCoinUt

# Helper functions

def FindLineInFile(fileName, searchString):
    try:
        with open(fileName, "r") as file:
            for line in file:
                if searchString in line:
                    print(f"\n---> Found file\n:{line.strip()}")
                    return line.strip()
        print(f"---> String '{searchString}' not found in file.")
        return None
    except FileNotFoundError:
        print(f"---> File {fileName} does not exist.")
        return None

    
def Run(failuresFileName, excludedFailuresFileName, nonExcludedFailuresFileName): 
    
    # Read failures DataFrame
    failures = pd.read_csv(failuresFileName)
    # Set up new empty DataFrames with the same headers as "failures"
    excludedFailures = pd.DataFrame(columns=failures.columns)
    nonExcludedFailures = pd.DataFrame(columns=failures.columns)

    # MDC2020ae file list
    # fileList = "../Txt/FileLists/MDC2020ae.txt" # pnfs
    fileList = "../Txt/FileLists/MDC2020aeOnExpData.txt" # /exp/data

    # Iterate through tags and files
    tagColumnName = "Start run/subrun"
    
    for tag in failures[tagColumnName]: # "Start run/subrun"]:

        # This guy has multiple, so it's good for testing
        # if tag!="001205_00000064": 
        #     continue
        
        print(f"\n---> Analysing {tag}")
        
        # Get file
        sourceFileName = FindLineInFile(fileList, tag)

        # Input file
        # sourceFileName = "/exp/mu2e/data/users/sgrant/CRVSim/CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.000/11946817/00/00033/nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000080.root"

        arrays = ak.Array([])
        with uproot.open(sourceFileName+":TrkAnaExt/trkana") as tree: 
            # Do some fancy read from /pnfs! 
            # For now use /data
            # For some reason this works with Jupyter, but not in a script? 
            # OK it works now? Is it caching the data or something? What's happening?
            # Don't understand this. Need more testing!
            arrays = tree.arrays(["evtinfo.", "crvcoincs", "crvcoincsmc", "kl", "klfit", "klkl"]) 
            
        # print(arrays)
        # Apply tracker cuts
        trkCoinUt.MarkCuts(arrays)
        withCRV = trkCoinUt.ApplyCuts(arrays, ["goodCRV", "goodTrk", "CRV1", "KLCRV1", "bestFit"])

        print(withCRV)

    #     break

    # return
        
        # Get failure ID(s)
        thisFailure = failures[failures[tagColumnName] == tag] # could be multiple!

        # There could be multiple failures, so iterate through the rows. 
        for index, failure in thisFailure.iterrows():
            run = failure[" evtinfo.run"]
            subrun = failure[" evtinfo.subrun"]
            event = failure[" evtinfo.event"]

            print(f"\n---> This failure: {run}:{subrun}:{event}.") 
            
            mask = ( (withCRV["evtinfo."]["evtinfo.run"] == run) 
                    & (withCRV["evtinfo."]["evtinfo.subrun"] == subrun) 
                    & (withCRV["evtinfo."]["evtinfo.event"] == event) )

            failureArray = withCRV[mask]

            # Check if the array is empty 
            if (len(ak.flatten(failureArray, axis=None)) == 0): 
                print(f"Failure excluded! {run}:{subrun}:{event}")
                excludedFailures.loc[len(excludedFailures)] = [tag, run, subrun, event]
                # excludedFailures = excludedFailures.append({
                #     'tag': tag,
                #     'run': run,
                #     'subrun': subrun,
                #     'event': event
                # }, ignore_index=True)
            else:
                print("Not excluded")
                nonExcludedFailures.loc[len(nonExcludedFailures)] = [tag, run, subrun, event]
                

            # print(f"---> Failure array:\n{failureArray}.")

            # Count number of events in the failure array, after enforcing tracker cuts. 
            
            
            # print(f"---> Failure array:\n{(len(ak.flatten(failureArray[, axis=None)))}.")
            # print(f"---> Failure array:\n{ak.num(failureArray["evtinfo."]["evtinfo.event"], axis=None)}") #, axis=None)))}.")
            # print(f"---> Failure array:\n{ak.num(ak.flatten(failureArray, axis=None))}.")
            
            # Apply tracker cuts, now or later? Probably do that first and see if this event passes? 
        
        # break

    # Write 
    excludedFailures.to_csv(excludedFailuresFileName, index=False)
    nonExcludedFailures.to_csv(nonExcludedFailuresFileName, index=False)
    print(f"---> Written:\n{excludedFailuresFileName}\n{nonExcludedFailuresFileName}")
    
    # Iterate through files 
    
    return    

def main():

    failuresFileName = "../Txt/MDC2020ae/JulyCRVMeeting/concatenated/failures_concise_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv"
    
    excludedFailuresFileName = "../Txt/MDC2020ae/AnaTrkFailuresWithTrk/trk_excluded_failures_concise_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv"
    nonExcludedFailuresFileName = "../Txt/MDC2020ae/AnaTrkFailuresWithTrk/trk_non_excluded_failures_concise_all_10PEs2Layers_one_coincidence_per_trigger_sector.csv"
    
    Run(failuresFileName=failuresFileName, excludedFailuresFileName=excludedFailuresFileName, nonExcludedFailuresFileName=nonExcludedFailuresFileName) 

    return

if __name__ == "__main__":
    main()