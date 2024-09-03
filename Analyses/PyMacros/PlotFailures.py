'''
Samuel Grant 2024

Read failure events and make plots. 

''' 
# External libraries
import numpy as np
import pandas as pd
import awkward as ak
from concurrent.futures import ThreadPoolExecutor, as_completed
import awkward as ak
from multiprocessing import Manager

# Internal libraries
import Utils as ut
from Mu2eEAF import ReadData as rd 

# Retreive the failures from the dataset and append them to a master array
def GetFailures(fileList_, failureFileName, max_workers=381):
    
    print("\n---> Getting failures...\n")

    print(f"\n---> Using info file {failureFileName}.")
    failureInfo_ = pd.read_csv(failureFileName)

    # Collect file list
    tags_ = list(set(failureInfo_["Tag"]))
    
    # Extract the tag from the file name
    def ExtractTag(fileName):
        parts = fileName.split('.')
        if len(parts) > 1:
            return parts[-2]
        return None

    # Filter and sort files based on tags
    fileList_ = sorted(
        [file for file in fileList_ if ExtractTag(file) in tags_]
        , key=lambda file: tags_.index(ExtractTag(file))
    )

    # Bug check
    if len(fileList_) == len(tags_):
        if False: print("\n---> Collected and sorted failure file names.")
    else:
        raise Exception("\n---> len(fileList_) != len(tags_)")

    # Testing
    tags_ = tags_[:2]
    fileList_ = fileList_[:2]

    completedFiles = 0
    totalFiles = len(fileList_)

    # Array shared among threads 
    # This needs to be an awkward array...
    # Use Manager to share data between threads
    with Manager() as manager:
        # Thread-safe list to collect arrays
        globalArray_ = manager.list()

        def processFunction(tag, fileName):
            nonlocal globalArray_
    
            # Find failure events
            thisFailureInfo_ = failureInfo_[failureInfo_["Tag"] == tag]
            outputStr = ( 
                f"\n--->\n" 
                f"fileName: {fileName}\n"
                f"tag: {tag}\n"
                f"failures:\n{thisFailureInfo_}\n"
                f"---"
            )
            print(outputStr)
        
            # Read the file
            file = ut.ReadFile(fileName, quiet=True)
            
            # Get array
            localArray_ = ut.GetData(file)
            
            if False: print(f"\n---> Loaded corresponding data.\n{thisData_}")
    
            if False: print(f"\n---> Applying masks.")
        
            # Extract unique values from DataFrame
            runs_ = set(thisFailureInfo_["evtinfo.run"])
            subruns_ = set(thisFailureInfo_[" evtinfo.subrun"])
            events_ = set(thisFailureInfo_[" evtinfo.event"])
    
            # Construct masks
            runCondition = ak.any([thisData_["evt"]["evtinfo.run"] == value for value in runs_], axis=0)
            subrunCondition = ak.any([thisData_["evt"]["evtinfo.subrun"] == value for value in subruns_], axis=0)
            eventCondition = ak.any([thisData_["evt"]["evtinfo.event"] == value for value in events_], axis=0)
    
            # Apply masks
            localArray_ = localArray_[runCondition & subrunCondition & eventCondition]
    
            # Append result to the managed list
            print(f"\n---> Appending data to globalArray_. Current length: {len(globalArray_)}")
            globalArray_.append(localArray_)
            
            return
    
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Prepare a list of futures and map them to (tag, fileName) tuples
            futures_ = {executor.submit(processFunction, tag, fileName): (tag, fileName) for tag, fileName in zip(tags_, fileList_)}
            
            # Process results as they complete
            for future in as_completed(futures_):
                tag, fileName = futures_[future]  # Get the (tag, fileName) tuple associated with this future
                try:
                    future.result()  # Retrieve the result or raise an exception if occurred
                    completedFiles += 1
                    percentComplete = (completedFiles / totalFiles) * 100
                    print(f'\n---> File {fileName} with tag {tag} processed successfully! ({percentComplete:.1f}% complete)')
                except Exception as exc:
                    print(f'\n---> File {fileName} with tag {tag} generated an exception!\n{exc}')

        print("\n---> Multithreading completed!")
    
        # Collect results
        array_ = ak.Array([])
    
        if len(globalArray_) > 0:
            array_ = ak.concatenate(list(globalArray_), axis=0)
        else: 
            raise Exception("--->  len(globalArray_) == 0!")
    
        inputEventList = failureInfo_[" evtinfo.event"]
        outputEventList = ak.flatten(array_["evt"]["evtinfo.event"])
        if len(inputEventList) == len(outputEventList):
            print("---> Output array contains the correct number of events!")
            print(f"Input array: {inputEventList}")
            print(f"Output array: {outputEventList}")
        else:
            raise Exception("---> Output array contains the wrong number of events!")
        
        print("\nDone!")
        return array_


# This is good. TODO: add this to Mu2eEAF. 
# def Multithread2(processFunction, tags_, fileList_, max_workers=381):
    
#     print("\n---> Starting multithreading...\n")

#     completedFiles = 0
#     totalFiles = len(fileList_)
    
#     with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
#         # Prepare a list of futures and map them to file names
#         futures_ = {executor.submit(processFunction, tag, fileName): (tag, fileName) for tag, fileName in zip(tags_, fileList_)}
        
#         # Process results as they complete
#         for future in as_completed(futures_):
#             fileName = futures_[future]  # Get the file name associated with this future
#             try:
#                 future.result()  # Retrieve the result or raise an exception if occurred
#                 completedFiles += 1
#                 percentComplete = (completedFiles / (totalFiles)) * 100
#                 print(f'\n---> {fileName} processed successfully! ({percentComplete:.1f}% complete)')
#             except Exception as exc:
#                 print(f'\n---> {fileName} generated an exception!\n{exc}')
                
#     print("\n---> Multithreading completed!")
#     return

# def Multithread(processFunction, fileList_, max_workers=381):
    
#     print("\n---> Starting multithreading...\n")

#     completedFiles = 0
#     totalFiles = len(fileList_)
    
#     with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
#         # Prepare a list of futures and map them to file names
#         futures_ = {executor.submit(processFunction, fileName): (fileName) for fileName in fileList_}
        
#         # Process results as they complete
#         for future in as_completed(futures_):
#             fileName = futures_[future]  # Get the file name associated with this future
#             try:
#                 future.result()  # Retrieve the result or raise an exception if occurred
#                 completedFiles += 1
#                 percentComplete = (completedFiles / (totalFiles)) * 100
#                 print(f'\n---> {fileName} processed successfully! ({percentComplete:.1f}% complete)')
#             except Exception as exc:
#                 print(f'\n---> {fileName} generated an exception!\n{exc}')
                
#     print("\n---> Multithreading completed!")
#     return
    

# ------------------------------------------------
#                      Run
# ------------------------------------------------

def Run(fileList_, PEs=10, layers="3", particle="all", cut="singles"): 

    '''
    Get failure info
    For each file, collect the events. 
    '''
    
    # Collect failure info
    failureFileName = f"../Txt/MDC2020ae/concatenated/failures_concise/failures_concise_all_{PEs}PEs{layers}Layers_{cut}.csv"

    failures_ = GetFailures(fileList_, failureFileName)
    

        
    # Iterate through the file list.
    # Collected failed events appending them into a master array. 
    
    # print("\n---> Iterating through file list.") 
    
    # # Initialise a master array
    # data_ = ak.Array([])

    # def processFunction(tag, fileName):
        
    #         # Find failure events
    #         thisFailureInfo_ = failureInfo_[failureInfo_["Tag"] == tag]
    #         outputStr = ( 
    #             f"\n--->\n" 
    #             f"fileName: {fileName}\n"
    #             f"tag: {tag}\n"
    #             f"failures:\n{thisFailureInfo_}\n"
    #             f"---"
    #         )
    #         print(outputStr)
    
    #         # Read the file
    #         file = ut.ReadFile(fileName, quiet=True)
    #         # Get array
    #         thisData_ = ut.GetData(file)
            
    #         if False: print(f"\n---> Loaded corresponding data.\n{thisData_}")
    
    #         if False: print(f"\n---> Applying masks.")
    
    #         # Extract unique values from DataFrame
    #         runs_ = set(thisFailureInfo_["evtinfo.run"])
    #         subruns_ = set(thisFailureInfo_[" evtinfo.subrun"])
    #         events_ = set(thisFailureInfo_[" evtinfo.event"])
    
    #         # Construct masks
    #         runCondition = ak.any([thisData_["evt"]["evtinfo.run"] == value for value in runs_], axis=0)
    #         subrunCondition = ak.any([thisData_["evt"]["evtinfo.subrun"] == value for value in subruns_], axis=0)
    #         eventCondition = ak.any([thisData_["evt"]["evtinfo.event"] == value for value in events_], axis=0)
    
    #         # Apply masks
    #         thisData_ = thisData_[runCondition & subrunCondition & eventCondition]
            
    #         # Append to master array
    #         if False: print(f"\n---> Appending failures to master array.")
    #         data_ = ak.concatenate([data_, thisData_], axis=0)

    # Multithread2(processFunction, tags_, fileList_)
        
    # # This takes quite a while. Can we multithread?
    # for tag, fileName in zip(tags_, fileList_): 

    #     # if (tag != "001205_00000064"):
    #     #     continue

    #     # Find failure events
    #     thisFailureInfo_ = failureInfo_[failureInfo_["Tag"] == tag]
    #     outputStr = ( 
    #         f"\n--->\n" 
    #         f"fileName: {fileName}\n"
    #         f"tag: {tag}\n"
    #         f"failures:\n{thisFailureInfo_}\n"
    #         f"---"
    #     )
    #     print(outputStr)

    #     # Read the file
    #     file = ut.ReadFile(fileName, quiet=True)
    #     # Get array
    #     thisData_ = ut.GetData(file)
        
    #     print(f"\n---> Loaded corresponding data.\n{thisData_}")

    #     print(f"\n---> Applying masks.")

    #     # Extract unique values from DataFrame
    #     runs_ = set(thisFailureInfo_["evtinfo.run"])
    #     subruns_ = set(thisFailureInfo_[" evtinfo.subrun"])
    #     events_ = set(thisFailureInfo_[" evtinfo.event"])

    #     # Construct masks
    #     runCondition = ak.any([thisData_["evt"]["evtinfo.run"] == value for value in runs_], axis=0)
    #     subrunCondition = ak.any([thisData_["evt"]["evtinfo.subrun"] == value for value in subruns_], axis=0)
    #     eventCondition = ak.any([thisData_["evt"]["evtinfo.event"] == value for value in events_], axis=0)

    #     # Apply masks
    #     thisData_ = thisData_[runCondition & subrunCondition & eventCondition]
        
    #     # Append to master array
    #     print(f"\n---> Appending failures to master array.")
    #     data_ = ak.concatenate([data_, thisData_], axis=0)
        

    # ut.PrintNEvents(data_, len(failureInfo_["Tag"]))


    
        # print(data_["evt"]["evtinfo.event"])
        # print()
        # print(data_)
        
        
        # return
        

    return

# ------------------------------------------------
#                      Main
# ------------------------------------------------

def main():

    defname = "nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.root"
    # fileList_ = rd.GetFileList(defname)

    fileList_ = ['nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000000.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000001.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000002.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000005.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000006.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000007.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000009.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000011.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000012.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000013.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000014.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000015.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000016.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000017.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000018.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000019.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000020.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000024.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000025.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000026.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000027.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000028.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000029.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000030.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000031.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000032.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000034.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000035.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000036.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000037.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000039.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000040.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000042.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000044.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000048.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000049.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000050.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000051.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000052.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000053.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000054.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000056.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000057.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000059.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000062.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000064.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000066.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000069.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000072.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000073.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000074.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000076.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000078.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000080.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000082.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000083.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000089.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000094.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000097.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000098.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000099.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000101.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000103.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000108.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000110.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000116.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000117.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000118.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000120.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000121.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000129.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000130.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000136.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000144.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000145.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000150.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000151.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000157.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000160.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000171.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000176.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000178.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000180.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000181.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000184.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000191.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000210.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000214.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000216.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000228.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000230.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000231.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000243.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000321.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000544.root', 'nts.sgrant.CosmicCRYExtractedCatTriggered.MDC2020ae_best_v1_3.001205_00000554.root']
    # print(fileList_)
    Run(fileList_) 

    return

if __name__ == "__main__":
    main()
