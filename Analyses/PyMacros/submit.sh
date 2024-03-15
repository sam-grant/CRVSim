
# nCores=3
# for filterCondition in "no_filter" "one_coincidence_per_sector" "one_coincidence_per_trigger_sector"; do 
#     for file in `mu2eDatasetFileList nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.tka`; do
#         echo $filterCondition $file
#         sleep 1
#     done | xargs -P $nCores -I {} bash -c ". run.sh {}"
#     sleep 1
# done

NCORES=4
for FILE in $(cat ../Txt/Lists/ReprocessedFiles.txt); do 
    for PARTICLE in "all" "muons" "non_muons"; do
        # for PE in 10 12 14 16 18 20; do
        # for PE in 10; do # 22 24 26 28 30 32; do
        for PE in 10; do #  12 14 16 18 20 22 24 26 28 30 32; do
        
            echo $FILE $PARTICLE $PE 
            # break 3 # UNCOMMENT FOR TESTING
        done 
    done
done | xargs -P $NCORES -I {} bash -c ". run.sh {}"