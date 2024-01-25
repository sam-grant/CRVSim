
nCores=3
for filterCondition in "no_filter" "one_coincidence_per_sector" "one_coincidence_per_trigger_sector"; do 
    for file in `mu2eDatasetFileList nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.tka`; do
        echo $filterCondition $file
        sleep 1
    done | xargs -P $nCores -I {} bash -c ". run.sh {}"
    sleep 1
done