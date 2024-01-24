for file in `mu2eDatasetFileList nts.mu2e.CosmicCRYExtractedTrk.MDC2020z1_best_v1_1_std_v04_01_00.tka`; do
    echo $file
done | xargs -P $nCores -I {} bash -c ". run.sh {}"