#!/bin/bash

# tags_=("001205_00000048" "001205_00000059" "001205_00000108" "001205_00000073" "001205_00000228" "001205_00000006" "001205_00000035" "001205_00000015" "001205_00000024" "001205_00000121" "001205_00000042" "001205_00000129" "001205_00000072" "001205_00000076" "001205_00000009" "001205_00000101" "001205_00000017")
# for tag in ${tags_[@]}; do
#     echo $tag; 
#     ls -ltrh ../Txt/MDC2020ae/results/${tag}
# done

recon="MDC2020ae" # "original"
# particle_=("all" "muons" "non_muons")
# layers_=(2 3) 
# particle_=("all" "muons" "non_muons")
# layers_=(2 3) 
# PEs_=($(seq 10 5 130.0))

particle_=("all" "non_muons") # ("all" "muons" "non_muons")
layers_=(3) #2 3)
PEs_=($(seq 15 5 130)) #($(seq 10 5 130))
# triggerModes_=("crv_trigger" "trk_crv_trigger" "trk_trigger" "trk_crv2_trigger" "crv_2layers_trigger" "trk_crv_2layers_trigger" "trk_crv2_2layers_trigger") #  ("crv_2layers_trigger" "crv_3layers_trigger" "trk_crv_2layers_trigger" "trk_crv_3layers_trigger") # "crv_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv2_trigger" "trk_crv3_trigger" "trk_crv2_2layers_trigger" "trk_crv2_2layers_trigger")
triggerModes_=("trk_trigger") # ("crv_2layers_trigger" "trk_crv_2layers_trigger" "trk_crv2_2layers_trigger")


# cut="track_crv12"
# triggerModes_=("crv_2layers_trigger" "crv_3layers_trigger" "trk_crv_2layers_trigger" "trk_crv_3layers_trigger")
# triggerModes_=("crv_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv2_trigger" "trk_crv3_trigger" "trk_crv2_2layers_trigger" "trk_crv2_3layers_trigger")

    # 1. CRV-DS and CRV-L trigger, crv_trigger X
    # 4. Tracker trigger, trk_triggerX
    # 5. CRV and tracker trigger, trk_crv_trigger X
    # 6. CRV-DS and tracker trigger, trk_crv2_trigger X
    # 7. CRV-L-end and tracker trigger, trk_crv3_trigger X
    # 8. CRV-DS (3 layers active) and tracker trigger, trk_crv2_3layers_trigger X
    # 9. CRV-DS (2 layers active) and tracker trigger, trk_crv2_2layers_trigger X
    # 10. CRV-DS and CRV-L trigger (2 layers active), crv_2layers_trigger X
    # 11. CRV-DS and CRV-L trigger (3 layers active), crv_3layers_trigger X
    # 12. CRV and tracker trigger (2 layers active), trk_crv_2layers_triggerX
    # 13. CRV and tracker trigger (3 layers active), trk_crv_3layers_trigger 
    
# triggerModes_=("crv_trigger" "crv_2layers_trigger" "crv_3layers_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv_2layers_trigger" "trk_crv_3layers_trigger" "trk_crv2_trigger" "trk_crv2_2layers_trigger" "trk_crv2_3layers_trigger" "trk_crv3_trigger")
# triggerModes_=("crv_trigger" "trk_crv_trigger" "trk_trigger" "trk_crv2_trigger" "crv_2layers_trigger" "trk_crv_2layers_trigger" "trk_crv2_2layers_trigger") #  ("crv_2layers_trigger"

output_file="../Txt/${recon}/FailedJobs/failures.csv"
previous_output_file="../Txt/${recon}/FailedJobs/previous_failures.csv"

if [ -f $output_file ]; then
    echo "---> Output file exists, making a copy and deleting original."
    cp $output_file ${previous_output_file}
    rm $output_file
fi

echo "Tag,PEs,Layer,Particle,Trigger" >> $output_file

# Find failed jobs 
baseDir=../Txt/${recon}/results
# i=0
for dir in `ls ${baseDir}`; do 

    count=$(ls ${baseDir}/${dir} | wc -l)
    
    # if (( count == 1050 )); then
    #     continue
    # fi

    # ((i++))

    # echo "*************************"
    # echo $dir
    
    # Now look and check which config failed
    for PE in ${PEs_[@]}; do
        for particle in ${particle_[@]}; do
            for layer in ${layers_[@]}; do
                for triggerMode in ${triggerModes_[@]}; do
                    # results_all_100PEs2Layers_singles.csv
                    file="results_${particle}_${PE}PEs${layer}Layers_${triggerMode}.csv"
                    filePath="${baseDir}/${dir}/${file}"
                    # echo $filePath 
                    if [ ! -f $filePath ]; then
                        # echo $filePath 
                        echo "$dir,$PE,$layer,$particle,$triggerMode" >> $output_file
                    fi
                done
            done
        done
    done
    
done

printf "%d failed configurations" $(($(cat $output_file | wc -l) - 1))
echo ", written configurations to ${output_file}"


