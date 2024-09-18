#!/bin/bash

# tags_=("001205_00000048" "001205_00000059" "001205_00000108" "001205_00000073" "001205_00000228" "001205_00000006" "001205_00000035" "001205_00000015" "001205_00000024" "001205_00000121" "001205_00000042" "001205_00000129" "001205_00000072" "001205_00000076" "001205_00000009" "001205_00000101" "001205_00000017")
# for tag in ${tags_[@]}; do
#     echo $tag; 
#     ls -ltrh ../Txt/MDC2020ae/results/${tag}
# done

recon="MDC2020ae" # "original"
# particle_=("all" "muons" "non_muons")
# layers_=(2 3) 
particle_=("all" "muons" "non_muons")
layers_=(2 3) 
PEs_=($(seq 10 5 130.0))
# cut="track_crv12"
triggerModes_=("crv_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv2_trigger" "trk_crv3_trigger" "trk_crv2_2layers_trigger" "trk_crv2_3layers_trigger")

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

    echo "*************************"
    echo $dir
    
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


