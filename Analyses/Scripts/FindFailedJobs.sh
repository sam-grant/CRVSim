# tags_=("001205_00000048" "001205_00000059" "001205_00000108" "001205_00000073" "001205_00000228" "001205_00000006" "001205_00000035" "001205_00000015" "001205_00000024" "001205_00000121" "001205_00000042" "001205_00000129" "001205_00000072" "001205_00000076" "001205_00000009" "001205_00000101" "001205_00000017")
# for tag in ${tags_[@]}; do
#     echo $tag; 
#     ls -ltrh ../Txt/MDC2020ae/results/${tag}
# done

recon="MDC2020ae" # "original"
particle_=("all" "muons" "non_muons")
filter_=("singles" "singles_track_cuts")
layers_=(2 3) 
PEs_=($(seq 10 5 130.0))

output_file="../Txt/${recon}/FailedJobs/failures.csv"

if [ -f $output_file ]; then
    echo "---> Output file exists, making a copy and deleting original."
    cp $output_file "${output_file}.BK"
    echo "${output_file}.BK"
    rm $output_file
fi

echo "Tag,PEs,Layer,Particle,Filter" >> $output_file

# Find failed jobs 
baseDir=../Txt/${recon}/results
i=0
for dir in `ls ${baseDir}`; do 

    count=$(ls ${baseDir}/${dir} | wc -l)
    
    if (( count == 300 )); then
        continue
    fi

    ((i++))

    echo "*************************"
    echo $dir
    
    # Now look and check which config failed
    for PE in ${PEs_[@]}; do
        for particle in ${particle_[@]}; do
            for filter in ${filter_[@]}; do
                for layer in ${layers_[@]}; do
                    # results_all_100PEs2Layers_singles.csv
                    file="results_${particle}_${PE}PEs${layer}Layers_${filter}.csv"
                    filePath="${baseDir}/${dir}/${file}"
                    # echo $filePath 
                    if [ ! -f $filePath ]; then
                        echo $filePath 
                        echo "$dir,$PE,$layer,$particle,$filter" >> $output_file
                    fi
                done
            done
        done
    done
    
done

printf "%d failed files" $i
echo ", written configurations to ${output_file}"


