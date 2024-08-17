#  Count completed tasks

recon="MDC2020ae" # "original"
particle_=("all" "muons" "non_muons")
filter_=("singles" "singles_track_cuts")
layers_=(2 3) 
PEs_=($(seq 10 5 130.0))

# Calculate lengths of arrays
len_particle_=${#particle_[@]}
len_filter_=${#filter_[@]}
len_layers_=${#layers_[@]}
len_PEs_=${#PEs_[@]}

# Find completed jobs 
baseDir="../Txt/${recon}/results"
numFiles=$(ls $baseDir | wc -l) 
# Calculate the total number of tasks
totalTasks=$(( numFiles * len_particle_ * len_filter_ * len_layers_ * len_PEs_ ))
completedTasks=0

for dir in `ls ${baseDir}`; do 

    # count=$(ls ${baseDir}/${dir} | wc -l) + count

    # ((i++))

    # echo "*************************"
    # echo $dir
    
    # Now look and check which config failed
    for PE in ${PEs_[@]}; do
        for particle in ${particle_[@]}; do
            for filter in ${filter_[@]}; do
                for layer in ${layers_[@]}; do
                    # results_all_100PEs2Layers_singles.csv
                    file="results_${particle}_${PE}PEs${layer}Layers_${filter}.csv"
                    filePath="${baseDir}/${dir}/${file}"
                    # echo $filePath 
                    if [ -f $filePath ]; then
                        ((completedTasks++))
                    fi
                    # break 4
                done
            done
        done
    done
    
done

echo "${completedTasks}/${totalTasks}"
# printf "%d failed files" $i
# echo ", written configurations to ${output_file}"


