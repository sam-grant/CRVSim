#!/bin/bash

recon="MDC2020ae"
particle_=("all" "muons" "non_muons")
layers_=(2 3)
PEs_=($(seq 10 5 130))

len_particle=${#particle_[@]}
len_layers=${#layers_[@]}
len_PEs=${#PEs_[@]}

baseDir="../Txt/${recon}/results"
numFiles=$(ls $baseDir | wc -l)

totalTasks=$((numFiles * len_particle * len_layers * len_PEs))

#logFile="../Txt/Monitoring/completion_log_" + (date +%s) +".csv"
logFile="../Txt/Monitoring/completion_log_$(date +%s).csv"

echo "time,completed_tasks" > $logFile

completedTasks=0
previousCompleted=-1

# Loop to continually monitor the tasks
while [ $completedTasks -lt $totalTasks ]; do
    completedTasks=0
    for dir in $(ls ${baseDir}); do
		for PE in ${PEs_[@]}; do
			for particle in ${particle_[@]}; do
				for layer in ${layers_[@]}; do
					file="results_${particle}_${PE}PEs${layer}Layers_track_cuts.csv"
					filePath="${baseDir}/${dir}/${file}"
					if [ -f $filePath ]; then
					((completedTasks++))
					fi
				done
			done
		done
    done

    # Log only if the number of completed tasks has increased
    if [ $completedTasks -gt $previousCompleted ]; then
	echo "$(date +%s),${completedTasks}" >> $logFile
	previousCompleted=$completedTasks
    fi

    #echo "${completedTasks}/${totalTasks} tasks completed."

    # Sleep for a specified duration before the next check
    sleep 300 # Check every 5 minutes; adjust this as needed
done

echo "---> All tasks completed."


