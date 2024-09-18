#!/bin/bash

recon="MDC2020ae"
particle_=("all" "muons" "non_muons")
layers_=(2 3)
PEs_=($(seq 10 5 130))
triggerModes_=("crv_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv2_trigger" "trk_crv3_trigger" "trk_crv2_2layers_trigger" "trk_crv2_2layers_trigger")
len_particle=${#particle_[@]}
len_layers=${#layers_[@]}
len_PEs=${#PEs_[@]}
len_triggers=${#triggerModes_[@]}

baseDir="../Txt/${recon}/results"
numFiles=$(ls $baseDir | wc -l)

totalTasks=$((numFiles * len_particle * len_layers * len_PEs * len_triggers))

#logFile="../Txt/Monitoring/completion_log_" + (date +%s) +".csv"
logFile="../Txt/Monitoring/completion_log_$(date +%s).csv"
# logFile="../Txt/Monitoring/completion_log_1726625479.csv"
echo $logFile

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
					for trigger in ${triggerModes_[@]}; do
						file="results_${particle}_${PE}PEs${layer}Layers_${trigger}.csv"
						filePath="${baseDir}/${dir}/${file}"
						if [ -f $filePath ]; then
							((completedTasks++))
						fi
					done
				done
			done
		done
    done

    # Log only if the number of completed tasks has increased
    if [ $completedTasks -gt $previousCompleted ]; then
	echo "$(date +%s),${completedTasks}" >> $logFile
	previousCompleted=$completedTasks
    fi

    echo "${completedTasks}/${totalTasks} tasks completed."

    # Sleep for a specified duration before the next check
    sleep 60 # Check every minute; adjust this as needed
done

echo "---> All tasks completed."


