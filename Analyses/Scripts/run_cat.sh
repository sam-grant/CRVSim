#!/bin/bash

# Run concatenation job
DATASET=$1
TRIGGER=$2
PARTICLE=$3
LAYER=$4
PE=$5

# Need 1 decimal place to match Python
# PE=$(printf "%.1f" "$PE")
            
coincidenceConditions="${PE}PEs${LAYER}Layers"

config="${PARTICLE}_${coincidenceConditions}_${TRIGGER}"

# Set up directory

baseDir="../Txt/${DATASET}/concatenated"

if [ ! -d $baseDir ]; then
    mkdir $baseDir
fi

# RESULTS 

i="results"

thisDir="${baseDir}/${i}"
if [ ! -d "${thisDir}" ]; then
    mkdir $thisDir
fi

fout="${thisDir}/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

# Set header
header=$(awk 'NR==1' $(ls ../Txt/${DATASET}/${i}/001205_00000000/results_${config}.csv -1))
echo "Tag,${header}" >> $fout 
# Reset header manually
# echo "Tag,Total,Successes,Failures,Efficiency [%],Inefficiency [%]" >> $fout

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
    results=$(awk 'NR==2' "$fin")
    echo "${id},${results}" >> $fout
done

echo "---> Written ${fout}"

# FAILURES_CONCISE 

i="failures_concise"

thisDir="${baseDir}/${i}"
if [ ! -d "${thisDir}" ]; then
    mkdir $thisDir
fi

fout="${thisDir}/${i}_${config}.csv"

# fout="../Txt/${DATASET}/concatenated/${i}/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

header=$(awk 'NR==1' $(ls ../Txt/${DATASET}/${i}/001205_00000000/${i}_${config}.csv -1))
echo "tag,${header}" >> $fout 
# echo "tag,evtinfo.run,evtinfo.subrun,evtinfo.event" >> $fout 

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
    awk 'NR>1 {print "'"$id"', " $0}' "$fin" >> "$fout"
done

echo "---> Written ${fout}"

# FAILURES_VERBOSE 

i="failures_verbose"

thisDir="${baseDir}/${i}"
if [ ! -d "${thisDir}" ]; then
    mkdir $thisDir
fi

fout="${thisDir}/${i}_${config}.txt"

# fout="../Txt/${DATASET}/concatenated/${i}/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.txt | sort -V); do
    cat $fin >> $fout
done

echo "---> Written ${fout}"