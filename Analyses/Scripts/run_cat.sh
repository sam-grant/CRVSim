# Run concatenation job
DATASET=$1
FILTER=$2
PARTICLE=$3
LAYER=$4
PE=$5

coincidenceConditions="${PE}PEs${LAYER}Layers"

config="${PARTICLE}_${coincidenceConditions}_${FILTER}"

i="failures_ntuple"

fout="../Txt/${DATASET}/concatenated/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

header=$(awk 'NR==1' $(ls ../Txt/${DATASET}/${i}/001205_00000000/${i}_${config}.csv -1))

echo "${header}" >> $fout 

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    # id=${fin%%*/${i}}
    # id=${id##*/${i}}
    id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')

    # Everything except the header
    results=$(awk 'NR>1' "$fin")
    # echo "${id}, ${results}" >> $fout
    echo "${results}" >> $fout

done

echo "---> Written ${fout}"

i="results"

fout="../Txt/${DATASET}/concatenated/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

header=$(awk 'NR==1' $(ls ../Txt/${DATASET}/${i}/001205_00000000/results_${config}.csv -1))

echo "Start run/subrun, ${header}" >> $fout 

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
    results=$(awk 'NR==2' "$fin")
    echo "${id}, ${results}" >> $fout
done

echo "---> Written ${fout}"

i="failures_concise"

fout="../Txt/${DATASET}/concatenated/${i}_${config}.csv"

if [ -f $fout ]; then
    rm $fout && touch $fout
else 
    touch $fout
fi

header=$(awk 'NR==1' $(ls ../Txt/${DATASET}/${i}/001205_00000000/${i}_${config}.csv -1))
echo "Start run/subrun, ${header}" >> $fout 

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
    awk 'NR>1 {print "'"$id"', " $0}' "$fin" >> "$fout"
done

echo "---> Written ${fout}"

i="failures_verbose"
fout="../Txt/${DATASET}/concatenated/${i}_${config}.csv"

for fin in $(ls ../Txt/${DATASET}/${i}/0*/${i}_${config}.csv | sort -V); do
    cat $fin >> $fout
done

echo "---> Written ${fout}"