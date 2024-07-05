# Concatenate all csv and txt result files 
# Script needs some work...

reproc="MDC2020ae" # "original"
particle_=("all" "muons" "non_muons")
# coincidenceConditions="ana1"
coincidenceFilter="one_coincidence_per_trigger_sector"

# for particle in "${particle_[@]}"; do
#     for layer in 2 3; do 
#         for PE in 10 12 14 16 18 20 22 24 26 28 30 32; do # 

# for FILE in $(cat ../Txt/Lists/MDC2020ae.txt); do 
# for PARTICLE in "all"; do 
for particle in "all"; do # "muons" "non_muons"; do
    # for LAYERS in 2; do
    for layer in 2 3; do
        for PE in 10; do #  12 14 16 18 20 22 24 26 28 30 32; do

            coincidenceConditions="${PE}PEs${layer}Layers"

            config="${particle}_${coincidenceConditions}_${coincidenceFilter}"

            i="failures_ntuple"

            fout="../Txt/${reproc}/concatenated/${i}_${config}.csv"

            if [ -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi

            header=$(awk 'NR==1' $(ls ../Txt/${reproc}/${i}/001205_00000000/${i}_${config}.csv -1))

            echo "${header}" >> $fout 

            for fin in $(ls ../Txt/${reproc}/${i}/0*/${i}_${config}.csv | sort -V); do
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

            fout="../Txt/${reproc}/concatenated/${i}_${config}.csv"

            if [ -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi

            header=$(awk 'NR==1' $(ls ../Txt/${reproc}/${i}/001205_00000000/results_${config}.csv -1))

            echo "Start run/subrun, ${header}" >> $fout 

            for fin in $(ls ../Txt/${reproc}/${i}/0*/${i}_${config}.csv | sort -V); do
                # id=${fin%%*/${i}}
                # id=${id##*/${i}}
                id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')

                results=$(awk 'NR==2' "$fin")
                echo "${id}, ${results}" >> $fout
            done

            echo "---> Written ${fout}"

            i="failures_concise"

            fout="../Txt/${reproc}/concatenated/${i}_${config}.csv"

            if [ -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi

            header=$(awk 'NR==1' $(ls ../Txt/${reproc}/${i}/001205_00000000/${i}_${config}.csv -1))
            echo "Start run/subrun, ${header}" >> $fout 

            for fin in $(ls ../Txt/${reproc}/${i}/0*/${i}_${config}.csv | sort -V); do
                # id=${fin%%${i}}
                # id=${id##*${i}}
                id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
                awk 'NR>1 {print "'"$id"', " $0}' "$fin" >> "$fout"
            done

            echo "---> Written ${fout}"

            i="failures_verbose"
            fout="../Txt/${reproc}/concatenated/${i}_${config}.csv"

            for fin in $(ls ../Txt/${reproc}/${i}/0*/${i}_${config}.csv | sort -V); do
                cat $fin >> $fout
            done

            echo "---> Written ${fout}"

        done

    done

done
