# Concatenate all csv and txt result files 
# Script needs some work...

recon="MDC2020ae" # "original"
particle_=("all" "muons" "non_muons")
coincidenceFilter="one_coincidence_per_trigger_sector"

PEs_=($(seq 10.0 5.0 130.0))

for particle in "all" "muons" "non_muons"; do
    for layer in 2 3; do
        for PE in ${PEs_[@]}; do 
        
            # Need 1 decimal place to match Python
            # PE=$(printf "%.1f" "$PE")
            
            coincidenceConditions="${PE}PEs${layer}Layers"

            config="${particle}_${coincidenceConditions}_${coincidenceFilter}"

            i="results"

            fout="../Txt/${recon}/concatenated2/${i}_${config}.csv"

            if [ ! -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi
                
            header=$(awk 'NR==1' $(ls ../Txt/${recon}/${i}/001205_00000000/results_${config}.csv -1))

            echo "Start run/subrun, ${header}" >> $fout 

            for fin in $(ls ../Txt/${recon}/${i}/0*/${i}_${config}.csv | sort -V); do
                # id=${fin%%*/${i}}
                # id=${id##*/${i}}
                id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')

                results=$(awk 'NR==2' "$fin")
                echo "${id}, ${results}" >> $fout
            done

            echo "---> Written ${fout}"


            i="failures_concise"

            fout="../Txt/${recon}/concatenated2/${i}_${config}.csv"

            if [ ! -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi
            
            header=$(awk 'NR==1' $(ls ../Txt/${recon}/${i}/001205_00000000/${i}_${config}.csv -1))
            echo "Start run/subrun, ${header}" >> $fout 

            for fin in $(ls ../Txt/${recon}/${i}/0*/${i}_${config}.csv | sort -V); do
                # id=${fin%%${i}}
                # id=${id##*${i}}
                id=$(echo "$fin" | awk -F'/' '{print $(NF-1)}')
                awk 'NR>1 {print "'"$id"', " $0}' "$fin" >> "$fout"
            done

            echo "---> Written ${fout}"
        
            i="failures_verbose"
            fout="../Txt/${recon}/concatenated2/${i}_${config}.csv"

            if [ ! -f $fout ]; then
                rm $fout && touch $fout
            else 
                touch $fout
            fi
            
            for fin in $(ls ../Txt/${recon}/${i}/0*/${i}_${config}.csv | sort -V); do
                cat $fin >> $fout
            done

            echo "---> Written ${fout}"

        done

    done

done
