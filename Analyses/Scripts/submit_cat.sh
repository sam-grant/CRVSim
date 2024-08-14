# Submit concatenation jobs
NCORES=8
DATASET="MDC2020ae"
PARTICLE_=("all" "muons" "non_muons")
FILTER="one_coincidence_per_trigger_sector"
PE_=($(seq 10 5 130))
LAYER_=(2 3) # ($(seq 2 1 3))

for PARTICLE in ${PARTICLE_[@]}; do 
    for LAYER in ${LAYER_[@]}; do
        for PE in ${PE_[@]}; do 
            echo "${DATASET} ${FILTER} ${PARTICLE} ${LAYER} ${PE}"
            # echo "${DATASET} ${FILTER} ${PARTICLE} ${LAYER} ${PE}"
            # break 4 # UNCOMMENT FOR TESTING
        done
    done
done | xargs -P $NCORES -I {} bash -c ". run_cat.sh {}"


# for LAYER in ${LAYER_[@]}; do 
#     echo $LAYER
# done