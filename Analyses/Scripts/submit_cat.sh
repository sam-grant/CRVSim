#!/bin/bash

# Submit concatenation jobs
NCORES=20
DATASET="MDC2020ae"
# PARTICLE_=("all") # "muons" "non_muons")
# # TRIGGER_=("crv_trigger" "crv_2layers_trigger" "crv_3layers_trigger" "trk_trigger" "trk_crv_trigger" "trk_crv_2layers_trigger" "trk_crv_3layers_trigger" "trk_crv2_trigger" "trk_crv2_2layers_trigger" "trk_crv2_3layers_trigger" "trk_crv3_trigger")
# # TRIGGER_=("crv_trigger" "trk_crv_trigger" "trk_trigger" "trk_crv2_trigger" "crv_2layers_trigger" "trk_crv_2layers_trigger" "trk_crv2_2layers_trigger")
# TRIGGER_=("crv_2layers_trigger" "trk_crv_2layers_trigger" "trk_crv2_2layers_trigger")
# PE_=($(seq 10 5 130)) # (10) # $(seq 10 5 130))
# LAYER_=(3) # 2 3)


PARTICLE_=("all" "non_muons") 
LAYER_=(3) 
PE_=($(seq 10 5 130))
TRIGGER_=("trk_trigger")

for PE in ${PE_[@]}; do 
    for PARTICLE in ${PARTICLE_[@]}; do 
        for LAYER in ${LAYER_[@]}; do
            for TRIGGER in ${TRIGGER_[@]}; do
                echo "${DATASET} ${TRIGGER} ${PARTICLE} ${LAYER} ${PE}"
                # break 4 # UNCOMMENT FOR TESTING
            done
        done
    done
done | xargs -P $NCORES -I {} bash -c ". run_cat.sh {}"