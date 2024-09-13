# Submit concatenation jobs
# NCORES=8
# DATASET="MDC2020ae"
# PARTICLE_=("all" "muons" "non_muons")
# FILTER_=("no_track_cuts" "track_cuts") 
# # PE_=($(seq 10 5 90))
# PE_=($(seq 10 5 130))
# LAYER_=(2 3)

NCORES=8
DATASET="MDC2020ae"
PARTICLE_=("all") #"muons" "non_muons")
FILTER_=("track_crv12") # "track_crv1_only") #track_crv12") # "no_track_cuts" "track_cuts") 
# PE_=($(seq 10 5 90))
PE_=($(seq 10 5 130))
LAYER_=(2 3)

# PE_=(10)
# LAYER_=(2)
# PARTICLE_=("all") #  "muons" "non_muons")

for PE in ${PE_[@]}; do 
    for PARTICLE in ${PARTICLE_[@]}; do 
        for LAYER in ${LAYER_[@]}; do
            for FILTER in ${FILTER_[@]}; do
                echo "${DATASET} ${FILTER} ${PARTICLE} ${LAYER} ${PE}"
                # break 4 # UNCOMMENT FOR TESTING
            done
        done
    done
done | xargs -P $NCORES -I {} bash -c ". run_cat.sh {}"


# for LAYER in ${LAYER_[@]}; do 
#     echo $LAYER
# done