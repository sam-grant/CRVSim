FILE=$1
PARTICLE=$2
LAYERS=$3
PE=$4 

COINC="${PE}PEs${LAYERS}Layers"

# FILENAME=$(basename "$FILE" ".root") 

ID=${FILE%%.root}
ID=${ID##*.}

LOGDIR="../Logs/$ID"

if [ -d $LOGDIR ]; then
    # rm $LOGDIR/*.log 
    echo "Using log directory ${LOGDIR}"
else 
    mkdir $LOGDIR
    # echo "Created log directory ${LOGDIR}"
fi

LOGFILE="${LOGDIR}/${PARTICLE}_${COINC}.log"

if [ -f $LOGFILE ]; then
    rm $LOGFILE && touch $LOGFILE
else 
    touch $LOGFILE
fi

echo "Log: ${LOGFILE}"

# Check if we have already done this one
if [ -f "../Txt/MDC2020ae/results/${ID}/results_${PARTICLE}_${PE}PEs${LAYERS}Layers_one_coincidence_per_trigger_sector.csv" ]; then
    echo "Already completed ${FILE} ${PARTICLE} ${COINC}... skipping"
else 
    echo "Running ${FILE} ${PARTICLE} ${COINC}"
    python Analyse.py $FILE $PARTICLE $COINC >> $LOGFILE 
fi

# done
# if [ -f $LOGFILE ]; then
#     # rm $LOGFILE && touch $LOGFILE
#     echo "Skipping" 
# else 
#     touch $LOGFILE
# fi

# echo "Log: ${LOGFILE}"
# python Analyse.py $FILE $PARTICLE $COINC >> $LOGFILE 