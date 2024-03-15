FILE=$1
PARTICLE=$2
PE=$3 

COINC="${PE}PEs2Layers"

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
python Analyse.py $FILE $PARTICLE $COINC >> $LOGFILE 

# if [ -f $LOGFILE ]; then
#     # rm $LOGFILE && touch $LOGFILE
#     echo "Skipping" 
# else 
#     touch $LOGFILE
# fi

# echo "Log: ${LOGFILE}"
# python Analyse.py $FILE $PARTICLE $COINC >> $LOGFILE 