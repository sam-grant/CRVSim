for FILE in $(cat ../Txt/Lists/MDC2020ae.txt); do 

    ID=${FILE%%.root}
    ID=${ID##*.}

    # mkdir ../results/$ID
    #mkdir failures_concise/$ID
    #mkdir failures_verbose/$ID
    # mkdir ../Txt/reprocessed/PEsPerLayer/$ID
    mkdir -p ../Txt/MDC2020ae/results/$ID
    mkdir -p ../Txt/MDC2020ae/failures_concise/$ID
    mkdir -p ../Txt/MDC2020ae/failures_verbose/$ID
    mkdir -p ../Txt/MDC2020ae/failures_ntuple/$ID

done
