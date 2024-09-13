for FILE in $(cat ../Txt/FileLists/MDC2020ae.txt); do 

    ID=${FILE%%.root}
    ID=${ID##*.}
    
    mkdir -p ../Txt/MDC2020ae/results/$ID
    mkdir -p ../Txt/MDC2020ae/failures_concise/$ID
    mkdir -p ../Txt/MDC2020ae/failures_verbose/$ID
    # mkdir -p  ../Txt/MDC2020ae/PEsPerLayer/$ID

done
