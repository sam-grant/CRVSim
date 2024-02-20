# Samuel Grant Feb 2024
# Container for coincidence conditions 

# Underlying cuts: 
# PEthreshold : 20  //PEs
# maxTimeDifferenceAdjacentPulses : 10  //ns
# maxTimeDifference : 20  //ns
# minOverlapTimeAdjacentPulses : 30  //ns
# minOverlapTime : 30  //ns
# minSlope :-11
# maxSlope : 11  //width direction over thickness direction
# maxSlopeDifference : 4
# coincidenceLayers : 3
# minClusterPEs : 0

coincidenceConditions_ = { 

    "default" : {
        "PEthreshold" : 20 
        ,"nLayers" : 3
        ,"minSlope" : -11
        ,"maxSlope" : 11
    },

    "strict" : {
        "PEthreshold" : 20
        ,"nLayers" : 3
        ,"minSlope" : -5
        ,"maxSlope" : 2
    }

}