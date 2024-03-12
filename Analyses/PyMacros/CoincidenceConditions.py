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

    # "original" : {
    #     "PEthreshold" : 20 
    #     ,"nLayers" : 3
    #     ,"minSlope" : -11
    #     ,"maxSlope" : 11
    # },

    "default" : {
        "PEthreshold" : 1 
        ,"nLayers" : 1
        ,"minSlope" : -100
        ,"maxSlope" : 100
    },

    "ana1" : {
        "PEthreshold" : 10
        ,"nLayers" : 3
        ,"minSlope" : -11
        ,"maxSlope" : 11
    },

    # Inefficient, but I think in this context it's perfectly OK. 

    "10PEs3Layers": {
        "PEthreshold": 10,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "11PEs3Layers": {
        "PEthreshold": 11,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "12PEs3Layers": {
        "PEthreshold": 12,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "13PEs3Layers": {
        "PEthreshold": 13,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "14PEs3Layers": {
        "PEthreshold": 14,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "15PEs3Layers": {
        "PEthreshold": 15,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "16PEs3Layers": {
        "PEthreshold": 16,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "17PEs3Layers": {
        "PEthreshold": 17,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "18PEs3Layers": {
        "PEthreshold": 18,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "19PEs3Layers": {
        "PEthreshold": 19,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    },
    "20PEs3Layers": {
        "PEthreshold": 20,
        "nLayers": 3,
        "minSlope": -11,
        "maxSlope": 11
    }

}