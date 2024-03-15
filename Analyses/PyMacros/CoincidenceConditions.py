# Samuel Grant Feb 2024
# Container for coincidence conditions 

# This needs some updates! 

# You should have a function to produce the dict 



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

    # Inefficient and lazy, oh well
    # "10PEs2Layers": {
    #     "PEthreshold": 10
    #     ,"nLayers": 2
    #     ,"minSlope": -11
    #     ,"maxSlope": 11
    # },
    "10PEs2Layers": {
        "PEthreshold": 10
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "11PEs2Layers": {
        "PEthreshold": 11
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "12PEs2Layers": {
        "PEthreshold": 12
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "13PEs2Layers": {
        "PEthreshold": 13
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "14PEs2Layers": {
        "PEthreshold": 14
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "15PEs2Layers": {
        "PEthreshold": 15
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "16PEs2Layers": {
        "PEthreshold": 16
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "17PEs2Layers": {
        "PEthreshold": 17
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "18PEs2Layers": {
        "PEthreshold": 18
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "19PEs2Layers": {
        "PEthreshold": 19
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "20PEs2Layers": {
        "PEthreshold": 20
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "22PEs2Layers": {
        "PEthreshold": 22
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "24PEs2Layers": {
        "PEthreshold": 24
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },       
    "26PEs2Layers": {
        "PEthreshold": 26
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "28PEs2Layers": {
        "PEthreshold": 28
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "30PEs2Layers": {
        "PEthreshold": 30
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "32PEs2Layers": {
        "PEthreshold": 32
        ,"nLayers": 2
        ,"minSlope": -11
        ,"maxSlope": 11
    },

    "10PEs3Layers": {
        "PEthreshold": 10
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "11PEs3Layers": {
        "PEthreshold": 11
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "12PEs3Layers": {
        "PEthreshold": 12
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "13PEs3Layers": {
        "PEthreshold": 13
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "14PEs3Layers": {
        "PEthreshold": 14
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "15PEs3Layers": {
        "PEthreshold": 15
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "16PEs3Layers": {
        "PEthreshold": 16
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "17PEs3Layers": {
        "PEthreshold": 17
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "18PEs3Layers": {
        "PEthreshold": 18
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "19PEs3Layers": {
        "PEthreshold": 19
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "20PEs3Layers": {
        "PEthreshold": 20
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "22PEs3Layers": {
        "PEthreshold": 22
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "24PEs3Layers": {
        "PEthreshold": 24
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },       
    "26PEs3Layers": {
        "PEthreshold": 26
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "28PEs3Layers": {
        "PEthreshold": 28
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "30PEs3Layers": {
        "PEthreshold": 30
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    },
    "32PEs3Layers": {
        "PEthreshold": 32
        ,"nLayers": 3
        ,"minSlope": -11
        ,"maxSlope": 11
    }
}