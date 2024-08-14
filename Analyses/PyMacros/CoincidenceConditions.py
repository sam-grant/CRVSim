# Samuel Grant Feb 2024
# Container for coincidence conditions 

import numpy as np

coincidenceConditions_ = {}
layers_ = [2, 3]
PEs_ = np.arange(10, 135, 5)

for layer in layers_:
    for PE in PEs_:
        key = f"{PE}PEs{layer}Layers"
        coincidenceConditions_[key] = {
            "PEs": PE
            ,"nLayers": layer
            ,"minSlope": -11
            ,"maxSlope": 11
        }