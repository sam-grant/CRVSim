# Could be interesting! 

import awkward as ak
from sklearn.cluster import DBSCAN
import numpy as np

def MarkTracks(arr_):
    print("\n---> Getting cosmic ray tracks")

    # Select only the coincidences for clustering
    coincidences = arr_[arr_["is_coincidence"]]

    # Extract relevant features for clustering
    features = ak.to_numpy(ak.zip(
        x=coincidences["crvhit.pos.fCoordinates.x"],
        y=coincidences["crvhit.pos.fCoordinates.y"],
        z=coincidences["crvhit.pos.fCoordinates.z"],
        time=coincidences["crvhit.time"]
    ))

    # Reshape the features array for clustering
    features = features.reshape(-1, 4)

    # Perform clustering using DBSCAN
    clustering = DBSCAN(eps=100, min_samples=2).fit(features)

    # Add a new field 'track_id' to mark tracks
    arr_["track_id"] = -1  # Initialize track_id to -1 for non-coincidence hits

    # Assign track IDs to coincidences
    arr_["track_id"][arr_["is_coincidence"]] = clustering.labels_

    print("Done!")

    return arr_

# Apply the MarkTracks function to your dataset
result = MarkTracks(your_dataset)

