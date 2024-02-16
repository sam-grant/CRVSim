# Samuel Grant 
# Oct 2023
# Common functions and utilities for CRV KPP analysis
# Utils.py should ONLY contain functions that are actually in-use, unused functions belong in ExtensiveUtils.py

# Colours 

colours = [
    (0., 0., 0.),                                                   # Black
    (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),  # Red
    (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),  # Blue
    (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), # Green
    (1.0, 0.4980392156862745, 0.054901960784313725),                # Orange
    (0.5803921568627451, 0.403921568627451, 0.7411764705882353),    # Purple
    (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),  # Cyan
    (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),   # Pink
    (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), # Brown
    (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),   # Gray 
    (0.7372549019607844, 0.7411764705882353, 0.13333333333333333)   # Yellow
]
    
# ---------------
# TTree wrangling 
# ---------------

branchNamesTrkAna_ = [

    # ---> evtinfo
    "evtinfo.runid" # run ID 
    ,"evtinfo.subrunid" # sub-run ID 
    ,"evtinfo.eventid" # event ID 
    
    # ---> crvhit (reco)
    ,"crvhit.sectorType" # CRV sector hit
    , "crvhit.pos.fCoordinates.fX" # Reconstructed position of the cluster in X 
    , "crvhit.pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
    , "crvhit.pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
    , "crvhit.timeStart" # Earliest time recorded at either end of all bars in the hit
    , "crvhit.timeEnd" # Latest time recorded at either end of all bars in the hit
    , "crvhit.time" # average reconstructed hit time of the cluster.
    , "crvhit.PEs" # total number of photoelectrons in this cluser
    , "crvhit.nHits" # Number of individual bar hits combined in this hit
    , "crvhit.nLayers" # Number of CRV layers that are part of this cluster
    , "crvhit.angle" # slope (in the plane perpendicular to the bar axis of a sector) of the track assumed to be responsible for the cluster (=change in the "layer direction" / change in the "thickness direction")

    # ---> crvhitmc (truth)
    , "crvhitmc.valid" # Records if there is a valid MC match to this CRV reco hit
    , "crvhitmc.pdgId" # PDG ID of the track mostly likely responsible for this cluster
    # , "crvhitmc.primaryPdgId" # PDG ID of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.primaryE" # energy of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.primary.pos.fCoordinates.fX" # start position of the primary particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parentPdgId" # PDG ID of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parentE" # start energy of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parent.fCoordinates.fX" # X start position of the parent particle of the track mostly likely responsible for this cluster
    # , "crvhitmc.parent.fCoordinates.fY" # Y
    # , "crvhitmc.parent.fCoordinates.fZ" # Z
    # , "crvhitmc.gparentPdgId" # grandparent info...
    # , "crvhitmc.gparentE" # "
    # , "crvhitmc.gparent" # "
    # , "crvhitmc.depositedEnergy" # total deposited energy of the cluster based on the CrvSteps
    
]

import uproot
import awkward as ak

# Awkward arrays are more suited to the nested tree structure of TrkAna
def TTreeToAwkwardArray(finName, treeName, branchNames):

    print("---> Converting TTree to awkward array")

    # Open the ROOT file and access the TTree
    file = uproot.open(finName)
    tree = file[treeName]

    # Read the data into Awkward Arrays, if you exclude branchNames it reads all of them 
    arrays = tree.arrays(branchNames, library="ak")

    # Open the ROOT file and access the TTree
    with uproot.open(finName) as file:
        tree = file[treeName]

        # Read the data into Awkward Arrays
        arrays = tree.arrays(branchNames, library="ak")
        
    print("Done!")

    return arrays

# ---------------------------------
# PDGid wrangling
# ---------------------------------

particle_dict = {
    2212: 'proton',
    211: 'pi+',
    -211: 'pi-',
    -13: 'mu+',
    13: 'mu-',
    -11: 'e+',
    11: 'e-',
    321: "kaon+",
    -321: "kaon-",
    311: "kaon0",
    130: "kaon0L",
    310: "kaon0S"
    # Add more particle entries as needed
    }

# --------
# Plotting
# --------

import matplotlib.pyplot as plt

# graph_ a list of xy pairs of lists, graphs_ [ (x0_, y0_), (x1_, y1_)]
def PlotGraphOverlay(graphs_, title=None, xlabel=None, ylabel=None, labels_=[], fout="scatter.png", NDPI=300):
    
    # Create figure and axes
    fig, ax = plt.subplots()

    # Iterate over each pair of xy lists
    for i, (x, y) in enumerate(graphs_):
        # Scatter plot each pair
        ax.scatter(x, y, s=0.5, color=colours[i+1], edgecolor=colours[i+1], marker='o', linestyle='None', label=labels_[i])

    # Set title, xlabel, and ylabel
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)  
    ax.tick_params(axis='y', labelsize=14)  

    # Check if x or y values exceed 9999 for scientific notation
    if any(max(x_) > 999 for x_, _ in graphs_):
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        ax.xaxis.offsetText.set_fontsize(14)
    if any(max(y_) > 999 for _, y_ in graphs_):
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.offsetText.set_fontsize(14)

    ax.legend(loc="best", frameon=False, fontsize=14, markerscale=5)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

def Plot1DOverlay(hists_, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", label_=None, legPos="best", NDPI=300, includeBlack=False, logY=False, legFontSize=12):

    # Create figure and axes
    fig, ax = plt.subplots()

    # Define a colormap
    # cmap = cm.get_cmap('tab10') # !!deprecated!!

    # Create the colormap
    cmap = ListedColormap(colours)

    # cmap = cm.get_cmap('tab10')

    # Iterate over the hists and plot each one
    for i, hist in enumerate(hists_):
        colour = colours[i]
        if not includeBlack: colour = colours[i+1]
        counts, bin_edges, _ = ax.hist(hist, bins=nbins, range=(xmin, xmax), histtype='step', edgecolor=colour, linewidth=1.0, fill=False, density=False, color=colour, label=label_[i], log=logY)

    # Set x-axis limits
    ax.set_xlim(xmin, xmax)

    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=14)  # Set y-axis tick label font size
    
    # Scientific notation
    if ax.get_xlim()[1] > 9999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(14)
    if ax.get_ylim()[1] > 9999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    # Add legend to the plot
    ax.legend(loc=legPos, frameon=False, fontsize=legFontSize)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

from matplotlib.ticker import ScalarFormatter

def BarChart(data_, label_dict, title=None, xlabel=None, ylabel=None, fout="bar_chart.png", percentage=False, bar_alpha=1.0, bar_color='black', NDPI=300):
    
    # This came from ChatGPT
    # it matches the key of the dict with row in the data array and returns the element as the label
    labels = [label_dict.get(p, 'other') for p in data_]

    # Count occurrences of each label
    unique_labels, label_counts = np.unique(labels, return_counts=True)

    # Only works for particles 

    # Sort labels and counts in descending order
    sorted_indices = np.argsort(label_counts)[::-1]
    unique_labels = unique_labels[sorted_indices]
    label_counts = label_counts[sorted_indices]

    if percentage: 
        label_counts = (label_counts / np.sum(label_counts))*100

    # Create figure and axes
    fig, ax = plt.subplots()

    # print(unique_labels)

    # Plot the bar chart
    indices = np.arange(len(unique_labels))

    # print(indices)
    # for i, index in enumerate(indices):
    #     indices[i] = GetLatexParticleName(index)

    # TODO: handle this better
    n_bars = len(indices)
    bar_width = 3.0 / n_bars
    if(n_bars == 3.0): 
        bar_width = 2.0 / n_bars
    elif(n_bars == 2.0):
        bar_width = 1.0 / n_bars


    ax.bar(indices, label_counts, align='center', alpha=bar_alpha, color=bar_color, width=bar_width, fill=False, hatch='/', linewidth=1, edgecolor='black')

    # Set x-axis labels
    ax.set_xticks(indices)
    ax.set_xticklabels(unique_labels, rotation=0) # 45)

    # Set labels for the chart
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=14)  # Set y-axis tick label font size

    # Scientific notation
    # if ax.get_xlim()[1] > 999:
    #     ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #     ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #     ax.xaxis.offsetText.set_fontsize(14)
    if ax.get_ylim()[1] > 999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    # ax.legend(loc="best", frameon=False, fontsize=14)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

def BarChartOverlay(data_, label_dict, title=None, xlabel=None, ylabel=None, label_ = [], fout="bar_chart.png", percentage=False, bar_alpha=1.0, NDPI=300):

    # Initialize figure and axis
    fig, ax = plt.subplots()

    # Initialize variables for bar width calculation
    n_data_sets = len(data_)
    bar_width = 0.5 / n_data_sets

    # Get unique labels from the label dictionary
    unique_labels = list(label_dict.values())

    # Loop through each dataset
    for i, dataset in enumerate(data_):
    
        labels = [label_dict.get(p, 'other') for p in dataset]
        unique_labels_data, label_counts_data = np.unique(labels, return_counts=True)

        # Reorder label counts based on unique_labels
        label_counts = [label_counts_data[unique_labels_data == label][0] if label in unique_labels_data else 0 for label in unique_labels]

        if percentage:
            label_counts = (np.array(label_counts) / sum(label_counts)) * 100

        # Calculate the position of bars
        indices = np.arange(len(unique_labels)) + i * bar_width

        # Plot the bar chart for the current dataset
        # ax.bar(indices, label_counts, align='center', alpha=bar_alpha, color=colours[i+1], edgecolor=colours[i+1], width=bar_width, fill=False, hatch='/', linewidth=1, label=label_[i])
        ax.bar(indices, label_counts, align='center', alpha=bar_alpha, color=colours[i+1], edgecolor=colours[i+1], width=bar_width, linewidth=1, label=label_[i])

    # Set x-axis labels
    ax.set_xticks(indices - bar_width * (n_data_sets - 1) / 2)
    ax.set_xticklabels(unique_labels, rotation=0)

    # Set labels for the chart
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Set scientific notation for y-axis if necessary
    if ax.get_ylim()[1] > 999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.offsetText.set_fontsize(14)

    # Add legend
    ax.legend(loc="best", frameon=False, fontsize=14)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()