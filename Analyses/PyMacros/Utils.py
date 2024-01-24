# Samuel Grant 
# Oct 2023
# Common functions and utilities for CRV KPP analysis

# ------------------------------------------------
# Calculations, value formatting, array operations
# ------------------------------------------------

import math
import numpy as np
from scipy import stats

def Round(value, sf):

    if value == 0.00:
        return "0"
    elif math.isnan(value):
        return "NaN"
    else:

        # Determine the order of magnitude
        magnitude = math.floor(math.log10(abs(value))) + 1

        # Calculate the scale factor
        scale_factor = sf - magnitude

        # Truncate the float to the desired number of significant figures
        truncated_value = math.trunc(value * 10 ** scale_factor) / 10 ** scale_factor

        # Convert the truncated value to a string
        truncated_str = str(truncated_value).rstrip('0').rstrip('.')

        return truncated_str

# Stats for histograms tends to assume a normal distribution
# ROOT does the same thing with TH1
def GetBasicStats(data, xmin, xmax):

    filtered_data = data[(data >= xmin) & (data <= xmax)]  # Filter data within range

    N = len(filtered_data)                      
    mean = np.mean(filtered_data)  
    meanErr = stats.sem(filtered_data) # Mean error (standard error of the mean from scipy)
    stdDev = np.std(filtered_data) # Standard deviation
    stdDevErr = np.sqrt(stdDev**2 / (2*N)) # Standard deviation error assuming normal distribution
    underflows = len(data[data < xmin]) # Number of underflows
    overflows = len(data[data > xmax])

    return N, mean, meanErr, stdDev, stdDevErr, underflows, overflows

# def Flatten(x): 
#     return np.concatenate([np.array(xi) for xi in x], axis=0)

# ---------------
# TTree wrangling 
# ---------------

import uproot
import pandas as pd

branchNamesTrkAna_ = [
    # ---> evtinfo
    "evtinfo.eventid" # the event number
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
    # , "crvhitmc.valid" # Records if there is a valid MC match to this CRV reco hit
    # , "crvhitmc.pdgId" # PDG ID of the track mostly likely responsible for this cluster
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

def TTreeToDataFrame(finName, treeName, branchNames):
    
    print("---> Reading", treeName, "in", finName)

    # Open the file
    fin = uproot.open(finName)

    # Get the tree
    tree = fin[treeName]

    if len(tree) == 0:
        return

    # Create an empty dictionary to store the selected columns as NumPy arrays
    branchData = {}

    # Iterate over the specified branch names
    for branchName in branchNames:
        # Check if the branch name exists in the TTree
        if branchName in tree:
            # Load values into an array
            branchData[branchName] = tree[branchName].array(library="np")

    # Create the DataFrame directly from the dictionary of column data
    df = pd.DataFrame(branchData)

    # Close the ROOT file
    fin.file.close()

    # Return the DataFrame
    return df

import awkward as ak

# Awkward arrays are more suited to the nested tree structure of TrkAna
def TTreeToAwkwardArray(finName, treeName, branchNames):

    print("---> Converting TTree to awkward array")

    # Open the ROOT file and access the TTree
    file = uproot.open(finName)
    tree = file[treeName]

    # Read the data into Awkward Arrays, if you exclude branchNames it reads all of them 
    arrays = tree.arrays(branchNames, library="ak")

    print("Done!")

    return arrays

# --------
# Plotting
# --------

import matplotlib.pyplot as plt

def ProfileX(x, y, nBinsX=100, xmin=-1.0, xmax=1.0, nBinsY=100, ymin=-1.0, ymax=1.0): 
   
    # Create 2D histogram with one bin on the y-axis 
    hist, xEdge_, yEdge_ = np.histogram2d(x, y, bins=[nBinsX, nBinsY], range=[[xmin, xmax], [ymin, ymax]])

    # hist, xEdge_ = np.histogram(x, bins=nBinsX, range=[xmin, xmax]) # , [ymin, ymax]])

    # bin widths
    xBinWidths = xEdge_[1]-xEdge_[0]

    # Calculate the mean and RMS values of each vertical slice of the 2D distribution
    xSlice_, xSliceErr_, ySlice_, ySliceErr_, ySliceRMS_ = [], [], [], [], []

    for i in range(len(xEdge_) - 1):

        # Average x-value
        xSlice = x[ (xEdge_[i] < x) & (x <= xEdge_[i+1]) ]

        # Get y-slice within current x-bin
        ySlice = y[ (xEdge_[i] < x) & (x <= xEdge_[i+1]) ]

        # Avoid empty slices
        if len(xSlice) == 0 or len(ySlice) == 0:
            continue

        # Central values are means and errors are standard errors on the mean
        xSlice_.append(np.mean(xSlice))
        xSliceErr_.append(stats.sem(xSlice)) # RMS/sqrt(n)
        ySlice_.append(ySlice.mean()) 
        ySliceErr_.append(stats.sem(ySlice)) 
        ySliceRMS_.append(np.std(ySlice))

    return np.array(xSlice_), np.array(xSliceErr_), np.array(ySlice_), np.array(ySliceErr_), np.array(ySliceRMS_)

import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
# Enable LaTeX rendering
# plt.rcParams["text.usetex"] = True

def Plot1D(data, nBins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", legPos="best", stats=True, peak=False, underOver=False, errors=False, NDPI=300):
    
    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot the histogram with outline
    counts, bin_edges, _ = ax.hist(data, bins=nBins, range=(xmin, xmax), histtype='step', edgecolor='black', linewidth=1.0, fill=False, density=False)

    # Set x-axis limits
    ax.set_xlim(xmin, xmax)

    # Calculate statistics
    N, mean, meanErr, stdDev, stdDevErr, underflows, overflows = GetBasicStats(data, xmin, xmax)
    # peak = np.max(counts)
    # peak_bin_edges = bin_edges[i_peak:i_peak + 2]
    # peak = counts[i_peak]
    # peakErr = (bin_edges[1] - bin_edges[0]) / 2
    # N, mean, meanErr, stdDev, stdDevErr = str(N), Round(mean, 3), Round(mean, 3), Round(meanErr, 1), Round(stdDev, 3), Round(stdDevErr, 1) 

    # Create legend text
    legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}\nStd Dev: {Round(stdDev, 3)}"
    # if errors: legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDevErr, 1)}"
    if errors: legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDevErr, 1)}"
    if peak and not errors: legend_text += f"\nPeak: {Round(GetMode(data, nBins / (xmax - xmin))[0], 3)}"
    if peak and errors: legend_text += f"\nPeak: {Round(GetMode(data, nBins / (xmax - xmin))[0], 3)}$\pm${Round(GetMode(data, nBins / (xmax - xmin))[1], 1)}"
    if underOver: legend_text += f"\nUnderflows: {underflows}\nOverflows: {overflows}"

    # legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDev, 1)}"

    # Add legend to the plot
    if stats: ax.legend([legend_text], loc=legPos, frameon=False, fontsize=14)

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
    # if ax.get_ylim()[1] > 999:
    #     ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #     ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #     ax.yaxis.offsetText.set_fontsize(14)

    if ax.get_xlim()[1] > 9999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(14)
    if ax.get_ylim()[1] > 9999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

def Plot2D(x, y, nBinsX=100, xmin=-1.0, xmax=1.0, nBinsY=100, ymin=-1.0, ymax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", log=False, cb=True, NDPI=300):

    # Filter out empty entries from x and y
    valid_indices = [i for i in range(len(x)) if np.any(x[i]) and np.any(y[i])]

    # Extract valid data points based on the indices
    x = [x[i] for i in valid_indices]
    y = [y[i] for i in valid_indices]

    # Check if the input arrays are not empty and have the same length
    if len(x) == 0 or len(y) == 0:
        print("Input arrays are empty.")
        return
    if len(x) != len(y):
        print("Input arrays x and y have different lengths.")
        return

    # Create 2D histogram
    hist, x_edges, y_edges = np.histogram2d(x, y, bins=[nBinsX, nBinsY], range=[[xmin, xmax], [ymin, ymax]])

    # Set up the plot
    fig, ax = plt.subplots()

    norm = colors.Normalize(vmin=0, vmax=np.max(hist))  
    if log: norm = colors.LogNorm(vmin=1, vmax=np.max(hist)) 

    # Plot the 2D histogram
    im = ax.imshow(hist.T, cmap='inferno', extent=[xmin, xmax, ymin, ymax], aspect='auto', origin='lower', norm=norm)  # , vmax=np.max(hist), norm=colors.LogNorm())
    # im = ax.imshow(hist.T, extent=[xmin, xmax, ymin, ymax], aspect='auto', origin='lower', vmax=np.max(hist))

    # Add colourbar
    if cb: plt.colorbar(im)

    plt.title(title, fontsize=16, pad=10)
    plt.xlabel(xlabel, fontsize=14, labelpad=10)
    plt.ylabel(ylabel, fontsize=14, labelpad=10)

    # Scientific notation
    if ax.get_xlim()[1] > 99999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(14)
    if ax.get_ylim()[1] > 99999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()


def PlotGraph(x, y, title=None, xlabel=None, ylabel=None, fout="scatter.png", NDPI=300):

    # Create a scatter plot with error bars using NumPy arrays 

    # Create figure and axes
    fig, ax = plt.subplots()

    # Fine graphs for CRV visualisation
    ax.scatter(x, y, color='black', s=0.25, edgecolor='black', marker='o', linestyle='None')

    # Set title, xlabel, and ylabel
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=14)  # Set y-axis tick label font size

    # Check if x or y values exceed 9999 for scientific notation
    if max(x) > 999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        ax.xaxis.offsetText.set_fontsize(14)
    if max(y) > 999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.offsetText.set_fontsize(14)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

from mpl_toolkits.mplot3d import Axes3D

def PlotGraph3D(x, y, z, title=None, xlabel=None, ylabel=None, zlabel=None, fout="scatter3D.png", NDPI=300):

    # Create a 3D scatter plot with error bars using NumPy arrays

    # Create figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Fine graphs for CRV visualisation
    ax.scatter(x, y, z) # , s=0.25, marker='o')

    # Set the color of the points to black
    # sc.set_facecolor('black')
    # sc.set_edgecolor('black')
    
    # Set title, labels, and font sizes
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)
    ax.set_zlabel(zlabel, fontsize=14, labelpad=10)

    # Set font size of tick labels on x, y, and z axes
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='z', labelsize=14)

    # Check if x, y, or z values exceed 999 for scientific notation
    if max(x) > 999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        ax.xaxis.offsetText.set_fontsize(14)
    if max(y) > 999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.offsetText.set_fontsize(14)
    if max(z) > 999:
        ax.zaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
        ax.zaxis.offsetText.set_fontsize(14)

    # Save the 3D scatter plot
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()


def PlotGraphErrors(x, y, xerr=np.array([]), yerr=np.array([]), title=None, xlabel=None, ylabel=None, fout="scatter.png", NDPI=300):

   # Create a scatter plot with error bars using NumPy arrays 

    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot scatter with error bars
    if len(xerr)==0: xerr = [0] * len(x) # Sometimes we only use yerr
    if len(yerr)==0: yerr = [0] * len(y) # Sometimes we only use xerr

    if len(x) != len(y): print("Warning: x has length", len(x),", while y has length", len(y))

    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='black', markersize=4, ecolor='black', capsize=2, elinewidth=1, linestyle='None')

    # Set title, xlabel, and ylabel
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


    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

# def PlotGraph3D(x, y, z, xerr=None, yerr=None, zerr=None, title=None, xlabel=None, ylabel=None, zlabel=None, fout="scatter3d.png", NDPI=300):
def PlotGraph3D(x, y, z, title=None, xlabel=None, ylabel=None, zlabel=None, fout="scatter3d.png", NDPI=300):

    # Create a 3D scatter plot with error bars using NumPy arrays

    # Create figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Ensure x, y, z are numpy arrays
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    # Loop through the data and add individual points
    # for i in range(len(x)):
    #     ax.scatter(x[i], y[i], z[i])
        
    ax.scatter(x, z) # , fmt='o', color='black', markersize=4, ecolor='black', capsize=2, elinewidth=1, linestyle='None')

    # # Plot 3D scatter with error bars
    # if xerr==None: xerr = np.zeros(len(x)) 
    # if yerr==None: yerr = np.zeros(len(y)) 
    # if zerr==None: zerr = np.zeros(len(z)) 

    # assert len(x) == len(y) == len(z) == len(xerr) == len(yerr) == len(zerr), "Array dimensions do not match."


    # if len(xerr) == 0:
    #     xerr = [0] * len(x)  # Sometimes we only use yerr
    # if len(yerr) == 0:
    #     yerr = [0] * len(y)  # Sometimes we only use yerr
    # if len(zerr) == 0:
    #     zerr = [0] * len(z)  # Sometimes we only use zerr


    # Set title, labels, and font sizes
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10)
    ax.set_zlabel(zlabel, fontsize=14, labelpad=10)

    # Set font size of tick labels on x, y, and z axes
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='z', labelsize=14)

    # Scientific notation
    if ax.get_xlim()[1] > 9999:
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        ax.xaxis.offsetText.set_fontsize(14)
    if ax.get_ylim()[1] > 9999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax.yaxis.offsetText.set_fontsize(14)
    if ax.get_zlim()[1] > 9999:
        ax.zaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='z', scilimits=(0, 0))
        ax.zaxis.offsetText.set_fontsize(14)

    # Save the 3D scatter plot
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

def Plot1DOverlay(hists, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", labels=None, legPos="upper right", NDPI=300, includeBlack=False, logY=False, legFontSize=12):

    # Create figure and axes
    fig, ax = plt.subplots()

    # Define a colormap
    # cmap = cm.get_cmap('tab10') # !!deprecated!!

    # Define the colourmap colours
    colours = [
        (0., 0., 0.),                                                   # Black
        (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),  # Blue
        (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),  # Red
        (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), # Green
        (1.0, 0.4980392156862745, 0.054901960784313725),                # Orange
        (0.5803921568627451, 0.403921568627451, 0.7411764705882353),    # Purple
        (0.09019607843137255, 0.7450980392156863, 0.8117647058823529),   # Cyan
        (0.8901960784313725, 0.4666666666666667, 0.7607843137254902),   # Pink
        (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), # Brown
        (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),   # Gray 
        (0.7372549019607844, 0.7411764705882353, 0.13333333333333333)  # Yellow
    ]

    # Create the colormap
    cmap = ListedColormap(colours)

    # cmap = cm.get_cmap('tab10')

    # Iterate over the hists and plot each one
    for i, hist in enumerate(hists):
        colour = cmap(i)
        if not includeBlack: colour = cmap(i+1)
        counts, bin_edges, _ = ax.hist(hist, bins=nbins, range=(xmin, xmax), histtype='step', edgecolor=colour, linewidth=1.0, fill=False, density=False, color=colour, label=labels[i], log=logY)

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

def BarChart(data, label_dict, title=None, xlabel=None, ylabel=None, fout="bar_chart.png", percentage=False, bar_alpha=1.0, bar_color='black', NDPI=300):
    
    # This came from ChatGPT
    # it matches the key of the dict with row in the data array and returns the element as the label
    labels = [label_dict.get(p, 'other') for p in data]

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

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

def BarChart2(data_dict, title=None, xlabel=None, ylabel=None, fout="bar_chart.png", percentage=False, bar_alpha=1.0, bar_color='black', NDPI=300):
    
    # Extract labels and counts from the data_dict
    labels = list(data_dict.keys())
    label_counts = list(data_dict.values())

    # Sort labels and counts in descending order
    sorted_indices = np.argsort(label_counts)[::-1]
    labels = np.array(labels)[sorted_indices]
    label_counts = np.array(label_counts)[sorted_indices]

    if percentage: 
        label_counts = (label_counts / np.sum(label_counts))*100

    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot the bar chart
    indices = np.arange(len(labels))

    # TODO: handle this better
    n_bars = len(indices)
    bar_width = 3.0 / n_bars
    if n_bars == 3:
        bar_width = 2.0 / n_bars
    elif n_bars == 2:
        bar_width = 1.0 / n_bars

    ax.bar(indices, label_counts, align='center', alpha=bar_alpha, color=bar_color, width=bar_width, fill=False, hatch='/', linewidth=1, edgecolor='black')

    # Set x-axis labels
    ax.set_xticks(indices)
    ax.set_xticklabels(labels, rotation=0)

    # Set labels for the chart
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlabel, fontsize=14, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=14, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Scientific notation
    if ax.get_ylim()[1] > 999:
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()

def Plot1DRatio(hists, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", labels=None, legPos="best", stats=False, errors=False, NDPI=300, peak=False, invertRatio=False, limitRatio=False, ratioMin=0, ratioMax=1):
    
    if len(hists) > 2: 
        print("!!! ERROR: Plot1DRatio must take two histograms as input !!!")
        return

    # Create figure and axes
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6))

    # Define a colormap
    colours = [
        (0.12156862745098039, 0.4666666666666667, 0.7058823529411765),  # Blue
        (0.8392156862745098, 0.15294117647058825, 0.1568627450980392),  # Red
    ]

    # Create the colormap
    cmap = ListedColormap(colours)

    counts_ = []

    # Iterate over the histograms and plot each one in the top frame
    for i, hist in enumerate(hists):

        colour = cmap(i)

        # Calculate statistics for the current histogram
        N, mean, meanErr, stdDev, stdDevErr, underflows, overflows = GetBasicStats(hist, xmin, xmax)

        # Create legend text
        legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}\nStd Dev: {Round(stdDev, 3)}"
        if errors:
            legend_text = f"Entries: {N}\nMean: {Round(mean, 4)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDevErr, 1)}"
        if peak and not errors:
            legend_text += f"\nPeak: {Round(GetMode(hist, nbins / (xmax - xmin))[0], 3)}"
        if peak and errors:
            legend_text += f"\nPeak: {Round(GetMode(hist, nbins / (xmax - xmin))[0], 3)}$\pm${Round(GetMode(hist, nbins / (xmax - xmin))[1], 1)}"

        if stats:
            label = r"$\bf{"+labels[i]+"}$"+"\n"+legend_text
        else:
            label = labels[i]

        # Plot the current histogram in the top frame
        counts, bin_edges, _ = ax1.hist(hist, bins=nbins, range=(xmin, xmax), histtype='step', edgecolor=colour, linewidth=1.0, fill=False, density=False, color=colour, label=label) 

        # Plot the current histogram in the top frame with error bars
        # hist_err = np.sqrt(hist)
        # bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        # ax1.bar(bin_centers, hist, width=0.1, align='center', alpha=0.7, label=label)
        # ax1.errorbar(bin_centers, hist, yerr=hist_err, fmt='none', color=colour, capsize=2)


        counts_.append(counts)

    # Calculate the ratio of the histograms with a check for division by zero
    ratio = np.divide(counts_[0], counts_[1], out=np.full_like(counts_[0], np.nan), where=(counts_[1] != 0))

    # Calculate the statistical uncertainty for the ratio
    ratio_err = np.divide(np.sqrt(counts_[0]), counts_[1], out=np.full_like(counts_[0], np.nan), where=(counts_[1] != 0))

    # Create a second y-axis for the ratio
    # ax2 = ax1.twinx() # This overlays them
    # Create a separate figure and axis for the ratio plot
    # fig2, ax2 = plt.subplots(figsize=(8, 2))  # Adjust the height as needed

    # Add line at 1.0 
    ax2.axhline(y=1.0, color='gray', linestyle='--', linewidth=1)


    if invertRatio: ratio = np.divide(1, ratio)

    # Plot the ratio in the lower frame with error bars
    ax2.errorbar(bin_edges[:-1], ratio, yerr=ratio_err, color='black', fmt='o', markersize=4, linewidth=1)

    # # Plot the ratio in the lower frame
    # ax2.plot(bin_edges[:-1], ratio, color='black', marker='o', markersize=4, linewidth=0)

    # Format 

    # Set x-axis limits for the top frame
    ax1.set_xlim(xmin, xmax)


    # Remove markers for main x-axis
    ax1.set_xticks([])

    ax2.set_xlabel(xlabel, fontsize=14, labelpad=10)
    ax2.set_ylabel("Ratio", fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.tick_params(axis='y', labelsize=14)

    # Create a second y-axis for the ratio
    ax2.yaxis.tick_left()
    ax2.xaxis.tick_bottom()
    ax2.xaxis.set_tick_params(width=0.5)
    ax2.yaxis.set_tick_params(width=0.5)

    # Set x-axis limits for the ratio plot to match the top frame
    if limitRatio:
        ax2.set_ylim(ratioMin, ratioMax)

    ax2.set_xlim(xmin, xmax)

    # Set titles and labels for both frames
    ax1.set_title(title, fontsize=16, pad=10)
    # ax1.set_xlabel("", fontsize=0, labelpad=10)
    ax1.set_ylabel(ylabel, fontsize=14, labelpad=10)

    # Set font size of tick labels on x and y axes for both frames
    # ax1.tick_params(axis='x', labelsize=14)
    ax1.tick_params(axis='y', labelsize=14)

    # Scientific notation for top frame
    if ax2.get_xlim()[1] > 9999:
        ax2.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax2.xaxis.offsetText.set_fontsize(14)
    if ax1.get_ylim()[1] > 9999:
        ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax1.yaxis.offsetText.set_fontsize(14)

    # Add legend to the top frame
    ax1.legend(loc=legPos, frameon=False, fontsize=14)

    # Adjust the spacing between subplots
    plt.tight_layout()

    # Adjust the spacing between subplots to remove space between them
    plt.subplots_adjust(hspace=0.0)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()