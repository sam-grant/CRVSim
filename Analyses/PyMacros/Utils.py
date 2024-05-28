# Samuel Grant 
# Oct 2023
# Common functions and utilities for CRV KPP analysis
# Utils.py should ONLY contain functions that are actually in-use, unused functions belong in ExtensiveUtils.py

coincsBranchName = "crvcoincs" # crvhit

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

# ---------------
# TTree wrangling 
# ---------------

branchNamesTrkAna_ = [

    # ---> evtinfo
    "evtinfo.runid" # run ID 
    ,"evtinfo.subrunid" # sub-run ID 
    ,"evtinfo.eventid" # event ID 
    
    # ---> crvhit (reco)
    , coincsBranchName+".sectorType" # CRV sector hit
    , coincsBranchName+".pos.fCoordinates.fX" # Reconstructed position of the cluster in X 
    , coincsBranchName+".pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
    , coincsBranchName+".pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
    , coincsBranchName+".timeStart" # Earliest time recorded at either end of all bars in the hit
    , coincsBranchName+".timeEnd" # Latest time recorded at either end of all bars in the hit
    , coincsBranchName+".time" # average reconstructed hit time of the cluster.
    , coincsBranchName+".PEs" # total number of photoelectrons in this cluser
    , coincsBranchName+".nHits" # Number of individual bar hits combined in this hit
    , coincsBranchName+".nLayers" # Number of CRV layers that are part of this cluster
    , coincsBranchName+".PEsPerLayer[4]"
    , coincsBranchName+".angle" # slope (in the plane perpendicular to the bar axis of a sector) of the track assumed to be responsible for the cluster (=change in the "layer direction" / change in the "thickness direction")

    # ---> crvhitmc (truth)
    , coincsBranchName+"mc.valid" # Records if there is a valid MC match to this CRV reco hit
    , coincsBranchName+"mc.pdgId" # PDG ID of the track mostly likely responsible for this cluster
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

extendedBranchNamesTrkAna_ = [

    # ---> evtinfo
    "evtinfo.runid" # run ID 
    ,"evtinfo.subrunid" # sub-run ID 
    ,"evtinfo.eventid" # event ID 
    
    # ---> crvhit (reco)
    , coincsBranchName+".sectorType" # CRV sector hit
    , coincsBranchName+".pos.fCoordinates.fX" # Reconstructed position of the cluster in X 
    , coincsBranchName+".pos.fCoordinates.fY" # Reconstructed position of the cluster in Y
    , coincsBranchName+".pos.fCoordinates.fZ" # Reconstructed position of the cluster in Z
    , coincsBranchName+".timeStart" # Earliest time recorded at either end of all bars in the hit
    , coincsBranchName+".timeEnd" # Latest time recorded at either end of all bars in the hit
    , coincsBranchName+".time" # average reconstructed hit time of the cluster.
    , coincsBranchName+".PEs" # total number of photoelectrons in this cluser
    , coincsBranchName+".nHits" # Number of individual bar hits combined in this hit
    , coincsBranchName+".nLayers" # Number of CRV layers that are part of this cluster
    , coincsBranchName+".PEsPerLayer[4]"
    , coincsBranchName+".angle" # slope (in the plane perpendicular to the bar axis of a sector) of the track assumed to be responsible for the cluster (=change in the "layer direction" / change in the "thickness direction")

    # ---> crvhitmc (truth)
    , coincsBranchName+"mc.valid" # Records if there is a valid MC match to this CRV reco hit
    , coincsBranchName+"mc.pdgId" # PDG ID of the track mostly likely responsible for this cluster
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

    # ----> kl (tracks)
    , "kl.status"
    , "kl.nactive"
    , "kl.nhits"
    , "kl.nplanes"
    , "kl.nnullambig"
    , "kl.ndof"
    # , "klfit" # arrays of structs
    # , "klkl" # arrays of structs
    # , "kl.fitcon"
    # , "klfit.sid"
    # , "klfit.sindex"
    # , "klfit.pos.X"
    # , "klfit.pos.Y"
    # , "klfit.pos.Z"
    # , "klfit.time"
    # , "klkl.z0err"
    # , "klkl.d0err"
    # , "klkl.thetaerr"
    # , "klkl.phi0err"
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

particleDict = {
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

def GetLatexParticleName(particle):

    if particle == "proton": return "$p$"
    elif particle == "pi+-": return "$\pi^{\pm}$"
    elif particle == "pi+": return "$\pi^{+}$"
    elif particle == "pi-": return "$\pi^{-}$"
    elif particle == "mu+-": return "$\mu^{\pm}$"
    elif particle == "mu+": return "$\mu^{+}$"
    elif particle == "mu-": return "$\mu^{-}$"
    elif particle == "e+": return "$e^{+}$"
    elif particle == "e-": return "$e^{-}$"
    # elif particle == "e+ >10 MeV": return "$e^{+} > 10$ MeV"
    # elif particle == "e- >10 MeV": return "$e^{-} > 10$ MeV"
    elif particle == "kaon+": return "$K^{+}$"
    elif particle == "kaon-": return "$K^{-}$"
    elif particle == "kaon0": return "$K^{0}$"
    elif particle == "kaon0L": return "$K^{0}_{L}$"
    elif particle == "kaon0S": return "$K^{0}_{S}$"
    # elif particle == "no_proton": return "No protons"
    # elif particle == "pi-_and_mu-": return "$\pi^{-}$ & $\mu^{-}$"
    # elif particle == "pi+_and_mu+": return "$\pi^{+}$ & $\mu^{+}$"
    # Add more as required
    else: return "other"

# get the latex names of the particles in the particle dictionary 
latexParticleDict = {}
for key, value in  particleDict.items():
    latexParticleDict[key] = GetLatexParticleName(value)

# --------
# Plotting
# --------

import matplotlib.pyplot as plt

def PlotGraph(x, y, title=None, xlabel=None, ylabel=None, fout="scatter.png", NDPI=300):

    # Create a scatter plot with error bars using NumPy arrays 

    # Create figure and axes
    fig, ax = plt.subplots()

    # Fine graphs for CRV visualisation
    ax.scatter(x, y, color='black', s=16, edgecolor='black', marker='o', linestyle='None')

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
    plt.close()

def PlotGraphErrors(x, xerr, y, yerr, title=None, xlabel=None, ylabel=None, fout="scatter.png", NDPI=300):

   # Create a scatter plot with error bars using NumPy arrays 

    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot scatter with error bars
    if len(xerr)==0: xerr = [0] * len(x) # Sometimes we only use yerr
    if len(yerr)==0: yerr = [0] * len(y) # Sometimes we only use yerr

    if len(x) != len(y): print("Warning: x has length", len(x),", while y has length", len(y))

    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='black', markersize=4, ecolor='black', capsize=2, elinewidth=1, linestyle='None')

    # Set title, xlabel, and ylabel
    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel(xlabel, fontsize=13, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=13, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=13)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=13)  # Set y-axis tick label font size

    if (ax.get_xlim()[1] > 9.999e3) or (ax.get_xlim()[1] < 9.999e-3) :
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(13)
    if (ax.get_ylim()[1] > 9.999e3) or (ax.get_ylim()[1] < 9.999e-3) :
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(13)


    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()


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

# graphs_ = { "label" : [ x_, xerr_, y_, yerr] }
def PlotGraphOverlay2(graphs_, title=None, xlabel=None, ylabel=None, labels_=[], fout="scatter.png", log=False, NDPI=300):
    
    # Create figure and axes
    fig, ax = plt.subplots()

    # Iterate over each pair of xy lists
    for i, (label, data_) in enumerate(graphs_.items()):

        x = data_[0] 
        xerr = data_[1]
        y = data_[2] 
        yerr = data_[3]

         # Plot scatter with error bars
        if len(xerr)==0: xerr = [0] * len(x) # Sometimes we only use yerr
        if len(yerr)==0: yerr = [0] * len(y) # Sometimes we only use yerr
        # Scatter plot each pair
        # ax.scatter(x, y, s=16, color=colours[i+1], edgecolor=colours[i+1], marker='o', linestyle='None', label=label)
        # ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=colours[i+1], markersize=4, ecolor=colours[i+1], capsize=2, elinewidth=1, linestyle='None',label=label)
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=colours[i+1], markersize=4, ecolor=colours[i+1], capsize=0, elinewidth=0, linestyle='None',label=label)

    if log: 
        # ax.set_xscale("log")
        ax.set_yscale("log")
        # ax.xaxis.set_major_formatter(ScalarFormatter())
        # ax.yaxis.set_major_formatter(ScalarFormatter())
        # ax.ticklabel_format(axis='y', style='plain')

    # Set title, xlabel, and ylabel
    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel(xlabel, fontsize=13, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=13, labelpad=10)

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=13)  
    ax.tick_params(axis='y', labelsize=13)  

    # Disable scientific notation
    # ax.ticklabel_format(useOffset=False)

    # Check if x or y values exceed 9999 for scientific notation
    # if any(max(data_[0]) > 9999 for _, data_ in graphs_.items()) or any(max(data_[0]) < 9.99e-4 for _, data_ in graphs_.items()):
    #     ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #     ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    #     ax.xaxis.offsetText.set_fontsize(13)
    # if any(max(data_[2]) > 9999 for _, data_ in graphs_.items()) or any(max(data_[2]) < 9.99e-4 for _, data_ in graphs_.items()):
    #     ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #     ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    #     ax.yaxis.offsetText.set_fontsize(13)

    ax.legend(loc="best", frameon=False, fontsize=13) # , markerscale=5)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

def Plot1D(data, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", legPos="best", stats=True, underOver=False, errors=False, NDPI=300):
    
    data = np.array(data)
    # data = np.flatten(data)

    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot the histogram with outline
    counts, bin_edges, _ = ax.hist(data, bins=nbins, range=(xmin, xmax), histtype='step', edgecolor='black', linewidth=1.0, fill=False, density=False)

    # Set x-axis limits
    ax.set_xlim(xmin, xmax)

    # # # Calculate statistics
    N, mean, meanErr, stdDev, stdDevErr, underflows, overflows = GetBasicStats(data, xmin, xmax)
    # # N, mean, meanErr, stdDev, stdDevErr = str(N), Round(mean, 3), Round(meanErr, 1), Round(stdDev, 3), Round(stdDevErr, 1) 

    # # Create legend text
    legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}\nStd Dev: {Round(stdDev, 3)}"
    # if errors: legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDevErr, 1)}"
    # if underOver: legend_text += f"\nUnderflows: {underflows}\nOverflows: {overflows}"

    # # legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDev, 1)}"

    # Add legend to the plot
    if stats: ax.legend([legend_text], loc=legPos, frameon=False, fontsize=13)

    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel(xlabel, fontsize=13, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=13, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=13)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=13)  # Set y-axis tick label font size

    if (ax.get_xlim()[1] > 9.999e3) or (ax.get_xlim()[1] < 9.999e-3):
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(13)
    if (ax.get_ylim()[1] > 9.999e3) or (ax.get_ylim()[1] < 9.999e-3):
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(13)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.close()

from scipy.optimize import curve_fit

# --------------------
# Fit function
# --------------------

# The Gaussian function
def gaussian(x, norm, mu, sigma):
    return norm * np.exp(-((x - mu) / (2 * sigma)) ** 2)

# Under development!
def Plot1DWithGaussFit(data, nbins=100, xmin=-1.0, xmax=1.0, norm=1.0, mu=0.0, sigma=1.0, fitMin=-1.0, fitMax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", legPos="best", stats=True, peak=False, underOver=False, errors=False, NDPI=300):
    
    data = np.array(data)
    
    # Create figure and axes
    fig, ax = plt.subplots()

    # Plot the histogram with outline
    counts, bin_edges, _ = ax.hist(data, bins=nbins, range=(xmin, xmax), histtype='step', edgecolor='black', linewidth=1.0, fill=False, density=False)

    # Set x-axis limits
    ax.set_xlim(xmin, xmax)

    # Fit gaussian

    # Calculate bin centers
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    # Fit the Gaussian function to the histogram data
    params, covariance = curve_fit(gaussian, bin_centres[(bin_centres >= fitMin) & (bin_centres <= fitMax)], counts, p0=[norm, mu, sigma])
    # Extract parameters from the fitting
    norm, mu, sigma = params
    # Plot the Gaussian curve
    ax.plot(bin_centres, gaussian(bin_centres, norm, mu, sigma), color="red", label=f"Norm: {norm}\n$\mu$: {mu}\n$sigma {sigma}")

    # Calculate statistics
    N, mean, meanErr, stdDev, stdDevErr, underflows, overflows = GetBasicStats(data, xmin, xmax)

    # Create legend text
    legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}\nStd Dev: {Round(stdDev, 3)}"
    # if errors: legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDevErr, 1)}"
    # if errors: legend_text = f"Entries: {N}\nMean: {Round(mean, 4)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 4)}$\pm${Round(stdDevErr, 1)}"
    # # if peak: legend_text += f"\nPeak: {Round(peak, 4)}$\pm${Round(peakErr, 1)}"
    # if underOver: legend_text += f"\nUnderflows: {underflows}\nOverflows: {overflows}"

    # legend_text = f"Entries: {N}\nMean: {Round(mean, 3)}$\pm${Round(meanErr, 1)}\nStd Dev: {Round(stdDev, 3)}$\pm${Round(stdDev, 1)}"

    # Add legend to the plot
    if stats: ax.legend([legend_text], loc=legPos, frameon=False, fontsize=13)

    ax.set_title(title, fontsize=15, pad=10)
    ax.set_xlabel(xlabel, fontsize=13, labelpad=10) 
    ax.set_ylabel(ylabel, fontsize=13, labelpad=10) 

    # Set font size of tick labels on x and y axes
    ax.tick_params(axis='x', labelsize=13)  # Set x-axis tick label font size
    ax.tick_params(axis='y', labelsize=13)  # Set y-axis tick label font size

    if (ax.get_xlim()[1] > 9.999e3) or (ax.get_xlim()[1] < 9.999e-3):
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(13)
    if (ax.get_ylim()[1] > 9.999e3) or (ax.get_ylim()[1] < 9.999e-3):
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(13)

    # Save the figure
    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("---> Written", fout)

    # Clear memory
    plt.clf()
    plt.close()


import matplotlib.colors as colors

def Plot2D(x, y, nbinsX=100, xmin=-1.0, xmax=1.0, nbinsY=100, ymin=-1.0, ymax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", log=False, cb=True, NDPI=300):

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
    hist, x_edges, y_edges = np.histogram2d(x, y, bins=[nbinsX, nbinsY], range=[[xmin, xmax], [ymin, ymax]])

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
    if (ax.get_xlim()[1] > 9.999e3) or (ax.get_xlim()[1] < 9.999e-3):
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.offsetText.set_fontsize(14)
    if (ax.get_ylim()[1] > 9.999e3) or (ax.get_ylim()[1] < 9.999e-3):
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.offsetText.set_fontsize(14)

    plt.savefig(fout, dpi=NDPI, bbox_inches="tight")
    print("\n---> Written:\n\t", fout)

    # Clear memory
    plt.close()

def Plot1DOverlay(hists_, nbins=100, xmin=-1.0, xmax=1.0, title=None, xlabel=None, ylabel=None, fout="hist.png", label_=None, legPos="best", NDPI=300, includeBlack=False, logY=False, legFontSize=12):

    # Create figure and axes
    fig, ax = plt.subplots()

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

    # # Accumulate label counts over all datasets
    # total_label_counts = np.zeros(len(unique_labels))
    # for dataset in data_:
    #     labels = [label_dict.get(p, 'other') for p in dataset]
    #     unique_labels_data, label_counts_data = np.unique(labels, return_counts=True)
    #     for i, label in enumerate(unique_labels_data):
    #         index = np.where(unique_labels == label)[0][0]
    #         total_label_counts[index] += label_counts_data[i]

    # # Loop through each dataset
    # for i, dataset in enumerate(data_):

    #     labels = [label_dict.get(p, 'other') for p in dataset]
    #     unique_labels_data, label_counts_data = np.unique(labels, return_counts=True)

    #     # Reorder label counts based on unique_labels
    #     label_counts = [label_counts_data[np.where(unique_labels_data == label)[0][0]] if label in unique_labels_data else 0 for label in unique_labels]

    #     if percentage:
    #         label_counts_percentage = [(count / sum(total_label_counts)) * 100 for count in label_counts]  # Calculate percentage for each label using the total count over all datasets
    #         label_counts = label_counts_percentage

    # Loop through each dataset
    for i, dataset in enumerate(data_):
    
        labels = [label_dict.get(p, 'other') for p in dataset]
        unique_labels_data, label_counts_data = np.unique(labels, return_counts=True)

        # Reorder label counts based on unique_labels
        label_counts = [label_counts_data[unique_labels_data == label][0] if label in unique_labels_data else 0 for label in unique_labels]

        # if percentage:
        #     total_count = sum(label_counts_data)  # Total count of all labels
        #     label_counts_percentage = [(count / total_count) * 100 for count in label_counts]  # Calculate percentage for each label
        #     label_counts = label_counts_percentage

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


