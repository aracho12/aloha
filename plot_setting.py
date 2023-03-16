from cycler import cycler
import matplotlib

linewidth=0.4
majorticksize=3.5
minorticksize=2.0
labelsize=6
fontsize=7

colors=[
    '#D32F2F', #dark red
    '#303F9F', #indigo
    '#eecc16', #yellow
    '#008176', #teal
    "#7847a2", #purple
    '#FF4081', #pink
    '#455A64', #blue grey
    '#b3b3b3' #grey
]

preamble  = r"""
              \usepackage{color}
              \usepackage[tx]{sfmath}
              \usepackage{helvet}
              \usepackage{sansmath}
           """

options={
    'font.size' : fontsize,

    'axes.labelsize': labelsize, #'medium'
    'axes.linewidth': linewidth, #0.8 ; 1.2
    'axes.prop_cycle': cycler('color',colors ),
    #'default': cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']),
    #'nature':['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2FF', '#DC0000', '#7E6148', '#B09C85']
    'axes.titlesize': fontsize, #large

    'boxplot.boxprops.color': '9ae1f9', #'black'
    'boxplot.boxprops.linewidth': 0, #1.0
    'boxplot.capprops.color': 'C0', #'black'
    'boxplot.flierprops.color': 'C0', #'black'
    'boxplot.flierprops.markeredgecolor': 'white', #'black'
    'boxplot.flierprops.markerfacecolor': 'auto', #'none'
    'boxplot.flierprops.markersize': 7, #6.0
    'boxplot.meanprops.color': 'C1', #'C2'
    'boxplot.meanprops.markeredgecolor': 'C1', #'C2'
    'boxplot.meanprops.markerfacecolor': 'C1', #'C2'
    'boxplot.meanprops.markersize': 7, #6.0
    'boxplot.medianprops.color': '9ae1f9',
    'boxplot.patchartist': True, #False
    'boxplot.whiskerprops.color': 'C0', #'black'

    'savefig.bbox':'tight',
    'savefig.pad_inches':0.01,

    'figure.figsize' : (3.5, 2.8),
    'figure.dpi': 300, #250
    'figure.edgecolor': 'white', #(1, 1, 1, 0)
    'figure.facecolor': 'white', #(1, 1, 1, 0)
    'figure.frameon': True, #True

    'font.family': ['sans-serif'], #['sans-serif']
    'font.sans-serif': ['Avenir','Helvetica','Arial'], 
    #'Helvetica, Computer Modern Sans Serif, DejaVu Sans, 'Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Arial, Avant Garde, sans-serif', 

    'text.latex.preamble': preamble,

    'grid.alpha': 0.2, #1.0
    'grid.color': '93939c', # #b0b0b0
    'grid.linestyle': '--', #'-'
    'grid.linewidth':linewidth,

    'legend.edgecolor': '0.8', 
    'legend.fancybox': False, #True
    'legend.frameon': False, #True
    'legend.fontsize':fontsize,

    'lines.markeredgecolor': 'black', #'auto'
    'lines.markersize':3,
    'lines.linewidth':1,

    'patch.linewidth':linewidth,
    'patch.facecolor': '000000', #'C0'
    'patch.force_edgecolor': True, #False
    'scatter.edgecolors': '000000', #'face'

    'xtick.labelsize': labelsize, #'medium'
    'xtick.direction': 'out',
    'xtick.minor.visible': False,
    'xtick.major.size': majorticksize, #3.5
    'xtick.minor.size': minorticksize,
    'xtick.major.top': False, #True
    'xtick.major.width': linewidth, #0.8
    'xtick.minor.width': linewidth,
    'xtick.minor.bottom': True, #True
    'xtick.minor.top': True, #True
    'xtick.alignment': 'center',

    'ytick.labelsize': labelsize, #'medium'
    'ytick.minor.visible' : False,
    'ytick.major.right': True, #True
    'ytick.major.size': majorticksize, #3.5  distance to the minor tick label in points
    'ytick.minor.size': minorticksize,
    'ytick.major.width': linewidth, #0.8
    'ytick.minor.width':linewidth,
    'ytick.minor.pad':3.4, #3.4
    'ytick.minor.right': True #True
} 

matplotlib.rcParams.update(options)

from matplotlib.ticker import AutoMinorLocator

def pretty_plot(width=None, height=None, plt=None, dpi=None):
    if plt is None:
        plt = matplotlib.pyplot
        if width is None:
            width = matplotlib.rcParams["figure.figsize"][0]
        if height is None:
            height = matplotlib.rcParams["figure.figsize"][1]
        if dpi is not None:
            matplotlib.rcParams["figure.dpi"] = dpi

        fig = plt.figure(figsize=(width, height))
        ax = fig.add_subplot(1, 1, 1)
        matplotlib.rcParams["xtick.minor.visible"]=True
        matplotlib.rcParams["ytick.minor.visible"]=True
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2)) 
    return plt