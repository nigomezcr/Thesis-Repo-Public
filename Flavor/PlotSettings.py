import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd

#Global Plot Setting

sns.set_style('ticks') # darkgrid, white grid, dark, white and ticks
plt.rc('axes', titlesize=18)     # fontsize of the axes title
plt.rc('axes', labelsize=16)     # fontsize of the x and y labels
plt.rc('axes', linewidth=1.65 )   # width of the frame
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize
plt.rc('font', size=14)          # controls default text sizes

#Color palette for the thesis

RegionAlpha = 0.5
MainColor1 = (0.580,0.706,0.231) # Green (148, 180, 59)
MainColor2 = (0.650,0.110,0.192) # Red (166, 28, 49)
BackgroundColor1 = (0.275,0.419,0.247) # Dark Green (70, 106, 63)
BackgroundColor2 = (0.467,0.137,0.184) # Dark Red (119, 35, 47)
Gray1 = (0.337,0.352,0.360) # (85, 90, 92)
Gray2 = (0.694,0.698,0.690) # (177, 178, 176)