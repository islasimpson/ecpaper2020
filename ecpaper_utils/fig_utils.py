import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap ## used to create custom colormaps

## defines colors and plotting order that are used throughout analysis

ols_col = [124/255, 79/255, 9/255, 1]
tls_col = [20/255, 120/255, 112/255, 1]
bhm_col = [73/255, 18/255, 133/255, 1]

cmip_col = [[254/255, 10/255, 2/255, 1], [0/255, 88/255, 153/255, 1]]

le_order = ["CanESM2", "CESM1-CAM5", "CSIRO-Mk3-6-0", "GFDL-CM3", "MPI-ESM"]
le_col = [[80/255, 179/255, 101/255, 1], 
         [136/255, 133/255, 190/255, 1], 
         [78/255, 154/255, 203/255, 1], 
         [240/255, 71/255, 52/255, 1], 
         [203/255, 58/255, 148/255, 1]]

## defines functions used for figure formatting

def combine_cmaps(n, white = 2, neg = "Blues", pos = "YlOrRd", start = 0.1):
    """ combine two existing cmaps to create a diverging cmap with white in the middle
    n = number of colors on each side
    white = number of white spots
    neg = name of cmap for negative side
    pos = name of cmap for positive side
    start = minimum value of existing cmap to use"""
    
    new_col = ListedColormap(np.concatenate(
    [plt.cm.get_cmap(neg)(np.flip(np.linspace(start = start, stop = 1, num = n))), 
    np.ones((white, 4)), 
    plt.cm.get_cmap(pos)(np.linspace(start = start, stop = 1, num = n))]))
    
    return new_col

def jlatscattersetup(axs,titlestr):
    """ Setting up axes for the jet latitude scatter plots """
    axs.set_xlabel("Jet latitude, $\phi_{o}$ ($^{\circ}$N)", fontsize=14)
    axs.set_ylabel("Jet shift, $\Delta \phi$, ($^{\circ}$N)", fontsize=14)
    axs.set_title(titlestr, fontsize=16)
    axs.set_xlim(-57,-33)
    axs.set_ylim(-9,3)
    axs.set_xticks([-55,-50,-45,-40,-35])
    axs.set_xticklabels(['-55','-50','-45','-40','-35'], fontsize=14)
    axs.set_yticks([-8,-6,-4,-2,0,2])
    axs.set_yticklabels(['-8','-6','-4','-2','0','2'], fontsize=14)
    return axs
