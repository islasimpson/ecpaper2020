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

def vwindscattersetup(axs,titlestr):
    """ Setting up axes for the delta v_sw vs psi scatter plots """
    axs.set_xlabel("$|\psi|$ (ms$^{-1}$)", fontsize=14)
    axs.set_ylabel("$\Delta V_{SW}$ (ms$^{-1}$)", fontsize=14)
    axs.set_title(titlestr, fontsize=16)
    axs.set_xlim(1,6)
    axs.set_ylim(-10,2)
    axs.set_xticks([1,2,3,4,5,6])
    axs.set_xticklabels(['1','2','3','4','5','6'], fontsize=14)
    axs.set_yticks([-10,-8,-6,-4,-2,0,2])
    axs.set_yticklabels(['-10','-8','-6','-4','-2','0','2'], fontsize=14)
    return axs

def calpscattersetup(axs, titlestr):
    """ Setting up axes for the CALP statter plots """
    axs.set_xlabel("r(Ni\~{n}o3.4, Cal pr)")
    axs.set_ylabel("Cal pr trend (mm/day/100y)")
    axs.set_title(titlestr, fontsize=16)
    axs.set_xlim(-0.5,1)
    axs.set_ylim(-1.5,3)
    axs.set_xticks([-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
    axs.set_xticklabels(['-0.4','-0.2','0','0.2','0.4','0.6','0.8','1'])
    axs.set_yticks([-1,0,1,2,3])
    axs.set_yticklabels(['-1','0','1','2','3'])
    return axs


def plotconstraintinfo(axs, olsdata, tlsdata, bhmdata, ylim):
    """Plot the constraint and associated text.  
    """
    axs.text(6, ylim[0] + 0.4*(ylim[1]-ylim[0]) ,'Forced', horizontalalignment='center', verticalalignment='center', fontsize=14, color='black')
    axs.text(13, ylim[0] + 0.4*(ylim[1]-ylim[0]) ,'Forced', horizontalalignment='center', verticalalignment='center', fontsize=14, color='black')
    axs.text(13,ylim[0] + 0.35*(ylim[1]-ylim[0]),'+', horizontalalignment='center', verticalalignment='center', fontsize=14, color='black')
    axs.text(13,ylim[0] + 0.3*(ylim[1]-ylim[0]),'Internal', horizontalalignment='center', verticalalignment='center', fontsize=14, color='black')

    meanforced = (olsdata['meanforced'] + tlsdata['meanforced'] + bhmdata['meanforced'])/3.
    meanwithiv = (olsdata['meanwithiv'] + tlsdata['meanwithiv'] + bhmdata['meanwithiv'])/3.
    min95forced = (olsdata['min95forced'] + tlsdata['min95forced'] + bhmdata['min95forced'])/3.
    max95forced = (olsdata['max95forced'] + tlsdata['max95forced'] + bhmdata['max95forced'])/3.
    min95withiv = (olsdata['min95withiv'] + tlsdata['min95withiv'] + bhmdata['min95withiv'])/3.
    max95withiv = (olsdata['max95withiv'] + tlsdata['max95withiv'] + bhmdata['max95withiv'])/3.
    meangtforced = (olsdata['gtymean_forced'] + tlsdata['gtymean_forced'] + bhmdata['gtymean_forced'])/3.
    meangtwithiv = (olsdata['gtymean_withiv'] + tlsdata['gtymean_withiv'] + bhmdata['gtymean_withiv'])/3.

    meanforcedstr = '{0:6.2f}'.format(np.array(meanforced)).strip()
    meanwithivstr = '{0:6.2f}'.format(np.array(meanwithiv)).strip()

    axs.text(6,ylim[0] + 0.16*(ylim[1]-ylim[0]),'$\Delta \phi$='+meanforcedstr, horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
    axs.text(13,ylim[0] + 0.16*(ylim[1]-ylim[0]),'$\Delta \phi$='+meanwithivstr, horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
    min95forcedstr = '{0:6.2f}'.format(min95forced).strip()
    max95forcedstr = '{0:6.2f}'.format(max95forced).strip()
    min95withivstr = '{0:6.2f}'.format(min95withiv).strip()
    max95withivstr = '{0:6.2f}'.format(max95withiv).strip()
    axs.text(6,ylim[0] + 0.1*(ylim[1]-ylim[0]),'ci=('+min95forcedstr+','+max95forcedstr+')', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
    axs.text(13,ylim[0] + 0.1*(ylim[1]-ylim[0]),'ci=('+min95withivstr+','+max95withivstr+')', horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
    gtforcedstr = '{0:6.2f}'.format(meangtforced).strip()
    gtwithivstr = '{0:6.2f}'.format(meangtwithiv).strip()
    axs.text(6,ylim[0] + 0.04*(ylim[1]-ylim[0]),gtforcedstr+"%>CMIP5", horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')
    axs.text(13,ylim[0] + 0.04*(ylim[1]-ylim[0]),gtwithivstr+"%>CMIP5", horizontalalignment='center', verticalalignment='center', fontsize=10, color='black')

    axs.bar(4,olsdata['max95forced']-olsdata['min95forced'],bottom=olsdata['min95forced'], color='saddlebrown', width=1, zorder=2)
    axs.plot([4-0.5,4+0.5],[olsdata['min66forced'],olsdata['min66forced']], color='white')
    axs.plot([4-0.5,4+0.5],[olsdata['max66forced'],olsdata['max66forced']], color='white')
    axs.plot([4-0.5,4+0.5],[olsdata['meanforced'],olsdata['meanforced']], color='black')
    
    axs.bar(6,tlsdata['max95forced']-tlsdata['min95forced'],bottom=tlsdata['min95forced'], color='forestgreen', width=1, zorder=2)
    axs.plot([6-0.5,6+0.5],[tlsdata['min66forced'],tlsdata['min66forced']], color='white')
    axs.plot([6-0.5,6+0.5],[tlsdata['max66forced'],tlsdata['max66forced']], color='white')
    axs.plot([6-0.5,6+0.5],[tlsdata['meanforced'],tlsdata['meanforced']], color='black')
    
    axs.bar(8,bhmdata['max95forced']-bhmdata['min95forced'],bottom=bhmdata['min95forced'], color='blueviolet', width=1, zorder=2)
    axs.plot([8-0.5,8+0.5],[bhmdata['min66forced'],bhmdata['min66forced']], color='white')
    axs.plot([8-0.5,8+0.5],[bhmdata['max66forced'],bhmdata['max66forced']], color='white')
    axs.plot([8-0.5,8+0.5],[bhmdata['meanforced'],bhmdata['meanforced']], color='black')
    
    axs.plot([9.5,9.5],[ylim[0],ylim[1]], color='black')
    
    axs.bar(11,olsdata['max95withiv']-olsdata['min95withiv'],bottom=olsdata['min95withiv'], color='saddlebrown', width=1, zorder=2)
    axs.plot([11-0.5,11+0.5],[olsdata['min66withiv'],olsdata['min66withiv']], color='white')
    axs.plot([11-0.5,11+0.5],[olsdata['max66withiv'],olsdata['max66withiv']], color='white')
    axs.plot([11-0.5,11+0.5],[olsdata['meanwithiv'],olsdata['meanwithiv']], color='black')
    
    axs.bar(13,tlsdata['max95withiv']-tlsdata['min95withiv'],bottom=tlsdata['min95withiv'], color='forestgreen', width=1, zorder=2)
    axs.plot([13-0.5,13+0.5],[tlsdata['min66withiv'],tlsdata['min66withiv']], color='white')
    axs.plot([13-0.5,13+0.5],[tlsdata['max66withiv'],tlsdata['max66withiv']], color='white')
    axs.plot([13-0.5,13+0.5],[tlsdata['meanwithiv'],tlsdata['meanwithiv']], color='black')
    
    axs.bar(15,bhmdata['max95withiv']-bhmdata['min95withiv'],bottom=bhmdata['min95withiv'], color='blueviolet', width=1, zorder=2)
    axs.plot([15-0.5,15+0.5],[bhmdata['min66withiv'],bhmdata['min66withiv']], color='white')
    axs.plot([15-0.5,15+0.5],[bhmdata['max66withiv'],bhmdata['max66withiv']], color='white')
    axs.plot([15-0.5,15+0.5],[bhmdata['meanwithiv'],bhmdata['meanwithiv']], color='black')












