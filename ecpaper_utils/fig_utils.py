import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap ## used to create custom colormaps

## defines colors and plotting order that are used throughout analysis

ols_col = [124/255, 79/255, 9/255, 1]
tls_col = [20/255, 120/255, 112/255, 1]
bhm_col = [73/255, 18/255, 133/255, 1]

#cmip_col = [[254/255, 10/255, 2/255, 1], [0/255, 88/255, 153/255, 1]]
cmip_col = ["firebrick", "royalblue"]

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
    axs.set_xlabel("r(Ni$\~{n}$o3.4, Cal pr)", fontsize=14)
    axs.set_ylabel("Cal pr trend (mm/day/100y)", fontsize=14)
    axs.set_title(titlestr, fontsize=16)
    axs.set_xlim(-0.5,1)
    axs.set_ylim(-1.5,3)
    axs.set_xticks([-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
    axs.set_xticklabels(['-0.4','-0.2','0','0.2','0.4','0.6','0.8','1'], fontsize=14)
    axs.set_yticks([-1,0,1,2,3])
    axs.set_yticklabels(['-1','0','1','2','3'], fontsize=14)
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

    return axs

def plotconstraint_sensitivity(axs, ols, ols_noxobsvar, ols_nocoefvar, ols_nodelta,
                               tls, tls_noxobsvar, tls_nocoefvar, tls_nodelta, 
                               bhm, bhm_noxobsvar, bhm_nocoefvar, bhm_nodelta):
    """ Plotting the ranges that demonstrate the sensitivity of the constraint to the
        ommition of each uncertainty.  Used in supplemental figure 5 """

    axs.bar(2,ols['max95withiv']-ols['min95withiv'], bottom=ols['min95withiv'], color='saddlebrown', width=1)
    axs.plot([1.5,2.5],[ols['min66withiv'],ols['min66withiv']], color='white', linewidth=2)
    axs.plot([1.5,2.5],[ols['max66withiv'],ols['max66withiv']], color='white', linewidth=2)
    axs.plot([1.5,2.5],[ols['meanwithiv'],ols['meanwithiv']], color='black', linewidth=3)

    axs.bar(4,ols['max95forced']-ols['min95forced'], bottom=ols['min95forced'], color='saddlebrown', width=1)
    axs.plot([3.5,4.5],[ols['min66forced'],ols['min66forced']], color='white', linewidth=2)
    axs.plot([3.5,4.5],[ols['max66forced'],ols['max66forced']], color='white', linewidth=2)
    axs.plot([3.5,4.5],[ols['meanforced'],ols['meanforced']], color='black', linewidth=3)

    axs.bar(6,ols_noxobsvar['max95withiv']-ols_noxobsvar['min95withiv'], bottom=ols_noxobsvar['min95withiv'], color='saddlebrown', width=1)
    axs.plot([5.5,6.5],[ols_noxobsvar['min66withiv'],ols_noxobsvar['min66withiv']], color='white', linewidth=2)
    axs.plot([5.5,6.5],[ols_noxobsvar['max66withiv'],ols_noxobsvar['max66withiv']], color='white', linewidth=2)
    axs.plot([5.5,6.5],[ols_noxobsvar['meanwithiv'],ols_noxobsvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(8,ols_nocoefvar['max95withiv']-ols_nocoefvar['min95withiv'], bottom=ols_nocoefvar['min95withiv'], color='saddlebrown', width=1)
    axs.plot([7.5,8.5],[ols_nocoefvar['min66withiv'],ols_nocoefvar['min66withiv']], color='white', linewidth=2)
    axs.plot([7.5,8.5],[ols_nocoefvar['max66withiv'],ols_nocoefvar['max66withiv']], color='white', linewidth=2)
    axs.plot([7.5,8.5],[ols_nocoefvar['meanwithiv'],ols_nocoefvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(10,ols_nodelta['max95withiv']-ols_nodelta['min95withiv'], bottom=ols_nodelta['min95withiv'], color='saddlebrown', width=1)
    axs.plot([9.5,10.5],[ols_nodelta['min66withiv'],ols_nodelta['min66withiv']], color='white', linewidth=2)
    axs.plot([9.5,10.5],[ols_nodelta['max66withiv'],ols_nodelta['max66withiv']], color='white', linewidth=2)
    axs.plot([9.5,10.5],[ols_nodelta['meanwithiv'],ols_nodelta['meanwithiv']], color='black', linewidth=3)


    axs.bar(13,tls['max95withiv']-tls['min95withiv'], bottom=tls['min95withiv'], color='forestgreen', width=1)
    axs.plot([12.5,13.5],[tls['min66withiv'],tls['min66withiv']], color='white', linewidth=2)
    axs.plot([12.5,13.5],[tls['max66withiv'],tls['max66withiv']], color='white', linewidth=2)
    axs.plot([12.5,13.5],[tls['meanwithiv'],tls['meanwithiv']], color='black', linewidth=3)

    axs.bar(15,tls['max95forced']-tls['min95forced'], bottom=tls['min95forced'], color='forestgreen', width=1)
    axs.plot([14.5,15.5],[tls['min66forced'],tls['min66forced']], color='white', linewidth=2)
    axs.plot([14.5,15.5],[tls['max66forced'],tls['max66forced']], color='white', linewidth=2)
    axs.plot([14.5,15.5],[tls['meanforced'],tls['meanforced']], color='black', linewidth=3)

    axs.bar(17,tls_noxobsvar['max95withiv']-tls_noxobsvar['min95withiv'], bottom=tls_noxobsvar['min95withiv'], color='forestgreen', width=1)
    axs.plot([16.5,17.5],[tls_noxobsvar['min66withiv'],tls_noxobsvar['min66withiv']], color='white', linewidth=2)
    axs.plot([16.5,17.5],[tls_noxobsvar['max66withiv'],tls_noxobsvar['max66withiv']], color='white', linewidth=2)
    axs.plot([16.5,17.5],[tls_noxobsvar['meanwithiv'],tls_noxobsvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(19,tls_nocoefvar['max95withiv']-tls_nocoefvar['min95withiv'], bottom=tls_nocoefvar['min95withiv'], color='forestgreen', width=1)
    axs.plot([18.5,19.5],[tls_nocoefvar['min66withiv'],tls_nocoefvar['min66withiv']], color='white', linewidth=2)
    axs.plot([18.5,19.5],[tls_nocoefvar['max66withiv'],tls_nocoefvar['max66withiv']], color='white', linewidth=2)
    axs.plot([18.5,19.5],[tls_nocoefvar['meanwithiv'],tls_nocoefvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(21,tls_nodelta['max95withiv']-tls_nodelta['min95withiv'], bottom=tls_nodelta['min95withiv'], color='forestgreen', width=1)
    axs.plot([20.5,21.5],[tls_nodelta['min66withiv'],tls_nodelta['min66withiv']], color='white', linewidth=2)
    axs.plot([20.5,21.5],[tls_nodelta['max66withiv'],tls_nodelta['max66withiv']], color='white', linewidth=2)
    axs.plot([20.5,21.5],[tls_nodelta['meanwithiv'],tls_nodelta['meanwithiv']], color='black', linewidth=3)

    axs.bar(24,bhm['max95withiv']-bhm['min95withiv'], bottom=bhm['min95withiv'], color='blueviolet', width=1)
    axs.plot([23.5,24.5],[bhm['min66withiv'],bhm['min66withiv']], color='white', linewidth=2)
    axs.plot([23.5,24.5],[bhm['max66withiv'],bhm['max66withiv']], color='white', linewidth=2)
    axs.plot([23.5,24.5],[bhm['meanwithiv'],bhm['meanwithiv']], color='black', linewidth=3)

    axs.bar(26,bhm['max95forced']-bhm['min95forced'], bottom=bhm['min95forced'], color='blueviolet', width=1)
    axs.plot([25.5,26.5],[bhm['min66forced'],bhm['min66forced']], color='white', linewidth=2)
    axs.plot([25.5,26.5],[bhm['max66forced'],bhm['max66forced']], color='white', linewidth=2)
    axs.plot([25.5,26.5],[bhm['meanforced'],bhm['meanforced']], color='black', linewidth=3)

    axs.bar(28,bhm_noxobsvar['max95withiv']-bhm_noxobsvar['min95withiv'], bottom=bhm_noxobsvar['min95withiv'], color='blueviolet', width=1)
    axs.plot([27.5,28.5],[bhm_noxobsvar['min66withiv'],bhm_noxobsvar['min66withiv']], color='white', linewidth=2)
    axs.plot([27.5,28.5],[bhm_noxobsvar['max66withiv'],bhm_noxobsvar['max66withiv']], color='white', linewidth=2)
    axs.plot([27.5,28.5],[bhm_noxobsvar['meanwithiv'],bhm_noxobsvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(30,bhm_nocoefvar['max95withiv']-bhm_nocoefvar['min95withiv'], bottom=bhm_nocoefvar['min95withiv'], color='blueviolet', width=1)
    axs.plot([29.5,30.5],[bhm_nocoefvar['min66withiv'],bhm_nocoefvar['min66withiv']], color='white', linewidth=2)
    axs.plot([29.5,30.5],[bhm_nocoefvar['max66withiv'],bhm_nocoefvar['max66withiv']], color='white', linewidth=2)
    axs.plot([29.5,30.5],[bhm_nocoefvar['meanwithiv'],bhm_nocoefvar['meanwithiv']], color='black', linewidth=3)

    axs.bar(32,bhm_nodelta['max95withiv']-bhm_nodelta['min95withiv'], bottom=bhm_nodelta['min95withiv'], color='blueviolet', width=1)
    axs.plot([31.5,32.5],[bhm_nodelta['min66withiv'],bhm_nodelta['min66withiv']], color='white', linewidth=2)
    axs.plot([31.5,32.5],[bhm_nodelta['max66withiv'],bhm_nodelta['max66withiv']], color='white', linewidth=2)
    axs.plot([31.5,32.5],[bhm_nodelta['meanwithiv'],bhm_nodelta['meanwithiv']], color='black', linewidth=3)
















    
    return axs









