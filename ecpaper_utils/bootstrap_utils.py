import xarray as xr
import numpy as np
import sys
from ecpaper_utils import linfit_utils as linfit

def bootgen_multimem(darray, nboots, nmems, seed=None):
    """ Generate nboots bootstrap samples from darray with nmems members for each sample

    Input: darray = an xarray data array with the sampling being performed on the first dimension
           nboots = the number of bootstrap samples
           nmems = the number of members in each bootstrap sample

    Output: bootdatxr = an xarray data array containing the bootstrap samples
            with dimensions (nboots, all the other coords of darray except the first)

    Option: a seed if you want to specify the seed for the random number generator 
    """

    # set up the dimensions and coordinates of the bootstrap array
    dims = darray.dims
    dimboot = [nmems*nboots]
    dimboot2d = [nmems, nboots]
    bootcoords = [("iboot", np.arange(0,nboots,1))]
   
    for icoord in range(1,len(dims)):
        dimboot.append(darray[dims[icoord]].size)
        dimboot2d.append(darray[dims[icoord]].size)
        bootcoords.append( (dims[icoord], darray[dims[icoord]] ))

    # generate random numbers for bootstrapping
    if (seed):
        np.random.seed(seed)

    nmemin = darray[dims[0]].size
    ranu = np.random.uniform(0,nmemin,nboots*nmems)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros(dimboot)
    bootdat = np.array(darray[ranu])
    bootdat = bootdat.reshape(dimboot2d)
    bootdatm = np.mean(bootdat, axis=0)
    bootdatxr = xr.DataArray(bootdatm, coords=bootcoords)


    return bootdatxr


def boot_corr_multimem(darrays, nboots, nmems, seed=None):
    """ Generate nboots bootstrap samples, with nmems members for each sample, from two darrays of the same coords. 
        Calculate correlation coefficient of the two variables across the nmems members for each sample.   

    Input: darrays = list of xarray data arrays with the sampling being performed on the first dimension
           nboots = the number of bootstrap samples
           nmems = the number of members in each bootstrap sample

    Output: bootdat = an xarray data array containing the bootstrap samples
            with dimensions nboots

    Option: a seed if you want to specify the seed for the random number generator 
    """

    dims = darrays[0].dims
    bootcoords = [("iboot", np.arange(0,nboots,1))]

    # generate random numbers for bootstrapping
    if (seed):
        np.random.seed(seed)

    nmemin = darrays[0][dims[0]].size
    ranu = np.random.uniform(0,nmemin,nboots*nmems)
    ranu = np.floor(ranu).astype(int)

    ## select ensemble members and reshape into long time series for each bootstrap sample
    bootdat_var1 = np.array(darrays[0][ranu]).reshape(nboots, -1)
    bootdat_var2 = np.array(darrays[1][ranu]).reshape(nboots, -1)

    ## calculate correlation for each bootstrap
    bootdat_corr = xr.corr(xr.DataArray(bootdat_var1), xr.DataArray(bootdat_var2), dim = "dim_1")
    
    bootdat = xr.DataArray(bootdat_corr, coords=bootcoords)

    return bootdat

def boot_corr_ci(a1,a2,conf, nboots=1000):
    """ Output the conf% confidence interval on correlation between 
    two 1 dimensional arrays by bootstrapping with replacement

    Input:
        a1 = first array
        a2 = second array
        nboots = the number of bootstrap samples used to generate the ci
        conf = the confidence interval you want e.g., 95 for 95% ci

    Output:
        minci = the minimum range of the confidence interval
        maxci = the maximum range of the confidence interval
 
    This assumes a two sided test.
    """

    ptilemin = (100.-conf)/2.
    ptilemax = conf + (100-conf)/2.

    if (a1.size != a2.size):
        print("The two arrays must have the same size")
        sys.exit()

    samplesize = a1.size
    ranu = np.random.uniform(0,samplesize,nboots*samplesize)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros([samplesize,nboots])
    bootdat1 = np.array(a1[ranu])
    bootdat2 = np.array(a2[ranu])
    bootdat1 = bootdat1.reshape([samplesize,nboots])
    bootdat2 = bootdat2.reshape([samplesize,nboots])

    bootcor = xr.corr(xr.DataArray(bootdat1),xr.DataArray(bootdat2), dim='dim_0')
    minci = np.percentile(bootcor,ptilemin)
    maxci = np.percentile(bootcor,ptilemax)

    return minci, maxci

def boot_regcoef_ci(a1,a2,conf,sigx=None,sigy=None,nboots=1000):
    """ Output the conf% confidence interval on regression coefficients between 
    two 1 dimensional arrays by bootstrapping with replacement

    Input:
        a1 = first array
        a2 = second array
        conf = the confidence interval you want e.g., 95 for 95% ci (2-sided)
    Optional input:
        nboots = the number of bootstrap samples used to generate the ci
        sigx = the standard deviation on the predictor points
        sigy = the standard deviation on the predictand points

    Output:
        aminci = the minimum range of the confidence interval on a
        amaxci = the maximum range of the confidence interval on a
        bminci = the minimum range of the confidence interval on b 
        bmaxci = the maximum range of the confidence interval on b
    
    where y = a + bx
 
    This assumes a two sided test.  Different regression methods are used 
    depending on the availability of sigx or sigy
    if no sigx then ordinary least squares regression
    if sigx and sigy then total least squares regression
    """

    ptilemin = (100.-conf)/2.
    ptilemax = conf + (100-conf)/2.

    if (a1.size != a2.size):
        print("The two arrays must have the same size")
        sys.exit()

    samplesize = a1.size
    ranu = np.random.uniform(0,samplesize,nboots*samplesize)
    ranu = np.floor(ranu).astype(int)

    bootdat = np.zeros([samplesize,nboots])
    bootdat1 = np.array(a1[ranu])
    bootdat2 = np.array(a2[ranu])
    bootdat1 = bootdat1.reshape([samplesize,nboots])
    bootdat2 = bootdat2.reshape([samplesize,nboots])

    if sigx is not None:
        bootdatsigx = np.array(sigx[ranu])
        bootdatsigx = bootdatsigx.reshape([samplesize,nboots])
    if sigy is not None:
        bootdatsigy = np.array(sigy[ranu])
        bootdatsigy = bootdatsigy.reshape([samplesize,nboots])

    acoef = np.zeros(nboots) ; bcoef=np.zeros(nboots)

    if sigx is not None:
        for iboot in range(0,nboots,1):
            acoef[iboot], bcoef[iboot] = linfit.tls(bootdat1[:,iboot],bootdat2[:,iboot],
                              bootdatsigx[:,iboot],bootdatsigy[:,iboot])
        
    else:
        for iboot in range(0,nboots,1):
            acoef[iboot], bcoef[iboot] = linfit.linfit_xy(bootdat1[:,iboot],bootdat2[:,iboot],
                                                         sigma=bootdatsigy[:,iboot])
    
    aminci = np.percentile(acoef,ptilemin)
    amaxci = np.percentile(acoef,ptilemax)
    bminci = np.percentile(bcoef,ptilemin)
    bmaxci = np.percentile(bcoef,ptilemax)

    arange=[aminci, amaxci]
    brange=[bminci, bmaxci]

    return arange, brange








