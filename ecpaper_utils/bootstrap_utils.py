import xarray as xr
import numpy as np

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
