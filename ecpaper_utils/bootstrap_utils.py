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

