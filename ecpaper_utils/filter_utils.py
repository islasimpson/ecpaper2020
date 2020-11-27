import xarray as xr
import numpy as np
from scipy.fft import fft, ifft

def filterk(darray,kmin,kmax):
    """filter an [:,lat,lon] array to retain only zonal wavenumbers greater than or 
    equal to kmin and less than or equal to kmax
    """
    darray_np = np.array(darray)
    darray_np2 = darray_np.reshape(-1, darray_np.shape[-1])
    nk = kmax - kmin + 1
    nlon = darray_np2[0,:].size

    tempft = fft(darray_np2, axis=1)
    tempft2 = np.zeros_like(tempft)
    tempft2[:,kmin:kmax+1]=tempft[:,kmin:kmax+1]

    if (kmin == 0):
        tempft2[:,nlon-nk:nlon] = tempft[:,nlon-nk:nlon]
    else:
        tempft2[:,nlon-kmin-nk+1:nlon-kmin+1] = tempft[:,nlon-kmin-nk+1:nlon-kmin+1]

    darray_filtered = np.real(ifft(tempft2,axis=1))
    darray_filtered = darray_filtered.reshape(darray_np.shape)
    darray_filtered_xr = xr.DataArray(darray_filtered,coords=darray.coords)

    return darray_filtered_xr

