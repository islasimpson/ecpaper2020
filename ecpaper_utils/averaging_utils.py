import xarray as xr
import numpy as np

def cosweightlonlat(darray,lon1,lon2,lat1,lat2):
    """Calculate the weighted average for an [:,lat,lon] array over the region
    lon1 to lon2, and lat1 to lat2
    """
    # flip latitudes if they are decreasing
    if (darray.lat[0] > darray.lat[darray.lat.size -1]):
        print("flipping latitudes")
        darray = darray.sortby('lat')


    region=darray.sel(lon=slice(lon1,lon2),lat=slice(lat1,lat2))
    weights = np.cos(np.deg2rad(region.lat))
    regionw = region.weighted(weights)
    regionm = regionw.mean(("lon","lat"))

    return regionm
