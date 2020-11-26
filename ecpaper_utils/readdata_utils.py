# routines for reading in data in various forms

import xarray as xr

def read_zonalmean_1lev(filepath, datestart, dateend, plev):
    """Read in a time slice for one pressure level from datestart to dateend
    calculate the zonal mean.  Try using datetime64 and if that doesn't work
    decode times manually.
    Args:
        filepath (string) = path to files e.g., "/path/to/files/*.nc"
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
        plev (string) = pressure level to select
    """

    # First try opening and doing the select assuming everything
    # is working ok with the time axis
    try:
        dat=xr.open_mfdataset(filepath, coords="minimal", join="override",
                             decode_times=True, use_cftime=True).\
                             sel(time=slice(datestart, dateend)).\
                             sel(plev=plev, method="nearest")
        try:
            datzm=dat.mean(dim="lon")
        except:
            # deal with problematic coordinate names
            dat=dat.rename({"longitude":"lon", "latitude":"lat"})
            datzm=dat.mean(dim="lon")
    except:
        print("Something's wierd about thte time axis, decoding manually")
        dat=xr.open_mfdataset(filepath, coords="minimal", join="override",
                             decode_times=False).sel(plev=plev, method="nearest")
        try:
            datzm=dat.mean(dim="lon")
        except:
            # deal with problematic coordinate names
            dat=dat.rename({"longitude":"lon", "latitude":"lat"})
            datzm=dat.mean(dim="lon")
   
        datzm=xr.decode_cf(datzm, use_cftime=True)
        datzm=datzm.sel(time=slice(datestart, dateend))
        datetimeindex=datzm.indexes['time'].to_datetimeindex()
        datzm['time'] = datetimeindex
   
    return datzm

