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
                             sel(time=slice(datestart, dateend))

        try:
            dat=dat.sel(plev=plev, method="nearest")
        except:
            # deal with different coordinate names
            try:
                dat=dat.rename({"pre":"plev"})
            except:
                dat=dat.rename({"lev":"plev"})
            dat=dat.sel(plev=plev, method="nearest")

        try:
            datzm=dat.mean(dim="lon")
        except:
            # deal with problematic coordinate names
            dat=dat.rename({"longitude":"lon", "latitude":"lat"})
            datzm=dat.mean(dim="lon")
    except:
        print("Something's wierd about thte time axis, decoding manually")
        dat=xr.open_mfdataset(filepath, coords="minimal", join="override",
                             decode_times=False)

        try:
            dat=dat.sel(plev=plev, method="nearest")
        except:
            try:
                dat=dat.rename({"pre":"plev"})
            except:
                dat=dat.rename({"lev":"plev"})
            dat=dat.sel(plev=plev, method="nearest")

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


def read_1lev(filepath, datestart, dateend, plev):
    """Read in a time slice for one pressure level from datestart to dateend.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
        plev (string) = pressure level to select
    """

    #First try opening and doing the select assuming everything is working ok with the time axis
    try:
        dat = \
        xr.open_mfdataset\
        (filepath, coords="minimal", join="override", decode_times=True, use_cftime=True).\
        sel(time=slice(datestart, dateend))

        try:
            dat = dat.sel(plev=plev, method="nearest")
        except:
            try:
                dat = dat.rename({"pre":"plev"})
            except:
                dat = dat.rename({"lev":"plev"})
            dat = dat.sel(plev=plev, method="nearest")

        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = False)

        try:
            dat = dat.sel(plev=plev, method="nearest")
        except:
            try:
                dat = dat.rename({"pre":"plev"})
            except:
                dat = dat.rename({"pre":"plev"})
            dat = dat.sel(plev=plev, method="nearest")

        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass

        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    return dat


def read_sfc(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    #First try opening and doing the select assuming everything is working ok with the time axis
    try:
        dat = \
        xr.open_mfdataset\
        (filepath, coords="minimal", join="override", decode_times=True, use_cftime=True).\
        sel(time=slice(datestart, dateend))
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = False)
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass

        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    return dat



