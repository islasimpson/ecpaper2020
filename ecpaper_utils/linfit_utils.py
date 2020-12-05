import xarray as xr
import numpy as np
import scipy.optimize as opt 

def linfit_lonlat(darray,dimname):
    """Calculate the linear trend coefficients for an [:, lat, lon] array"""
    darray_2d=darray.stack(allpoints=['lat','lon'])
    result=np.polyfit(darray_2d[dimname], darray_2d, 1)
    b=result[0,:] ; a=result[1,:]
    a_xr=xr.DataArray(a, coords=[darray_2d.allpoints], name='a')
    b_xr=xr.DataArray(b, coords=[darray_2d.allpoints], name='b')
    a_xr=a_xr.unstack('allpoints')
    b_xr=b_xr.unstack('allpoints')

    return a_xr, b_xr

def linfit_xy(x,y, sigma=None):
    """Calculate a weighted least squares linear fit between 1d arrays x and y
    Input: x (the predictor)
           y (the predictand)
           sigma (optional) the standard deviation of the y values
    Output: a and b and the residuals
    """

    if (sigma.any() != None):
        w = 1./(sigma)
        coefs = np.polyfit(x,y,1,w=w)
    else:
        coefs = np.polyfit(x,y,1)

    a = coefs[1]
    b = coefs[0]
    return a, b

def tls(xin ,yin, sigmaxin, sigmayin):
    """Calculate the total least squares regression with errors in both x and y.
       by minimizing the chi-sq ( y - a -bx)/(sigmay^2 + b^2*sigmax^2)
    """
    x = np.array(xin)
    y = np.array(yin)
    sigmax = np.array(sigmaxin)
    sigmay = np.array(sigmayin)

    # normalize
    xm = np.mean(x) ; ym =np.mean(y)
    xsd = np.std(x) ; ysd = np.std(y)
    x = (x[:] - xm)/xsd
    y = (y[:] - ym)/ysd
    sigmay = sigmay[:]/ysd
    sigmax = sigmax[:]/xsd

    # compute first guess for B slope using standard 1-D least squares fit
    # where the errors in x are ignored 
    a, b = linfit_xy(x,y,sigma=sigmay)
    firstfit = (b,a)

    # Find a and b that minimize chisq 
    def chisq_func(coefs):
        denom = (sigmay[:]**2 + (coefs[1]**2)*(sigmax[:]**2))
        num = ( y[:] - (coefs[0] + coefs[1]*x[:]))**2
        chisq = np.sum(num/denom)
        return chisq

    coefs = np.array([a, b])
    coefsin=coefs
    result = opt.minimize(chisq_func, coefs)
    coefs=result.x
    a = coefs[0]
    b = coefs[1]

    # un-normalize
    bout = b*(ysd/xsd)
    aout = ym + a*ysd - b*(ysd/xsd)*xm 
    return aout, bout
