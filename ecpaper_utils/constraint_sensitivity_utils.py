import xarray as xr
import numpy as np
from ecpaper_utils import linfit_utils as linfit
from ecpaper_utils import bootstrap_utils as boot
import sys

def dotheconstraint_onlyxvar(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS'):

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit()

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS")
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)
        a = np.mean(aboots)
        b = np.mean(bboots)

    # sampling the uncertainty on the observed predictor.
    # 250 values for each observational value.
    nobs=obsx.size
    obstrue = np.zeros([nobs*250])
    obspdf = np.zeros([nobs*250])
    obssample = np.random.normal(0,sigx1mem,250)
    for iobs in range(0, obsx.size,1):
        obspdf[iobs*250:(iobs+1)*250] = obsx[iobs] + obssample[:]
        if (method == "BHM"):
            obstrue[iobs*250:(iobs+1)*250] = obsx[iobs]
  
    y = np.zeros([nobs*250])
    y[:] = a + b*obspdf[:] 

    vary = np.var(y) 

    return vary

def dotheconstraint_onlycoefs(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS'):

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit()

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS")
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)
        a = np.mean(aboots)
        b = np.mean(bboots)

    # sampling the uncertainty on the observed predictor.
    # 250 values for each observational value.
    nobs=obsx.size
    obsmean = np.mean(obsx)
 
    if (method != "BHM"):
        if (method == 'OLS'):
            sigxin=None
        if (method == 'TLS'):
            sigxin=sigxem
        aboots, bboots = boot.boot_regcoefs(xem, yem, sigx=sigxin, sigy=sigyem)
    else:
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3) 


 
    y = np.zeros([nobs*250])
    y[:] = aboots[:] + bboots[:]*obsmean

    vary = np.var(y) 

    return vary


def dotheconstraint_onlydelta(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS'):

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit()

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS")
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)
        a = np.mean(aboots)
        b = np.mean(bboots)

    # sampling the uncertainty on the observed predictor.
    # 250 values for each observational value.
    nobs=obsx.size
    obsmean = np.mean(obsx)
 
    if (method != "BHM"):
        if (method == 'OLS'):
            sigxin=None
        if (method == 'TLS'):
            sigxin=sigxem
        aboots, bboots = boot.boot_regcoefs(xem, yem, sigx=sigxin, sigy=sigyem)
    else:
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3) 

    if (method == "BHM"):
        sigforced = np.sqrt(del2[:])
    else:
        eps = y1mem[:] - (a+b*x1mem[:])
        sigeps = np.std(eps)
        sigwithiv = np.sqrt(sigeps**2 - (b**2)*sigx1mem)
        sigyiv = sigy1mem
        sigforced = np.sqrt(sigwithiv**2. -sigyiv**2.)


    randomvals = np.random.normal(0,1,nboots)
    if (method != "BHM"):
        noise_forced = randomvals*np.array(sigforced)
        yforced = np.zeros([nboots])
        yforced = a +b*obsmean + noise_forced[:]
    else:
        noise_forced = randomvals[:]*np.array(sigforced[:])
        yforced = np.zeros([nboots])
        yforced = a + b*obsmean + noise_forced[:] 

 
    vary = np.var(yforced) 

    return vary

def dotheconstraint_onlyiv(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS'):

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit()

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS")
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)
        a = np.mean(aboots)
        b = np.mean(bboots)

    # sampling the uncertainty on the observed predictor.
    # 250 values for each observational value.
    nobs=obsx.size
    obsmean = np.mean(obsx)
 
    if (method != "BHM"):
        if (method == 'OLS'):
            sigxin=None
        if (method == 'TLS'):
            sigxin=sigxem
        aboots, bboots = boot.boot_regcoefs(xem, yem, sigx=sigxin, sigy=sigyem)
    else:
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3) 

    if (method == "BHM"):
        sigforced = np.sqrt(del2[:])
        sigyiv = np.sqrt( (sigy1mem**2)*(1.-rxy1mem**2.))
    else:
        eps = y1mem[:] - (a+b*x1mem[:])
        sigeps = np.std(eps)
        sigwithiv = np.sqrt(sigeps**2 - (b**2)*sigx1mem)
        sigyiv = sigy1mem
        sigforced = np.sqrt(sigwithiv**2. -sigyiv**2.)


    randomvals = np.random.normal(0,1,nboots)
    if (method != "BHM"):
        noise_forced = randomvals*np.array(sigyiv)
        y = np.zeros([nboots])
        y = a +b*obsmean + noise_forced[:]
    else:
        noise_forced = randomvals[:]*np.array(sigyiv)
        y = np.zeros([nboots])
        y = a + b*obsmean + noise_forced[:] 

 
    vary = np.var(y) 

    return vary






def dotheconstraint_all(xem, yem, x1mem, y1mem, obsx, sigxem=None, sigyem=None,sigx1mem=None,
                    sigy1mem=None, rxyem=None, rxy1mem=None, seed=None, nboots=1000, method='OLS'):

    # check that all the needed things are there for the method
    if sigyem is None:
        print("You need to specify the sigma_y's for the ensemble mean for any constraint")
        sys.exit()
    if sigx1mem is None:
        print("You need to specify the sigma_x for 1 member for any constraint")
        sys.exit()
    if sigy1mem is None:
        print("You need to specify the sigma_y for 1 member for any constraint")
        sys.exit()

    if ((method == "TLS") or (method == "BHM")):
        if sigxem is None:
            print("You need to specify the sigma_x for the ensemble mean for TLS or BHM")
            sys.exit()

    if ((method != "OLS") and (method !="TLS") and (method != "BHM")):
        print("You have chosen an invalid method.  Choose OLS, TLS or BHM")
        sys.exit()

    if (method == "BHM"):
        if rxyem is None:
            print("You need to specify rxy for the ensemble mean to use the BHM")
            sys.exit()
        if rxy1mem is None:
            print("You need to specify rxy for 1 member to use the BHM")
            sys.exit()

    # calculate the regression coefficients using the ensemble mean
    if (method == 'OLS'):
        print("Constraining using OLS")
        a, b = linfit.linfit_xy(xem, yem, sigma=sigyem)
    if (method == 'TLS'):
        print("Constraining using TLS")
        a, b = linfit.tls(xem, yem, sigxem, sigyem)
    if (method == 'BHM'):
        print("Constraining using the BHM")
        aboots, bboots, del2, mux, delx2 = linfit.bhm(xem, yem, sigxem, sigyem, rxyem, iseed=3)

    if (method == "BHM"):
        # an array of standard deviations for the forced noise
        sigforced = np.sqrt(del2[:])
        # standard deviation for the internal variability noise component
        sigyiv=np.sqrt( (sigy1mem**2)*(1.-rxy1mem**2.))
    else:
        # calculate the single member residuals from the linear regression fit
        # their standard deviation and the standard deviation of the noise term
        # both with (sigwithiv) and without (sigforced) internal variability
        eps = y1mem[:] - (a + b*x1mem[:])
        sigeps = np.std(eps)
        sigwithiv = np.sqrt(sigeps**2 - (b**2)*sigx1mem)
        sigyiv = sigy1mem
        sigforced = np.sqrt(sigwithiv**2. - sigyiv**2.)

    # sampling the uncertainty on the observed predictor
    # 250 values for each observational value
    nobs=obsx.size
    obspdf = np.zeros([nobs*250])
    obstrue = np.zeros([nobs*250])
    obssample = np.random.normal(0,sigx1mem,250)
    for iobs in range(0,obsx.size,1):
        obspdf[iobs*250:(iobs+1)*250] = obsx[iobs] + obssample[:]
        if (method == "BHM"):
            obstrue[iobs*250:(iobs+1)*250] = obsx[iobs]


    # combine all the sampling
    # OLS and TLS
    if (method != "BHM"):

        # sample the noise terms and regression coefficients
        randomvals = np.random.normal(0,1,nboots)
        # forced + internal
        noise_withiv = randomvals*np.array(sigwithiv)
        # forced
        noise_forced = randomvals*np.array(sigforced)
        if (method == "OLS"):
            sigxin=None
        if (method == "TLS"):
            sigxin=sigxem
        aboots, bboots = boot.boot_regcoefs(xem, yem, sigx=sigxin, sigy=sigyem)

        # first, regression coefficient uncertainty with observational predictor uncertainty
        y = np.zeros([nobs*250*nboots])
        for iboot in range(0,nboots,1):
            y[iboot*nobs*250:(iboot+1)*nobs*250] = aboots[iboot] + bboots[iboot]*obspdf[:]

        # now adding on the noise terms
        yplusiv = np.zeros([nobs*250*nboots*nboots])
        yforced = np.zeros([nobs*250*nboots*nboots])
        for iboot in range(0,nboots,1):
            yplusiv[iboot*(nobs*250*nboots):(iboot+1)*(nobs*250*nboots)] = y[:] + noise_withiv[iboot]
            yforced[iboot*(nobs*250*nboots):(iboot+1)*(nobs*250*nboots)] = y[:] + noise_forced[iboot]

    else:
        #sample the noise terms
        # only do internal variability here because forced noise is 
        # dependend on delxdelx which is paired with alpha and beta's
        randomvals = np.random.normal(0,1,nboots*nboots)
        randomvals2 = np.random.normal(0,1,nboots*nboots)
        noiseiv = np.array(sigyiv)*randomvals2[:]

        yplusiv = np.zeros([nobs*250*nboots*nboots])
        yforced = np.zeros([nobs*250*nboots*nboots])
        for iboot in range(0,nboots,1):
            y=np.zeros([nobs*250])
            y[:] = aboots[iboot] + bboots[iboot]*obspdf[:]
            noise_forced = randomvals[:]*sigforced[iboot]

            yplusivt = np.zeros([nobs*250*nboots])
            yforcedt = np.zeros([nobs*250*nboots])
            for inoise in range(0,nboots,1):
                yforcedt[inoise*nboots:(inoise+1)*nboots]=\
                y[inoise]+noise_forced[inoise*nboots:(inoise+1)*nboots]

                yplusivt[inoise*nboots:(inoise+1)*nboots]=\
                yforcedt[inoise*nboots:(inoise+1)*nboots]+\
                noiseiv[inoise] + \
                np.array(rxy1mem)*(np.array(sigy1mem)/np.array(sigx1mem))*(obspdf[:]-obstrue[:])

            yforced[iboot*nobs*250*nboots:(iboot+1)*nobs*250*nboots]=yforcedt[:]
            yplusiv[iboot*nobs*250*nboots:(iboot+1)*nobs*250*nboots]=yplusivt[:]

    varforced = np.var(yforced)
    varplusiv = np.var(yplusiv)

    return varforced,varplusiv






