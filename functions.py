# In this File, all necessary functions are defined. This first includes the Metropolis algorithm used for sampling. Then the functions returning the EDSD, GGD and Photogeometric priors and the respective posteriors are defined. For the EDSD posterior, there is a function returning the mode, which will later be used as an initialisation fo the MCMC sampling of the GGD and EDSD posterior. Finally, the functions are defined that return the quantiles of the posteriors after being sampled.

import packages
from packages import *

# functions

# Metropolis (MCMC) algorithm

# Samples from function func. The function func() has to return a two-element array, the logPrior and logLike (log base 10), the sum of which is taken to be the log of the density function (i.e. unnormalized posterior) (we will later directly use the unnormalized posterior as the first element and simply add 0 by setting the second element 0)

# thetaInit: initial guess
# Nburnin : Number of Burn-In's
# Nsamp: Samplenumber
# sampleCov: Covariance of samples
# returns a Nsamp * (2+Ntheta) array, where the columns are
# 1:  log10 prior PDF
# 2:  log10 likelihood
# 3+: Ntheta parameters



def metrop(func,thetaInit,Nburnin,Nsamp,sampleCov,**kwargs):
    
    if not np.isscalar(thetaInit): 
        Ntheta = len(thetaInit)
    else:
        Ntheta = 1
        
    thetaCur = thetaInit
    funcCur = func(thetaInit,**kwargs)
    funcSamp = np.empty((Nsamp,2+Ntheta))
    funcSamp[:] = np.nan
    
    nAccept = 0
    acceptRate = 0
    
    for n in np.arange(1,Nburnin+Nsamp+1):
        
        if np.isscalar(sampleCov):
            thetaProp = np.random.normal(loc=thetaCur, scale=np.sqrt(sampleCov), size=1)
        else:
            thetaProp = np.random.multivariate_normal(mean=thetaCur, cov=sampleCov, size = 1)
            
        funcProp = func(thetaProp,**kwargs)
        logMR = np.sum(funcProp) - np.sum(funcCur) 
        
        if logMR >= 0 or logMR > np.log10(np.random.uniform(low=0, high=1, size=1)):
            thetaCur = thetaProp
            funcCur = funcProp
            nAccept = nAccept + 1
            acceptRate = nAccept/n
        
        if n > Nburnin:
            funcSamp[n-Nburnin-1,0:2] = funcCur
            funcSamp[n-Nburnin-1,2:(2+Ntheta)] = thetaCur
    
    return funcSamp


# function returns 1 if rlo <= r <= rhi

def lim(r,rlo,rhi):
    return 1*np.bitwise_and(r>=rlo,r<=rhi) 
    
# Normalized Gaussian likelihood in w

def d_like(w,r,wsd):
    return norm.pdf(w,loc=1/r,scale=wsd)

# EDSD prior---------------------------------------------------------------------------------------------------------------------------

# Unnormalized posterior in distance using exponentially decreasing space density prior (EDSD) (with length scale rlen)

def ud_distpost3(r,w,wsd,rlen):
    
    return d_like(w=w,r=r,wsd=wsd)*lim(r=r,rlo=0,rhi=np.inf)*np.exp(-r/rlen)*r**2

# logarithm of unnormalized posterior using the EDSD prior

def func_post3(r,w,wsd,rlen):
    return np.array([np.log10(ud_distpost3(r=r,w=w,wsd=wsd,rlen=rlen)),0],dtype=object)

# Find quantiles at specified probs for ud_distpost3 for given data
# using MCMC with specified initialization and step size.

# Returns a tuple of a dictionary and the MCMC samples. The entries in the dictionary are defined as follows:

# 1. code:     0=failed due to incorrect input or results which should not occur
#             +1=quantiles found
# 2. val:     if code=+1 : array containing the specified quantiles
#             if code = 0 : nan
# 3. message: if code=0 : string error message 

def quantile_distpost3(w,wsd,rlen,rInit,rStep,Nburnin,Nsamp,probs):
    if np.any(np.array([w,wsd,rlen,rInit,probs],dtype=object))==np.nan:
        return {'code':0,'val':np.nan,'message':'some inputs NA'},np.nan
    if w == np.inf:
        return {'code':0,'val':np.nan,'message':'parallax not finite'},np.nan
    if not (wsd > 0 and wsd != np.inf):
        return {'code':0,'val':np.nan,'message':'parallax uncertainty not (finite and positive)'},np.nan
    if not (rlen > 0 and rlen != np.inf):
        return {'code':0,'val':np.nan,'message':'rlen not (finite and positive)'},np.nan
    if not (rInit > 0 and rInit != np.inf):
        return {'code':0,'val':np.nan,'message':'rInit not (finite and positive)'},np.nan
    if not (rStep > 0 and rStep != np.inf):
        return {'code':0,'val':np.nan,'message':'rStep not (finite and positive)'},np.nan
    if not (Nsamp > 0 and Nsamp != np.inf):
        return {'code':0,'val':np.nan,'message':'Nsamp not (finite and positive)'},np.nan
    if not (Nburnin > 0 and Nburnin != np.inf):
        return {'code':0,'val':np.nan,'message':'Nburnin not (finite and positive)'},np.nan
    #if Nburnin >= Nsamp:
    #    return {'code':0,'val':np.nan,'message':'Nburnin >= Nsamp'},np.nan
    if np.any(probs<0) or np.any(probs>1):
        return {'code':0,'val':np.nan,'message':'probs not in range 0-1'},np.nan
   
    if ud_distpost3(r=rInit,w=w,wsd=wsd,rlen=rlen) <= 0:
        return {'code':0,'val':np.nan,'message':'metrop fails as posterior=zero at rInit'},np.nan
    
    samp = metrop(func=func_post3,thetaInit = rInit,Nburnin=Nburnin,Nsamp= Nsamp,sampleCov=rStep**2,w=w,wsd=wsd,rlen=rlen)
    samp = np.array([i[2] for i in samp])
    
    return {'code':1,'val': np.quantile(samp,probs),'message':np.nan}, samp

# Returns the mode of the EDSD posterior or nan if the inputs are other than defined below. 
# retall=True : return all roots of the derivative of the EDSD posterior sorted in increasing order of real part
# 
# Inputs:
# w    - parallax,             unrestricted
# wsd  - parallax uncertainty, must be >0
# rlen - prior length scale,   must be >0 and <Inf

#retall=False (default): return the mode of the EDSD posterior. This is selected as follows:

# There are two types of solutions of the cubic equation: 
# 1. only one root is real => one maxima. Take this as the mode.
# 2. all three roots are real => two maxima
#    if w>=0, take the smallest.
#    if w<0   take the positive one (there should be only one).
# The other two possibilities, 0 or 2 real roots, should not occur.
# If they do, return nan. 

# This is later used to get the initial guess rInit for quantile_distpost3

def mode_post3(w,wsd,rlen,retall = False):
    # special cases:
    # w<=0 works normally, except for w=-Inf
    if np.isnan(np.array([w,wsd,rlen])).any(): 
        return np.nan
    if w == -np.inf:
        return np.inf
    if w == np.inf: 
        return 0
    if wsd == 0: 
        return 1/w
    if wsd == np.inf:
        return 2*rlen
    if wsd<0:
        return np.nan
    if rlen <=0:
        return np.nan
    
    coeff = [1/rlen,-2,w/wsd**2,-1/wsd**2] 
    roots = np.roots(coeff)    
    
    r = sum(np.isreal(roots)) # gives number of real modes
    
    if r == 0:
        rMode = np.nan
    if r == 1:
        rMode = np.extract(np.isreal(roots),roots)
    if r == 2: 
        rMode = np.nan
    if r == 3: 
        rMode = np.extract(np.isreal(roots),roots)
        if w>0:
            rMode = min(roots)
        else:
            rMode = roots[roots>0]
    
    if retall is True:
        return(roots)
    else:
        return(rMode.real)
    
# GGD-Prior---------------------------------------------------------------------------------------------------------------------------

# generalized gamma distribution prior (GGD)

def prior4(r,rlen,alpha,beta):
    r = np.where(r > 0,r,0)
    return 1/gamma((beta+1)/alpha)*alpha/(rlen**(beta+1))*r**beta*np.exp(-((r/rlen)**alpha))

#Unnormalized Posterior using generalized gamma distribution (GGD) as prior

def ud_distpost4(r,w,wsd,rlen,alpha,beta):
    return d_like(w=w,r=r,wsd=wsd)*lim(r=r,rlo=0,rhi=np.inf)*prior4(r=r,rlen=rlen,alpha=alpha,beta=beta)

# logartithm of unnormalized posterior using the GGD prior

def func_post4(r,w,wsd,rlen,alpha,beta): 
    return np.array([np.log10(ud_distpost4(r=r,w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta)),0],dtype=object)

# Same as quantile_distpost4, only for ud_distpost4
        
def quantile_distpost4(w,wsd,rlen,alpha,beta,rInit,rStep,Nburnin,Nsamp,probs):
    
    if not (alpha > 0 and alpha != np.inf):
        return {'code':0,'val':np.nan,'message':'alpha not (finite and positive)'},np.nan
    if not (beta > 0 and beta != np.inf):
        return {'code':0,'val':np.nan,'message':'beta not (finite and positive)'},np.nan
    if np.any(np.array([w,wsd,rlen,alpha,beta,rInit,probs],dtype=object))==np.nan:
        return {'code':0,'val':np.nan,'message':'some inputs NA'},np.nan
    if w == np.inf:
        return {'code':0,'val':np.nan,'message':'parallax not finite'},np.nan
    if not (wsd > 0 and wsd != np.inf):
        return {'code':0,'val':np.nan,'message':'parallax uncertainty not (finite and positive)'},np.nan
    if not (rlen > 0 and rlen != np.inf):
        return {'code':0,'val':np.nan,'message':'rlen not (finite and positive)'},np.nan
    if not (rInit > 0 and rInit != np.inf):
        return {'code':0,'val':np.nan,'message':'rInit not (finite and positive)'},np.nan
    if not (rStep > 0 and rStep != np.inf):
        return {'code':0,'val':np.nan,'message':'rStep not (finite and positive)'},np.nan
    if not (Nsamp > 0 and Nsamp != np.inf):
        return {'code':0,'val':np.nan,'message':'Nsamp not (finite and positive)'},np.nan
    if not (Nburnin > 0 and Nburnin != np.inf):
        return {'code':0,'val':np.nan,'message':'Nburnin not (finite and positive)'},np.nan
    #if Nburnin >= Nsamp:
    #    return {'code':0,'val':np.nan,'message':'Nburnin >= Nsamp'},np.nan
    if np.any(probs<0) or np.any(probs>1):
        return {'code':0,'val':np.nan,'message':'probs not in range 0-1'},np.nan 
    
    if ud_distpost4(r=rInit,w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta) <= 0.:
        return {'code':0,'val':np.nan,'message':'metrop fails as posterior=zero at rInit'},np.nan
    
    samp = metrop(func=func_post4,thetaInit = rInit,Nburnin=Nburnin,Nsamp=Nsamp,sampleCov=rStep**2,w=w,wsd=wsd,rlen=rlen,alpha=alpha,beta=beta)
    samp = np.array([i[2] for i in samp])
    
    return {'code':1,'val': np.quantile(samp,probs),'message':np.nan},samp


# Photogeometric prior ---------------------------------------------------------------------------------------------------------------------------    
    
# Extract the photogeometric functions from the R-code. The file that is sourced is photogeo_dist_est.r, which contains the foloowing functions:
# r_photogeo_samples: returns the MCMC samples from sampleing the photogeometric distance posterior
# r_ud_distpost5_photogeo: photogeometric distance posterior
# r_photogeo_dist_prior: photogeometric distance prior


r_source = robjects.r['source']
r_source('./Rfiles/photogeo_dist_est.r')
r_photogeo_samples  = robjects.globalenv['distance.photgeo.samples']
r_ud_distpost5_photogeo  = robjects.globalenv['photgeo.dist.posterior']
r_photogeo_dist_prior = robjects.globalenv['photgeo.dist.prior']    
    
# returns the quantiles of the photogeometric distance posterior, like quantile_distpost3 and quantile_distpost4, additional input: magnitude phot_g_mean_mag, colour bp_rp, pseudocolour.
    
def quantile_distpost5(w,wsd,hp,phot_g_mean_mag,bp_rp,pseudocolour,rInit,Nsamp,Nburnin,probs):
    
    if w == np.inf:
        return {'code':0,'val':np.nan,'message':'parallax not finite'},np.nan
    if not (wsd > 0 and wsd != np.inf):
        return {'code':0,'val':np.nan,'message':'parallax uncertainty not (finite and positive)'},np.nan
    if not (rInit > 0 and rInit != np.inf):
        return {'code':0,'val':np.nan,'message':'rInit not (finite and positive)'},np.nan
    if not (Nsamp > 0 and Nsamp != np.inf):
        return {'code':0,'val':np.nan,'message':'Nsamp not (finite and positive)'},np.nan
    if not (Nburnin > 0 and Nburnin != np.inf):
        return {'code':0,'val':np.nan,'message':'Nburnin not (finite and positive)'},np.nan
    #if Nburnin >= Nsamp:
    #    return {'code':0,'val':np.nan,'message':'Nburnin >= Nsamp'},np.nan
    if np.any(probs<0) or np.any(probs>1):
        return {'code':0,'val':np.nan,'message':'probs not in range 0-1'},np.nan   
    if not np.isfinite(phot_g_mean_mag):
        return {'code':0,'val':np.nan,'message':'G not finite'},np.nan  
    if not np.isfinite(bp_rp): 
        return {'code':0,'val':np.nan,'message':'colours not finite'},np.nan
    #if hp <= 0 or hp >= 199:
     #   return {'code':0,'val':np.nan,'message':'HEALpixel out of range, only 0-199 available'},np.nan
    if w == np.nan or wsd == np.nan or hp == np.nan or rInit == np.nan or probs.any() == np.nan or phot_g_mean_mag == np.nan or bp_rp == np.nan:
        return {'code':0,'val':np.nan,'message':'some inputs NA'},np.nan
    
    samp0 = r_photogeo_samples(parallax=w,parallax_error=wsd,p=hp,phot_g_mean_mag=phot_g_mean_mag,bp_rp=bp_rp,pseudocolour = pseudocolour,rEDSDModeGeo=rInit,Nsamp=Nsamp,Nburnin=Nburnin)
    samp = samp0[0]
    flag = samp0[1]
    
    return {'code':1,'val': np.quantile(samp,probs),'message':np.nan,},samp,flag    

# function that resolves a simbad name to a gaia name by first determining the ra and deg in degree using simbad and then searching for the source of respective ra and deg in gaia (take sources in a circle of 0.01 deg closest to center). if the name does not exist or there is no corresponding source_id, an error-string is returned.

def resolve_simbad_to_gaia(simbad_name):
    result_table = Simbad.query_object(simbad_name)
    if result_table is not None:
        ra = Angle(result_table['RA'][0],unit='hour').degree
        dec = Angle(result_table['DEC'][0], unit='deg').degree 
        gaia_query = "SELECT source_id FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',"+str(ra)+","+str(dec)+",0.01))=1;"
        job = Gaia.launch_job(gaia_query)
        gaia_table = job.get_results()
        
        if len(gaia_table) > 0:
            
            return gaia_table['source_id'][0]
        
        else:
            
            return 'Error: source_id of simbad_name not found!'

    else:
        return 'Error: simbad_name not found!'