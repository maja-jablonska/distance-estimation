library(data.table) # provides fread
library(bit64)      # enables fread to read source_id as long integer
#library(stableGR)   # provides stable.GR for MCMC Gelman-Rubin test    (not usually needed)
#library(boa)        # provides boa.geweke for MCMC convergence test    (not usually needed)
#library(diptest)    # provides dip.test for Hartigan's unimodalty test (only for flag)
source("./Rfiles/metropolis.R")
source("./Rfiles/functions.R")

# Function returns photogeometric samples

distance.photgeo.samples <- function(parallax,parallax_error,p,phot_g_mean_mag,bp_rp,pseudocolour,rEDSDModeGeo,Nsamp,Nburnin){
rcharQGmodelMax <- 1e4
verbose <- 1000    
bfac <- 1/2 
fpuMeas <- parallax_error/parallax


# G-magnitude correction
    
if (is.finite(pseudocolour) & is.finite(phot_g_mean_mag) & is.finite(bp_rp)){
    phot_g_mean_mag <- gmag.corrected(phot_g_mean_mag, bp_rp)   
    }
    
load(url(paste0("https://www.mpia.de/homes/calj/gedr3_distances/prior_fits_final/Robj/prior_",p, ".Robj")))
    
# provides
# - GGDfit, a dataframe of parameters specifying the geometric prior
# - CQDfit, list of bprpnbins objects (each is of class smooth.spline, 
#   two-gaussian, one-gaussian, or NULL)
# - bprpbincent, vector of bprpnbins bp_rp bin centres (for each CQGfit object)
# - bprpbinwidth, scalar float (common to all CQGfit objects)
# - bprprange, 2-element vector: range of applicability of CQDfit
# - priormaglim, scalar float: faintest mag used in fit    
    
flag <- matrix(nrow=1, ncol=4) # 3xinteger, 1xcharacter
colnames(flag) <- c("gmaglim", "diptest_geo", "diptest_photogeo", "QGmodel")
# Set Gmaglim flag, and other flags to default values
  flag[1,"gmaglim"] <- 0L
  if(is.na(phot_g_mean_mag)) {
    flag[1,"gmaglim"] <- 0L
  } else {
    flag[1,"gmaglim"] <- ifelse(phot_g_mean_mag <= priormaglim, 1, 2)
  }
  flag[1,"diptest_geo"]      <- 0L
  flag[1,"diptest_photogeo"] <- 0L
  flag[1,"QGmodel"]          <- "88"

if(is.na(phot_g_mean_mag) || is.na(bp_rp) ) {
    flag[1,"QGmodel"] <- "99"
  }

QGmodel1 <- NULL
QGmodel2 <- NULL
bprpofQGmodels <- NULL
if(bp_rp>=bprprange[1] && bp_rp<=bprprange[2]) {
  modSel <- which(bp_rp - bprpbincent < 0)[1] # position of model that makes upper bracket
  if(is.na(modSel)) { # colour is above upper end, so set only upper end model
    QGmodel2 <- CQDfit[[length(bprpbincent)]]
  } else {
    if(modSel==1) { # colour is below lower end, so set only lower end model
      QGmodel1 <- CQDfit[[modSel]]
    } else {
      QGmodel1 <- CQDfit[[modSel-1]] # lower bracketing model
      QGmodel2 <- CQDfit[[modSel]]   # upper bracketing model
      bprpofQGmodels <- c(bprpbincent[c(modSel-1, modSel)])
    }
  }
}
flag[1,"QGmodel"] = paste0(switch(class(QGmodel1), "NULL"=0, "one-gaussian"=1, "two-gaussian"=2, "smooth.spline"=3),
                           switch(class(QGmodel2), "NULL"=0, "one-gaussian"=1, "two-gaussian"=2, "smooth.spline"=3))
# Initialization
# Compute the characteristic distance for the QG model from its characteristic QG value.
# If we have two QG models (i.e. for interpolating), compute it from the weighted 
# combination of the two characteristic QG values.
# Then truncate rcharQGmodel to rcharQGmodelMax

if(!is.null(QGmodel1)) {
    if(!is.null(QGmodel2)) {
      modfrac <- (bp_rp - bprpofQGmodels[1]) / (bprpofQGmodels[2]-bprpofQGmodels[1])
      rcharQGmodel <- 10^( (phot_g_mean_mag - 
                              (modfrac*QGmodel1$QGcharval + (1-modfrac)*QGmodel2$QGcharval) + 5 ) / 5 )
    } else {
      rcharQGmodel <- 10^((phot_g_mean_mag-QGmodel1$QGcharval+5)/5)
    }
  } else {
    rcharQGmodel <- 10^((phot_g_mean_mag-QGmodel2$QGcharval+5)/5)
  }
  rcharQGmodel <- ifelse(rcharQGmodel<rcharQGmodelMax, rcharQGmodel, rcharQGmodelMax)
  
  if(fpuMeas <0) { # form simple average
    rInitPhotogeo <- (rEDSDModeGeo + rcharQGmodel)/2
  } else {
    if(fpuMeas <=0.01) { # fpuMeas is [0,0.01]
      rInitPhotogeo <- rEDSDModeGeo
    } else { # fpuMeas>0.01 => form variable-weighted average
      rInitPhotogeo <- ((bfac/fpuMeas)*rEDSDModeGeo + rcharQGmodel)/(bfac/fpuMeas + 1)
    }
  }
 # Only try to compute quantiles if MCMC initialization is finite

if(all(is.finite(func.post15(r=rInitPhotogeo, w=parallax, 
                                wsd=parallax_error, rlen=GGDfit$rlen, 
                                alpha=GGDfit$alpha, beta=GGDfit$beta,
                                gmag=phot_g_mean_mag, bprp=bp_rp,
                                QGmodel1=QGmodel1, QGmodel2=QGmodel2, bprpofQGmodels=bprpofQGmodels)))) {
    
    rStepPhotogeo <- 0.75*rInitPhotogeo*min(1/3, abs(fpuMeas))
    
    if (is.finite(rStepPhotogeo) && rStepPhotogeo>0) {
        
        sampPhotogeo  <- metrop(func=func.post15, thetaInit=rInitPhotogeo, Nburnin=Nburnin, 
                          Nsamp=Nsamp, sampleCov=rStepPhotogeo^2, verbose=verbose, 
                          w=parallax, wsd=parallax_error, 
                          rlen=GGDfit$rlen, alpha=GGDfit$alpha, beta=GGDfit$beta, 
                          gmag=phot_g_mean_mag, bprp=bp_rp,
                          QGmodel1=QGmodel1, QGmodel2=QGmodel2, bprpofQGmodels=bprpofQGmodels)[,3]
        
        #diptest <- try(dip.test(sampPhotogeo, simulate.p.value=FALSE)$p.value, silent=TRUE)
        #if(is.finite(diptest)) {
        #    flag[1,"diptest_photogeo"] <- ifelse(diptest<1e-3, 1, 0)
        #    }
        return(list(sampPhotogeo,flag))

        }
    }
}

#Returns photogeometric distance posterior; needed for plotting

photgeo.dist.posterior <- function(r,parallax,parallax_error,p,phot_g_mean_mag,bp_rp,pseudocolour){
    
    if (is.finite(pseudocolour) & is.finite(phot_g_mean_mag) & is.finite(bp_rp)){
    phot_g_mean_mag <- gmag.corrected(phot_g_mean_mag, bp_rp)   
    }
    
    load(url(paste0("https://www.mpia.de/homes/calj/gedr3_distances/prior_fits_final/Robj/prior_",p, ".Robj")))
    
    
    QGmodel1 <- NULL
    QGmodel2 <- NULL
    bprpofQGmodels <- NULL
    if(bp_rp>=bprprange[1] && bp_rp<=bprprange[2]) {
      modSel <- which(bp_rp - bprpbincent < 0)[1] # position of model that makes upper bracket
      if(is.na(modSel)) { # colour is above upper end, so set only upper end model
        QGmodel2 <- CQDfit[[length(bprpbincent)]]
      } else {
        if(modSel==1) { # colour is below lower end, so set only lower end model
          QGmodel1 <- CQDfit[[modSel]]
        } else {
          QGmodel1 <- CQDfit[[modSel-1]] # lower bracketing model
          QGmodel2 <- CQDfit[[modSel]]   # upper bracketing model
          bprpofQGmodels <- c(bprpbincent[c(modSel-1, modSel)])
        }
      }
    }
    
    
    
postPhotogeo <- ud.post15(r=r, w=parallax, wsd=parallax_error, 
                              rlen=GGDfit$rlen, alpha=GGDfit$alpha, beta=GGDfit$beta,
                              gmag=phot_g_mean_mag, bprp=bp_rp,
                              QGmodel1=QGmodel1, QGmodel2=QGmodel2, bprpofQGmodels=bprpofQGmodels)
    
}

# Photogeo Prior; needed for plotting

photgeo.dist.prior <- function(r,rlen,beta,alpha,p,phot_g_mean_mag,bp_rp,pseudocolour){
    
    if (is.finite(pseudocolour) & is.finite(phot_g_mean_mag) & is.finite(bp_rp)){
    phot_g_mean_mag <- gmag.corrected(phot_g_mean_mag, bp_rp)   
    }
    
    load(url(paste0("https://www.mpia.de/homes/calj/gedr3_distances/prior_fits_final/Robj/prior_",p, ".Robj")))
   
    
    QGmodel1 <- NULL
    QGmodel2 <- NULL
    bprpofQGmodels <- NULL
    
    if(bp_rp>=bprprange[1] && bp_rp<=bprprange[2]) {
      modSel <- which(bp_rp - bprpbincent < 0)[1] # position of model that makes upper bracket
      if(is.na(modSel)) { # colour is above upper end, so set only upper end model
        QGmodel2 <- CQDfit[[length(bprpbincent)]]
      } else {
        if(modSel==1) { # colour is below lower end, so set only lower end model
          QGmodel1 <- CQDfit[[modSel]]
        } else {
          QGmodel1 <- CQDfit[[modSel-1]] # lower bracketing model
          QGmodel2 <- CQDfit[[modSel]]   # upper bracketing model
          bprpofQGmodels <- c(bprpbincent[c(modSel-1, modSel)])
        }
      }
    }  
    
    prior <- d.prior5(r, rlen, alpha, beta) * Vectorize(interpolate.QGmodel.density, "r")(r=r, gmag=phot_g_mean_mag, bprp=bp_rp, QGmodel1=QGmodel1, QGmodel2=QGmodel2, bprpofQGmodels=bprpofQGmodels)

    return(prior)
}


