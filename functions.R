# CBJ Jan. & Oct. 2020
# Functions for fitting priors and inferring (photo)geometric distances

library(fields)     # provides image.plot
library("PolynomF") # required by mode.post3()


############ EDR3 parallax zeropoint correction

# Return parallax zeropoint in arcseconds for a number of sources.
# Units: phot_g_mean_mag [mag], nu_eff_used_in_astrometry and pseudocolour [micron^-1], ecl_lat [deg]
# phot_g_mean_mag must be catalogue G magnitude, i.e. before applying gmag.corrected().
# All quantities passed in should be vectors of same length (not checked).
# If any of the following are true:
# phot_g_mean_mag is NULL (can happen), ecl_lat is NULL (should not happen), 
# both nu_eff_used_in_astrometry and pseudocolour are NULL (can happen: if no parallax),
# then return default value.
# Otherwise:
# If nu_eff_used_in_astrometry is finite, apply 5p solution
# If pseudocolour is finite, apply 6p solution 
# (in that order, i.e. nu_eff_used_in_astrometry is used if both are defined).
# Formulae for the correction are taken from EDR3 release paper of Lindegren et al. (2020)
# (implementation is a translation of their MATLAB code) and is computed in microas.
edr3.par.zp <- function(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat) {
  
  Nsource <- length(ecl_lat)
  sinbeta <- sin(0.01745329 * ecl_lat)
  sel <- which(is.finite(phot_g_mean_mag) & is.finite(sinbeta))
  sel5p <- intersect(sel, which(is.finite(nu_eff_used_in_astrometry)))
  sel6p <- setdiff(intersect(sel, which(is.finite(pseudocolour))), sel5p)
  cat("Computing EDR3 parallax zeropoint:\n")
  cat("  Number with 5p solutions:", length(sel5p), "\n")
  cat("  Number with 6p solutions:", length(sel6p), "\n")
  cat("  Number with default zp (as G, nueff, or ecliptic latitude missing):", 
      Nsource-length(sel5p)-length(sel6p), "\n")
  rm(sel, sel5p, sel6p)
  
  tab <- vector("list", 2) # 5p, 6p
  
  # 5p solution
  tab[[1]]$tab <- matrix(nrow=13, byrow=TRUE, c(
    -26.98,-9.62,27.40,-25.1,-0.0,-1257,0,0,
    -27.23,-3.07,23.04,35.3,15.7,-1257,0,0,
    -30.33,-9.23,9.08,-88.4,-11.8,-1257,0,0,
    -33.54,-10.08,13.28,-126.7,11.6,-1257,0,0,
    -13.65,-0.07,9.35,-111.4,40.6,-1257,0,0,
    -19.53,-1.64,15.86,-66.8,20.6,-1257,0,0,
    -37.99,+2.63,+16.14,-5.7,+14.0,-1257,+107.9,+104.3,
    -38.33,+5.61,+15.42,0,+18.7,-1189,+243.8,+155.2,
    -31.05,+2.83,+8.59,0,+15.5,-1404,+105.5,+170.7,
    -29.18,-0.09,+2.41,0,+24.5,-1165,+189.7,+325.0,
    -18.40,+5.98,-6.46,0,+5.5,0,0,+276.6,
    -12.65,-4.57,-7.46,0,+97.9,0,0,0,
    -18.22,-15.24,-18.54,0,+128.2,0,0,0) )
  tab[[1]]$gmag  <- c(6.0,10.8,11.2,11.8,12.2,12.9,13.1,15.9,16.1,17.5,19.0,20.0,21.0) # length=Ngmag
  tab[[1]]$j     <- c(0,0,0,1,1,2,3,4) # length=Ncoef
  tab[[1]]$k     <- c(0,1,2,0,1,0,0,0) # length=Ncoef
  tab[[1]]$Ngmag <- nrow(tab[[1]]$tab)
  tab[[1]]$Ncoef <- ncol(tab[[1]]$tab)
  
  # 6p solution
  tab[[2]]$tab <- matrix(nrow=13, byrow=TRUE, c(
    -27.85,-7.78,27.47,-32.1,14.4,9.5,-67,
    -28.91,-3.57,22.92,7.7,12.6,1.6,-572,
    -26.72,-8.74,9.36,-30.3,5.6,17.2,-1104,
    -29.04,-9.69,13.63,-49.4,36.3,17.7,-1129,
    -12.39,-2.16,10.23,-92.6,19.8,27.6,-365,
    -18.99,-1.93,15.90,-57.2,-8.0,19.9,-554,
    -38.29,2.59,16.20,-10.5,1.4,0.4,-960,
    -36.83,4.20,15.76,22.3,11.1,10.0,-1367,
    -28.37,1.99,9.28,50.4,17.2,13.7,-1351,
    -24.68,-1.37,3.52,86.8,19.8,21.3,-1380,
    -15.32,4.01,-6.03,29.2,14.1,0.4,-563,
    -13.73,-10.92,-8.30,-74.4,196.4,-42.0,536,
    -29.53,-20.34,-18.74,-39.5,326.8,-262.3,1598) )
  tab[[2]]$gmag  <- c(6.0,10.8,11.2,11.8,12.2,12.9,13.1,15.9,16.1,17.5,19.0,20.0,21.0) # length=Ngmag
  tab[[2]]$j     <- c(0,0,0,1,1,1,2) # length=Ncoef
  tab[[2]]$k     <- c(0,1,2,0,1,2,0) # length=Ncoef
  tab[[2]]$Ngmag <- nrow(tab[[2]]$tab)
  tab[[2]]$Ncoef <- ncol(tab[[2]]$tab)
  
  left.index <- function(phot_g_mean_mag, gmag) { 
    # Return element of gmag that is closest less than or equal to phot_g_mean_mag.
    # Return 0 if below gmag range and length(gmag) is above.
    # gmag vector of non-decreasing values.
    sel <- which(gmag>phot_g_mean_mag)[1]-1
    return(ifelse(is.na(sel), length(gmag), sel))
  }
  func.b <- function(sinbeta) { # 3-element vector k=0,1,2
    c(1, sinbeta, sinbeta^2-1/3)
  }
  func.c <- function(nueff) { # 5-element vector j=0,1,2,3,4
    c(1, max(-0.24, min(0.24, nueff-1.48)), 
      min(0.24, max(0, 1.48-nueff))^3, min(0, nueff-1.24), max(0, nueff-1.72) )
  }
  
  parzp <- rep(-17e-6, Nsource) # [arcsec] default zp
  
  for(n in 1:Nsource) {
    
    if(!is.finite(phot_g_mean_mag[n]) || !is.finite(sinbeta[n])) {
      next
    }
    if(is.finite(nu_eff_used_in_astrometry[n])) {
      nueff <- nu_eff_used_in_astrometry[n]
      sol <- 1 # 5p
    } else {
      if(is.finite(pseudocolour[n])) {
        nueff <- pseudocolour[n]
        sol <- 2 # 6p
      } else {
        next
      }
    }
    
    termb <- func.b(sinbeta=sinbeta[n])
    termc <- func.c(nueff=nueff)
    left <- max(1, min(tab[[sol]]$Ngmag-1, left.index(phot_g_mean_mag[n], tab[[sol]]$gmag)))
    hfac <- max(0, min(1, (phot_g_mean_mag[n] - tab[[sol]]$gmag[left])/(tab[[sol]]$gmag[left+1] - tab[[sol]]$gmag[left])))
    zpt <- 0
    for(i in 1:tab[[sol]]$Ncoef) {
      coef <- (1-hfac) * tab[[sol]]$tab[left,i] + hfac*tab[[sol]]$tab[left+1,i]
      zpt  <- zpt + coef * termc[tab[[sol]]$j[i]+1] * termb[tab[[sol]]$k[i]+1]
    }
    parzp[n] <- 1e-6*zpt # convert to arcsec
    
  }
  
  return(parzp)
  
}



############ Cycle 3 G-magnitude correction

# Given vector of G magnitudes [mag] and corresponding BP-RP colours [mag]
# from EDR3, return vector of corrected G magnitudes. All vectors have same size.
# This should only be applied to 6p solutions, and only to
# finite values of gmag and bprp (none of these are checked here).
gmag.corrected <- function(gmag, bprp) {
  wc1 <- c(1.0087583646311324, -0.0254043532944348,   0.017466488186085014, -0.0027693464181843207)
  wc2 <- c(1.0052497473798703, -0.023230818732958947, 0.017399605901929606, -0.002533056831155478)
  Nsource <- length(gmag)
  clcol <- ifelse(bprp<0.25, 0.25, bprp) 
  clcol <- ifelse(clcol>3.0, 3.0, clcol) # these two lines do: min(3.0, max(0.25, bprp))
  polyclcol <- cbind(1, poly(clcol, degree=3, raw=TRUE))
  
  corfac <- rep(1, length=Nsource) # default, for G<13
  sel1 <- which(gmag>=13 & gmag<=16)
  corfac[sel1] <- drop(polyclcol[sel1,] %*% wc1)
  sel2 <- which(gmag>16)
  corfac[sel2] <- drop(polyclcol[sel2,] %*% wc2)
  
  return(gmag - 2.5*log10(corfac))
}


############ Likelihood

# Gaussian likelihood in w. Vectorized in all parameters
d.like <- function(w, r, wsd) dnorm(x=w, mean=1/r, sd=wsd)



############ Exponentially decreasing space density prior 
############ and corresponding geometric distance posterior

# Return normalized prior density. Vectorized in r and rlen.
d.prior3 <- function(r, rlen) ifelse(r>0, (1/(2*rlen^3))*r^2*exp(-r/rlen), 0) # normalized. Vectorized in r and rlen

# Return unnormalized posterior density.
ud.post3 <- function(r, w, wsd, rlen) {
  return( d.prior3(r, rlen)*d.like(w, r, wsd) )
}

# Return posterior mode - a single real number (or NA if something is wrong).
# If retall=TRUE, return all three roots of derivative of PDF (default is FALSE).
# Inputs:
# w    - parallax,             unrestricted
# wsd  - parallax uncertainty, must be >0
# rlen - prior length scale,   must be >0 and <Inf
# I think there are only two types of solutions of the cubic equation: 
# 1. only one root is real => one maxima. Take this as the mode.
# 2. all three roots are real => two maxima. 
#    if w>=0, take the smallest.
#    if w<0   take the positive one (there should be only one).
#   Is that correct/sufficient?
# The other two possibilities, 0 or 2 real roots, should not occur.
# If they do, return NA. But maybe 2 can, i.e. identical roots?
mode.post3 <- function(w, wsd, rlen, retall=FALSE) {
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- r^3/rlen - 2*r^2 + (w/wsd^2)*r - 1/wsd^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = NA,
                  "3" = ifelse(w>0, min(roots), roots[roots>0]) # should be real and unique
  )
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}


############ Generalized gamma distribution (GGD) prior
############ P(r) ~ r^{beta}exp[-(r/rlen)^{alpha}]
############ where rlen>0, alpha>0, beta>-1

# 2020-11-22: The correct definition should have r>=0 instead of r>0.
# This will have no impact in practice (MCMC will almost never return exactly zero),
# and I may anyway be using this as a (likely redundant) check to prevent zero
# distances, which can be a problem. 
# Note that in EDR3 all priors have a mode greater than zero.

# Return normalized density. Vectorized in r and rlen.
d.prior5 <- function(r, rlen, alpha, beta) {
  ifelse(r>0 & rlen>0 & alpha>0 & beta> -1, ( alpha/((rlen^(beta+1)) * gamma((beta+1)/alpha)) ) * 
           r^beta*exp(-(r/rlen)^alpha), 0)
}
# rlen <- 1e3 ; alpha <- 2 ; beta <- 2
# r=seq(from=1,to=(10*rlen),length.out=1e3)
# plot(r, d.prior5(r, rlen, alpha, beta), type="l")


############ Geometric distance posterior: 
############ GGD distance prior, Gaussian likelihood

# Return unnormalized posterior density. Vectorized in r.
ud.post5 <- function(r, w, wsd, rlen, alpha, beta) {
  return( ifelse(r<=0, 0, d.like(w, r, wsd) * 
                   d.prior5(r, rlen, alpha, beta)) )
                   # (1/r)) ) # For P(r) ~ 1/r
}

# Define func() required by metrop() in posterior sampling.
# This returns x = (log10(distance prior), log10(likelihood)).
# Not vectorized in anything.
# direct=TRUE uses direct coding of log posterior. This is faster, but
# omits the constants, so the posterior is not normalized.
# direct=FALSE calls normalized functions. This is slower (as it involves 
# taking logs of exponentials) but the posterior is normalized.
func.post5 <- function(r, w, wsd, rlen, alpha, beta, direct=TRUE) {
  if(direct) {
    if(!is.finite(r) || r<=0) {
      lnprior <- -Inf
      lnlike  <- -Inf
    } else {
      lnprior <- beta*log(r) - (r/rlen)^alpha
      # lnprior <- -log(r) # For P(r) ~ 1/r
      lnlike  <- -(w-1/r)^2/(2*wsd^2)
    }
    return(0.4342944819*c(lnprior, lnlike))
  } else {
    if(!is.finite(r) || r<=0) {
      like  <- 0
      prior <- 0
    } else {
      like  <- d.like(w, r, wsd)
      prior <- d.prior5(r, rlen, alpha, beta)
    }
    return( c(log10(prior), log10(like)) ) 
  }
}


############ Photogeometric distance posterior: 
############ GGD distance+CQD prior, Gaussian likelihood

# QGmodel is either NULL or is an element of CQDfit, 
# which is a model for P*(Q_G | bp_rp). This can be
# smooth.spline, or two- or one-component Gaussian,
# with additional fields mincount and QGrange .

# Evaluate density of QGmodel at a given r and gmag.
# If model is NULL, return NA.
# If smooth.spline is extrapolated beyond QGrange or 
# if any model returns a value below QGmodel$mincount, 
# return QGmodel$mincount.
density.QGmodel <- function(r, gmag, QGmodel) {
  if(r<=0) { # Should never happen
    return(NA) 
  }
  if(is.null(QGmodel)) {
    return(NA) # Can happen
  }
  QGdensity <- NA # default value
  if(class(QGmodel)=="smooth.spline") {
    q <- gmag-5*log10(r)+5
    if(q<QGmodel$QGrange[1] || q>QGmodel$QGrange[2]) {
      QGdensity <- QGmodel$mincount
    } else {
      QGdensity <- predict(QGmodel, x=q)$y
    }
  }
  if(class(QGmodel)=="two-gaussian") {
    q <- gmag-5*log10(r)+5
    QGdensity <- QGmodel$wt[1]*dnorm(x=q, mean=QGmodel$mean[1], sd=QGmodel$sd[1]) +
                 QGmodel$wt[2]*dnorm(x=q, mean=QGmodel$mean[2], sd=QGmodel$sd[2])
  }
  if(class(QGmodel)=="one-gaussian") {
    q <- gmag-5*log10(r)+5
    QGdensity <- dnorm(x=q, mean=QGmodel$mean, sd=QGmodel$sd)
  }
  QGdensity <- max(QGdensity, QGmodel$mincount) # NA if QGmodel is NULL
  return(QGdensity)
}

# Compute the value (density) using the two QGmodels and return.
# If neither model is NULL (so the density is NA), interpolate densities in bprp.
# If one model is NULL, return the density from the other.
# If both models are NULL, return a density of 1.
# This value is chosen so that it does not modify the prior, 
# and so has no impact on the inference (see func.post15).
interpolate.QGmodel.density <- function(r, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels) {
  den1 <- density.QGmodel(r, gmag, QGmodel1) 
  den2 <- density.QGmodel(r, gmag, QGmodel2)
  if(!is.na(den1) && !is.na(den2)) { 
    return( den1 + (bprp-bprpofQGmodels[1])*(den2-den1)/(bprpofQGmodels[2]-bprpofQGmodels[1]) )
  }
  if(!is.na(den1) && is.na(den2)) { 
    return(den1)
  } 
  if(is.na(den1) && !is.na(den2)) { 
    return(den2)
  } 
  return(1) # if both densities are NA
}
  
# Return unnormalized posterior density. Vectorized in r.
ud.post15 <- function(r, w, wsd, rlen, alpha, beta, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels) {
  return( ifelse(r<=0, 0, d.like(w, r, wsd) * 
                   d.prior5(r, rlen, alpha, beta) *
                   # (1/r) * # For P(r) ~ 1/r
                   Vectorize(interpolate.QGmodel.density, "r")(r, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels)) )
}

# Define func() required by metrop() in posterior sampling.
# This returns x = (log10(distance prior), log10(likelihood)).
# Not vectorized in anything.
# direct=TRUE uses direct coding of log posterior. This is faster, but
# omits the constants, so the posterior is not normalized.
# direct=FALSE calls normalized functions. This is slower (as it involves 
# taking logs of exponentials) but the posterior is normalized.
# The CQD prior is linearly interpolated by colour over the two QGmodels provided
# (which have colours given by bprpofQGmodels)
func.post15 <- function(r, w, wsd, rlen, alpha, beta, gmag, bprp,
                        QGmodel1=NULL, QGmodel2=NULL, bprpofQGmodels=NULL, direct=TRUE) {
  if(direct) {
    if(!is.finite(r) || r<=0) {
      lnprior <- -Inf
      lnlike  <- -Inf
    } else {
      lnprior <- beta*log(r) - (r/rlen)^alpha + 
        log(interpolate.QGmodel.density(r, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels))
      # For P(r) ~ 1/r
      #lnprior <- -log(r) + 
      #  log(interpolate.QGmodel.density(r, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels))
      lnlike  <- -(w-1/r)^2/(2*wsd^2)
    }
    return(0.4342944819*c(lnprior, lnlike))
  } else {
    if(!is.finite(r) || r<=0) {
      prior <- 0
      like  <- 0
    } else {
      prior <- d.prior5(r, rlen, alpha, beta) * 
        interpolate.QGmodel.density(r, gmag, bprp, QGmodel1, QGmodel2, bprpofQGmodels)
      like  <- d.like(w, r, wsd)
    }
    return( c(log10(prior), log10(like)) ) 
  }
}

# An interface to func.post15 that takes logr (base e) rather than r
func.post15.logr <- function(logr, ...) {
  func.post15(r=exp(logr), ...)
}




############ Plotting and analysis functions

# Compute a histogram with Nbreaks over range xlim.
# If xlim is not provided, the full range of the data defines xlim.
# If scale=TRUE, divide counts in each bin by the largest counts in a bin, i.e.
# the highest peak is 1.0
# If add=TRUE, add to an existing plot (in which case xlab is not used).
plot.hist <- function(dat, Nbreaks, xlim=NULL, xlab="x", 
                      col="black", lty=1, lwd=1.5, add=FALSE, scale=TRUE) {
  # If limits defined and 10+ samples lie within those limits, only plot within those limits,
  # Otherwise set range to data (previously I had this from 0 to maximum of data)
  if(all(!is.null(xlim))) { 
    sel <- which(dat>xlim[1] & dat<xlim[2])
    if(length(sel)>=10) {
      dat <- dat[sel]
    } 
  }
  else {
    xlim <- range(dat)
  }
  f <- hist(dat, breaks=Nbreaks, plot=FALSE)
  scalefac <- ifelse(scale, max(f$counts), 1)
  if(!add) {
    plot(f$breaks, c(f$counts,0)/scalefac, type="n", xlim=xlim, 
         ylim=c(0,1.05)*ifelse(scale, 1, max(f$counts)),
         xlab=xlab, ylab="", yaxt="n", yaxs="i", bty="n")
  }
  lines(f$breaks, c(f$counts,0)/scalefac, type="s", col=col, lty=lty, lwd=lwd)
  lines(x=f$breaks[1]*c(1,1), y=c(0,f$counts[1])/scalefac, col=col, lty=lty, lwd=lwd) # connect left of first bin to zero
}

# Compute two histograms (one from dat1 in col1, one from dat2 in col2) with Nbreaks over
# range xlim and plot both over range xlim with ylim computed to include peak of both histograms.
# If xlim is not provided, the full range of the data defines xlim.
plot.2hist <- function(dat1, dat2, Nbreaks, xlim=NULL, xlab="x",
                       col1="blue", col2="orange2", lty=1, lwd=1.5) {
  # If limits defined and 10+ samples lie within those limits, only plot within those limits,
  # Otherwise set range to data (previously I had this from 0 to maximum of data)
  if(all(!is.null(xlim))) { 
    sel1 <- which(dat1>xlim[1] & dat1<xlim[2])
    sel2 <- which(dat2>xlim[1] & dat2<xlim[2])
    if(length(sel1) + length(sel2) >= 20) {
      dat1 <- dat1[sel1]
      dat2 <- dat2[sel2]
    } 
  }
  else {
    xlim <- range(c(dat1,dat2))
  }
  f1 <- hist(dat1, breaks=Nbreaks, plot=FALSE)
  f2 <- hist(dat2, breaks=Nbreaks, plot=FALSE)
  ylim <- c(0, 1.05*max(c(f1$counts,f2$counts)))
  plot(c(f1$breaks, f2$breaks), c(f1$counts, 0, f2$counts, 0), type="n", xlim=xlim, ylim=ylim, 
       xlab=xlab, ylab="", yaxt="n", yaxs="i", bty="n")
  lines(f1$breaks, c(f1$counts,0), type="s", col=col1, lty=lty, lwd=lwd)
  lines(x=f1$breaks[1]*c(1,1), y=c(0,f1$counts[1]), col=col1, lty=lty, lwd=lwd) # connect left of first bin to zero
  lines(f2$breaks, c(f2$counts,0), type="s", col=col2, lty=lty, lwd=lwd)
  lines(x=f2$breaks[1]*c(1,1), y=c(0,f2$counts[1]), col=col2, lty=lty, lwd=lwd) # connect left of first bin to zero
}

# Compute a normalized histogram with Nbreaks over range xnormlim, and plot over range xlim.
# Either of these limits defaults to the range of data if they are NULL.
# Function also returns the maximum of the density.
plot.hist.normalized <- function(dat, Nbreaks, xnormlim=NULL, xlim=NULL, xlab="x", cex=1) {
  # If limits defined and 10+ samples lie within those limits, compute over those limits.
  # Otherwise set range to data.
  if(all(!is.null(xnormlim))) { 
    sel <- which(dat>xnormlim[1] & dat<xnormlim[2])
    if(length(sel)>=10) {
      dat <- dat[sel]
    } 
  }
  else {
    xnormlim <- range(dat)
  }
  f <- hist(dat, breaks=Nbreaks, plot=FALSE)
  # Plotting with type="s" puts vertical lines at points specified. As these are the
  # edges of the bins, I need to shift them left by half a bin width. I ste the edge
  # of the right-most bin to zero. If you instead used type="l" or "p", need deltax=0.
  deltax <- (f$mids[2]-f$mids[1])/2 
  plot(f$mids-deltax, f$density, type="s", lwd=1.5, xlim=xlim, yaxs="i",
       xlab=xlab, ylab="", yaxt="n", bty="n", cex.lab=cex, cex.axis=cex)
  lines(x=(f$mids[1]-deltax)*c(1,1), y=c(0,f$density[1]), lwd=2) # connect left of first bin to zero
  return(max(f$density))
}

# Given equal-length vectors of distances (estimate, lower CI, upper CI, true)
# compute various statistics and print to screen. Returns nothing.
print.stats <- function(rEst, rLo, rHi, rTrue) { 
  residual=rEst-rTrue 
  uncertainty=0.5*(rHi-rLo)
  cat("mean(residual) =", mean(residual, na.rm=TRUE), "\n")
  cat("  sd(residual) =", sd(residual, na.rm=TRUE),   "\n")
  cat("mean(residual/uncertainty) =", mean(residual/uncertainty, na.rm=TRUE), "\n")
  cat("  sd(residual/uncertainty) =", sd(residual/uncertainty, na.rm=TRUE),   "\n")
  sel <- which(residual>0)
  cat("Number of positive residuals:", length(sel), "\n")
  cat("residual>0: mean(residual/upperCI) =", 
      mean(residual[sel]/(rHi[sel]-rEst[sel]), na.rm=TRUE), "\n")
  cat("              sd(residual/upperCI) =", 
      sd(residual[sel]/(rHi[sel]-rEst[sel]), na.rm=TRUE), "\n")
  sel <- which(residual<0)
  cat("Number of negative residuals:", length(sel), "\n")
  cat("residual<0: mean(residual/lowerCI) =", 
      mean(residual[sel]/(rEst[sel]-rLo[sel]), na.rm=TRUE), "\n")
  cat("              sd(residual/lowerCI) =", 
      sd(residual[sel]/(rEst[sel]-rLo[sel]), na.rm=TRUE), "\n")
}

# Plot CQD from data.frame dat using colour vector mycols.
# dat must have named columns "bp_rp" and "QG", all of which are finite.
# If xlim!=NULL, then it specifies range of bp_rp for plot.
# If ylim!=NULL, then it specifies range of QG for plot.
# If zlim!=NULL, then it specifies range of density for plot.
# If logdensity=TRUE, plot density scale on logarithmic scale.
# If logdensity=TRUE and zlim=NULL, upper limit of density scale is
# maximum and lower limit is dynamicrange orders of magnitude smaller.
# (dynamicrange is only used for this case)
plot.CQD <- function(dat, mycols, bandwidth=c(0.02,0.02), 
                     gridsize=c(500,500), logdensity=FALSE,
                     dynamicrange=NULL, xlim=NULL, ylim=NULL, zlim=NULL,
                     legend.width=2) {
  if(any(is.null(xlim))) {
    xlim <- range(dat$bp_rp, na.rm=TRUE)
  }
  if(any(is.null(ylim))) {
    ylim <- range(dat$QG, na.rm=TRUE)
  }
  im2 <- bkde2D(x=dat[,c("bp_rp", "QG")], bandwidth=bandwidth,
                gridsize=gridsize, 
                range.x=list(xlim, ylim))
  im2 <- list(x=im2$x1, y=im2$x2, z=im2$fhat)
  if(logdensity) {
    im2$z <- log10(im2$z)
    if(any(is.null(zlim))) {
      zlim <- c(-1*dynamicrange,0) + max(im2$z, na.rm=TRUE)
    } else {
      zlim <- log10(zlim)
    }
    image.plot(im2, xlim=xlim, ylim=rev(ylim), 
               zlim=zlim, col=mycols, axis.args=list(cex.axis=1), 
               legend.lab=expression(paste(log[10],"(number density)")), 
               legend.width=legend.width, legend.line=3.5, legend.cex=1, xlab="BP-RP [mag]", 
               ylab=expression(paste(Q[G], " = ", M[G]+A[G], " [mag]")))
    
  } else {
    if(any(is.null(zlim))) {
      zlim <- c(0, max(im2$z)) 
    }
    image.plot(im2, xlim=xlim, ylim=rev(ylim), 
               zlim=zlim, col=mycols, axis.args=list(cex.axis=1), 
               legend.lab="number density",
               legend.width=legend.width, legend.line=3.5, legend.cex=1, xlab="BP-RP [mag]", 
               ylab=expression(paste(Q[G], " = ", M[G]+A[G], " [mag]")))
  }
}


############ Misc. functions

# Given the list compgrid, which specify vectors bprpcent and QGcent, 
# which make a regular 2D grid, and maxpercell, find those sources
# (elements of the vectors bprp, QG that give the data on those source)
# that fall within each cell of the grid (cell size equal to grid spacing).
# If there are more found than N=compgrid$maxpercell, then take the first N found.
# You might want to replace this with a random selection of N.
# Return single vector of all the elements that fall in a cell.
# This works even if bprp and QG are numeric(0): it returns NULL
find.sources.for.compgrid <- function(compgrid, bprp, QG) {
  
  if(length(QG)!=length(bprp)) {
    cat("find.sources.for.compgrid: QG and bprp have different lengths. Ignoring compgrid\n")
    return(NULL)
  }
  bprphcs <- (compgrid$bprpcent[2]-compgrid$bprpcent[1])/2 # half cell size in bprp
  QGhcs   <- (compgrid$QGcent[2]-compgrid$QGcent[1])/2 # half cell size in QG
  retain <- NULL
  for(i in 1:length(compgrid$bprpcent)) {
    for(j in 1:length(compgrid$QGcent)) {
      sel <- which(bprp>=compgrid$bprpcent[i]-bprphcs & bprp<compgrid$bprpcent[i]+bprphcs &
                   QG>=compgrid$QGcent[j]-QGhcs & QG<compgrid$QGcent[j]+QGhcs)
      if(length(sel)>0) {
        retain <- c(retain, sel[1:(min(length(sel), compgrid$maxpercell))])
      }
    }
  }
  return(sort(retain))
  
}

# Like find.sources.for.compgrid(), but now compute mean and sd of (inf-true)
# for all the points in a cell and return as 3D array [bprp, QG, statistic].
# bprp, QG, true, inf must all be the same length (only partially checked below).
# Many cells will be empty, and thus gridstats will be NaN (mean) or NA (sd).
# compgrid$maxpercell is not used.
# This is very slow if bprp, QG, true, inf are very large.
compute.residuals.for.compgrid <- function(compgrid, bprp, QG, true, inf) {
  
  if(length(true)!=length(inf)) {
    cat("compute.residuals.for.compgrid: true and inf have different lengths. Ignoring compgrid\n")
    return(NULL)
  }
  bprphcs <- (compgrid$bprpcent[2]-compgrid$bprpcent[1])/2 # half cell size in bprp
  QGhcs   <- (compgrid$QGcent[2]-compgrid$QGcent[1])/2 # half cell size in QG
  gridstats <- array(dim=c(length(compgrid$bprpcent), length(compgrid$QGcent), 2)) # mean, SD
  for(i in 1:length(compgrid$bprpcent)) {
    for(j in 1:length(compgrid$QGcent)) {
      sel <- which(bprp>=compgrid$bprpcent[i]-bprphcs & bprp<compgrid$bprpcent[i]+bprphcs &
                     QG>=compgrid$QGcent[j]-QGhcs & QG<compgrid$QGcent[j]+QGhcs)
      resid <- inf[sel]-true[sel]
      gridstats[i,j,1] <- mean(resid, na.rm=TRUE)
      gridstats[i,j,2] <-   sd(resid, na.rm=TRUE)
    }
  }
  return(gridstats)
  
}



