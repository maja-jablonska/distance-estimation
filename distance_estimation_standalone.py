#!/usr/bin/env python3
"""
Standalone command-line distance estimation tool for stellar parallax data.
Supports EDSD, GGD, and Photogeometric models without R dependencies.
"""

import click
import pandas as pd
import numpy as np
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import time
import traceback
from scipy.stats import norm
from scipy.special import gamma
import csv
import math
import multiprocessing as mp
import os
import json
from functools import partial, lru_cache

# For HEALPix calculations
try:
    from astropy_healpix import HEALPix
    import astropy.units as u
    HEALPIX_AVAILABLE = True
except ImportError:
    HEALPIX_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('distance_estimation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Log the HEALPix availability after logger is configured
if not HEALPIX_AVAILABLE:
    logger.warning("astropy_healpix not available. HEALPix calculations from coordinates will not work.")


class PhotogeometricModel:
    """Class to handle photogeometric model data and computations."""
    
    def __init__(self):
        self.models = {}  # Cache for loaded models
        
    def gmag_corrected(self, gmag: float, bp_rp: float) -> float:
        """
        Apply G-magnitude correction for 6p solutions.
        
        Args:
            gmag: G magnitude
            bp_rp: BP-RP color
            
        Returns:
            Corrected G magnitude
        """
        wc1 = np.array([1.0087583646311324, -0.0254043532944348, 
                       0.017466488186085014, -0.0027693464181843207])
        wc2 = np.array([1.0052497473798703, -0.023230818732958947, 
                       0.017399605901929606, -0.002533056831155478])
        
        clcol = np.clip(bp_rp, 0.25, 3.0)
        poly_terms = np.array([1, clcol, clcol**2, clcol**3])
        
        if gmag < 13:
            corfac = 1.0
        elif gmag <= 16:
            corfac = np.dot(wc1, poly_terms)
        else:
            corfac = np.dot(wc2, poly_terms)
        
        return gmag - 2.5 * np.log10(corfac)
    
    def load_prior_data(self, healpix: int) -> Dict:
        """
        Load prior data for a given HEALPix pixel.
        This is a simplified version - in practice you'd load from files or URLs.
        
        Args:
            healpix: HEALPix pixel number
            
        Returns:
            Dictionary with prior parameters
        """
        # For now, return mock data - in real implementation, load from files
        logger.warning(f"Using mock photogeometric prior data for HEALPix {healpix}")
        
        return {
            'GGDfit': {
                'rlen': 1250.0,
                'alpha': 1.2,
                'beta': 2.8
            },
            'CQDfit': [],  # Empty for mock
            'bprpbincent': np.array([]),
            'bprpbinwidth': 0.1,
            'bprprange': np.array([0.0, 3.0]),
            'priormaglim': 20.0
        }
    
    def density_qgmodel(self, r: float, gmag: float, qgmodel: Dict) -> float:
        """
        Evaluate density of QG model at given r and gmag.
        
        Args:
            r: Distance in parsecs
            gmag: G magnitude
            qgmodel: QG model dictionary
            
        Returns:
            Model density
        """
        if r <= 0 or qgmodel is None:
            return np.nan
            
        q = gmag - 5 * np.log10(r) + 5  # Calculate absolute magnitude
        model_type = qgmodel.get('type', 'null')
        
        if model_type == 'smooth_spline':
            qg_range = qgmodel.get('QGrange', [-np.inf, np.inf])
            if q < qg_range[0] or q > qg_range[1]:
                return qgmodel.get('mincount', 1e-10)
            # Simplified spline evaluation
            return qgmodel.get('mincount', 1e-10)
            
        elif model_type == 'two_gaussian':
            wt = qgmodel.get('wt', [0.5, 0.5])
            mean = qgmodel.get('mean', [0.0, 0.0])
            sd = qgmodel.get('sd', [1.0, 1.0])
            density = (wt[0] * norm.pdf(q, mean[0], sd[0]) + 
                      wt[1] * norm.pdf(q, mean[1], sd[1]))
            return max(density, qgmodel.get('mincount', 1e-10))
            
        elif model_type == 'one_gaussian':
            mean = qgmodel.get('mean', 0.0)
            sd = qgmodel.get('sd', 1.0)
            density = norm.pdf(q, mean, sd)
            return max(density, qgmodel.get('mincount', 1e-10))
            
        else:
            # For null models, use a simple exponential function of absolute magnitude
            # This provides a more realistic density than constant 1.0
            return np.exp(-0.5 * (q - 5.0)**2)  # Peak at M_G = 5.0
    
    def interpolate_qgmodel_density(self, r: float, gmag: float, bp_rp: float,
                                   qgmodel1: Optional[Dict], qgmodel2: Optional[Dict],
                                   bp_rp_of_qgmodels: Optional[np.ndarray]) -> float:
        """
        Interpolate QG model density between two models.
        
        Args:
            r: Distance in parsecs
            gmag: G magnitude  
            bp_rp: BP-RP color
            qgmodel1: First QG model
            qgmodel2: Second QG model
            bp_rp_of_qgmodels: BP-RP values for the two models
            
        Returns:
            Interpolated density
        """
        den1 = self.density_qgmodel(r, gmag, qgmodel1) if qgmodel1 else np.nan
        den2 = self.density_qgmodel(r, gmag, qgmodel2) if qgmodel2 else np.nan
        
        if not np.isnan(den1) and not np.isnan(den2) and bp_rp_of_qgmodels is not None:
            # Linear interpolation
            t = ((bp_rp - bp_rp_of_qgmodels[0]) / 
                 (bp_rp_of_qgmodels[1] - bp_rp_of_qgmodels[0]))
            return den1 + t * (den2 - den1)
        elif not np.isnan(den1):
            return den1
        elif not np.isnan(den2):
            return den2
        else:
            # If both models are NaN, use a simple exponential function
            # This provides a more realistic density than constant 1.0
            q = gmag - 5 * np.log10(r) + 5  # Calculate absolute magnitude
            return np.exp(-0.5 * (q - 5.0)**2)  # Peak at M_G = 5.0


# Optimized photogeometric model implementation for better parallelization

# Global worker state to avoid repeated initialization
_worker_photogeo_model = None
_worker_prior_cache = {}

class OptimizedPhotogeometricModel:
    """Optimized photogeometric model with minimal overhead."""
    
    def __init__(self):
        # Pre-compute correction coefficients
        self.wc1 = np.array([1.0087583646311324, -0.0254043532944348, 
                            0.017466488186085014, -0.0027693464181843207])
        self.wc2 = np.array([1.0052497473798703, -0.023230818732958947, 
                            0.017399605901929606, -0.002533056831155478])
    
    @lru_cache(maxsize=1000)  # Cache results for repeated calculations
    def gmag_corrected_cached(self, gmag: float, bp_rp: float) -> float:
        """Cached G-magnitude correction."""
        clcol = np.clip(bp_rp, 0.25, 3.0)
        poly_terms = np.array([1, clcol, clcol**2, clcol**3])
        
        if gmag < 13:
            corfac = 1.0
        elif gmag <= 16:
            corfac = np.dot(self.wc1, poly_terms)
        else:
            corfac = np.dot(self.wc2, poly_terms)
        
        return gmag - 2.5 * np.log10(corfac)
    
    @staticmethod
    def fast_qg_density(r: float, gmag: float, bp_rp: float) -> float:
        """Fast QG density calculation without object overhead."""
        if r <= 0:
            return 1e-10
        
        q = gmag - 5 * np.log10(r) + 5  # Absolute magnitude
        # Realistic but fast density model
        return np.exp(-0.5 * (q - 5.0)**2)


def init_worker_optimized(seed_base: int, log_level: str, prior_data: List):
    """Optimized worker initialization for photogeometric model."""
    global _worker_photogeo_model, _worker_prior_cache
    
    # Standard worker setup
    worker_seed = seed_base + mp.current_process().pid
    np.random.seed(worker_seed)
    np.seterr(divide='ignore')
    
    # Initialize model once per worker
    _worker_photogeo_model = OptimizedPhotogeometricModel()
    
    # Pre-cache all prior parameters to avoid repeated lookups
    _worker_prior_cache = {}
    for row in prior_data[1:]:  # Skip header
        try:
            hp = int(row[0])
            _worker_prior_cache[hp] = {
                'glon': float(row[1]),
                'glat': float(row[2]),
                'ggd_rlen': float(row[5]),
                'alpha': float(row[6]),
                'beta': float(row[7]),
                'edsd_rlen': float(row[10])
            }
        except (ValueError, IndexError):
            continue
    
    logging.basicConfig(level=getattr(logging, log_level))


def func_post5_optimized(r: float, params: Dict) -> np.ndarray:
    """Optimized photogeometric posterior for MCMC.
    
    Args:
        r: Distance parameter
        params: Pre-computed parameters dictionary
    
    Returns:
        Log posterior array [prior, likelihood]
    """
    if not np.isfinite(r) or r <= 0:
        return np.array([-np.inf, -np.inf])
    
    try:
        # Extract pre-computed parameters
        w = params['w']
        wsd = params['wsd']
        rlen = params['rlen']
        alpha = params['alpha']
        beta = params['beta']
        gmag = params['gmag']
        bp_rp = params['bp_rp']
        
        # Fast density calculation
        qg_density = OptimizedPhotogeometricModel.fast_qg_density(r, gmag, bp_rp)
        
        # Direct computation without intermediate objects
        ln_prior = (beta * np.log(r) - (r/rlen)**alpha + np.log(qg_density))
        ln_like = -(w - 1/r)**2 / (2 * wsd**2)
        
        # Convert to log10
        return 0.4342944819 * np.array([ln_prior, ln_like])
    
    except Exception:
        return np.array([-np.inf, -np.inf])


def metrop_optimized(func, params: Dict, r_init: float, nburnin: int, 
                    nsamp: int, r_step: float) -> np.ndarray:
    """Optimized Metropolis sampler with pre-computed parameters."""
    r_cur = r_init
    func_cur = func(r_cur, params)
    samples = np.full((nsamp, 3), np.nan)
    
    n_accept = 0
    
    for n in range(1, nburnin + nsamp + 1):
        # Propose new state
        r_prop = np.random.normal(r_cur, r_step)
        
        if r_prop > 0:  # Only evaluate if positive
            func_prop = func(r_prop, params)
            
            # Metropolis acceptance
            log_mr = np.sum(func_prop) - np.sum(func_cur)
            
            if log_mr >= 0 or log_mr > np.log10(np.random.uniform()):
                r_cur = r_prop
                func_cur = func_prop
                n_accept += 1
        
        # Store samples after burn-in
        if n > nburnin:
            idx = n - nburnin - 1
            samples[idx, :2] = func_cur
            samples[idx, 2] = r_cur
    
    return samples


def process_source_photogeo_optimized(row_data: Dict) -> Dict:
    """Optimized source processing using global worker state."""
    global _worker_photogeo_model, _worker_prior_cache
    
    source_id = row_data['gaia_dr3_source_id']
    
    try:
        # Extract data
        w = float(row_data['parallax_gaia']) * 1e-3  # mas to arcsec
        wsd = float(row_data['parallax_error_gaia']) * 1e-3
        phot_g_mean_mag = float(row_data['g_mag'])
        bp_rp = float(row_data['bp_rp'])
        
        # Validate
        if not (wsd > 0 and np.isfinite(w) and np.isfinite(phot_g_mean_mag) and np.isfinite(bp_rp)):
            return {'source_id': source_id, 'status': 'failed', 'message': 'Invalid input data'}
        
        # Calculate HEALPix (simplified version)
        hp = int(source_id) // (2**35 * 4**7)  # Fast integer division
        
        # Get cached prior parameters
        if hp not in _worker_prior_cache:
            return {'source_id': source_id, 'status': 'failed', 'message': f'HEALPix {hp} not in cache'}
        
        prior_params = _worker_prior_cache[hp]
        
        # Apply G-magnitude correction using cached model
        gmag = _worker_photogeo_model.gmag_corrected_cached(phot_g_mean_mag, bp_rp)
        
        # MCMC parameters
        rlen = prior_params['ggd_rlen']
        alpha = prior_params['alpha']
        beta = prior_params['beta']
        
        # Initialize
        r_init = 1.0 / w if w > 0 else 1000.0
        fpu_meas = wsd / w if w != 0 else 1.0
        r_step = 0.75 * r_init * min(1/3, abs(fpu_meas))
        
        # Pre-compute parameters dictionary
        mcmc_params = {
            'w': w,
            'wsd': wsd,
            'rlen': rlen,
            'alpha': alpha,
            'beta': beta,
            'gmag': gmag,
            'bp_rp': bp_rp
        }
        
        # Run optimized MCMC
        samples_full = metrop_optimized(
            func_post5_optimized, mcmc_params, r_init, 
            500, 5000, r_step  # nburnin, nsamp
        )
        
        samples = samples_full[:, 2]  # Extract distance samples
        quantiles = np.quantile(samples, [0.5, 0.159, 0.841])
        
        return {
            'source_id': source_id,
            'status': 'success',
            'message': 'Success',
            'samples': samples,
            'quantiles': quantiles,
            'hp': hp,
            'prior_params': prior_params
        }
        
    except Exception as e:
        return {
            'source_id': source_id,
            'status': 'failed',
            'message': f'Error: {str(e)}',
            'samples': None
        }


def process_batch_photogeo(sources_batch: List[Dict]) -> List[Dict]:
    """Process a batch of sources in a single worker."""
    return [process_source_photogeo_optimized(source) for source in sources_batch]


class OptimizedDistanceEstimator:
    """Distance estimator with optimized photogeometric parallelization."""
    
    def __init__(self, model: str, prior_data: List, n_processes: int = 4, batch_size: int = 50):
        self.model = model
        self.prior_data = prior_data
        self.n_processes = n_processes
        self.batch_size = batch_size
    
    def process_dataframe_photogeo_optimized(self, df, output_file: str):
        """Process dataframe with optimized photogeometric parallelization."""
        
        # Convert DataFrame to list of dictionaries
        sources_data = [row.to_dict() for _, row in df.iterrows()]
        
        # Create batches to reduce process creation overhead
        batches = [sources_data[i:i + self.batch_size] 
                  for i in range(0, len(sources_data), self.batch_size)]
        
        print(f"Processing {len(sources_data)} sources in {len(batches)} batches")
        print(f"Using {self.n_processes} processes with batch size {self.batch_size}")
        
        # Process with optimized multiprocessing
        with mp.Pool(
            processes=self.n_processes,
            initializer=init_worker_optimized,
            initargs=(42, 'INFO', self.prior_data)
        ) as pool:
            
            batch_results = pool.map(process_batch_photogeo, batches)
        
        # Flatten results
        all_results = [result for batch in batch_results for result in batch]
        
        # Process and save results (same as before)
        successful = sum(1 for r in all_results if r['status'] == 'success')
        print(f"Completed: {successful}/{len(all_results)} successful")
        
        return all_results


# Essential statistical functions
def lim(r: Union[float, np.ndarray], rlo: float, rhi: float) -> Union[float, np.ndarray]:
    """Return 1 if rlo <= r <= rhi, else 0."""
    return np.where((r >= rlo) & (r <= rhi), 1.0, 0.0)

def d_like(w: float, r: Union[float, np.ndarray], wsd: float) -> Union[float, np.ndarray]:
    """Normalized Gaussian likelihood in w."""
    return norm.pdf(w, loc=1/r, scale=wsd)

def metrop(func, theta_init: Union[float, np.ndarray], nburnin: int, nsamp: int, 
          sample_cov: Union[float, np.ndarray], verbose: int = 1000, **kwargs) -> np.ndarray:
    """
    Metropolis MCMC algorithm.
    
    Args:
        func: Function to sample from
        theta_init: Initial parameter values
        nburnin: Number of burn-in samples
        nsamp: Number of samples to keep
        sample_cov: Proposal covariance
        verbose: Print diagnostics every verbose samples
        **kwargs: Additional arguments for func
        
    Returns:
        Array of samples
    """
    if np.isscalar(theta_init):
        ntheta = 1
        theta_init = np.array([theta_init])
    else:
        ntheta = len(theta_init)
        
    theta_cur = theta_init.copy()
    func_cur = func(theta_cur, **kwargs)
    func_samp = np.full((nsamp, 2 + ntheta), np.nan)
    
    n_accept = 0
    
    for n in range(1, nburnin + nsamp + 1):
        # Propose new state
        if np.isscalar(sample_cov):
            theta_prop = np.random.normal(theta_cur, np.sqrt(sample_cov))
        else:
            theta_prop = np.random.multivariate_normal(theta_cur, sample_cov)
            
        func_prop = func(theta_prop, **kwargs)
        
        # Metropolis acceptance
        log_mr = np.sum(func_prop) - np.sum(func_cur)
        
        if log_mr >= 0 or log_mr > np.log10(np.random.uniform()):
            theta_cur = theta_prop.copy()
            func_cur = func_prop.copy()
            n_accept += 1
            
        # Store samples after burn-in
        if n > nburnin:
            idx = n - nburnin - 1
            func_samp[idx, :2] = func_cur
            func_samp[idx, 2:] = theta_cur
            
        # Diagnostics
        if verbose > 0 and (n % verbose == 0 or n == nburnin + nsamp):
            accept_rate = n_accept / n
            logger.debug(f"MCMC step {n}/{nburnin + nsamp}, accept rate: {accept_rate:.4f}")
    
    return func_samp

# EDSD functions
def d_prior3(r: Union[float, np.ndarray], rlen: float) -> Union[float, np.ndarray]:
    """Normalized EDSD prior density."""
    return np.where(r > 0, (1/(2*rlen**3)) * r**2 * np.exp(-r/rlen), 0.0)

def ud_distpost3(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float) -> Union[float, np.ndarray]:
    """Unnormalized EDSD posterior."""
    return d_like(w, r, wsd) * lim(r, 0, np.inf) * d_prior3(r, rlen)

def func_post3(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float) -> np.ndarray:
    """Log unnormalized EDSD posterior for MCMC."""
    if np.isscalar(r):
        result = ud_distpost3(r, w, wsd, rlen)
        if result <= 0:
            return np.array([-np.inf, 0.0])
        return np.array([np.log10(result), 0.0])
    else:
        # Handle array case
        results = ud_distpost3(r, w, wsd, rlen)
        log_results = np.where(results > 0, np.log10(results), -np.inf)
        return np.column_stack([log_results, np.zeros_like(log_results)])

def mode_post3(w: float, wsd: float, rlen: float, retall: bool = False) -> float:
    """Calculate mode of EDSD posterior."""
    # Handle special cases
    if np.isnan([w, wsd, rlen]).any():
        return np.nan
    if w == -np.inf:
        return np.inf
    if w == np.inf:
        return 0.0
    if wsd == 0:
        return 1.0 / w if w != 0 else np.nan
    if wsd == np.inf:
        return 2 * rlen
    if wsd < 0 or rlen <= 0:
        return np.nan
        
    # Solve cubic equation for mode
    coeff = [1/rlen, -2, w/wsd**2, -1/wsd**2]
    roots = np.roots(coeff)
    
    real_roots = roots[np.isreal(roots)].real
    n_real = len(real_roots)
    
    if n_real == 0:
        return np.nan
    elif n_real == 1:
        return real_roots[0]
    elif n_real == 2:
        return np.nan
    elif n_real == 3:
        if w > 0:
            return np.min(real_roots)
        else:
            positive_roots = real_roots[real_roots > 0]
            return positive_roots[0] if len(positive_roots) > 0 else np.nan
    
    return np.nan

def quantile_distpost3(w: float, wsd: float, rlen: float, r_init: float, r_step: float,
                      nburnin: int, nsamp: int, probs: np.ndarray) -> Tuple[Dict, np.ndarray]:
    """Calculate quantiles for EDSD posterior."""
    # Input validation
    if np.isnan([w, wsd, rlen, r_init]).any():
        return {'code': 0, 'val': np.nan, 'message': 'some inputs NA'}, np.array([np.nan])
    if w == np.inf:
        return {'code': 0, 'val': np.nan, 'message': 'parallax not finite'}, np.array([np.nan])
    if not (wsd > 0 and wsd != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'parallax uncertainty not (finite and positive)'}, np.array([np.nan])
    if not (rlen > 0 and rlen != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rlen not (finite and positive)'}, np.array([np.nan])
    if not (r_init > 0 and r_init != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rInit not (finite and positive)'}, np.array([np.nan])
    if not (r_step > 0 and r_step != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rStep not (finite and positive)'}, np.array([np.nan])
    if np.any(probs < 0) or np.any(probs > 1):
        return {'code': 0, 'val': np.nan, 'message': 'probs not in range 0-1'}, np.array([np.nan])
   
    if ud_distpost3(r_init, w, wsd, rlen) <= 0:
        return {'code': 0, 'val': np.nan, 'message': 'metrop fails as posterior=zero at rInit'}, np.array([np.nan])
    
    try:
        samp = metrop(func_post3, r_init, nburnin, nsamp, r_step**2, w=w, wsd=wsd, rlen=rlen)
        samples = samp[:, 2]  # Extract distance samples
        
        return {'code': 1, 'val': np.quantile(samples, probs), 'message': None}, samples
    except Exception as e:
        logger.error(f"MCMC sampling failed: {e}")
        return {'code': 0, 'val': np.nan, 'message': f'MCMC failed: {str(e)}'}, np.array([np.nan])

# GGD functions
def d_prior5(r: Union[float, np.ndarray], rlen: float, alpha: float, beta: float) -> Union[float, np.ndarray]:
    """Normalized GGD prior density."""
    condition = (r > 0) & (rlen > 0) & (alpha > 0) & (beta > -1)
    result = np.where(condition, 
                     (alpha / (rlen**(beta+1) * gamma((beta+1)/alpha))) * 
                     r**beta * np.exp(-(r/rlen)**alpha), 
                     0.0)
    return result

def ud_distpost4(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float, 
                alpha: float, beta: float) -> Union[float, np.ndarray]:
    """Unnormalized GGD posterior."""
    return d_like(w, r, wsd) * lim(r, 0, np.inf) * d_prior5(r, rlen, alpha, beta)

def func_post4(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float, 
              alpha: float, beta: float) -> np.ndarray:
    """Log unnormalized GGD posterior for MCMC."""
    if np.isscalar(r):
        result = ud_distpost4(r, w, wsd, rlen, alpha, beta)
        if result <= 0:
            return np.array([-np.inf, 0.0])
        return np.array([np.log10(result), 0.0])
    else:
        # Handle array case
        results = ud_distpost4(r, w, wsd, rlen, alpha, beta)
        log_results = np.where(results > 0, np.log10(results), -np.inf)
        return np.column_stack([log_results, np.zeros_like(log_results)])

def quantile_distpost4(w: float, wsd: float, rlen: float, alpha: float, beta: float,
                      r_init: float, r_step: float, nburnin: int, nsamp: int, 
                      probs: np.ndarray) -> Tuple[Dict, np.ndarray]:
    """Calculate quantiles for GGD posterior."""
    # Input validation
    if not (alpha > 0 and alpha != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'alpha not (finite and positive)'}, np.array([np.nan])
    if not (beta > 0 and beta != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'beta not (finite and positive)'}, np.array([np.nan])
    if np.isnan([w, wsd, rlen, alpha, beta, r_init]).any():
        return {'code': 0, 'val': np.nan, 'message': 'some inputs NA'}, np.array([np.nan])
    if w == np.inf:
        return {'code': 0, 'val': np.nan, 'message': 'parallax not finite'}, np.array([np.nan])
    if not (wsd > 0 and wsd != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'parallax uncertainty not (finite and positive)'}, np.array([np.nan])
    if not (rlen > 0 and rlen != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rlen not (finite and positive)'}, np.array([np.nan])
    if not (r_init > 0 and r_init != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rInit not (finite and positive)'}, np.array([np.nan])
    if not (r_step > 0 and r_step != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rStep not (finite and positive)'}, np.array([np.nan])
    if np.any(probs < 0) or np.any(probs > 1):
        return {'code': 0, 'val': np.nan, 'message': 'probs not in range 0-1'}, np.array([np.nan])
    
    if ud_distpost4(r_init, w, wsd, rlen, alpha, beta) <= 0:
        return {'code': 0, 'val': np.nan, 'message': 'metrop fails as posterior=zero at rInit'}, np.array([np.nan])
    
    try:
        samp = metrop(func_post4, r_init, nburnin, nsamp, r_step**2, 
                     w=w, wsd=wsd, rlen=rlen, alpha=alpha, beta=beta)
        samples = samp[:, 2]  # Extract distance samples
        
        return {'code': 1, 'val': np.quantile(samples, probs), 'message': None}, samples
    except Exception as e:
        logger.error(f"MCMC sampling failed: {e}")
        return {'code': 0, 'val': np.nan, 'message': f'MCMC failed: {str(e)}'}, np.array([np.nan])

# Photogeometric functions
def ud_distpost5(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float, 
                alpha: float, beta: float, gmag: float, bp_rp: float,
                qgmodel1: Optional[Dict], qgmodel2: Optional[Dict], 
                bp_rp_of_qgmodels: Optional[np.ndarray], 
                photogeo_model: PhotogeometricModel) -> Union[float, np.ndarray]:
    """Unnormalized photogeometric posterior."""
    prior_density = d_prior5(r, rlen, alpha, beta)
    likelihood = d_like(w, r, wsd)
    
    # Vectorized QG model density
    if np.isscalar(r):
        qg_density = photogeo_model.interpolate_qgmodel_density(
            r, gmag, bp_rp, qgmodel1, qgmodel2, bp_rp_of_qgmodels)
    else:
        qg_density = np.array([photogeo_model.interpolate_qgmodel_density(
            ri, gmag, bp_rp, qgmodel1, qgmodel2, bp_rp_of_qgmodels) for ri in r])
    
    return np.where(r > 0, likelihood * prior_density * qg_density, 0.0)

def func_post5(r: Union[float, np.ndarray], w: float, wsd: float, rlen: float, 
              alpha: float, beta: float, gmag: float, bp_rp: float,
              qgmodel1: Optional[Dict], qgmodel2: Optional[Dict], 
              bp_rp_of_qgmodels: Optional[np.ndarray],
              photogeo_model: PhotogeometricModel) -> np.ndarray:
    """Log unnormalized photogeometric posterior for MCMC."""
    if np.isscalar(r):
        if not np.isfinite(r) or r <= 0:
            return np.array([-np.inf, -np.inf])
        
        # Direct computation for speed
        try:
            ln_prior = (beta * np.log(r) - (r/rlen)**alpha + 
                       np.log(photogeo_model.interpolate_qgmodel_density(
                           r, gmag, bp_rp, qgmodel1, qgmodel2, bp_rp_of_qgmodels)))
            ln_like = -(w - 1/r)**2 / (2 * wsd**2)
            
            # Convert natural log to log10 and ensure 1D array
            result = 0.4342944819 * np.array([ln_prior, ln_like])
            return result.flatten()  # Ensure it's a 1D array
        except Exception:
            return np.array([-np.inf, -np.inf])
    else:
        # Handle array case - this shouldn't be called in normal MCMC operation
        results = []
        for ri in np.atleast_1d(r):
            if not np.isfinite(ri) or ri <= 0:
                results.append([-np.inf, -np.inf])
            else:
                try:
                    ln_prior = (beta * np.log(ri) - (ri/rlen)**alpha + 
                               np.log(photogeo_model.interpolate_qgmodel_density(
                                   ri, gmag, bp_rp, qgmodel1, qgmodel2, bp_rp_of_qgmodels)))
                    ln_like = -(w - 1/ri)**2 / (2 * wsd**2)
                    result = 0.4342944819 * np.array([ln_prior, ln_like])
                    results.append(result)
                except Exception:
                    results.append([-np.inf, -np.inf])
        return np.array(results)

def quantile_distpost5(w: float, wsd: float, hp: int, phot_g_mean_mag: float, 
                      bp_rp: float, pseudocolour: float, r_init: float,
                      nsamp: int, nburnin: int, probs: np.ndarray, 
                      prior_params: Dict[str, float]) -> Tuple[Dict, np.ndarray, Dict]:
    """Calculate quantiles for photogeometric posterior."""
    # Input validation
    if w == np.inf:
        return {'code': 0, 'val': np.nan, 'message': 'parallax not finite'}, np.array([np.nan]), {}
    if not (wsd > 0 and wsd != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'parallax uncertainty not (finite and positive)'}, np.array([np.nan]), {}
    if not (r_init > 0 and r_init != np.inf):
        return {'code': 0, 'val': np.nan, 'message': 'rInit not (finite and positive)'}, np.array([np.nan]), {}
    if not np.isfinite(phot_g_mean_mag):
        return {'code': 0, 'val': np.nan, 'message': 'G not finite'}, np.array([np.nan]), {}
    if not np.isfinite(bp_rp):
        return {'code': 0, 'val': np.nan, 'message': 'colours not finite'}, np.array([np.nan]), {}
    
    try:
        photogeo_model = PhotogeometricModel()
        
        # Apply G-magnitude correction if needed
        if np.isfinite(pseudocolour) and np.isfinite(phot_g_mean_mag) and np.isfinite(bp_rp):
            gmag = photogeo_model.gmag_corrected(phot_g_mean_mag, bp_rp)
        else:
            gmag = phot_g_mean_mag
            
        # Use prior parameters from DistanceEstimator instead of mock data
        rlen = prior_params['ggd_rlen']  # Use GGD rlen for photogeometric
        alpha = prior_params['alpha'] 
        beta = prior_params['beta']
        
        # Setup QG models (simplified - returns constant value of 1.0)
        qgmodel1 = None
        qgmodel2 = None
        bp_rp_of_qgmodels = None
        
        # Flag for compatibility
        flag = {
            'gmaglim': 1 if gmag <= 20.0 else 2,  # Default magnitude limit
            'diptest_geo': 0,
            'diptest_photogeo': 0,
            'QGmodel': '00'  # Mock value
        }
        
        # MCMC setup
        fpu_meas = wsd / w
        r_step = 0.75 * r_init * min(1/3, abs(fpu_meas))
        
        # Check if initialization is valid
        test_post = func_post5(r_init, w, wsd, rlen, alpha, beta, gmag, bp_rp,
                              qgmodel1, qgmodel2, bp_rp_of_qgmodels, photogeo_model)
        
        if not np.all(np.isfinite(test_post)):
            return {'code': 0, 'val': np.nan, 'message': 'invalid initialization'}, np.array([np.nan]), flag
        
        # Run MCMC
        samp = metrop(func_post5, r_init, nburnin, nsamp, r_step**2,
                     w=w, wsd=wsd, rlen=rlen, alpha=alpha, beta=beta,
                     gmag=gmag, bp_rp=bp_rp, qgmodel1=qgmodel1, qgmodel2=qgmodel2,
                     bp_rp_of_qgmodels=bp_rp_of_qgmodels, photogeo_model=photogeo_model)
        
        samples = samp[:, 2]  # Extract distance samples
        
        return {'code': 1, 'val': np.quantile(samples, probs), 'message': None}, samples, flag
        
    except Exception as e:
        logger.error(f"Photogeometric sampling failed: {e}")
        return {'code': 0, 'val': np.nan, 'message': f'Photogeometric failed: {str(e)}'}, np.array([np.nan]), {}


class DistanceEstimator:
    """Main class for distance estimation processing."""
    
    def __init__(self, model: str, seed: int = 42, prior_file: str = 'prior_summary.csv',
                 n_processes: int = 1, chunk_size: int = 1000, log_level: str = 'INFO'):
        """
        Initialize the distance estimator.
        
        Args:
            model: Model type ('EDSD', 'GGD', or 'Photogeometric')
            seed: Random seed for reproducibility
            prior_file: Path to prior summary CSV file (None for worker processes)
            n_processes: Number of processes for multiprocessing (1 = no multiprocessing)
            chunk_size: Number of sources to process before saving intermediate results
            log_level: Logging level string
        """
        self.model = model
        self.seed = seed
        self.prior_file = prior_file
        self.n_processes = n_processes
        self.chunk_size = chunk_size
        self.log_level = log_level
        self.use_optimized = False  # Default to standard processing
        self.setup_random_seeds()
        
        # Load prior data (skip for worker processes)
        if prior_file is not None:
            self.prior_data = self.load_prior_data()
        else:
            self.prior_data = None  # Will be set by worker process
        
        # MCMC parameters
        self.nsamp = int(5e3)
        self.nburnin = int(5e2)
        
        # Processing statistics
        self.total_sources = 0
        self.successful_sources = 0
        self.failed_sources = 0
        self.start_time = None
        
        if prior_file is not None:  # Only log for main process
            logger.info(f"Initialized DistanceEstimator with model: {model}, seed: {seed}")
            logger.info(f"Multiprocessing: {'Enabled' if n_processes > 1 else 'Disabled'} "
                       f"({n_processes} processes)")
            logger.info(f"Chunk size for periodic saving: {chunk_size}")
    
    def setup_random_seeds(self):
        """Set random seeds for reproducibility."""
        np.random.seed(self.seed)
        np.seterr(divide='ignore')  # Ignore divide by 0 errors
        logger.info(f"Set random seed to {self.seed}")
    
    def validate_input_data(self, df: pd.DataFrame) -> bool:
        """
        Validate that the input dataframe has all required columns.
        
        Args:
            df: Input dataframe
            
        Returns:
            bool: True if valid, False otherwise
        """
        # Base columns required for all models
        required_base_cols = ['gaia_dr3_source_id', 'parallax_gaia', 'parallax_error_gaia']
        
        # Check for optional coordinate columns for HEALPix calculation
        has_coordinates = 'ra' in df.columns and 'dec' in df.columns
        
        if self.model == 'EDSD':
            required_cols = required_base_cols
        elif self.model == 'GGD':
            required_cols = required_base_cols
        elif self.model == 'Photogeometric':
            # Photogeometric needs photometric data
            required_cols = required_base_cols + ['g_mag', 'bp_rp', 'pseudocolour']
        else:
            logger.error(f"Unknown model: {self.model}")
            return False
        
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return False
        
        # Check that we can calculate HEALPix either from source_id or coordinates
        if not has_coordinates and 'gaia_dr3_source_id' not in df.columns:
            logger.error("Need either 'gaia_dr3_source_id' or both 'ra' and 'dec' columns to calculate HEALPix")
            return False
        
        logger.info(f"Input validation passed for {len(df)} sources")
        logger.info(f"Available columns: {list(df.columns)}")
        if has_coordinates:
            logger.info("Will calculate HEALPix from RA/Dec coordinates")
        else:
            logger.info("Will calculate HEALPix from source_id")
        return True
    
    def calculate_mcmc_parameters(self, w: float, wsd: float, rlen: float) -> Tuple[float, float]:
        """
        Calculate MCMC initialization and step size parameters.
        
        Args:
            w: Parallax in arcsec
            wsd: Parallax uncertainty in arcsec
            rlen: Length scale parameter
            
        Returns:
            Tuple of (rInit, rStep)
        """
        try:
            # Get initial guess from EDSD mode
            r_init = mode_post3(w, wsd, rlen, retall=False)
            if np.isnan(r_init) or r_init <= 0:
                r_init = 1.0 / w if w > 0 else 1000.0  # Fallback
                logger.warning(f"Using fallback rInit: {r_init}")
            
            r_step = 0.75 * r_init * min(1/3, abs(wsd/w)) if w != 0 else 0.1 * r_init
            
            return r_init, r_step
            
        except Exception as e:
            logger.warning(f"Error calculating MCMC parameters: {e}. Using defaults.")
            r_init = 1000.0  # Default 1 kpc
            r_step = 100.0   # Default step
            return r_init, r_step
    
    def process_single_source(self, row: pd.Series) -> Dict:
        """
        Process a single source to get posterior samples.
        
        Args:
            row: Pandas series with source data
            
        Returns:
            Dictionary with processing results
        """
        source_id = row['gaia_dr3_source_id']
        logger.debug(f"Processing source {source_id}")
        
        try:
            # Extract and convert data
            w = float(row['parallax_gaia']) * 1e-3  # Convert mas to arcsec
            wsd = float(row['parallax_error_gaia']) * 1e-3  # Convert mas to arcsec
            
            # Validate parallax data
            if not (wsd > 0 and wsd != np.inf):
                return {
                    'source_id': source_id,
                    'status': 'failed',
                    'message': 'Invalid parallax uncertainty',
                    'samples': None
                }
            
            if w == np.inf:
                return {
                    'source_id': source_id,
                    'status': 'failed', 
                    'message': 'Infinite parallax',
                    'samples': None
                }
            
            # Calculate HEALPix pixel
            if 'ra' in row and 'dec' in row:
                # Calculate from coordinates if available
                ra = float(row['ra'])
                dec = float(row['dec'])
                hp = self.calculate_healpix_from_coordinates(ra, dec)
            else:
                # Calculate from source_id
                hp = self.calculate_healpix_from_source_id(int(source_id))
            
            # Get prior parameters for this HEALPix pixel
            prior_params = self.get_prior_parameters(hp)
            
            # Get model-specific parameters
            if self.model == 'EDSD':
                rlen = prior_params['edsd_rlen']
            elif self.model == 'GGD':
                rlen = prior_params['ggd_rlen']
                alpha = prior_params['alpha']
                beta = prior_params['beta']
            elif self.model == 'Photogeometric':
                # For photogeometric, use EDSD rlen as initial parameters
                rlen = prior_params['edsd_rlen']  
                alpha = prior_params['alpha']
                beta = prior_params['beta']
                # Extract photometric data
                phot_g_mean_mag = float(row['g_mag'])
                bp_rp = float(row['bp_rp'])
                pseudocolour = float(row['pseudocolour'])
            
            # Calculate MCMC parameters using EDSD rlen for initialization
            r_init, r_step = self.calculate_mcmc_parameters(w, wsd, prior_params['edsd_rlen'])
            
            logger.debug(f"Source {source_id}: HP={hp}, w={w*1e3:.3f} mas, wsd={wsd*1e3:.3f} mas, rInit={r_init:.1f} pc")
            
            # Run MCMC sampling based on model
            if self.model == 'EDSD':
                result = quantile_distpost3(
                    w=w, wsd=wsd, rlen=rlen, r_init=r_init, r_step=r_step,
                    nburnin=self.nburnin, nsamp=self.nsamp,
                    probs=np.array([0.5, 0.159, 0.841])
                )
                
            elif self.model == 'GGD':
                result = quantile_distpost4(
                    w=w, wsd=wsd, rlen=rlen, alpha=alpha, beta=beta,
                    r_init=r_init, r_step=r_step, nburnin=self.nburnin, nsamp=self.nsamp,
                    probs=np.array([0.5, 0.159, 0.841])
                )
                
            elif self.model == 'Photogeometric':
                result = quantile_distpost5(
                    w=w, wsd=wsd, hp=hp, phot_g_mean_mag=phot_g_mean_mag,
                    bp_rp=bp_rp, pseudocolour=pseudocolour, r_init=r_init,
                    nsamp=self.nsamp, nburnin=self.nburnin,
                    probs=np.array([0.5, 0.159, 0.841]),
                    prior_params=prior_params
                )
            
            # Check if sampling was successful
            quant_dict = result[0]
            samples = result[1]
            
            if quant_dict['code'] == 1:
                logger.debug(f"Source {source_id}: Successfully obtained {len(samples)} samples")
                return {
                    'source_id': source_id,
                    'status': 'success',
                    'message': 'Sampling successful',
                    'samples': samples,
                    'quantiles': quant_dict['val'],
                    'hp': hp,
                    'prior_params': prior_params
                }
            else:
                logger.warning(f"Source {source_id}: Sampling failed - {quant_dict['message']}")
                return {
                    'source_id': source_id,
                    'status': 'failed',
                    'message': quant_dict['message'],
                    'samples': None,
                    'hp': hp,
                    'prior_params': prior_params
                }
                
        except Exception as e:
            logger.error(f"Source {source_id}: Unexpected error - {str(e)}")
            logger.debug(f"Source {source_id}: Traceback - {traceback.format_exc()}")
            return {
                'source_id': source_id,
                'status': 'failed',
                'message': f'Unexpected error: {str(e)}',
                'samples': None
            }
    
    def process_dataframe(self, df: pd.DataFrame, output_file: str) -> None:
        """
        Process entire dataframe and save all samples to output file.
        
        Args:
            df: Input dataframe with source data
            output_file: Path to output file for samples
        """
        self.total_sources = len(df)
        self.start_time = time.time()
        
        logger.info(f"Starting processing of {self.total_sources} sources")
        logger.info(f"Model: {self.model}")
        logger.info(f"MCMC parameters: {self.nsamp} samples, {self.nburnin} burn-in")
        
        # Prepare output data
        samples_by_source = {}  # Dictionary to group samples by source_id
        summary_data = []
        
        # Process each source
        for idx, row in df.iterrows():
            if idx % 100 == 0:
                elapsed = time.time() - self.start_time
                rate = idx / elapsed if elapsed > 0 else 0
                eta = (self.total_sources - idx) / rate if rate > 0 else float('inf')
                logger.info(f"Processed {idx}/{self.total_sources} sources "
                           f"({100*idx/self.total_sources:.1f}%) - "
                           f"Rate: {rate:.1f} sources/sec - "
                           f"ETA: {eta/60:.1f} min")
            
            result = self.process_single_source(row)
            
            if result['status'] == 'success':
                self.successful_sources += 1
                
                # Store samples grouped by source_id
                samples = result['samples']
                if samples is not None and len(samples) > 0:
                    samples_by_source[result['source_id']] = samples.tolist()
                else:
                    samples_by_source[result['source_id']] = []
                
                # Store summary data
                quantiles = result['quantiles']
                prior_params = result.get('prior_params', {})
                hp = result.get('hp', np.nan)
                
                summary_data.append({
                    'source_id': result['source_id'],
                    'status': 'success',
                    'message': result['message'],
                    'hp': hp,
                    'glon': prior_params.get('glon', np.nan),
                    'glat': prior_params.get('glat', np.nan),
                    'r_median': float(quantiles[0]),
                    'r_lo': float(quantiles[1]), 
                    'r_hi': float(quantiles[2]),
                    'rlen': prior_params.get('edsd_rlen' if self.model == 'EDSD' else 'ggd_rlen', np.nan),
                    'alpha': prior_params.get('alpha', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                    'beta': prior_params.get('beta', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                    'n_samples': len(samples) if samples is not None else 0
                })
                
            else:
                self.failed_sources += 1
                summary_data.append({
                    'source_id': result['source_id'],
                    'status': 'failed',
                    'message': result['message'],
                    'hp': result.get('hp', np.nan),
                    'glon': np.nan,
                    'glat': np.nan,
                    'r_median': np.nan,
                    'r_lo': np.nan,
                    'r_hi': np.nan,
                    'rlen': np.nan,
                    'alpha': np.nan,
                    'beta': np.nan,
                    'n_samples': 0
                })
        
        # Save results
        self.save_results(samples_by_source, summary_data, output_file)
        self.print_summary()
    
    def process_dataframe_chunk(self, df: pd.DataFrame, start_idx: int = 0, 
                               end_idx: Optional[int] = None) -> pd.DataFrame:
        """
        Process a chunk of the dataframe.
        
        Args:
            df: Input dataframe with source data
            start_idx: Starting index (inclusive)
            end_idx: Ending index (exclusive). If None, process to end of dataframe
            
        Returns:
            DataFrame with chunk subset
        """
        if end_idx is None:
            end_idx = len(df)
        
        end_idx = min(end_idx, len(df))
        
        if start_idx >= len(df):
            raise ValueError(f"Start index {start_idx} is beyond dataframe length {len(df)}")
        
        logger.info(f"Processing chunk: indices {start_idx} to {end_idx-1} "
                   f"({end_idx - start_idx} sources)")
        
        return df.iloc[start_idx:end_idx].copy()
    
    def process_dataframe_parallel(self, df: pd.DataFrame, output_file: str,
                                  start_idx: int = 0, end_idx: Optional[int] = None) -> None:
        """
        Process dataframe with multiprocessing and periodic saving.
        
        Args:
            df: Input dataframe with source data
            output_file: Path to output file for samples
            start_idx: Starting index for processing
            end_idx: Ending index for processing (None = process to end)
        """
        # Get the chunk to process
        if end_idx is not None or start_idx > 0:
            df_chunk = self.process_dataframe_chunk(df, start_idx, end_idx)
        else:
            df_chunk = df
            
        self.total_sources = len(df_chunk)
        self.start_time = time.time()
        
        logger.info(f"Starting parallel processing of {self.total_sources} sources")
        logger.info(f"Model: {self.model}")
        logger.info(f"MCMC parameters: {self.nsamp} samples, {self.nburnin} burn-in")
        logger.info(f"Processes: {self.n_processes}")
        logger.info(f"Chunk size for saving: {self.chunk_size}")
        
        # Use optimized processing for photogeometric model
        if self.model == 'Photogeometric' and self.n_processes > 1 and self.use_optimized:
            self._process_photogeo_optimized(df_chunk, output_file)
        else:
            # Prepare for parallel processing
            if self.n_processes > 1:
                self._process_parallel(df_chunk, output_file)
            else:
                self._process_sequential(df_chunk, output_file)
        
        self.print_summary()
    
    def _process_photogeo_optimized(self, df: pd.DataFrame, output_file: str) -> None:
        """Process dataframe using optimized photogeometric parallelization."""
        logger.info("Using optimized photogeometric parallelization")
        
        # Convert DataFrame to list of dictionaries
        sources_data = [row.to_dict() for _, row in df.iterrows()]
        
        # Create batches to reduce process creation overhead
        batch_size = max(1, self.chunk_size // 10)  # Smaller batches for better load balancing
        batches = [sources_data[i:i + batch_size] 
                  for i in range(0, len(sources_data), batch_size)]
        
        logger.info(f"Processing {len(sources_data)} sources in {len(batches)} batches")
        logger.info(f"Using {self.n_processes} processes with batch size {batch_size}")
        
        # Setup output files
        output_path = Path(output_file)
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize storage for results
        all_samples = {}
        all_summary = []
        processed_count = 0
        
        # Process with optimized multiprocessing
        with mp.Pool(
            processes=self.n_processes,
            initializer=init_worker_optimized,
            initargs=(self.seed, self.log_level, self.prior_data)
        ) as pool:
            
            # Process batches
            for batch_idx, batch in enumerate(batches):
                logger.info(f"Processing batch {batch_idx + 1}/{len(batches)}: "
                           f"{len(batch)} sources")
                
                # Process this batch
                batch_results = pool.map(process_source_photogeo_optimized, batch)
                
                # Process results
                batch_samples = {}
                batch_summary = []
                
                for result in batch_results:
                    processed_count += 1
                    
                    if result['status'] == 'success':
                        self.successful_sources += 1
                        
                        # Store samples
                        samples = result['samples']
                        if samples is not None and len(samples) > 0:
                            batch_samples[result['source_id']] = samples.tolist()
                            all_samples[result['source_id']] = samples.tolist()
                        else:
                            batch_samples[result['source_id']] = []
                            all_samples[result['source_id']] = []
                        
                        # Store summary data
                        quantiles = result['quantiles']
                        prior_params = result.get('prior_params', {})
                        hp = result.get('hp', np.nan)
                        
                        summary_row = {
                            'source_id': result['source_id'],
                            'status': 'success',
                            'message': result['message'],
                            'hp': hp,
                            'glon': prior_params.get('glon', np.nan),
                            'glat': prior_params.get('glat', np.nan),
                            'r_median': float(quantiles[0]),
                            'r_lo': float(quantiles[1]), 
                            'r_hi': float(quantiles[2]),
                            'rlen': prior_params.get('ggd_rlen', np.nan),
                            'alpha': prior_params.get('alpha', np.nan),
                            'beta': prior_params.get('beta', np.nan),
                            'n_samples': len(samples) if samples is not None else 0
                        }
                        
                    else:
                        self.failed_sources += 1
                        summary_row = {
                            'source_id': result['source_id'],
                            'status': 'failed',
                            'message': result['message'],
                            'hp': result.get('hp', np.nan),
                            'glon': np.nan,
                            'glat': np.nan,
                            'r_median': np.nan,
                            'r_lo': np.nan,
                            'r_hi': np.nan,
                            'rlen': np.nan,
                            'alpha': np.nan,
                            'beta': np.nan,
                            'n_samples': 0
                        }
                    
                    batch_summary.append(summary_row)
                    all_summary.append(summary_row)
                
                # Save intermediate results for this batch
                batch_num = batch_idx + 1
                self._save_chunk_results(batch_samples, batch_summary, output_file, batch_num)
                
                # Progress update
                elapsed = time.time() - self.start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (self.total_sources - processed_count) / rate if rate > 0 else float('inf')
                
                logger.info(f"Completed batch {batch_num}: {processed_count}/{self.total_sources} sources "
                           f"({100*processed_count/self.total_sources:.1f}%) - "
                           f"Rate: {rate:.1f} sources/sec - ETA: {eta/60:.1f} min")
        
        # Save final consolidated results
        self.save_results(all_samples, all_summary, output_file)
    
    def _process_parallel(self, df: pd.DataFrame, output_file: str) -> None:
        """Process dataframe using multiprocessing with periodic saving."""
        # Prepare arguments for worker processes
        args_list = []
        for idx, row in df.iterrows():
            row_data = row.to_dict()
            args_list.append((row_data, self.model, self.nsamp, self.nburnin, self.prior_data))
        
        # Setup output files
        output_path = Path(output_file)
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize storage for results
        all_samples = {}
        all_summary = []
        processed_count = 0
        
        # Process in chunks with multiprocessing
        with mp.Pool(processes=self.n_processes, 
                     initializer=init_worker, 
                     initargs=(self.seed, self.log_level)) as pool:
            
            # Process chunks of chunk_size
            for chunk_start in range(0, len(args_list), self.chunk_size):
                chunk_end = min(chunk_start + self.chunk_size, len(args_list))
                chunk_args = args_list[chunk_start:chunk_end]
                
                logger.info(f"Processing chunk {chunk_start//self.chunk_size + 1}: "
                           f"sources {chunk_start} to {chunk_end-1}")
                
                # Process this chunk
                chunk_results = pool.map(process_source_worker, chunk_args)
                
                # Process results
                chunk_samples = {}
                chunk_summary = []
                
                for result in chunk_results:
                    processed_count += 1
                    
                    if result['status'] == 'success':
                        self.successful_sources += 1
                        
                        # Store samples
                        samples = result['samples']
                        if samples is not None and len(samples) > 0:
                            chunk_samples[result['source_id']] = samples.tolist()
                            all_samples[result['source_id']] = samples.tolist()
                        else:
                            chunk_samples[result['source_id']] = []
                            all_samples[result['source_id']] = []
                        
                        # Store summary data
                        quantiles = result['quantiles']
                        prior_params = result.get('prior_params', {})
                        hp = result.get('hp', np.nan)
                        
                        summary_row = {
                            'source_id': result['source_id'],
                            'status': 'success',
                            'message': result['message'],
                            'hp': hp,
                            'glon': prior_params.get('glon', np.nan),
                            'glat': prior_params.get('glat', np.nan),
                            'r_median': float(quantiles[0]),
                            'r_lo': float(quantiles[1]), 
                            'r_hi': float(quantiles[2]),
                            'rlen': prior_params.get('edsd_rlen' if self.model == 'EDSD' else 'ggd_rlen', np.nan),
                            'alpha': prior_params.get('alpha', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                            'beta': prior_params.get('beta', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                            'n_samples': len(samples) if samples is not None else 0
                        }
                        
                    else:
                        self.failed_sources += 1
                        summary_row = {
                            'source_id': result['source_id'],
                            'status': 'failed',
                            'message': result['message'],
                            'hp': result.get('hp', np.nan),
                            'glon': np.nan,
                            'glat': np.nan,
                            'r_median': np.nan,
                            'r_lo': np.nan,
                            'r_hi': np.nan,
                            'rlen': np.nan,
                            'alpha': np.nan,
                            'beta': np.nan,
                            'n_samples': 0
                        }
                    
                    chunk_summary.append(summary_row)
                    all_summary.append(summary_row)
                
                # Save intermediate results for this chunk
                chunk_num = chunk_start // self.chunk_size + 1
                self._save_chunk_results(chunk_samples, chunk_summary, output_file, chunk_num)
                
                # Progress update
                elapsed = time.time() - self.start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (self.total_sources - processed_count) / rate if rate > 0 else float('inf')
                
                logger.info(f"Completed chunk {chunk_num}: {processed_count}/{self.total_sources} sources "
                           f"({100*processed_count/self.total_sources:.1f}%) - "
                           f"Rate: {rate:.1f} sources/sec - ETA: {eta/60:.1f} min")
        
        # Save final consolidated results
        self.save_results(all_samples, all_summary, output_file)
        
    def _process_sequential(self, df: pd.DataFrame, output_file: str) -> None:
        """Process dataframe sequentially with periodic saving."""
        # Prepare output data
        all_samples = {}
        all_summary = []
        
        # Track chunk progress
        chunk_samples = {}
        chunk_summary = []
        processed_count = 0
        chunk_num = 1
        
        # Process each source
        for idx, row in df.iterrows():
            result = self.process_single_source(row)
            processed_count += 1
            
            if result['status'] == 'success':
                self.successful_sources += 1
                
                # Store samples
                samples = result['samples']
                if samples is not None and len(samples) > 0:
                    chunk_samples[result['source_id']] = samples.tolist()
                    all_samples[result['source_id']] = samples.tolist()
                else:
                    chunk_samples[result['source_id']] = []
                    all_samples[result['source_id']] = []
                
                # Store summary data
                quantiles = result['quantiles']
                prior_params = result.get('prior_params', {})
                hp = result.get('hp', np.nan)
                
                summary_row = {
                    'source_id': result['source_id'],
                    'status': 'success',
                    'message': result['message'],
                    'hp': hp,
                    'glon': prior_params.get('glon', np.nan),
                    'glat': prior_params.get('glat', np.nan),
                    'r_median': float(quantiles[0]),
                    'r_lo': float(quantiles[1]), 
                    'r_hi': float(quantiles[2]),
                    'rlen': prior_params.get('edsd_rlen' if self.model == 'EDSD' else 'ggd_rlen', np.nan),
                    'alpha': prior_params.get('alpha', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                    'beta': prior_params.get('beta', np.nan) if self.model in ['GGD', 'Photogeometric'] else np.nan,
                    'n_samples': len(samples) if samples is not None else 0
                }
                
            else:
                self.failed_sources += 1
                summary_row = {
                    'source_id': result['source_id'],
                    'status': 'failed',
                    'message': result['message'],
                    'hp': result.get('hp', np.nan),
                    'glon': np.nan,
                    'glat': np.nan,
                    'r_median': np.nan,
                    'r_lo': np.nan,
                    'r_hi': np.nan,
                    'rlen': np.nan,
                    'alpha': np.nan,
                    'beta': np.nan,
                    'n_samples': 0
                }
            
            chunk_summary.append(summary_row)
            all_summary.append(summary_row)
            
            # Save chunk if we've reached chunk_size
            if len(chunk_summary) >= self.chunk_size:
                self._save_chunk_results(chunk_samples, chunk_summary, output_file, chunk_num)
                
                # Progress update
                elapsed = time.time() - self.start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (self.total_sources - processed_count) / rate if rate > 0 else float('inf')
                
                logger.info(f"Completed chunk {chunk_num}: {processed_count}/{self.total_sources} sources "
                           f"({100*processed_count/self.total_sources:.1f}%) - "
                           f"Rate: {rate:.1f} sources/sec - ETA: {eta/60:.1f} min")
                
                # Reset for next chunk
                chunk_samples = {}
                chunk_summary = []
                chunk_num += 1
            
            # Regular progress updates
            if processed_count % 100 == 0:
                elapsed = time.time() - self.start_time
                rate = processed_count / elapsed if elapsed > 0 else 0
                eta = (self.total_sources - processed_count) / rate if rate > 0 else float('inf')
                logger.info(f"Processed {processed_count}/{self.total_sources} sources "
                           f"({100*processed_count/self.total_sources:.1f}%) - "
                           f"Rate: {rate:.1f} sources/sec - ETA: {eta/60:.1f} min")
        
        # Save any remaining results in the last partial chunk
        if chunk_summary:
            self._save_chunk_results(chunk_samples, chunk_summary, output_file, chunk_num)
            logger.info(f"Completed final chunk {chunk_num}: {processed_count}/{self.total_sources} sources")
        
        # Save final consolidated results
        self.save_results(all_samples, all_summary, output_file)
    
    def _save_chunk_results(self, chunk_samples: Dict, chunk_summary: List[Dict], 
                           output_file: str, chunk_num: int) -> None:
        """Save intermediate chunk results."""
        try:
            output_path = Path(output_file)
            output_dir = output_path.parent
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Save chunk samples
            if chunk_samples:
                chunk_samples_file = output_dir / f"{output_path.stem}_chunk_{chunk_num:04d}_samples.parquet"
                chunk_samples_df = pd.DataFrame([
                    {'source_id': source_id, 'samples': samples} 
                    for source_id, samples in chunk_samples.items()
                ])
                chunk_samples_df.to_parquet(chunk_samples_file, index=False)
            
            # Save chunk summary
            if chunk_summary:
                chunk_summary_file = output_dir / f"{output_path.stem}_chunk_{chunk_num:04d}_summary.csv"
                chunk_summary_df = pd.DataFrame(chunk_summary)
                chunk_summary_df.to_csv(chunk_summary_file, index=False)
            
            logger.debug(f"Saved chunk {chunk_num} results: {len(chunk_samples)} samples, "
                        f"{len(chunk_summary)} summaries")
            
        except Exception as e:
            logger.error(f"Failed to save chunk {chunk_num} results: {e}")

    def save_results(self, samples_by_source: Dict, summary_data: List[Dict], output_file: str) -> None:
        """
        Save all samples and summary data to files.
        
        Args:
            samples_by_source: Dictionary with source_id as key and list of samples as value
            summary_data: List of summary statistics for each source
            output_file: Base path for output files
        """
        output_path = Path(output_file)
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save samples as parquet file with source_id: [samples] format
        if samples_by_source:
            samples_file = output_dir / f"{output_path.stem}_samples.parquet"
            samples_df = pd.DataFrame([
                {'source_id': source_id, 'samples': samples} 
                for source_id, samples in samples_by_source.items()
            ])
            samples_df.to_parquet(samples_file, index=False)
            total_samples = sum(len(samples) for samples in samples_by_source.values())
            logger.info(f"Saved {total_samples} samples from {len(samples_by_source)} sources to {samples_file}")
        
        # Save summary data
        summary_file = output_dir / f"{output_path.stem}_summary.csv"
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Saved summary for {len(summary_data)} sources to {summary_file}")
    
    def print_summary(self) -> None:
        """Print processing summary."""
        elapsed = time.time() - self.start_time
        rate = self.total_sources / elapsed if elapsed > 0 else 0
        
        logger.info("="*50)
        logger.info("PROCESSING SUMMARY")
        logger.info("="*50)
        logger.info(f"Total sources processed: {self.total_sources}")
        logger.info(f"Successful: {self.successful_sources}")
        logger.info(f"Failed: {self.failed_sources}")
        logger.info(f"Success rate: {100*self.successful_sources/self.total_sources:.1f}%")
        logger.info(f"Total time: {elapsed/60:.1f} minutes")
        logger.info(f"Processing rate: {rate:.1f} sources/second")
        logger.info("="*50)

    def load_input_data(self, input_file: str) -> pd.DataFrame:
        """
        Load input data from CSV or Parquet file.
        
        Args:
            input_file: Path to input file
            
        Returns:
            DataFrame with input data
        """
        input_path = Path(input_file)
        file_ext = input_path.suffix.lower()
        
        logger.info(f"Loading input data from {input_file}")
        logger.info(f"Detected file format: {file_ext}")
        
        try:
            if file_ext == '.csv':
                df = pd.read_csv(input_file)
                logger.info(f"Successfully loaded CSV file with {len(df)} rows")
            elif file_ext in ['.parquet', '.pq']:
                df = pd.read_parquet(input_file)
                logger.info(f"Successfully loaded Parquet file with {len(df)} rows")
            else:
                # Try to auto-detect based on content
                try:
                    df = pd.read_csv(input_file)
                    logger.info(f"Auto-detected as CSV file with {len(df)} rows")
                except Exception:
                    try:
                        df = pd.read_parquet(input_file)
                        logger.info(f"Auto-detected as Parquet file with {len(df)} rows")
                    except Exception as e:
                        raise ValueError(f"Unable to read file as CSV or Parquet: {str(e)}")
            
            if len(df) == 0:
                raise ValueError("Input file is empty")
                
            return df
            
        except Exception as e:
            logger.error(f"Error loading input file {input_file}: {str(e)}")
            raise

    def load_prior_data(self) -> List[List[str]]:
        """
        Load prior data from CSV file.
        
        Returns:
            List of rows from the prior summary CSV file
        """
        try:
            prior_path = Path(self.prior_file)
            if not prior_path.exists():
                raise FileNotFoundError(f"Prior summary file not found: {self.prior_file}")
            
            with open(self.prior_file, 'r', newline='') as f:
                reader = csv.reader(f)
                rows = list(reader)
            
            logger.info(f"Loaded prior data from {self.prior_file} with {len(rows)-1} HEALPix pixels")
            return rows
            
        except Exception as e:
            logger.error(f"Failed to load prior data: {e}")
            raise
    
    def calculate_healpix_from_source_id(self, source_id: int) -> int:
        """
        Calculate HEALPix pixel from Gaia DR3 source_id.
        
        Args:
            source_id: Gaia DR3 source_id
            
        Returns:
            HEALPix level 5 pixel number
        """
        # Formula from original code: hp = math.floor(source_id / (2**(35)*4**(12-5)))
        hp = math.floor(source_id / (2**35 * 4**7))
        return hp
    
    def calculate_healpix_from_coordinates(self, ra: float, dec: float) -> int:
        """
        Calculate HEALPix pixel from RA/Dec coordinates.
        
        Args:
            ra: Right ascension in degrees
            dec: Declination in degrees
            
        Returns:
            HEALPix level 5 pixel number
        """
        if not HEALPIX_AVAILABLE:
            raise ImportError("astropy_healpix is required for coordinate-based HEALPix calculations")
        
        HP = HEALPix(nside=2**5, order='nested')  # level 5
        hp = HP.lonlat_to_healpix(ra * u.deg, dec * u.deg)
        return int(hp)
    
    def get_prior_parameters(self, hp: int) -> Dict[str, float]:
        """
        Get prior parameters for a given HEALPix pixel.
        
        Args:
            hp: HEALPix pixel number
            
        Returns:
            Dictionary with prior parameters
        """
        try:
            # Index is hp+1 because first row is header
            if hp + 1 >= len(self.prior_data):
                raise ValueError(f"HEALPix pixel {hp} not found in prior data")
            
            row = self.prior_data[hp + 1]
            
            # Column indices based on the CSV structure:
            # 0: healpix, 1: glon, 2: glat, 3: Nstar, 4: Gmaglim, 
            # 5: GGDrlen, 6: GGDalpha, 7: GGDbeta, 8: GGDmedian, 9: GGDmode,
            # 10: EDSDrlen, 11: QGmin, 12: QGmax, 13: bprpmin, 14: bprpmax, 15: bprpnbins, 16: Nnofit
            
            return {
                'hp': int(row[0]),
                'glon': float(row[1]),
                'glat': float(row[2]),
                'ggd_rlen': float(row[5]),  # Column 5: GGDrlen
                'alpha': float(row[6]),     # Column 6: GGDalpha  
                'beta': float(row[7]),      # Column 7: GGDbeta
                'edsd_rlen': float(row[10]) # Column 10: EDSDrlen
            }
            
        except (IndexError, ValueError) as e:
            logger.error(f"Error getting prior parameters for HEALPix {hp}: {e}")
            raise


def init_worker(seed_base: int, log_level: str):
    """
    Initialize worker process for multiprocessing.
    
    Args:
        seed_base: Base seed for random number generation
        log_level: Logging level string
    """
    # Set unique seed for each worker
    worker_seed = seed_base + os.getpid()
    np.random.seed(worker_seed)
    
    # Configure logging for worker
    logging.basicConfig(
        level=getattr(logging, log_level),
        format=f'[Worker {os.getpid()}] %(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    
    # Suppress numpy errors
    np.seterr(divide='ignore')


def process_source_worker(args: Tuple) -> Dict:
    """
    Worker function to process a single source (for multiprocessing).
    
    Args:
        args: Tuple containing (row_data, model, nsamp, nburnin, prior_data)
        
    Returns:
        Dictionary with processing results
    """
    row_data, model, nsamp, nburnin, prior_data = args
    
    # Create a temporary estimator instance for this worker
    # Note: We pass an empty prior_file since we already have the prior_data
    estimator = DistanceEstimator(model, seed=0, prior_file=None)
    estimator.nsamp = nsamp
    estimator.nburnin = nburnin
    estimator.prior_data = prior_data
    
    # Convert row_data dict back to pandas Series
    row = pd.Series(row_data)
    
    return estimator.process_single_source(row)


@click.command()
@click.option('--input-file', '-i', required=True, type=click.Path(exists=True),
              help='Input CSV or Parquet file with source data')
@click.option('--output-file', '-o', required=True, type=click.Path(),
              help='Output file prefix (will create _samples.parquet and _summary.csv)')
@click.option('--model', '-m', required=True,
              type=click.Choice(['EDSD', 'GGD', 'Photogeometric'], case_sensitive=False),
              help='Distance model to use')
@click.option('--seed', '-s', default=42, type=int,
              help='Random seed for reproducibility (default: 42)')
@click.option('--nsamp', default=5000, type=int,
              help='Number of MCMC samples (default: 5000)')
@click.option('--nburnin', default=500, type=int,
              help='Number of MCMC burn-in samples (default: 500)')
@click.option('--log-level', default='INFO',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              help='Logging level (default: INFO)')
@click.option('--prior-file', default='prior_summary.csv', type=click.Path(exists=True),
              help='Path to prior summary CSV file (default: prior_summary.csv)')
@click.option('--validate-only', is_flag=True,
              help='Only validate input file, do not process')
@click.option('--n-processes', '-p', default=1, type=int,
              help='Number of parallel processes to use (default: 1, no multiprocessing)')
@click.option('--chunk-size', '-c', default=1000, type=int,
              help='Number of sources to process before saving intermediate results (default: 1000)')
@click.option('--start-index', default=0, type=int,
              help='Starting index for processing (inclusive, default: 0)')
@click.option('--end-index', default=None, type=int,
              help='Ending index for processing (exclusive, default: None = process to end)')
@click.option('--use-optimized', is_flag=True,
              help='Use optimized photogeometric processing for better parallelization')
def main(input_file, output_file, model, seed, nsamp, nburnin, log_level, prior_file, 
         validate_only, n_processes, chunk_size, start_index, end_index, use_optimized):
    """
    Standalone distance estimation tool for stellar parallax data.
    
    Processes a CSV or Parquet file containing stellar parallax measurements and estimates
    distances using EDSD, GGD, or Photogeometric models via MCMC sampling.
    
    NEW FEATURES:
    - Multiprocessing: Use --n-processes to enable parallel processing
    - Chunking: Use --start-index and --end-index to process specific data ranges
    - Periodic saving: Intermediate results saved every --chunk-size sources
    - Optimized photogeometric processing: Use --use-optimized for better parallelization performance
    
    Supported input formats:
    - CSV (.csv)
    - Parquet (.parquet, .pq)
    
    Required input columns:
    - source_id: Unique identifier for each source (Gaia DR3 source_id)
    - parallax: Parallax measurement in milliarcseconds
    - parallax_error: Parallax uncertainty in milliarcseconds
    
    Optional input columns:
    - ra: Right ascension in degrees (if provided, used with dec to calculate HEALPix)
    - dec: Declination in degrees (if provided, used with ra to calculate HEALPix)
    
    Model-specific additional required columns:
    - Photogeometric: phot_g_mean_mag, bp_rp, pseudocolour
    
    Note: HEALPix pixel numbers and prior parameters (rlen, alpha, beta) are calculated
    automatically from the prior summary file based on source_id or coordinates.
    
    Output files:
    - {output_file}_samples.parquet: All posterior samples grouped by source_id
      (format: source_id column with corresponding samples as list)
    - {output_file}_summary.csv: Summary statistics for each source
    - {output_file}_chunk_XXXX_samples.parquet: Intermediate chunk results
    - {output_file}_chunk_XXXX_summary.csv: Intermediate chunk summaries
    
    Examples:
    
    Basic usage:
      python distance_estimation_standalone.py -i data.csv -o results -m EDSD
    
    Parallel processing with 4 cores:
      python distance_estimation_standalone.py -i data.csv -o results -m GGD -p 4
    
    Process specific range (sources 1000-1999):
      python distance_estimation_standalone.py -i data.csv -o results -m EDSD --start-index 1000 --end-index 2000
    
    Large dataset with frequent saves every 500 sources:
      python distance_estimation_standalone.py -i data.csv -o results -m Photogeometric -p 8 -c 500
    
    Optimized photogeometric processing:
      python distance_estimation_standalone.py -i data.csv -o results -m Photogeometric -p 8 --use-optimized
    """
    # Set logging level
    logger.setLevel(getattr(logging, log_level))
    
    logger.info("Starting standalone distance estimation CLI")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output file: {output_file}")
    logger.info(f"Model: {model}")
    logger.info(f"Seed: {seed}")
    
    if n_processes > 1:
        logger.info(f"Multiprocessing enabled: {n_processes} processes")
    if start_index > 0 or end_index is not None:
        logger.info(f"Processing range: {start_index} to {end_index if end_index else 'end'}")
    logger.info(f"Chunk size for periodic saving: {chunk_size}")
    
    if use_optimized and model == 'Photogeometric':
        logger.info("Optimized photogeometric processing enabled")
    elif use_optimized and model != 'Photogeometric':
        logger.warning("--use-optimized flag ignored (only available for Photogeometric model)")
    
    try:
        # Validate parameters
        if n_processes < 1:
            raise ValueError("Number of processes must be >= 1")
        if chunk_size < 1:
            raise ValueError("Chunk size must be >= 1")
        if start_index < 0:
            raise ValueError("Start index must be >= 0")
        if end_index is not None and end_index <= start_index:
            raise ValueError("End index must be greater than start index")
        
        # Initialize estimator
        estimator = DistanceEstimator(
            model=model, 
            seed=seed, 
            prior_file=prior_file,
            n_processes=n_processes,
            chunk_size=chunk_size,
            log_level=log_level
        )
        estimator.nsamp = nsamp
        estimator.nburnin = nburnin
        
        # Set optimized processing flag
        if use_optimized and model == 'Photogeometric':
            estimator.use_optimized = True
        else:
            estimator.use_optimized = False
        
        # Load input data using the new method
        df = estimator.load_input_data(input_file)
        
        # Validate input data
        if not estimator.validate_input_data(df):
            logger.error("Input validation failed")
            sys.exit(1)
        
        # Validate index ranges
        if start_index >= len(df):
            raise ValueError(f"Start index {start_index} is beyond dataframe length {len(df)}")
        if end_index is not None and end_index > len(df):
            logger.warning(f"End index {end_index} is beyond dataframe length {len(df)}, "
                          f"will process to end ({len(df)})")
            end_index = len(df)
        
        if validate_only:
            logger.info("Validation successful. Exiting (--validate-only flag used).")
            return
        
        # Process data with new parallel/chunking functionality
        estimator.process_dataframe_parallel(df, output_file, start_index, end_index)
        
        logger.info("Distance estimation completed successfully")
        
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        logger.debug(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)


if __name__ == '__main__':
    main() 