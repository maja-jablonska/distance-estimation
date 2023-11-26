# Load all of the required packages 
 
import numpy as np
import csv 
import ipywidgets as widgets
from IPython.display import display
from scipy.stats import norm
import scipy.integrate
from scipy.special import gamma, factorial
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import math
import re

# Zeropoint-correction-package

from zero_point import zpt # need to install gaiadr3-zeropoint
zpt.load_tables()

# astropy-modules

from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy_healpix import HEALPix
from astropy.coordinates import Angle

# rpy2 modules: used to convert R to python code

import rpy2
from rpy2.robjects.packages import importr, data
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
from rpy2.robjects import r

numpy2ri.activate() #activate conversion from R-code to numpy arrays