# File: main.py
# Author: Justin Blythe
#
# Created on January 5, 2013, 3:57 pm

###############################################################################
######################## Spectrum of Neutron Star #############################
###############################################################################

from numpy import *
from ctypes import *
from pylab import *
from utilities import *
from extinction import *
from grids import *
from temperature import *
from radtransfer import *
from numerical import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

decadestauTh = 4
decadesE1    = 3
M            = 20
maxiteration = 1000
ppdE1        = 20
ppdtauTh     = 20
E1min        = 1.0e-3
g            = 1.0
tauThmin     = 1.0e-3
teff         = 1.0

D = decadestauTh * ppdtauTh
K = decadesE1 * ppdE1
I = K * M

chiR   = array([0.0] * D)
tauR   = array([0.0] * D)

E1     = array([0.0] * K)
h      = array([0.0] * K)
mu     = array([0.0] * M)
w      = array([0.0] * M)
rho    = array([0.0] * D)
tauTh  = array([0.0] * D) 
T6     = array([0.0] * D)
oldT6  = array([0.0] * D)
oldT61 = array([0.0] * D)
chi    = array([[0.0] * K] * D)

setup_tauTh_grid(decadestauTh, ppdtauTh, tauThmin, tauTh)

setup_E1_grid(decadesE1, ppdE1, E1min, E1, h)

setup_mu_grid(M, mu, w)

update_rossld_profile(D, teff, T6, tauTh)
update_density_profile(D, g, rho, T6, tauTh)
copy_array(D, oldT61, T6)
for i in range(0, 100):
    
    calc_rossld_extinc(D, K, chiR, E1, h, rho, T6)
    calc_rossld_depth(D, chiR, tauR, tauTh)
    update_rossld_profile(D, teff, T6, tauR)
    update_density_profile(D, g, rho, T6, tauTh)
    
    print residual(D, 0.0, 1e-4, oldT6, T6)

    if abs(residual(D, 0.0, 1e-4, oldT6, T6)) < 1.0:
        break;

    copy_array(D, oldT6, T6)

