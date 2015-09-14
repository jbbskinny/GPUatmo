# File: main.py
# Author: Justin Blythe
#
# Created on January 5, 2013, 3:57 pm

###############################################################################
######################## Spectrum of Neutron Star #############################
###############################################################################

from numpy import *
from scipy.interpolate import *
from ctypes import *
from pylab import *
from utilities import *
from extinction import *
from grids import *
from temperature import *
from radtransfer import *
from numerical import *
from spectral import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

g            = 2.4
teff         = 1.0
maxiteration = 50

decadestauTh = 8
ppdtauTh     = 40
tauThmin     = 1.0e-6

decadesE1    = 5
ppdE1        = 20
E1min        = 1.0e-4

D = decadestauTh * ppdtauTh
K = decadesE1 * ppdE1
M = 20
I = K * M

dT6    = array([0.0] * D)
j      = array([0.0] * D)
fbar   = array([0.0] * D)
b      = array([0.0] * D)
chiJ   = array([0.0] * D)
chiF   = array([0.0] * D)
chiP   = array([0.0] * D)
jnu    = array([[0.0] * K] * D)
fnu    = array([[0.0] * K] * D)

E1     = array([0.0] * K)
h      = array([0.0] * K)
mu     = array([0.0] * M)
w      = array([0.0] * M)
rho    = array([0.0] * D)
tauTh  = array([0.0] * D) 
T6     = array([0.0] * D)
dBdtau = array([0.0] * I)
iminus = array([[0.0] * I] * D)
iplus  = array([[0.0] * I] * D)
source = array([[0.0] * I] * D)
dtau   = array([[0.0] * I] * D)
dtaum  = array([[0.0] * I] * D)
dtaup  = array([[0.0] * I] * D)
chi    = array([[0.0] * K] * D)

setup_tauTh_grid(decadestauTh, ppdtauTh, tauThmin, tauTh)

setup_E1_grid(decadesE1, ppdE1, E1min, E1, h)

setup_mu_grid(M, mu, w)

rosseland_temp_profile(D, K, maxiteration, g, teff, E1, h, rho, T6, tauTh)

update_chi_grid(D, K, E1, rho, T6, chi)

update_dtau_grids(D, K, M, mu, tauTh, chi, dtau, dtaum, dtaup)

update_source_grid(D, I, K, M, E1, T6, source)

solve_rte(D, I, K, M, teff, E1, h, mu, T6, chi, dtau, dtaum, dtaup, iminus, 
          iplus, source)

update_jnu_grid(D, K, M, w, iminus, iplus, jnu)
update_fnu_grid(D, K, M, mu, w, fnu, iminus, iplus)
update_j_profile(D, K, h, j, jnu)
update_fbar_profile(D, K, fbar, h, fnu)
update_b_profile(D, b, T6)
update_chiJ_profile(D, K, chiJ, h, j, chi, jnu)
update_chiF_profile(D, K, chiF, fbar, h, chi, fnu)
update_chiP_profile(D, K, b, chiP, h, chi, source)

fbarc  = 4.51245e-4 * teff**4.0
delta  = (log10(tauTh[1]) - log10(tauTh[0]))

#for l in range(0, maxiteration):
for l in range(0, 10):
    intval = 0.0

    for d in range(0, D):

        if d == 0:
            y = 2.0 * (chiJ[0] / chiP[0]) * (fbarc - fbar[0])
            z = (chiJ[0] / chiP[0]) * j[0] - b[0]
            dT6[d]   = T6[d] / (4.0 * b[0]) * (y + z)
            T6[d]   += dT6[d]
            continue
        
        x = 3.0 * (chiJ[d] / chiP[d])
        y = 2.0 * (fbarc - fbar[0]) * (chiJ[d] / chiP[d])
        z = (chiJ[d] / chiP[d]) * j[d] - b[d]
        
        intval += delta * (tauTh[d] * chiF[d] * (fbarc - fbar[d]) + \
                 tauTh[d-1] * chiF[d-1] * (fbarc - fbar[d-1]))
        
        dT6[d]   = (T6[d] / (4.0 * b[d])) * (x * intval + y + z)

        T6[d] += dT6[d]

    update_density_profile(D, g, rho, T6, tauTh)
    update_chi_grid(D, K, E1, rho, T6, chi)
    update_source_grid(D, I, K, M, E1, T6, source)
    update_dtau_grids(D, K, M, mu, tauTh, chi, dtau, dtaum, dtaup)
    solve_rte(D, I, K, M, teff, E1, h, mu, T6, chi, dtau, dtaum, dtaup, 
              iminus, iplus, source)

    update_jnu_grid(D, K, M, w, iminus, iplus, jnu)
    update_fnu_grid(D, K, M, mu, w, fnu, iminus, iplus)
    update_j_profile(D, K, h, j, jnu)
    update_fbar_profile(D, K, fbar, h, fnu)
    update_b_profile(D, b, T6)
    update_chiJ_profile(D, K, chiJ, h, j, chi, jnu)
    update_chiF_profile(D, K, chiF, fbar, h, chi, fnu)
    update_chiP_profile(D, K, b, chiP, h, chi, source)
    
    print l, max(abs(fbar - fbarc) / fbarc)
#    if max(abs(fbar - fbarc) / fbarc) < 0.01:
#        break

if l == maxiteration - 1:
    print "ERROR, no convergence!"
    exit(-1);
