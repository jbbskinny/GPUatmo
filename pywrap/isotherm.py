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

chiJ   = array([0.0] * D)
j      = array([0.0] * D)
chiF   = array([0.0] * D)
fbar   = array([0.0] * D)
chiP   = array([0.0] * D)
b      = array([0.0] * D)
dT6    = array([0.0] * D)

chik1  = array([0.0] * K)
chik2  = array([0.0] * K)
chik3  = array([0.0] * K)
ipk1   = array([0.0] * K)
ipk2   = array([0.0] * K)
ipk3   = array([0.0] * K)
imk1   = array([0.0] * K)
imk2   = array([0.0] * K)
imk3   = array([0.0] * K)
bk1    = array([0.0] * K)
bk2    = array([0.0] * K)
bk3    = array([0.0] * K)
ipd1   = array([0.0] * D)
ipd2   = array([0.0] * D)
ipd3   = array([0.0] * D)
imd1   = array([0.0] * D)
imd2   = array([0.0] * D)
imd3   = array([0.0] * D)
bd1    = array([0.0] * D)
bd2    = array([0.0] * D)
bd3    = array([0.0] * D)

E1     = array([0.0] * K)
h      = array([0.0] * K)
mu     = array([0.0] * M)
w      = array([0.0] * M)
rho    = array([0.0] * D)
tauTh  = array([0.0] * D) 
T6     = array([0.0] * D)
oldT6  = array([0.0] * D)
flux   = array([0.0] * D)
iminus = array([[0.0] * I] * D)
iplus  = array([[0.0] * I] * D)
source = array([[0.0] * I] * D)
chi    = array([[0.0] * K] * D)
f1     = array([[0.0] * M] * K)
f2     = array([[0.0] * K] * D)

setup_tauTh_grid(decadestauTh, ppdtauTh, tauThmin, tauTh)

setup_E1_grid(decadesE1, ppdE1, E1min, E1, h)

setup_mu_grid(M, mu, w)

T6[:] = teff

update_density_profile(D, g, rho, T6, tauTh)

#update_chi_grid(D, K, E1, rho, T6, chi_c)

setup_source_func(D, I, K, E1, T6, source)

solve_rte(D, I, K, M, E1, mu, rho, T6, tauTh, iminus, iplus, source)

#for l in range(0, maxiteration):
##for l in range(0, 1):
#
#    calc_chiJ(D, K, M, chiJ, E1, h, j, mu, w, chi_c, iminus_c, iplus_c)
#    calc_chiF(D, K, M, chiF, E1, fbar, h, mu, w, chi_c, iminus_c, iplus_c)
#    calc_chiP(D, K, b, chiP, E1, h, T6, chi_c, source_c)
#
#    fbarc  = 4.51245e-4 * teff**4.0
#    delta  = 1.15130 * (log10(tauTh[1]) - log10(tauTh[0]))
#
#    for d in range(0, D):
#        intval = 0.0
#        if d == 0:
#            y = 0.5 * (fbarc - fbar[0]) * (chiJ[d] / chiP[d]) / b[d]
#            z = 0.25 * (chiJ[d] / chiP[d]) * (j[d] / b[d]) - 0.25
#            dT6[d]   = T6[d] * (y + z)
#            oldT6[d] = T6[d]
#            T6[d]   += dT6[d]
#            continue
#
#        x = 0.75 * (chiJ[d] / chiP[d]) / b[d]
#        y = 0.5 * (fbarc - fbar[0]) * (chiJ[d] / chiP[d]) / b[d]
#        z = 0.25 * (chiJ[d] / chiP[d]) * (j[d] / b[d]) - 0.25
#
#        for dd in range(1, d):
#            intval += delta * (tauTh[dd] * chiF[dd] * (fbarc - fbar[dd]) + \
#                      tauTh[dd-1] * chiF[dd-1] * (fbarc - fbar[dd-1]))
#        dT6[d]   = T6[d] * (x * intval + y + z)
#        oldT6[d] = T6[d]
#        T6[d]   += dT6[d]
#
#    update_density_profile(D, g, rho, T6, tauTh)
#    update_chi_grid(D, K, E1, rho, T6, chi_c)
#    setup_source_func(D, I, K, E1, T6, source_c)
#    solve_rte(D, I, K, M, E1, mu, rho, T6, tauTh, iminus_c, iplus_c, source_c)
#    
#    print l, max(abs(fbar-fbarc)/fbarc)
#    if max(abs(fbar - fbarc) / fbarc) < 0.01:
#        break;
#
#delta = 1.15130 * (log10(tauTh[1]) - log10(tauTh[0]))
#td1   = array([0.0] * D)
#td2   = array([0.0] * D)
#td3   = array([0.0] * D)
#taunu = array([[0.0] * K] * D)
#
#for d in range(0, D):
#    for i in range(0, I):
#        k = mod(i, K)
#        iplus[d][i]  = iplus_c[d][i]
#        iminus[d][i] = iminus_c[d][i]
#        source[d][i] = source_c[d][i]
#        chi[d][k]    = chi_c[d][k]
#
#    ipd1[d] = iplus[d][3*K+(K/3)]
#    ipd2[d] = iplus[d][3*K+(2*K/3)]
#    ipd3[d] = iplus[d][3*K+(K-1)]
#    imd1[d] = iminus[d][3*K+(K/3)]
#    imd2[d] = iminus[d][3*K+(2*K/3)]
#    imd3[d] = iminus[d][3*K+(K-1)]
#    bd1[d]  = source[d][3*K+(K/3)]
#    bd2[d]  = source[d][3*K+(2*K/3)]
#    bd3[d]  = source[d][3*K+(K-1)]
#
#for d in range(0, D):
#    for k in range(0, K):
#        for m in range(0, M):
#            f1[k][m] = 0.5 * mu[m] * (iplus[d][K*m+k] - iminus[d][K*m+k])
#
#        f2[d][k] = ext_trap_quad(M, f1[k], w)
#
#        chik1[k] = chi[D/3][k]
#        chik2[k] = chi[2*D/3][k]
#        chik3[k] = chi[D-2][k]
#        ipk1[k] = iplus[1][12*K+k]
#        ipk2[k] = iplus[D/2][12*K+k]
#        ipk3[k] = iplus[D-2][12*K+k]
#        imk1[k] = iminus[1][12*K+k]
#        imk2[k] = iminus[D/2][12*K+k]
#        imk3[k] = iminus[D-2][12*K+k]
#        bk1[k]  = source[0][3*K+k]
#        bk2[k]  = source[D/2][3*K+k]
#        bk3[k]  = source[D-2][3*K+k]
#
#    flux[d] = ext_trap_quad(K, f2[d], h)
#
#free_2d(D, chi_c)
#free_2d(D, iminus_c)
#free_2d(D, iplus_c)
#free_2d(D, source_c)
