# File: radtransfer.py
# Author: Justin Blythe
#
# Created on January 5, 2013, 6:56 pm

###############################################################################
######################### Radiative Transfer Equation #########################
###############################################################################

from numpy import *
from scipy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Solve Radiative Transfer Equation

nsatmo.solve_rte.restype = None
nsatmo.solve_rte.argtype = [c_int, c_int, c_int, c_int, c_double, \
                            POINTER(c_double), POINTER(c_double), \
                            POINTER(c_double), POINTER(c_double), \
                            POINTER(c_double), POINTER(c_double), \
                            POINTER(c_double), POINTER(c_double), \
                            POINTER(c_double)]

def solve_rte(D, I, K, M, teff, E1, h, mu, T6, chi, dtau, dtaum, dtaup, 
              iminus, iplus, source):
    
    D_c      = c_int(D)
    I_c      = c_int(I)
    K_c      = c_int(K)
    M_c      = c_int(M)
    teff_c   = c_double(teff)
    E1_c     = E1.ctypes.data_as(POINTER(c_double))
    h_c      = h.ctypes.data_as(POINTER(c_double))
    mu_c     = mu.ctypes.data_as(POINTER(c_double))
    T6_c     = T6.ctypes.data_as(POINTER(c_double))
    chi_c    = (POINTER(c_double) * len(chi))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])
    dtau_c   = (POINTER(c_double) * len(dtau))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtau])
    dtaum_c  = (POINTER(c_double) * len(dtaum))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtaum])
    dtaup_c  = (POINTER(c_double) * len(dtaup))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtaup])
    iminus_c = (POINTER(c_double) * len(iminus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iminus])
    iplus_c  = (POINTER(c_double) * len(iplus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iplus])
    source_c = (POINTER(c_double) * len(source))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in source])

    nsatmo.solve_rte(D_c, I_c, K_c, M_c, teff_c, E1_c, h_c, mu_c, T6_c, chi_c, 
                     dtau_c, dtaum_c, dtaup_c, iminus_c, iplus_c, source_c)

# Setup optical depth grid spacing

nsatmo.update_dtau_grids.restype = None
nsatmo.update_dtau_grids.argtype = [c_int, c_int, c_int, POINTER(c_double), \
                                    POINTER(c_double), POINTER(c_double), \
                                    POINTER(c_double), POINTER(c_double), \
                                    POINTER(c_double), POINTER(c_double), \
                                    POINTER(c_double), POINTER(c_double)]

def update_dtau_grids(D, K, M, mu, tauTh, chi, dtau, dtaum, dtaup):

    D_c     = c_int(D)
    K_c     = c_int(K)
    M_c     = c_int(M)
    mu_c    = mu.ctypes.data_as(POINTER(c_double))
    tauTh_c = tauTh.ctypes.data_as(POINTER(c_double))
    chi_c = (POINTER(c_double) * len(chi))\
            (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])
    dtau_c = (POINTER(c_double) * len(dtau))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtau])
    dtaum_c = (POINTER(c_double) * len(dtaum))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtaum])
    dtaup_c = (POINTER(c_double) * len(dtaup))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in dtaup])

    nsatmo.update_dtau_grids(D_c, K_c, M_c, mu_c, tauTh_c, chi_c, dtau_c, 
                             dtaum_c, dtaup_c)

# Initialize Source Function

nsatmo.update_source_grid.restype = None
nsatmo.update_source_grid.argtype = [c_int, c_int, c_int, \
                                     POINTER(c_double), POINTER(c_double), \
                                     POINTER(c_double)]
    
def update_source_grid(D, I, K, M, E1, T6, source):

    D_c      = c_int(D)
    I_c      = c_int(I)
    K_c      = c_int(K)
    M_c      = c_int(M)
    E1_c     = E1.ctypes.data_as(POINTER(c_double))
    T6_c     = T6.ctypes.data_as(POINTER(c_double))
    source_c = (POINTER(c_double) * len(source))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in source])

    nsatmo.update_source_grid(D_c, I_c, K_c, M_c, E1_c, T6_c, source_c)
            
# Planck Function

nsatmo.planck_func.restype = c_double
nsatmo.planck_func.argtype = [c_double, c_double]

def planck_func(E1, T6):

    E1_c = c_double(E1)
    T6_c = c_double(T6)

    nsatmo.planck_func(E1_c, T6_c)

# Calculate dBdtau grid

nsatmo.calc_dBdtau_grid.restype = None
nsatmo.calc_dBdtau_grid.argtype = [c_int, c_int, c_int, c_double, \
                                   POINTER(c_double), POINTER(c_double), \
                                   POINTER(c_double), POINTER(c_double), \
                                   POINTER(c_double), POINTER(c_double)]

def calc_dBdtau_grid(D, K, M, teff, dBdtau, E1, h, mu, T6, chi):

    D_c      = c_int(D)
    K_c      = c_int(K)
    M_c      = c_int(M)
    teff_c   = c_double(teff)
    dBdtau_c = dBdtau.ctypes.data_as(POINTER(c_double))
    E1_c     = E1.ctypes.data_as(POINTER(c_double))
    h_c      = h.ctypes.data_as(POINTER(c_double))
    mu_c     = mu.ctypes.data_as(POINTER(c_double))
    T6_c     = T6.ctypes.data_as(POINTER(c_double))
    chi_c    = (POINTER(c_double) * len(chi))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])

    nsatmo.calc_dBdtau_grid(D_c, K_c, M_c, teff_c, dBdtau_c, E1_c, h_c, mu_c,
                            T6_c, chi_c)
