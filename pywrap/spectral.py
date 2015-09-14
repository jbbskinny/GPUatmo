# File: spectral.py
# Author: Justin Blythe
#
# Created on March 11, 2013, 8:47 PM

###############################################################################
########################## Spectral Quantities  ###############################
###############################################################################

from numpy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Update jnu grid

nsatmo.update_jnu_grid.restype = None
nsatmo.update_jnu_grid.argtype = [c_int, c_int, c_int, POINTER(c_double), \
                                  POINTER(c_double), POINTER(c_double), \
                                  POINTER(c_double)]

def update_jnu_grid(D, K, M, w, iminus, iplus, jnu):

    D_c      = c_int(D)
    K_c      = c_int(K)
    M_c      = c_int(M)
    w_c      = w.ctypes.data_as(POINTER(c_double))
    iminus_c = (POINTER(c_double) * len(iminus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iminus])
    iplus_c  = (POINTER(c_double) * len(iplus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iplus])
    jnu_c    = (POINTER(c_double) * len(jnu))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in jnu])

    nsatmo.update_jnu_grid(D_c, K_c, M_c, w_c, iminus_c, iplus_c, jnu_c)

# Update j profile

nsatmo.update_j_profile.restype = None
nsatmo.update_j_profile.argtype = [c_int, c_int, POINTER(c_double), \
                                   POINTER(c_double), POINTER(c_double)]

def update_j_profile(D, K, h, j, jnu):

    D_c   = c_int(D)
    K_c   = c_int(K)
    h_c   = h.ctypes.data_as(POINTER(c_double))
    j_c   = j.ctypes.data_as(POINTER(c_double))
    jnu_c = (POINTER(c_double) * len(jnu))\
            (*[row.ctypes.data_as(POINTER(c_double)) for row in jnu])

    nsatmo.update_j_profile(D_c, K_c, h_c, j_c, jnu_c)

# Update chiJ profile

nsatmo.update_chiJ_profile.restype = None
nsatmo.update_chiJ_profile.argtype = [c_int, c_int, POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double), \
                                      POINTER(c_double)]

def update_chiJ_profile(D, K, chiJ, h, j, chi, jnu):

    D_c    = c_int(D)
    K_c    = c_int(K)
    chiJ_c = chiJ.ctypes.data_as(POINTER(c_double))
    h_c    = h.ctypes.data_as(POINTER(c_double))
    j_c    = j.ctypes.data_as(POINTER(c_double))
    chi_c  = (POINTER(c_double) * len(chi))\
             (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])
    jnu_c  = (POINTER(c_double) * len(jnu))\
             (*[row.ctypes.data_as(POINTER(c_double)) for row in jnu])

    nsatmo.update_chiJ_profile(D_c, K_c, chiJ_c, h_c, j_c, chi_c, jnu_c)

# Update fnu grid

nsatmo.update_fnu_grid.restype = None
nsatmo.update_fnu_grid.argtype = [c_int, c_int, c_int, POINTER(c_double), \
                                  POINTER(c_double), POINTER(c_double), \
                                  POINTER(c_double), POINTER(c_double)]

def update_fnu_grid(D, K, M, mu, w, fnu, iminus, iplus):

    D_c      = c_int(D)
    K_c      = c_int(K)
    M_c      = c_int(M)
    mu_c     = mu.ctypes.data_as(POINTER(c_double))
    w_c      = w.ctypes.data_as(POINTER(c_double))
    fnu_c    = (POINTER(c_double) * len(fnu))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in fnu])
    iminus_c = (POINTER(c_double) * len(iminus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iminus])
    iplus_c  = (POINTER(c_double) * len(iplus))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in iplus])

    nsatmo.update_fnu_grid(D_c, K_c, M_c, mu_c, w_c, fnu_c, iminus_c, \
                           iplus_c)

# Update fbar profile

nsatmo.update_fbar_profile.restype = None
nsatmo.update_fbar_profile.argtype = [c_int, c_int, POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double)]

def update_fbar_profile(D, K, fbar, h, fnu):

    D_c    = c_int(D)
    K_c    = c_int(K)
    fbar_c = fbar.ctypes.data_as(POINTER(c_double))
    h_c    = h.ctypes.data_as(POINTER(c_double))
    fnu_c  = (POINTER(c_double) * len(fnu))\
             (*[row.ctypes.data_as(POINTER(c_double)) for row in fnu])

    nsatmo.update_fbar_profile(D_c, K_c, fbar_c, h_c, fnu_c)

# Update chiF profile

nsatmo.update_chiF_profile.restype = None
nsatmo.update_chiF_profile.argtype = [c_int, c_int, POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double)]

def update_chiF_profile(D, K, chiF, fbar, h, chi, fnu):

    D_c    = c_int(D)
    K_c    = c_int(K)
    chiF_c = chiF.ctypes.data_as(POINTER(c_double))
    fbar_c = fbar.ctypes.data_as(POINTER(c_double))
    h_c    = h.ctypes.data_as(POINTER(c_double))
    chi_c  = (POINTER(c_double) * len(chi))\
             (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])
    fnu_c  = (POINTER(c_double) * len(fnu))\
             (*[row.ctypes.data_as(POINTER(c_double)) for row in fnu])

    nsatmo.update_chiF_profile(D_c, K_c, chiF_c, fbar_c, h_c, chi_c, fnu_c)

# Update b profile

nsatmo.update_b_profile.restype = None
nsatmo.update_b_profile.argtype = [c_int, POINTER(c_double), POINTER(c_double)]

def update_b_profile(D, b, T6):

    D_c  = c_int(D)
    b_c  = b.ctypes.data_as(POINTER(c_double))
    T6_c = T6.ctypes.data_as(POINTER(c_double))

    nsatmo.update_b_profile(D_c, b_c, T6_c)

# Update chiP profile

nsatmo.update_chiP_profile.restype = None
nsatmo.update_chiP_profile.argtype = [c_int, c_int, POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double), \
                                      POINTER(c_double), POINTER(c_double)]

def update_chiP_profile(D, K, b, chiP, h, chi, source):

    D_c      = c_int(D)
    K_c      = c_int(K)
    b_c      = b.ctypes.data_as(POINTER(c_double))
    chiP_c   = chiP.ctypes.data_as(POINTER(c_double))
    h_c      = h.ctypes.data_as(POINTER(c_double))
    chi_c    = (POINTER(c_double) * len(chi))\
               (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])
    source_c = (POINTER(c_double) * len(source))\
                (*[row.ctypes.data_as(POINTER(c_double)) for row in source])

    nsatmo.update_chiP_profile(D_c, K_c, b_c, chiP_c, h_c, chi_c, source_c)
