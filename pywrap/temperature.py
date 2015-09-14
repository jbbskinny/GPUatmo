# File: temperature.py
# Author: Justin Blythe
#
# Created on January 5, 2013, 6:36 pm

###############################################################################
######################## Rosseland Temperature Profile ########################
###############################################################################

from numpy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Initialize Temperature Profile

nsatmo.rosseland_temp_profile.restype = None
nsatmo.rosseland_temp_profile.argtype = [c_int, c_int, c_int, c_double, \
                                         c_double, POINTER(c_double), \
                                         POINTER(c_double), POINTER(c_double),\
                                         POINTER(c_double), POINTER(c_double)]

def rosseland_temp_profile(D, K, maxiteration, g, teff, E1, h, rho, T6, tauTh):
    
    D_c            = c_int(D)
    K_c            = c_int(K)
    maxiteration_c = c_int(maxiteration)
    g_c            = c_double(g)
    teff_c         = c_double(teff)
    E1_c           = E1.ctypes.data_as(POINTER(c_double))
    h_c            = h.ctypes.data_as(POINTER(c_double))
    rho_c          = rho.ctypes.data_as(POINTER(c_double))
    T6_c           = T6.ctypes.data_as(POINTER(c_double))
    tauTh_c        = tauTh.ctypes.data_as(POINTER(c_double))

    nsatmo.rosseland_temp_profile(D_c, K_c, maxiteration_c, g_c, teff_c, E1_c,\
                                  h_c, rho_c, T6_c, tauTh_c)

# Update Rosseland temperature profile

nsatmo.update_rossld_profile.restype = None
nsatmo.update_rossld_profile.argtype = [c_int, c_double, POINTER(c_double),\
                                        POINTER(c_double)]

def update_rossld_profile(D, teff, T6, tau):

    D_c    = c_int(D)
    teff_c = c_double(teff)
    T6_c   = T6.ctypes.data_as(POINTER(c_double))
    tau_c  = tau.ctypes.data_as(POINTER(c_double))

    nsatmo.update_rossld_profile(D_c, teff_c, T6_c, tau_c)

# Calculate Rosslend Depth

nsatmo.calc_rossld_depth.restype = None
nsatmo.calc_rossld_depth.argtype = [c_int, POINTER(c_double), \
                                    POINTER(c_double), POINTER(c_double)]

def calc_rossld_depth(D, chiR, tauR, tauTh):

    D_c     = c_int(D)
    chiR_c  = chiR.ctypes.data_as(POINTER(c_double))
    tauR_c  = tauR.ctypes.data_as(POINTER(c_double))
    tauTh_c = tauTh.ctypes.data_as(POINTER(c_double))

    nsatmo.calc_rossld_depth(D_c, chiR_c, tauR_c, tauTh_c)

# Calculate Rosseland Extinction

nsatmo.calc_rossld_extinc.restype = None
nsatmo.calc_rossld_extinc.argtype = [c_int, c_int, POINTER(c_double), \
                                     POINTER(c_double), POINTER(c_double), \
                                     POINTER(c_double)] 

def calc_rossld_extinc(D, K, chiR, E1, h, rho, T6):

    D_c    = c_int(D)
    K_c    = c_int(K)
    chiR_c = chiR.ctypes.data_as(POINTER(c_double))
    E1_c   = E1.ctypes.data_as(POINTER(c_double))
    h_c    = h.ctypes.data_as(POINTER(c_double))
    rho_c  = rho.ctypes.data_as(POINTER(c_double))
    T6_c   = T6.ctypes.data_as(POINTER(c_double))

    nsatmo.calc_rossld_extinc(D_c, K_c, chiR_c, E1_c, h_c, rho_c, T6_c)
