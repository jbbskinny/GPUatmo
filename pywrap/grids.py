# File: grids.py
# Author: Justin Blythe
#
# Created on January 5, 2013

###############################################################################
############################ Grid Initializer  ################################
###############################################################################

from numpy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Setup Thompson optical depth (tauTh) Grid

nsatmo.setup_tauTh_grid.restype = None
nsatmo.setup_tauTh_grid.argtype = [c_int, c_int, c_double, POINTER(c_double)]
def setup_tauTh_grid(decades, ppd, init, tauTh):

    decades_c = c_int(decades)
    ppd_c     = c_int(ppd)
    init_c    = c_double(init)
    tauTh_c   = tauTh.ctypes.data_as(POINTER(c_double))

    nsatmo.setup_tauTh_grid(decades_c, ppd_c, init_c, tauTh_c)

# Setup Energy (E1) grid

nsatmo.setup_E1_grid.restype = None
nsatmo.setup_E1_grid.argtype =  [c_int, c_int, c_double, POINTER(c_double), \
                                 POINTER(c_double)]
def setup_E1_grid(decades, ppd, init, E1, h):

    decades_c = c_int(decades)
    ppd_c     = c_int(ppd)
    init_c    = c_double(init)
    E1_c      = E1.ctypes.data_as(POINTER(c_double))
    h_c       = h.ctypes.data_as(POINTER(c_double))

    nsatmo.setup_E1_grid(decades_c, ppd_c, init_c, E1_c, h_c)

#Setup Angular (mu) grid

def setup_mu_grid(M, mu, w):

    M_c  = c_int(M)
    mu_c = mu.ctypes.data_as(POINTER(c_double))
    w_c  = w.ctypes.data_as(POINTER(c_double))

    nsatmo.setup_mu_grid(M_c, mu_c, w_c)

# Update Density Profile

nsatmo.update_density_profile.restype = None
nsatmo.update_density_profile.argtype = [c_int, c_double, POINTER(c_double),\
                                         POINTER(c_double), POINTER(c_double)]

def update_density_profile(D, g, rho, T6, tauTh):

    D_c     = c_int(D)
    g_c     = c_double(g)
    rho_c   = rho.ctypes.data_as(POINTER(c_double))
    T6_c    = T6.ctypes.data_as(POINTER(c_double))
    tauTh_c = tauTh.ctypes.data_as(POINTER(c_double))

    nsatmo.update_density_profile(D_c, g_c, rho_c, T6_c, tauTh_c)
