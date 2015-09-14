# File: extinction.py
# Author: Justin Blythe
#
# Created January 17, 2013, 4:55pm

from numpy import *
from ctypes import *
import os

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Update chi grid

nsatmo.update_chi_grid.restype = None
nsatmo.update_chi_grid.argtype = [c_int, c_int, POINTER(c_double), \
                                 POINTER(c_double), POINTER(c_double), \
                                 POINTER(c_double)]

def update_chi_grid(D, K, E1, rho, T6, chi):
    
    D_c = c_int(D)
    K_c = c_int(K)
    E1_c = E1.ctypes.data_as(POINTER(c_double))
    rho_c = rho.ctypes.data_as(POINTER(c_double))
    T6_c = T6.ctypes.data_as(POINTER(c_double))
    chi_c = (POINTER(c_double) * len(chi))\
            (*[row.ctypes.data_as(POINTER(c_double)) for row in chi])

    nsatmo.update_chi_grid(D_c, K_c, E1_c, rho_c, T6_c, chi_c)

# Calculate Extinction Coefficient

nsatmo.chi_ff.restype = c_double
nsatmo.chi_ff.argtype = [c_double, c_double, c_double]

def chi_ff(E1, rho, T6):

    E1_c  = c_double(E1)
    rho_c = c_double(rho)
    T6_c  = c_double(T6)

    nsatmo.chi_ff(E1_c, rho_c, T6_c)

nsatmo.g_ff.restype = c_double
nsatmo.g_ff.argtype = [c_double, c_double]

def g_ff(gamma2, u):

    gamma2_c = c_double(gamma2)
    u_c      = c_double(u)

    nsatmo.g_ff(gamma2_c, u_c)
