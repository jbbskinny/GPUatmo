# File: numerical.py
# Author: Matt van Adelsberg
#
# Created December of 2012

###############################################################################
############################ Numerical Routines ###############################
###############################################################################

from numpy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Clenshaw's algorithm

nsatmo.clenshaw.restype = c_double
nsatmo.clenshaw.argtype = [c_int, c_double, POINTER(c_double)]

def clenshaw(N, x, C):

    N_c = c_int(N)
    x_c = c_double(x)
    C_c = C.ctypes.data_as(POINTER(c_double))

    result = nsatmo.clenshaw(N_c, x_c, C_c)
    return result

# Gauss-Legendre quadrature

nsatmo.setup_gauss_legendre.restype = None
nsatmo.setup_gauss_legendre.argtype = [c_int, c_double, c_double, \
                                       POINTER(c_double), POINTER(c_double)]

def setup_gauss_legendre(L, init, final, absys, weight):

    L_c      = c_int(L)
    init_c   = c_double(init)
    final_c  = c_double(final)
    absys_c  = absys.ctypes.data_as(POINTER(c_double))
    weight_c = weight.ctypes.data_as(POINTER(c_double))

    nsatmo.setup_gauss_legendre(L_c, init_c, final_c, absys_c, weight_c)

# Extended trapezoidal rule

nsatmo.quadrature.restype = c_double
nsatmo.quadrature.argtype = [c_int, POINTER(c_double), POINTER(c_double)]

def quadrature(L, f, weight):

    L_c      = c_int(L)
    f_c      = f.ctypes.data_as(POINTER(c_double))
    weight_c = weight.ctypes.data_as(POINTER(c_double))

    result = nsatmo.quadrature(L_c, f_c, weight_c)
    return result

