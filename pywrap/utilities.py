# File: utilities.py
# Author: Justin Blythe
#
# Created January 17, 2013, 4:07pm

###############################################################################
############################### Utilities #####################################
###############################################################################

from numpy import *
from ctypes import *
import os.path

curr_dir = os.path.dirname(os.path.realpath(__file__))
nsatmo = ctypeslib.load_library('libnsatmo.so', curr_dir+'/../build')

# Allocate memory for 2d array

nsatmo.malloc_double_2d.restype = POINTER(POINTER(c_double))
nsatmo.malloc_double_2d.argtype = [c_int, c_int]

def malloc_double_2d(dim1, dim2):

    dim1_c = c_int(dim1)
    dim2_c = c_int(dim2)

    pointer = nsatmo.malloc_double_2d(dim1_c, dim2_c)
    return(pointer)

# Free 2d array

nsatmo.free_2d.restype = None
nsatmo.free_2d.argtype = [c_int, POINTER(POINTER(c_double))]

def free_2d(dim1, pointer_c):

    dim1_c = c_int(dim1)

    nsatmo.free_2d(dim1_c, pointer_c)

# Allocate memory for 1d array

nsatmo.malloc_double_1d.restype = POINTER(c_double)
nsatmo.malloc_double_1d.argtype = c_int

def malloc_double_1d(dim1):

    dim1_c = c_int(dim1)

    pointer = nsatmo.malloc_double_1d(dim1_c)

    return(pointer)

# Free 1d array

nsatmo.free_1d.restype = None
nsatmo.free_1d.argtype = POINTER(c_double)

def free_1d(pointer):

    pointer_c = pointer.ctypes.data_as(POINTER(c_double))

    nsatmo.free_1d(pointer_c)

# Maximum relative deviation calculator

nsatmo.residual.restype = c_double
nsatmo.residual.argtype = [c_int, c_double, c_double, POINTER(c_double), \
                           POINTER(c_double)]

def residual(D, relerr, abserr, old, new):

    D_c      = c_int(D)
    relerr_c = c_double(relerr)
    abserr_c = c_double(abserr)
    old_c    = old.ctypes.data_as(POINTER(c_double))
    new_c    = new.ctypes.data_as(POINTER(c_double))

    result   = nsatmo.residual(D_c, relerr_c, abserr_c, old_c, new_c)
    return result

# Array copying utility

nsatmo.copy_array.restype = None
nsatmo.copy_array.argtypes = [c_int, POINTER(c_double), POINTER(c_double)]

def copy_array(size, target, source):

    size_c   = c_int(size)
    target_c = target.ctypes.data_as(POINTER(c_double))
    source_c = source.ctypes.data_as(POINTER(c_double))

    nsatmo.copy_array(size_c, target_c, source_c)
