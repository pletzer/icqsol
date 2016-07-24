#!/usr/bin/env python

from __future__ import print_function
import pkg_resources
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long
import glob
import sys

PY_MAJOR_VERSION = sys.version_info[0]

# On some distributions the fully qualified shared library name 
# includes suffixes such as '.cpython-35m-darwin.so'
def getFullyQualifiedSharedLibraryName(libName):
    paths = glob.glob(libName + '*')
    # return the first path
    return paths[0]

def getSharedLibraryName(libName):
    if PY_MAJOR_VERSION < 3:
        res = pkg_resources.resource_filename('icqsol', libName + '.so')
    else:
        libName = pkg_resources.resource_filename('icqsol', libName)
        res = getFullyQualifiedSharedLibraryName(libName)
    return res




