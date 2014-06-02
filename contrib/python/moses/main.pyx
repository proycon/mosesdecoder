from libcpp.string cimport string
from libcpp.vector cimport vector
import os
import cython
cimport cmain


cdef class Decoder(object):
    """High-level interface to the Moses decoder, wrapping around moses-cmd TranslationTask"""

    def __init__(self, *args):
        PRECISION = 3
        cmain.fix(cmain.cout, 3)
        cmain.fix(cmain.cerr, 3)

        cdef int argc = len(args)
        cdef char* argv[argc+1]

        for i, arg in enumerate(args):
            argv[i] = arg.c_str()

        cdef Parameter params cmain.Parameter

        if not params.LoadParam(argv,argv):
            raise Exception("Unable to parse Moses parameters")










