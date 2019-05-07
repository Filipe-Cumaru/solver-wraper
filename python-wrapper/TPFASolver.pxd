from libcpp.vector cimport vector
from libcpp.string cimport string as cxx_string

cdef extern from "TPFASolver.cpp":
    pass

cdef extern from "TPFASolver.h":
    cdef cppclass TPFASolver:
        pass
