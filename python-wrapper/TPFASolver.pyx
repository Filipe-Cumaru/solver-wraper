cimport numpy as np
import numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string as cxx_string

from . cimport TPFASolver

# TODO: Check how to add compatibility between PyMOAB structures
# and C++ MOAB ones.
cdef class TPFA(object):

    # TODO: Add initialization with PyMOAB core instance.
    def __cinit__(self, moab_inst=None):
        self.inst = new TPFASolver.TPFASolver()

    def __del__(self):
        del self.inst

    def run(self):
        pass

    def load_file(self, cxx_string fname):
        pass

    def write_file(self, cxx_string fname):
        pass

    def calculate_centroid_dist(self, c1, c2):
        pass

    def calculate_equivalent_perm(self, k1, k2, u):
        pass

    def calculate_unit_vector(self, c1, c2):
        pass

    def setup_tags(self, tag_handles):
        pass

    def assembly_matrix(self, A, b, tag_handles):
        pass

    def set_pressure_tags(self, X, volumes):
        pass
