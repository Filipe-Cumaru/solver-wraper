from libcpp.vector cimport vector
from libcpp.string cimport string as cxx_string

cdef extern from "serial/TPFASolver.cpp":
    pass

cdef extern from "serial/TPFASolver.h":
    ctypedef enum TagsID:
        global_id
        permeability
        centroid
        dirichlet
        neumann

    cdef cppclass TPFASolver:
        TPFASolver ()
        TPFASolver (Interface* mb)
        void run ()
        void load_file (cxx_string fname)
        void write_file (cxx_string fname)
        double calculate_centroid_dist (vector[double] c1, vector[double] c2)
        double calculate_equivalent_perm (vector[double] k1, vector[double] k2, double u[3])
        double* calculate_unit_vector (vector[double] c1, vector[double] c2)
        void setup_tags (Tag tag_handles[5])
        void assembly_matrix (Epetra_CrsMatrix& A, Epetra_Vector& b, Range volumes, Tag* tag_handles)
        void set_pressure_tags (Epetra_Vector& X, Range& volumes)
