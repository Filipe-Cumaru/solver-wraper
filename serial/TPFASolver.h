#ifndef TPFASOLVER_H
#define TPFASOLVER_H

/* C++ STL includes */
#include <iostream>	/* std::cout, std::cin */
#include <numeric>	/* std::accumulate */
#include <cstdlib>	/* calloc, free */
#include <cstdio>	/* printf */
#include <cmath>	/* sqrt, pow */
#include <ctime>
#include <string>
#include <stdexcept>

/* MOAB includes */
#include "moab/Core.hpp"
#include "moab/MeshTopoUtil.hpp"

/* Trilinos includes */
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "AztecOO.h"

#define ALL_PROCS -1
#define ALL_DIM -1
#define GHOST_DIM 3
#define BRIDGE_DIM 2

using namespace std;
using namespace moab;

// Enumeration created to make the access to tags more readable.
enum TagsID {global_id, permeability, centroid, dirichlet, neumann};

class TPFASolver {
private:
    Interface *mb;
    MeshTopoUtil *topo_util;
    string perm_tag_name;
    string centroid_tag_name;
    string dirichlet_tag_name;
    string neumann_tag_name;
public:
    TPFASolver ();
    TPFASolver (Interface *moab_interface);
    void run ();
    void load_file (string fname);
    void write_file (string fname);
private:
    double calculate_centroid_dist (std::vector<double> c1, std::vector<double> c2);
    double calculate_equivalent_perm (std::vector<double> k1, std::vector<double> k2, double u[3]);
    double* calculate_unit_vector (std::vector<double> c1, std::vector<double> c2);
    void setup_tags (Tag tag_handles[5]);
    void assembly_matrix (Epetra_CrsMatrix& A, Epetra_Vector& b, Range volumes, Tag* tag_handles);
    void set_pressure_tags (Epetra_Vector& X, Range& volumes);
};

#endif
