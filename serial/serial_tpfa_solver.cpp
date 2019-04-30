/* C++ STL includes */
#include <iostream>	/* std::cout, std::cin */
#include <numeric>	/* std::accumulate */
#include <cstdlib>	/* calloc, free */
#include <cstdio>	/* printf */
#include <cmath>	/* sqrt, pow */
#include <ctime>
#include <string>

/* MOAB includes */
#include "moab/Core.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/Skinner.hpp"

/* Trilinos includes */
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"
#include "AztecOO.h"
// #include "ml_include.h"
// #include "ml_epetra_preconditioner.h"

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
    TPFASolver () : mb(new Core()),
                    topo_util(new MeshTopoUtil(mb)),
                    perm_tag_name("PERMEABILITY"),
                    centroid_tag_name("CENTROID"),
                    dirichlet_tag_name("DIRICHLET_BC"),
                    neumann_tag_name("NEUMANN_BC") {}
    TPFASolver (Interface *moab_interface) : mb(moab_interface),
                                            topo_util(new MeshTopoUtil(mb)),
                                            perm_tag_name("PERMEABILITY"),
                                            centroid_tag_name("CENTROID"),
                                            dirichlet_tag_name("DIRICHLET_BC"),
                                            neumann_tag_name("NEUMANN_BC") {}
    ErrorCode run () {
        ErrorCode rval;
        Range volumes;

        // Get all volumes in the mesh.
        cout << "Getting volumes" << endl;
        rval = this->mb->get_entities_by_dimension(0, 3, volumes, false);
        MB_CHK_SET_ERR(rval, "get_entitites_by_dimension failed");

        cout << "Setting tags" << endl;
        Tag tag_handles[5];
        int* gids = (int*) calloc(volumes.size(), sizeof(int));
        if (gids == NULL) {
            printf("Error: Null pointer\n");
            exit(EXIT_FAILURE);
        }
        this->setup_tags(tag_handles);
        rval = this->mb->tag_get_data(tag_handles[global_id], volumes, (void*) gids);
        MB_CHK_SET_ERR(rval, "tag_get_data for gids failed");

        Epetra_SerialComm comm;
        Epetra_Map row_map (volumes.size(), volumes.size(), gids, 0, comm);
        Epetra_CrsMatrix A (Copy, row_map, 7);
        Epetra_Vector b (row_map);
        Epetra_Vector X (row_map);

        cout << "Assembling matrix" << endl;
        clock_t ts = clock();
        this->assembly_matrix(A, b, volumes, tag_handles);
        ts = clock() - ts;
        printf("Done. Time elapsed: %f\n", ((double) ts)/CLOCKS_PER_SEC);

        Epetra_LinearProblem linear_problem (&A, &X, &b);
        AztecOO solver (linear_problem);

        // TODO: Add ML compatibility.
        solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
        solver.Iterate(1000, 1e-14);

        this->set_pressure_tags(X, volumes);

        free(gids);
        return MB_SUCCESS;
    }

    ErrorCode load_file (string fname) {
        ErrorCode rval;
        rval = this->mb->load_file(fname.c_str()); MB_CHK_ERR(rval);
        return MB_SUCCESS;
    }

    ErrorCode write_file (string fname) {
        ErrorCode rval;
        EntityHandle volumes_meshset;
        Range volumes;
        this->mb->create_meshset(0, volumes_meshset);
        rval = this->mb->get_entities_by_dimension(0, 3, volumes, false);
        MB_CHK_ERR(rval);
        rval = this->mb->add_entities(volumes_meshset, volumes);
        MB_CHK_ERR(rval);
        rval = this->mb->write_file(fname.c_str(), 0, 0, &volumes_meshset, 1);
        MB_CHK_ERR(rval);
        return MB_SUCCESS;
    }

private:
    double calculate_centroid_dist (double c1[3], double c2[3]) {
    	/*
    		Calculate the distance between two points.

    		Parameters
    		----------
    		c1: double*
    			An array with the coordinates of the first point.
    		c2: double*
    			An array with the coordinates of the second point.

    		Returns
    		-------
    		Double value of the distance between c1 and c2.
    	*/

        return sqrt(pow(c1[0] - c2[0], 2) + pow(c1[1] - c2[1], 2) + pow(c1[2] - c2[2], 2));
    }

    double calculate_equivalent_perm (double k1[9], double k2[9], double u[3]) {
	    /*
	    	Calculate the equivalent permeability between k1 and k2. To obtain the
	    	equivalent permeability for each element, the permeability tensor is
	    	multiplied (using the inner product) by the unit vector u twice, i.e,
	    	K1_eq = <<K1, u>, u>, where <,> denotes the inner product.

	    	Parameters
	    	----------
	    	k1: double*
	    		An array representing the permeability tensor of the first element.
	    	k2: double*
		    	An array representing the permeability tensor of the second element.
		    u: double*
		    	An array representing the coordinates of the unit vector.

		    Returns
		    -------
		    Double value of the equivalent permeability of the two elements.
	    */

        // REVIEW: Check BLAS library for dot product routines and verify
        // overhead to wrap it with Cython.
        double k1_pre[3] = {k1[0]*u[0] + k1[3]*u[0] + k1[6]*u[0],
  		                    k1[1]*u[1] + k1[4]*u[1] + k1[7]*u[1],
  		                    k1[2]*u[2] + k1[5]*u[2] + k1[8]*u[2]};
        double k2_pre[3] = {k2[0]*u[0] + k2[3]*u[0] + k2[6]*u[0],
  		                    k2[1]*u[1] + k2[4]*u[1] + k2[7]*u[1],
  		                    k2[2]*u[2] + k2[5]*u[2] + k2[8]*u[2]};
        double k1_eq = k1_pre[0]*u[0] + k1_pre[1]*u[1] + k1_pre[2]*u[2];
        double k2_eq = k2_pre[0]*u[0] + k2_pre[1]*u[1] + k2_pre[2]*u[2];
        return 2*k1_eq*k2_eq/(k1_eq + k2_eq);
    }

    double* calculate_unit_vector (double c1[3], double c2[3]) {
    	/*
    		Calculate the unit vector given two points.

    		Parameters
    		----------
    		c1: double*
    			An array with the coordinates of the first point.
    		c2: double*
    			An array with the coordinates of the second point.

    		Returns
    		-------
    		Array with the vector coordinates.
    	*/

        double d = this->calculate_centroid_dist(c1, c2);
        static double u[3] = {(c2[0] - c1[0])/d, (c2[1] - c1[1])/d, (c2[2] - c1[2])/d};
        return u;
    }

    ErrorCode setup_tags (Tag tag_handles[5]) {
        ErrorCode rval;

        rval = this->mb->tag_get_handle("GLOBAL_ID", tag_handles[global_id]); MB_CHK_ERR(rval);
        rval = this->mb->tag_get_handle(this->centroid_tag_name.c_str(), tag_handles[centroid]); MB_CHK_ERR(rval);
        rval = this->mb->tag_get_handle(this->perm_tag_name.c_str(), tag_handles[permeability]); MB_CHK_ERR(rval);
        rval = this->mb->tag_get_handle(this->dirichlet_tag_name.c_str(), tag_handles[dirichlet]); MB_CHK_ERR(rval);
        rval = this->mb->tag_get_handle(this->neumann_tag_name.c_str(), tag_handles[neumann]); MB_CHK_ERR(rval);

        return MB_SUCCESS;
    }

    ErrorCode assembly_matrix (Epetra_CrsMatrix& A, Epetra_Vector& b, Range volumes, Tag* tag_handles) {
        /*
    		Assembly the transmissibility matrix, a.k.a, the coeficient matrix A of the linear
    		system to be solved.

    		Parameters
    		----------
    		A: Epetra_CrsMatrix
    			Epetra matrix used to store the coeficients.
    		b: Epetra_Vector*
    			Epetra vector used to store the boundary values.
    		volumes: moab::Range
    			MOAB Range containing the entity handles for all volumes.
    		tag_handles: Tag*
    			Array of tag handles

    		Returns
    		-------
    		MOAB error code.
    	*/

        ErrorCode rval;
        Range adjacencies;
        std::vector<double> row_values;
        std::vector<int> row_indexes;
        double c1[3], c2[3], k1[9], k2[9], *u;
        // double equiv_perm = 0, centroid_dist = 0, pressure_bc = 0, flux_bc = 0, diag_coef = 0;
        double equiv_perm = 0, centroid_dist = 0, diag_coef = 0;
        int row_id = -1, n = 0;
        
        int num_vols = volumes.size();
        double *pressure_bc = (double*) calloc(volumes.size(), sizeof(double));
        double *flux_bc = (double*) calloc(volumes.size(), sizeof(double));
        double *centroids = (double*) calloc(volumes.size(), 3*sizeof(double));
        double *perms = (double*) calloc(volumes.size(), 9*sizeof(double));
        int* gids = (int*) calloc(volumes.size(), sizeof(int));
        rval = this->mb->tag_get_data(tag_handles[dirichlet], volumes, pressure_bc);
        rval = this->mb->tag_get_data(tag_handles[neumann], volumes, flux_bc);
        rval = this->mb->tag_get_data(tag_handles[centroid], volumes, centroids);
        rval = this->mb->tag_get_data(tag_handles[permeability], volumes, perms);
        rval = this->mb->tag_get_data(tag_handles[global_id], volumes, (void*) gids);

        for (int i = 0; i < num_vols; i++) {
            if (pressure_bc[i] != 0) {
                diag_coef = 1;
            }
            else {
                rval = this->topo_util->get_bridge_adjacencies(volumes[i], BRIDGE_DIM, 3, adjacencies); MB_CHK_ERR(rval);
                c1[0] = centroids[num_vols*i]; c1[1] = centroids[num_vols*i+1]; c1[2] = centroids[num_vols*i+2];
                k1[0] = perms[num_vols*i]; k1[1] = perms[num_vols*i+1]; k1[2] = perms[num_vols*i+2];
                k1[3] = perms[num_vols*i+3]; k1[4] = perms[num_vols*i+4]; k1[5] = perms[num_vols*i+5];
                k1[6] = perms[num_vols*i+6]; k1[7] = perms[num_vols*i+7]; k1[8] = perms[num_vols*i+8];
                for (Range::iterator itt = adjacencies.begin(); itt != adjacencies.end(); itt++) {
                    rval = this->mb->tag_get_data(tag_handles[centroid], &(*itt), 1, &c2); MB_CHK_ERR(rval);
                    rval = this->mb->tag_get_data(tag_handles[permeability], &(*itt), 1, &k2); MB_CHK_ERR(rval);
                    rval = this->mb->tag_get_data(tag_handles[global_id], &(*itt), 1, &row_id); MB_CHK_ERR(rval);
                    centroid_dist = this->calculate_centroid_dist(c1, c2);
                    u = this->calculate_unit_vector(c1, c2);
                    equiv_perm = this->calculate_equivalent_perm(k1, k2, u);
                    // TODO: Generalize to unstructured grids, i.e., calculate
                    // distance for each element.
                    row_values.push_back(-equiv_perm/centroid_dist);
                    row_indexes.push_back(row_id);
                }
                diag_coef = -accumulate(row_values.begin(), row_values.end(), 0.0);
            }
            row_values.push_back(diag_coef);
            row_indexes.push_back(gids[i]);
            A.InsertGlobalValues(row_id, row_values.size(), &row_values[0], &row_indexes[0]);
            if (flux_bc[i] != -1)
                b[n] = pressure_bc[i] + flux_bc[i];
            else
                b[n] = pressure_bc[i];
            n += 1;
            row_values.clear();
            row_indexes.clear();
            adjacencies.clear();
        }

        // for (Range::iterator it = volumes.begin(); it != volumes.end(); it++) {
            // rval = this->mb->tag_get_data(tag_handles[dirichlet], &(*it), 1, &pressure_bc); MB_CHK_ERR(rval);
            // rval = this->mb->tag_get_data(tag_handles[neumann], &(*it), 1, &flux_bc); MB_CHK_ERR(rval);
            // if (pressure_bc != 0) {
            //     diag_coef = 1;
            // }
            // else {
            //     rval = this->topo_util->get_bridge_adjacencies(*it, BRIDGE_DIM, 3, adjacencies); MB_CHK_ERR(rval);
            //     rval = this->mb->tag_get_data(tag_handles[centroid], &(*it), 1, &c1); MB_CHK_ERR(rval);
            //     rval = this->mb->tag_get_data(tag_handles[permeability], &(*it), 1, &k1); MB_CHK_ERR(rval);
            //     for (Range::iterator itt = adjacencies.begin(); itt != adjacencies.end(); itt++) {
            //         rval = this->mb->tag_get_data(tag_handles[centroid], &(*itt), 1, &c2); MB_CHK_ERR(rval);
            //         rval = this->mb->tag_get_data(tag_handles[permeability], &(*itt), 1, &k2); MB_CHK_ERR(rval);
            //         rval = this->mb->tag_get_data(tag_handles[global_id], &(*itt), 1, &row_id); MB_CHK_ERR(rval);
            //         centroid_dist = this->calculate_centroid_dist(c1, c2);
            //         u = this->calculate_unit_vector(c1, c2);
            //         equiv_perm = this->calculate_equivalent_perm(k1, k2, u);
            //         // TODO: Generalize to unstructured grids, i.e., calculate
            //         // distance for each element.
            //         row_values.push_back(-equiv_perm/centroid_dist);
            //         row_indexes.push_back(row_id);
            //     }
            //     diag_coef = -accumulate(row_values.begin(), row_values.end(), 0.0);
            // }
            // rval = this->mb->tag_get_data(tag_handles[global_id], &(*it), 1, &row_id); MB_CHK_ERR(rval);
            // row_values.push_back(diag_coef);
            // row_indexes.push_back(row_id);
            // A.InsertGlobalValues(row_id, row_values.size(), &row_values[0], &row_indexes[0]);
            // if (flux_bc != -1)
            //     b[n] = pressure_bc + flux_bc;
            // else
            //     b[n] = pressure_bc;
            // n += 1;
            // row_values.clear();
            // row_indexes.clear();
            // adjacencies.clear();
        // }
        A.FillComplete();
        return MB_SUCCESS;
    }

    ErrorCode set_pressure_tags (Epetra_Vector& X, Range& volumes) {
    	/*
    		Create and set values for the pressure on each volume.

    		Parameters
    		----------
    		X: Epetra_Vector
    			Epetra vector containing the pressure values for each volume.
    		volumes: moab::Range
    			MOAB Range containing the entity handles for all volumes.

			Returns
			-------
			MOAB error code.
    	*/

        Tag pressure_tag;
        ErrorCode rval;
        rval = this->mb->tag_get_handle("PRESSURE", 1, MB_TYPE_DOUBLE, pressure_tag, MB_TAG_DENSE | MB_TAG_CREAT);
        MB_CHK_SET_ERR(rval, "tag_get_handle for pressure tag failed");
        rval = this->mb->tag_set_data(pressure_tag, volumes, &X[0]);
        MB_CHK_SET_ERR(rval, "tag_set_data for pressure tag failed");
        return MB_SUCCESS;
    }
};

int main () {
    TPFASolver* solver = new TPFASolver();
    // string input_file = "/home/facsa/Documents/TPFA/mesh_files/tpfa_mesh.h5m";
    string input_file = "/home/facsa/Documents/TPFA/tpfa_mesh.h5m";
    string output_file = "solution_mesh.h5m";

    cout << "Loading file..." << endl;
    solver->load_file(input_file);
    cout << "Done." << endl;
    solver->run();
    cout << "Writing file..." << endl;
    solver->write_file(output_file);
    cout << "Done." << endl;

    return 0;
}
