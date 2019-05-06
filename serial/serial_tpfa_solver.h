#define ALL_PROCS -1
#define ALL_DIM -1
#define GHOST_DIM 3
#define BRIDGE_DIM 2

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
    ErrorCode run ();
    ErrorCode load_file ();
    ErrorCode write_file ();
private:
    double calculate_centroid_dist (std::vector<double> c1, std::vector<double> c2);
    double calculate_equivalent_perm (std::vector<double> k1, std::vector<double> k2, double u[3]);
    double* calculate_unit_vector (std::vector<double> c1, std::vector<double> c2);
    ErrorCode setup_tags (Tag tag_handles[5]);
    ErrorCode assembly_matrix (Epetra_CrsMatrix& A, Epetra_Vector& b, Range volumes, Tag* tag_handles);
    ErrorCode set_pressure_tags (Epetra_Vector& X, Range& volumes);
}
