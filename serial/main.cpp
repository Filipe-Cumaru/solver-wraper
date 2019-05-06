#include <iostream>
#include <string>

#include "serial_tpfa_solver.h"

using namespace std;

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
