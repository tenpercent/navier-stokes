#ifndef FILLSYSTEM
#define FILLSYSTEM

#include "structures.h"
#include "sparse_matrix.h"

/* fill the system of linear equations
 * which purpose is to find the values of sought functions
 * in the mesh nodes of next time layer
 */
void fill_system (
    Sparse_matrix * lh_side,
    double * rh_side,
    Grid const * grid,
    Node_status const * node_status,
    Gas_parameters const * parameters,
    double const * space_coordinates,
    unsigned const time_step,
    double const * G,
    double const * V);

#endif
