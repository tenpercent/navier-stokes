#ifndef CALCULATION
#define CALCULATION

#include "structures.h"
#include "initialize.h"

// find the values of sought functions
// in mesh nodes with coordinates (space_coordinates)
// with the information of node positions stored in (node_status)
// with the information about parameters of system of PDEs stored in (gas_parameters)
// with the information about grid stored in (grid)
// using the iterative method and preconditioner specified in (iterative_method_parameters)
// and write it to arrays (V) and (G)
void find_approximate_solution (
  double * G, 
  double * V, 
  Node_status * node_status, 
  double * space_coordinate, 
  Gas_parameters * gas_parameters,
  Grid * grid,
  Iterative_Method_parameters * iterative_method_parameters);

// fill the system of linear equations
// which purpose is to find the values of sought functions
// in the mesh nodes of next time layer
void fill_system (
    Sparse_matrix * lh_side,
    double * rh_side,
    Grid const * grid,
    Node_status const * node_status,
    Gas_parameters * parameters,
    double const * space_coordinates,
    unsigned const time_step,
    double const * G,
    double const * V);

// fill the arrays (V[]) and (G[])
// with the values of sought functions
// in the nodes of the initial time layer.
// These values are obtained by applying
// (g(x,t)) and (v(x,t)) to the points of the initial time layer
void fill_mesh_at_initial_time (
  double * G,
  double * V,
  double (*g) (double, double),
  double (*v) (double, double),
  double * space_coordinates,
  unsigned space_nodes);

// fill the arrays (V[]) and (G[])
// with the values obtained
// after solving the system of linear equations
// from (solutions)
void fill_approximation (double * G, double * V, double * solutions, unsigned total_values);

#endif
