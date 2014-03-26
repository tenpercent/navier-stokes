#include "structures.h"

// fill (parameters) with information
// hard code inside
void gas_parameters_Initialize (Gas_parameters * parameters);

// fill (grid) with information
// hard code inside
void grid_Initialize (
  Grid * grid, 
  Gas_parameters * parameters,
  unsigned time_steps_magnifier, 
  unsigned space_steps_magnifier);

// fill (node_status[]) and (space_coordinates[]) with information
// hard code inside
void mesh_Initialize (Node_status * node_statuses, double * space_coordinates, Grid * grid);

// fill (parameters) with information from command line
void initialize_iterative_algorithm_parameters (Iterative_Method_parameters * parameters, int argc, char ** argv);
