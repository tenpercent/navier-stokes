#include "structures.h"

void gas_parameters_Initialize (Gas_parameters * parameters);

void grid_Initialize (
  Grid * grid, 
  Gas_parameters * parameters,
  unsigned time_steps_magnifier, 
  unsigned space_steps_magnifier);

void mesh_Initialize (Node_status * node_status, double * space_coordinate, Grid * grid); // no need... actually, it will be executed once.

void initialize_iterative_algorithm_parameters (Preconditioner_type * type, Iterative_method * method, int argc, char ** argv);
