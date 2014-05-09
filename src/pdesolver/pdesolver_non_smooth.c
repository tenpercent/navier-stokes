#include <time.h>

#include "initialize.h"
#include "calculation.h"
#include "construction.h"
#include "start_conditions.h"
#include "export.h"
#include "print.h"

static const unsigned time_maxiter = 10;
static const unsigned MAX_FILENAME_SIZE = 64;

void process_iteration(Iterative_Method_parameters const * method_parameters,
                       Gas_parameters const * gas_parameters,
                       Grid const * grid,
                       Node_status const * node_statuses,
                       double const * space_coordinates,
                       double * G,
                       double * V,
                       unsigned time_iter) {

  char filename[MAX_FILENAME_SIZE];

  find_approximate_solution (G, V, node_statuses, space_coordinates, gas_parameters,
                             grid, method_parameters);

  generate_table_filename ("G", time_iter, filename);
  write_value_table (G, space_coordinates, grid->X_nodes, filename);

  generate_table_filename ("V", time_iter, filename);
  write_value_table (V, space_coordinates, grid->X_nodes, filename);

  return;
}

int main (int argc, char * argv[])
{
  Iterative_Method_parameters method_parameters;
  Gas_parameters gas_parameters;

  unsigned time_iter   = 0;

  /* Initialize global structures */
  initialize_iterative_algorithm_parameters (&method_parameters, argc, argv);
  print_iterative_algorithm_info (&method_parameters);
  gas_parameters_Initialize (&gas_parameters);

  Grid grid;
  Node_status * node_statuses;
  double * V;
  double * G;
  double * space_coordinates;

  grid_Initialize (&grid, &gas_parameters, 3, 3);

  value_arrays_Construct (&G, &V, grid.X_nodes);
  mesh_elements_Construct (&node_statuses, &space_coordinates, grid.X_nodes);

  mesh_Initialize (node_statuses, space_coordinates, &grid);
  
  fill_mesh_at_initial_time (G, V, g_start, u_start, space_coordinates, grid.X_nodes);

  /* The main loop */
  for (time_iter = 0; time_iter < time_maxiter; time_iter++) {
    process_iteration(&method_parameters,
                      &gas_parameters,
                      &grid,
                      node_statuses,
                      space_coordinates,
                      G,
                      V,
                      time_iter);
  }

  value_arrays_Destruct (G, V);
  mesh_elements_Destruct (node_statuses, space_coordinates);

  return 0;
}
