#include <stdio.h>

#include "initialize.h"
#include "calculation.h"
#include "construction.h"
#include "start_conditions.h"
#include "export.h"
#include "print.h"

#define MAX_FILENAME_SIZE 64

void time_iteration(Iterative_Method_parameters const * method_parameters,
                       Gas_parameters const * gas_parameters,
                       Grid const * grid,
                       Node_status const * node_statuses,
                       double const * space_coordinates,
                       double * G,
                       double * V,
                       unsigned time_iter,
                       unsigned global_iteration) {

  char filename[MAX_FILENAME_SIZE];
  double time_elapsed;
  const clock_t start_time = clock();

  find_approximate_solution (G, V, node_statuses, space_coordinates, gas_parameters,
                             grid, method_parameters);

  time_elapsed = (clock() - start_time) / (double) CLOCKS_PER_SEC;

  generate_table_filename ("G", time_iter, global_iteration, filename);
  write_value_table (G, space_coordinates, grid->X_nodes, filename);

  generate_table_filename ("V", time_iter, global_iteration, filename);
  write_value_table (V, space_coordinates, grid->X_nodes, filename);

  /* print_results_at_current_iteration (time_elapsed, global_iteration);
   */

  return;
}

void test_iteration (double mu, double p_rho, double eta, unsigned global_iteration, int argc, char ** argv) {

  Iterative_Method_parameters method_parameters;
  Gas_parameters gas_parameters;

  Grid grid;
  Node_status * node_statuses;
  double * V;
  double * G;
  double * space_coordinates;

  unsigned time_iter   = 0;
  const unsigned time_maxiter = 10;
  /* Initialize global structures */
  initialize_iterative_algorithm_parameters (&method_parameters, argc, argv);
  print_iterative_algorithm_info (&method_parameters);
  gas_parameters_Initialize (&gas_parameters);

  gas_parameters.viscosity = mu;
  gas_parameters.artificial_viscosity = eta;
  gas_parameters.p_ro = p_rho;

  grid_Initialize (&grid, &gas_parameters, 3, 3);

  value_arrays_Construct (&G, &V, grid.X_nodes);
  mesh_elements_Construct (&node_statuses, &space_coordinates, grid.X_nodes);

  mesh_Initialize (node_statuses, space_coordinates, &grid);
  
  fill_mesh_at_initial_time (G, V, g_start, u_start, space_coordinates, grid.X_nodes);

  /* The main loop */
  for (time_iter = 0; time_iter < time_maxiter; ++time_iter) {
    time_iteration(&method_parameters,
                      &gas_parameters,
                      &grid,
                      node_statuses,
                      space_coordinates,
                      G,
                      V,
                      time_iter,
                      global_iteration);
  }

  print_iteration_info (mu, p_rho, eta, global_iteration);

  print_nonfancy_results_at_current_iteration(global_iteration);

  value_arrays_Destruct (G, V);
  mesh_elements_Destruct (node_statuses, space_coordinates);
  return;
}

int main (int argc, char * argv[])
{
  double mu_values[] = {.1, .01, .001};
  const unsigned MU_VALUES_SIZE = 3;
  double eta_values[] = {.1, 1., 10.};
  const unsigned ETA_VALUES_SIZE = 3;
  double p_rho_values[] = {1., 10., 100.};
  const unsigned P_RHO_VALUES_SIZE = 3;

  unsigned i = 0, 
           j = 0, 
           k = 0, 
           global_iteration = 0;

  for (i = 0; i < MU_VALUES_SIZE; ++i) {
    for (j = 0; j < P_RHO_VALUES_SIZE; ++j) {
      for (k = 0; k < ETA_VALUES_SIZE; ++k) {
        test_iteration (mu_values[i], p_rho_values[j], eta_values[k], global_iteration, argc, argv);
        global_iteration += 1;
      }
    }
  }

  return 0;
}

#undef MAX_FILENAME_SIZE
