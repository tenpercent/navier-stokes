#include <time.h>

#include "initialize.h"
#include "calculation.h"
#include "construction.h"
#include "export.h"
#include "print.h"

static const unsigned space_maxiter = 4;
static const unsigned time_maxiter  = 4;

void process_iteration(Iterative_Method_parameters * method_parameters,
                       Gas_parameters * gas_parameters,
                       unsigned time_iter,
                       unsigned space_iter,
                       unsigned global_iter,
                       double ** results)
{
  Grid grid;
  Node_status * node_statuses;
  double * V;
  double * G;
  double * space_coordinates;
  clock_t start_time;
  clock_t finish_time;

  /* Initialize grid and mesh */
  grid_Initialize (&grid, gas_parameters, space_iter, time_iter);
  scheme_elements_Construct (&G, &V, grid.X_nodes);
  mesh_elements_Construct (&node_statuses, &space_coordinates, grid.X_nodes);
  mesh_Initialize (node_statuses, space_coordinates, &grid);

  /* Run the solver and calculate time */
  start_time = clock();
  find_approximate_solution (G, V, node_statuses, space_coordinates, gas_parameters,
                             &grid, method_parameters);
  finish_time = clock();

  write_current_results (results, global_iter, start_time, finish_time,
                         gas_parameters, method_parameters,
                         &grid, space_coordinates, V, G);

  /* Print results to file */
  print_results_at_current_iteration (results, global_iter);
  global_iter++;

  /* Clean up */
  scheme_elements_Destruct (G, V);
  mesh_elements_Destruct (node_statuses, space_coordinates);
}

int main (int argc, char * argv[])
{
  Iterative_Method_parameters method_parameters;
  Gas_parameters gas_parameters;

  double * results[RESULTS_SIZE];
  unsigned global_iter = 0,
           time_iter   = 0,
           space_iter  = 0;

  /* Initialize global structures */
  initialize_iterative_algorithm_parameters (&method_parameters, argc, argv);
  print_iterative_algorithm_info (&method_parameters);
  gas_parameters_Initialize (&gas_parameters);
  results_Construct (results, time_maxiter * space_maxiter);

  /* The main loop */
  for (time_iter = 0; time_iter < time_maxiter; time_iter++) {
    for (space_iter = 0; space_iter < space_maxiter; space_iter++) {
      process_iteration(&method_parameters, &gas_parameters,
                        time_iter, space_iter, global_iter, results);
      global_iter++;
    }
  }

  /* Export results and cleanup */
  export_results (results, time_maxiter * space_maxiter);
  results_Destruct (results);

  return 0;
}
