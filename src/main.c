#include <stdio.h>
#include <time.h>

#include "initialize.h"
#include "functions.h"
#include "norm.h"
#include "calculation.h"
#include "construction.h"
#include "export.h"
#include "print.h"

const unsigned RELAXATION_CONSTANT_STEPS = 50;
const double RELAXATION_CONSTANT_MIN = 0.;
const double RELAXATION_CONSTANT_MAX = 3.;
const double RELAXATION_CONSTANT_INCREMENT = 
  (RELAXATION_CONSTANT_MAX - RELAXATION_CONSTANT_MIN) / 
  (RELAXATION_CONSTANT_STEPS - 1);

int main (int argc, char ** argv) {

  Gas_parameters parameters;
  Grid grid;

  double * residual_norm_V_C;
  double * residual_norm_V_L2;
  double * residual_norm_G_C;
  double * residual_norm_G_L2;
  double * time_elapsed_at_iteration;
  double * results[RESULTS_SIZE];
  double * time_step_at_iteration;
  double * space_step_at_iteration;
  double * relaxation_constant_at_iteration;

  clock_t start_time,
          finish_time;

  double * space_coordinates = NULL;
  Node_status * node_statuses = NULL;
  double * V = NULL;
  double * G = NULL;

  Iterative_Method_parameters iterative_method_parameters;

  initialize_iterative_algorithm_parameters (&iterative_method_parameters, argc, argv);

  print_iterative_algorithm_info (&iterative_method_parameters);

  unsigned global_iteration = 0;
  unsigned const max_iteration_space = 3;
  unsigned const max_iteration_time = 3;
  unsigned const max_global_iteration = (max_iteration_space) * (max_iteration_time) * RELAXATION_CONSTANT_STEPS;

  unsigned iteration_time, iteration_space, iteration_relaxation_constant;

  gas_parameters_Initialize (&parameters);
  results_Construct (results, max_global_iteration);

  time_elapsed_at_iteration = results[0];
  residual_norm_V_C = results[1];
  residual_norm_V_L2 = results[2];
  residual_norm_G_C = results[3];
  residual_norm_G_L2 = results[4];
  time_step_at_iteration = results[5];
  space_step_at_iteration = results[6];
  relaxation_constant_at_iteration = results[7];

  for (iteration_time = 0; iteration_time < max_iteration_time; ++iteration_time) {
    for (iteration_space = 0; iteration_space < max_iteration_space; ++iteration_space) {
      grid_Initialize (&grid, &parameters, iteration_space, iteration_time);
      scheme_elements_Construct (&G, &V, grid.X_nodes);
      mesh_elements_Construct (&node_statuses, &space_coordinates, grid.X_nodes);
      mesh_Initialize (node_statuses, space_coordinates, &grid);
      for (iteration_relaxation_constant = 0;
           iteration_relaxation_constant < RELAXATION_CONSTANT_STEPS;
           ++iteration_relaxation_constant) {
        iterative_method_parameters.relaxation_constant = 
          RELAXATION_CONSTANT_MIN + iteration_relaxation_constant * RELAXATION_CONSTANT_INCREMENT;
        start_time = clock();
        find_approximate_solution (G, 
                                   V, 
                                   node_statuses, 
                                   space_coordinates, 
                                   &parameters, 
                                   &grid, 
                                   &iterative_method_parameters);
        finish_time = clock();

        time_elapsed_at_iteration[global_iteration] = (finish_time - start_time) / (double) CLOCKS_PER_SEC;
        residual_norm_V_C[global_iteration]         = residual_norm_C (V, 
                                                                      grid.X_nodes, 
                                                                      space_coordinates, 
                                                                      parameters.time_upper_boundary, 
                                                                      u_exact);
        residual_norm_V_L2[global_iteration]        = residual_norm_L2 (V,
                                                                       grid.X_nodes, 
                                                                       space_coordinates, 
                                                                       parameters.time_upper_boundary, 
                                                                       u_exact);
        residual_norm_G_C[global_iteration]         = residual_norm_C (G, 
                                                                      grid.X_nodes, 
                                                                      space_coordinates, 
                                                                      parameters.time_upper_boundary, 
                                                                      g_exact);
        residual_norm_G_L2[global_iteration]        = residual_norm_L2 (G, 
                                                                       grid.X_nodes, 
                                                                       space_coordinates, 
                                                                       parameters.time_upper_boundary, 
                                                                       g_exact);
        space_step_at_iteration[global_iteration] = grid.X_step;
        time_step_at_iteration[global_iteration] = grid.T_step;
        relaxation_constant_at_iteration[global_iteration] = iterative_method_parameters.relaxation_constant;

  #ifdef ALTERNATIVE_OUTPUT
        printf ("\r[ \x1b[32;01mok\x1b[00m ] Iteration %u finished in %.1lf seconds.\n", global_iteration + 1,
                time_elapsed_at_iteration[global_iteration]);
  #else
        printf ("\rFinished iteration %u in %.1lf seconds.\n", global_iteration + 1,
                time_elapsed_at_iteration[global_iteration]);
  #endif /* ALTERNATIVE_OUTPUT */
        ++global_iteration;
      }
      
      scheme_elements_Destruct (G, V);
      mesh_elements_Destruct (node_statuses, space_coordinates);
    }
  }
  export_results (results, max_global_iteration);
  results_Destruct (results);
  return 0;
}
