#include <stdlib.h>
#include <time.h>
#include "initialize.h"
#include "functions.h"

int main () {
  Gas_parameters parameters;
  Grid grid;

  double * residual_norm_V_C;
  double * residual_norm_V_L2;
  double * residual_norm_G_C;
  double * residual_norm_G_L2;

  double * time_elapsed_at_iteration;
  clock_t start_time;
  clock_t finish_time;

  // what about double*[5] results ?

  double * space_coordinates;
  Node_status * node_statuses;
  double * V;
  double * G;

  unsigned global_iteration = 0;
  unsigned max_iteration_space = 3;
  unsigned max_iteration_time = 3;
  unsigned max_global_iteration = (max_iteration_space + 1) * (max_iteration_time + 1);

  gas_parameters_Initialize (&parameters);
  residual_norm_V_C = (double *) malloc (max_global_iteration * sizeof (double));
  residual_norm_V_L2 = (double *) malloc (max_global_iteration * sizeof (double));
  residual_norm_G_C = (double *) malloc (max_global_iteration * sizeof (double));
  residual_norm_G_L2 = (double *) malloc (max_global_iteration * sizeof (double));
  time_elapsed_at_iteration = (double *) malloc (max_global_iteration * sizeof (double));

  for (unsigned iteration_time = 0; iteration_time < max_iteration_time; ++iteration_time) {
    for (unsigned iteration_space = 0; iteration_space < max_iteration_space; ++iteration_space) {
      grid_Initialize (&grid, &parameters, iteration_space, iteration_time);

      V = (double *) malloc (grid.X_nodes * sizeof (double));
      G = (double *) malloc (grid.X_nodes * sizeof (double));
      space_coordinates = (double *) malloc (grid.X_nodes * sizeof (double));
      node_statuses = (Node_status *) malloc (grid.X_nodes * sizeof (Node_status));

      mesh_Initialize (node_statuses, space_coordinates, &grid);

      start_time = clock();

      next_time_layer_Calculate (G, V, node_statuses, space_coordinates, &grid);

      finish_time = clock();

      time_elapsed_at_iteration[global_iteration] = (finish_time - start_time) / (double) CLOCKS_PER_SEC;
      residual_norm_V_C[global_iteration]         = residual_norm_C (V, 
                                                                    grid.X_nodes, 
                                                                    space_coordinates, 
                                                                    parameters.time_upper_boundary, 
                                                                    u_bound);
      residual_norm_V_L2[global_iteration]        = residual_norm_L2 (V,
                                                                     grid.X_nodes, 
                                                                     space_coordinates, 
                                                                     parameters.time_upper_boundary, 
                                                                     u_bound);
      residual_norm_G_C[global_iteration]         = residual_norm_C (G, 
                                                                    grid.X_nodes, 
                                                                    space_coordinates, 
                                                                    parameters.time_upper_boundary, 
                                                                    logp_bound);
      residual_norm_G_L2[global_iteration]        = residual_norm_L2 (G, 
                                                                     grid.X_nodes, 
                                                                     space_coordinates, 
                                                                     parameters.time_upper_boundary, 
                                                                     logp_bound);

      free (V);
      free (G);
      free (space_coordinates);
      free (node_statuses);

      ++global_iteration;
    }
  }

  free (residual_norm_V_C);
  free (residual_norm_G_C);
  free (residual_norm_V_L2);
  free (residual_norm_G_L2);
  free (time_elapsed_at_iteration);

  return 0;
}
// too much code
// I want to create something like constructors and destructors for all of this
