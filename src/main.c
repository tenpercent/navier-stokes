#include <time.h>

#include <xmmintrin.h>

#include "initialize.h"
#include "functions.h"
#include "norm.h"
#include "calculation.h"
#include "construction.h"
#include "export.h"

int main (int argc, char ** argv) {

  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

  Gas_parameters parameters;
  Grid grid;

  double * residual_norm_V_C;
  double * residual_norm_V_L2;
  double * residual_norm_G_C;
  double * residual_norm_G_L2;
  double * time_elapsed_at_iteration;
  double * results[RESULTS_SIZE];

  clock_t start_time;
  clock_t finish_time;

  double * space_coordinates = NULL;
  Node_status * node_statuses = NULL;
  double * V = NULL;
  double * G = NULL;

  unsigned global_iteration = 0;
  unsigned max_iteration_space = 3;
  unsigned max_iteration_time = 3;
  unsigned max_global_iteration = (max_iteration_space) * (max_iteration_time);

  unsigned iteration_time, iteration_space;

  gas_parameters_Initialize (&parameters);
  results_Construct (results, max_global_iteration);

  time_elapsed_at_iteration = results[0];
  residual_norm_V_C = results[1];
  residual_norm_V_L2 = results[2];
  residual_norm_G_C = results[3];
  residual_norm_G_L2 = results[4];

  for (iteration_time = 0; iteration_time < max_iteration_time; ++iteration_time) {
    for (iteration_space = 0; iteration_space < max_iteration_space; ++iteration_space) {
      grid_Initialize (&grid, &parameters, iteration_space, iteration_time);
      scheme_elements_Construct (&G, &V, grid.X_nodes);
      mesh_elements_Construct (&node_statuses, &space_coordinates, grid.X_nodes);
      mesh_Initialize (node_statuses, space_coordinates, &grid);

      start_time = clock();
      next_TimeLayer_Calculate (G, V, node_statuses, space_coordinates, &parameters, &grid);
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
      printf ("\rFinished iteration %u in %.1lf seconds.\n", global_iteration,
              time_elapsed_at_iteration[global_iteration]);
      scheme_elements_Destruct (G, V);
      mesh_elements_Destruct (node_statuses, space_coordinates);
      ++global_iteration;
    }
  }
  export_results (results, max_global_iteration);
  results_Destruct (results);
  return 0;
}
