#include <time.h>

#include "initialize.h"
#include "calculation.h"
#include "construction.h"
#include "export.h"
#include "print.h"

int main (int argc, char ** argv) {
// local variables declaration
// let's begin with user defined constants
  unsigned const RELAXATION_CONSTANT_STEPS = 50; // something that is greater than 1
  double const RELAXATION_CONSTANT_MIN = .1,
               RELAXATION_CONSTANT_MAX = 4,
               RELAXATION_CONSTANT_INCREMENT = 
                (RELAXATION_CONSTANT_STEPS > 1) ?
                  (RELAXATION_CONSTANT_MAX - RELAXATION_CONSTANT_MIN) / (RELAXATION_CONSTANT_STEPS - 1) : 0;
  unsigned long const max_iteration_space  = 3,
                      max_iteration_time   = 3,
                      max_global_iteration = (max_iteration_space) * (max_iteration_time) * RELAXATION_CONSTANT_STEPS;
// user defined problem parameters and grid parameters
// iterative method parameters may also fall in this cathegory                      
  Gas_parameters gas_parameters;
  Grid grid;
  Iterative_Method_parameters iterative_method_parameters;
// mesh parameters
  double * space_coordinates = NULL;
  Node_status * node_statuses = NULL;
// function values in mesh nodes
  double * V = NULL;
  double * G = NULL;
// results to be exported
  double * results[RESULTS_SIZE];
// tick-tock
  clock_t start_time,
          finish_time;
// counters
  unsigned global_iteration              = 0,
           iteration_time                = 0, 
           iteration_space               = 0, 
           iteration_relaxation_constant = 0;
  /**********************************/
  /*** program evaluation starts! ***/
  /**********************************/
  initialize_iterative_algorithm_parameters (&iterative_method_parameters, argc, argv);
  print_iterative_algorithm_info (&iterative_method_parameters);

  gas_parameters_Initialize (&gas_parameters);
  results_Construct (results, max_global_iteration);  

  for (iteration_time = 0; iteration_time < max_iteration_time; ++iteration_time) {
    for (iteration_space = 0; iteration_space < max_iteration_space; ++iteration_space) {

      grid_Initialize (&grid, &gas_parameters, iteration_space, iteration_time);
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
                                   &gas_parameters, 
                                   &grid, 
                                   &iterative_method_parameters);
        finish_time = clock();

        write_current_results (results,
                               global_iteration, 
                               start_time,
                               finish_time,
                               &gas_parameters,
                               &iterative_method_parameters,
                               &grid,
                               space_coordinates,
                               V,
                               G);
        print_results_at_current_iteration (results, global_iteration);
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
