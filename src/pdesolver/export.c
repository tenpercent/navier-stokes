#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "export.h"
#include "functions.h"
#include "norm.h"

#define MAX_BUFFER_SIZE (1U << 17)
#define MAX_FILENAME_SIZE (1U << 6)

void write_current_results (
    double ** results,
    unsigned global_iteration, 
    clock_t start_time,
    clock_t finish_time,
    Gas_parameters const * gas_parameters,
    Iterative_Method_parameters const * iterative_method_parameters,
    Grid const * grid,
    double const * space_coordinates,
    double const * V,
    double const * G) {

  double * time_elapsed_at_iteration        = results[0],
         * residual_norm_V_C                = results[1],
         * residual_norm_V_L2               = results[2],
         * residual_norm_G_C                = results[3],
         * residual_norm_G_L2               = results[4],
         * time_step_at_iteration           = results[5],
         * space_step_at_iteration          = results[6],
         * relaxation_constant_at_iteration = results[7];

  time_elapsed_at_iteration[global_iteration] = (finish_time - start_time) / (double) CLOCKS_PER_SEC;
  residual_norm_V_C[global_iteration]         = residual_norm_C (V, 
                                                                 grid->X_nodes, 
                                                                 space_coordinates, 
                                                                 gas_parameters->time_upper_boundary, 
                                                                 u_exact);
  residual_norm_V_L2[global_iteration]        = residual_norm_L2 (V,
                                                                  grid->X_nodes, 
                                                                  space_coordinates, 
                                                                  gas_parameters->time_upper_boundary, 
                                                                  u_exact);
  residual_norm_G_C[global_iteration]         = residual_norm_C (G, 
                                                                 grid->X_nodes, 
                                                                 space_coordinates, 
                                                                 gas_parameters->time_upper_boundary, 
                                                                 g_exact);
  residual_norm_G_L2[global_iteration]        = residual_norm_L2 (G, 
                                                                  grid->X_nodes, 
                                                                  space_coordinates, 
                                                                  gas_parameters->time_upper_boundary, 
                                                                  g_exact);
  space_step_at_iteration[global_iteration] = grid->X_step;
  time_step_at_iteration[global_iteration] = grid->T_step;
  relaxation_constant_at_iteration[global_iteration] = iterative_method_parameters->relaxation_constant;
  return;
}

void export_results (double *const * results, unsigned total_experiments) {
  FILE * csv_data; 
  unsigned experiments_step = 0;
  // should be careful with c-strings
  char result_to_string[MAX_BUFFER_SIZE];
  char current_experiment_to_string[MAX_BUFFER_SIZE];

  char filename[MAX_FILENAME_SIZE];
  char time_representation[MAX_FILENAME_SIZE];

  struct tm * time_info;
  time_t timer = time(NULL);
  time_info = localtime (&timer);
  strftime (time_representation, MAX_FILENAME_SIZE, "%G_%b_%d_%H-%M-%S", time_info);

  strncpy (result_to_string, 
    ",Time elapsed,omega,tau,h,V residual in C,V residual in L2,G residual in C,G residual in L2\n",
    MAX_BUFFER_SIZE);
  for (experiments_step = 0; experiments_step < total_experiments; ++experiments_step) {
    snprintf (
      current_experiment_to_string,
      MAX_BUFFER_SIZE,
      "%d,\
      %.3lf s,\
      %.3lf,\
      %le,\
      %le,\
      %le,\
      %le,\
      %le,\
      %le\n",
        experiments_step + 1,
        results[0][experiments_step],
        results[7][experiments_step],
        results[5][experiments_step],
        results[6][experiments_step],
        results[1][experiments_step],
        results[2][experiments_step],
        results[3][experiments_step],
        results[4][experiments_step]
    );
    // |^ looks ugly, should fix later
    strncat (result_to_string, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result_to_string) - 1);
  }

  // might be cross-platform issues
  mkdir ("results", S_IRWXU);
  snprintf (filename, MAX_FILENAME_SIZE, "results/%s%s", time_representation, " !! results.csv");

  csv_data = fopen (filename, "w");

  if (!csv_data) {
    printf ("Could not open file to write results\n");
    return;
  }

  fputs (result_to_string, csv_data);
  fclose (csv_data);
  return;
}

#undef MAX_BUFFER_SIZE
