#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "export.h"
#include "functions.h"
#include "norm.h"
#include "print.h"

#define MAX_BUFFER_SIZE (1U << 17)
#define SMALL_BUFFER_SIZE (1U << 5)
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

void export_results_to_string (double *const * results, unsigned total_experiments, char * result_to_string) {
  char current_experiment_to_string[MAX_BUFFER_SIZE] = "";
  const char delimiter = '&';
  unsigned experiments_step = 0;

  strncpy (result_to_string,
    ",Time elapsed,omega,tau,h,V residual in C,V residual in L2,G residual in C,G residual in L2\n",
    MAX_BUFFER_SIZE);
  for (experiments_step = 0; experiments_step < total_experiments; ++experiments_step) {
    snprintf (
      current_experiment_to_string,
      MAX_BUFFER_SIZE,
      "%d %c\
      %.3lf s %c\
      %.3lf %c\
      %le %c\
      %le %c\
      %le %c\
      %le %c\
      %le %c\
      %le\n",
        experiments_step + 1,
        delimiter,
        results[0][experiments_step],
        delimiter,
        results[7][experiments_step],
        delimiter,
        results[5][experiments_step],
        delimiter,
        results[6][experiments_step],
        delimiter,
        results[1][experiments_step],
        delimiter,
        results[2][experiments_step],
        delimiter,
        results[3][experiments_step],
        delimiter,
        results[4][experiments_step]
    );
    // |^ looks ugly, should fix later
    strncat (result_to_string, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result_to_string) - 1);
  }
  return;
}

void export_results (double *const * results, unsigned total_experiments) {
  char filename[MAX_FILENAME_SIZE] = "";
  char time_representation[MAX_FILENAME_SIZE] = "";
  char result_to_string[MAX_BUFFER_SIZE] = "";

  struct tm * time_info;
  time_t timer = time(NULL);

  time_info = localtime (&timer);
  strftime (time_representation, MAX_FILENAME_SIZE, "%G_%b_%d_%H-%M-%S", time_info);

  export_results_to_string (results, total_experiments, result_to_string);
  // might be cross-platform issues
  mkdir ("results", S_IRWXU);
  snprintf (filename, MAX_FILENAME_SIZE, "results/%s_results.csv", time_representation);

  rewrite_file (filename, result_to_string);

  printf(FANCY_INFO "Results saved to file `%s'.\n", filename);

  return;
}

void export_residual_table_to_string (double const * residuals,
                            double const * header,
                            double const * left_column,
                            unsigned const time_steps,
                            unsigned const space_steps,
                            char * result) {
  // result is output string

  char current_experiment_to_string[MAX_BUFFER_SIZE] = "";
  unsigned time_step, space_step;

  for (space_step = 0; space_step < space_steps; ++space_step) {
    snprintf (current_experiment_to_string,
            MAX_BUFFER_SIZE,
            "%.6lf",
            header[space_step]);
    strncat (result, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result) - 1);
  }

  for (space_step = 0; space_step < space_steps; ++space_step) {
    snprintf (current_experiment_to_string,
              MAX_BUFFER_SIZE,
              "%.6lf",
              left_column[space_step]);
    strncat (result, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result) - 1);

    for (time_step = 0; time_step < time_steps; ++time_step) {
      snprintf (current_experiment_to_string,
              MAX_BUFFER_SIZE,
              "%le",
              residuals[space_step * time_steps + time_step]);
      strncat (result, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result) - 1);
    }

    snprintf (current_experiment_to_string,
              MAX_BUFFER_SIZE,
              "%s",
              "\n");
    strncat (result, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result) - 1);
  }

  return;
}

void write_value_table (double const * values,
                        double const * space_coordinates,
                        unsigned space_nodes,
                        char * filename) {

  char values_to_string[MAX_BUFFER_SIZE] = "";

  char iteration_buffer[SMALL_BUFFER_SIZE] = "";

  unsigned space_step = 0;

  for (space_step = 0; space_step < space_nodes; ++space_step) {
    snprintf (iteration_buffer,
              SMALL_BUFFER_SIZE,
              "%.6lf %.6lf\n",
              space_coordinates[space_step],
              values[space_step]);
    strncat (values_to_string, iteration_buffer, MAX_BUFFER_SIZE - strlen(values_to_string) - 1);
  }

  rewrite_file (filename, values_to_string);
  return;
}

void rewrite_file (char const * filename,
                   char const * data) {

  FILE * file_handler = fopen (filename, "w");
  if (!file_handler) {
    printf ("Could not open file %s to write results\n", filename);
    return;
  }
  fputs (data, file_handler);
  fclose (file_handler);
  return;
}

void generate_table_filename (char const * sort,
                              unsigned time_iter,
                              unsigned global_iteration,
                              char * filename) {

  mkdir ("results", S_IRWXU);
  mkdir ("results/dat", S_IRWXU);

  char dirname[MAX_FILENAME_SIZE] = "";
  snprintf(dirname, MAX_FILENAME_SIZE, "results/dat/%u", global_iteration);
  mkdir (dirname, S_IRWXU);

  snprintf (filename,
            MAX_FILENAME_SIZE,
            "results/dat/%u/%s_at_%u_iteration.dat",
            global_iteration,
            sort,
            time_iter);
  return;
}

#undef MAX_BUFFER_SIZE
#undef SMALL_BUFFER_SIZE
#undef MAX_FILENAME_SIZE

