#include <stdio.h>
#include <string.h>
#include "export.h"

#define MAX_BUFFER_SIZE (1 << 20)

void export_results (double ** results, unsigned total_experiments) {
  FILE * csv_data; 
  unsigned experiments_step = 0;
  // dangerous shit begins
  char result_to_string[MAX_BUFFER_SIZE];
  strncpy (result_to_string, 
    ",Time elapsed,omega,tau,h,V residual in C,V residual in L2,G residual in C,G residual in L2\n",
    MAX_BUFFER_SIZE);
  char current_experiment_to_string[MAX_BUFFER_SIZE];
  for (experiments_step = 0; experiments_step < total_experiments; ++experiments_step) {
    snprintf (
      current_experiment_to_string,
      MAX_BUFFER_SIZE,
      "%d,%.3lf s,%le,%le,%.3lf,%le,%le,%le,%le\n",
        experiments_step + 1,
        results[0][experiments_step],
        results[5][experiments_step],
        results[6][experiments_step],
        results[7][experiments_step],
        results[1][experiments_step],
        results[2][experiments_step],
        results[3][experiments_step],
        results[4][experiments_step]
    );
    // |^ looks ugly, should fix later
    strncat (result_to_string, current_experiment_to_string,
             MAX_BUFFER_SIZE - strlen(result_to_string) - 1);
  }

  csv_data = fopen ("results.csv", "w");

  if (!csv_data) {
    printf ("Could not open file to write results\n");
    return;
  }

  fputs (result_to_string, csv_data);
  fclose (csv_data);
  return;
}

#undef MAX_BUFFER_SIZE
