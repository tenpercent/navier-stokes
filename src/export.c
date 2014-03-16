#include <stdio.h>
#include <string.h>
#include "export.h"

#define MAX_BUFFER_SIZE 2048

void export_results (double ** results, unsigned total_experiments) {
  FILE * csv_data; 
  unsigned experiments_step = 0;
  // dangerous shit begins
  char * result_to_string = 
    ",Time elapsed,V residual in C,V residual in L2,G residual in C,G residual in L2\n\r";
  char * current_experiment_to_string = "";
  char * temp = "";
  for (experiments_step = 0; experiments_step < total_experiments; ++experiments_step) {
    snprintf (
      current_experiment_to_string,
      MAX_BUFFER_SIZE,
      ",%.3lf s,%.6lf,%.6lf,%.6lf,%.6lf\n\r", 
        results[0][experiments_step], 
        results[1][experiments_step], 
        results[2][experiments_step], 
        results[3][experiments_step], 
        results[4][experiments_step]
    );
    snprintf (temp, MAX_BUFFER_SIZE, "%s%s", result_to_string, current_experiment_to_string);
    result_to_string = temp;
  }
  csv_data = fopen ("results.csv", "w");
  fputs (result_to_string, csv_data);
  fclose (csv_data);
  return;
}

#undef MAX_BUFFER_SIZE
