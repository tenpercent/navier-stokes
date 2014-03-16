#include <stdio.h>
#include <string.h>
#include "export.h"

#define MAX_BUFFER_SIZE 10000

#ifdef _WIN32
#define ENDLINE \n\r
#else
#define ENDLINE \n
#endif

void export_results (double ** results, unsigned total_experiments) {
  FILE * csv_data; 
  unsigned experiments_step = 0;
  // dangerous shit begins
  char * result_to_string = 
    ",Time elapsed,V residual in C,V residual in L2,G residual in C,G residual in L2 ENDLINE";
  char * current_experiment_to_string = "";
  char * combined = "";
  for (experiments_step = 0; experiments_step < total_experiments; ++experiments_step) {
    snprintf (
      current_experiment_to_string,
      MAX_BUFFER_SIZE,
      ",%.3lf s,%.6lf,%.6lf,%.6lf,%.6lf ENDLINE", 
        results[0][experiments_step], 
        results[1][experiments_step], 
        results[2][experiments_step], 
        results[3][experiments_step], 
        results[4][experiments_step]
    );
    // |^ looks ugly, should fix later
    snprintf (combined, MAX_BUFFER_SIZE, "%s%s", result_to_string, current_experiment_to_string);
    result_to_string = combined;
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
