#ifndef PRINT
#define PRINT

#include "iterative_method.h"
#include "structures.h"

#if defined (_WIN32) || defined (WIN32)
  #define FANCY(color, text) text
#else
  #define FANCY(color, text) "\x1b[3" #color ";01m" text "\x1b[00m"
#endif
#define FANCY_INFO "[" FANCY(5, "info") "] "
#define FANCY_OK "[ " FANCY(2, "ok") " ] "

// print the information about used iterative method
// and precontitioner type
// to stdio
void print_iterative_algorithm_info (Iterative_Method_parameters const * iterative_method_parameters);

// print the information about results at current iteration
// to stdio
// currently implemented: time elapsed only
// because shorter is more beautiful
void print_results_at_current_iteration (
    double const time_elapsed_at_iteration,
    unsigned global_iteration);

// print information about something during {time_step}
// of iterative algorithm
void print_info_about_current_iteration (unsigned time_step, Grid const * grid);

#endif
