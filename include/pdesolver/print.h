#ifndef PRINT 
#define PRINT

#include "iterative_method.h"
#include "structures.h"

// print the information about used iterative method
// and precontitioner type
// to stdio
void print_iterative_algorithm_info (Iterative_Method_parameters const * iterative_method_parameters);

// print the information about results at current iteration
// to stdio
// currently implemented: time elapsed only
// because shorter is more beautiful
void print_results_at_current_iteration (
    double *const *results, 
    unsigned global_iteration);

// print information about something during {time_step}
// of iterative algorithm
void print_info_about_current_iteration (unsigned time_step, Grid const * grid);

#endif
