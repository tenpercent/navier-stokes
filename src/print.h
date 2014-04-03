#include "iterative_method.h"

// print the information about used iterative method
// and precontitioner type
// to stdio
void print_iterative_algorithm_info (Iterative_Method_parameters * iterative_method_parameters);

// print the information about current iteration
// to stdio
// currently implemented: time elapsed only
// because shorter is more beautiful
void print_info_about_current_iteration (
    double *const *results, 
    unsigned global_iteration);