#ifndef EXPORT 
#define EXPORT

#include <time.h>
#include "structures.h"
#include "iterative_method.h"

// write results at current iteration
// to array results[][]
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
    double const * G);
// export the results of calculation
// to csv file "./results.csv"
void export_results (double *const * results, unsigned total_experiments);

#endif
