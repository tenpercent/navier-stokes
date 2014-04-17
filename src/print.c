#include "print.h"
#include <stdio.h>
#include <string.h>

#ifndef NO_LASPACK
#include "itersolv.h"
#include "precond.h"
#endif /* NO_LASPACK */

#if defined (_WIN32) || defined (WIN32)
#define FANCY(color, text) text
#else
#define FANCY(color, text) "\x1b[3" #color ";01m" text "\x1b[00m"
#endif

void print_iterative_algorithm_info (Iterative_Method_parameters const * parameters) {
  unsigned const BUFFER_SIZE = 128;
  char method_to_string[BUFFER_SIZE];
  char preconditioner_to_string[BUFFER_SIZE];

  if (parameters->implementation == Implementation_Native) {
    if (parameters->method == Iterative_method_BiCGSTAB) {
      sprintf (method_to_string, "BiCGSTAB");
    } else {
      /* nothing to do here */
      sprintf (method_to_string, "UNKNOWN");
    }

    if (parameters->preconditioner_type == Precond_Jacobi) {
      sprintf (preconditioner_to_string, "Jacobi");
    } else {
      /* nothing to do here */
      sprintf (preconditioner_to_string, "UNKNOWN");
    }
  }

#ifndef NO_LASPACK
  if (parameters->implementation == Implementation_Laspack) {
    if (parameters->method_laspack == CGNIter) {
      sprintf (method_to_string, "CGN (Laspack)");
    } else if (parameters->method_laspack == CGSIter) {
      sprintf (method_to_string, "CGS (Laspack)");
    } else if (parameters->method_laspack == BiCGSTABIter) {
      sprintf (method_to_string, "BiCGStab (Laspack)");
    } else if (parameters->method_laspack == QMRIter) {
      sprintf (method_to_string, "QMR (Laspack)");
    } else if (parameters->method_laspack == GMRESIter) {
      sprintf (method_to_string, "GMRES (Laspack)");
    } else {
      /* nothing to do here */
      sprintf (method_to_string, "UNKNOWN (Laspack)");
    }

    if (parameters->preconditioner_type_laspack == JacobiPrecond) {
      sprintf (preconditioner_to_string, "Jacobi");
    } else if (parameters->preconditioner_type_laspack == SSORPrecond) {
      sprintf (preconditioner_to_string, "SSOR");
    } else {
      /* nothing to do here */
      sprintf (preconditioner_to_string, "UNKNOWN");
    }
  }
#endif /* NO_LASPACK */

  printf ("["FANCY(5, "info")"] Method: %s.\n", method_to_string);
  printf ("["FANCY(5, "info")"] Preconditioner: %s.\n", preconditioner_to_string);
  printf ("["FANCY(5, "info")"] Relaxation constant: %lg.\n", parameters->relaxation_constant);

  return;
}

void print_results_at_current_iteration (
    double *const *results,
    unsigned global_iteration) {

  double const *const time_elapsed_at_iteration = results[0];
  printf ("\r[ "FANCY(2, "ok")" ] Iteration %u finished in %.1lf seconds.\n",
          global_iteration + 1,
          time_elapsed_at_iteration[global_iteration]);

  return;
}

void print_info_about_current_iteration (unsigned time_step, Grid const * grid) {
  printf ("\r[....] Time step is %u of %u.", time_step + 1, grid->T_nodes);
  fflush (stdout);
  return;
}
