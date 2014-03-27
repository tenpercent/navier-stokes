#include "print.h"
#include <stdio.h>
#include <string.h>

#include "itersolv.h"
#include "precond.h"


void print_iterative_algorithm_info (Iterative_Method_parameters * parameters) {
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

  printf ("Using %s iterative method with %s preconditioner\n", method_to_string, preconditioner_to_string);

  return;
}
