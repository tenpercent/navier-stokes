#include "print.h"
#include <stdio.h>
#include <string.h>

#include "itersolv.h"
#include "precond.h"


void print_iterative_algorithm_info (Iterative_Method_parameters * iterative_method_parameters) {
/*
  unsigned const BUFFER_SIZE = 128;
  char method_to_string[BUFFER_SIZE];
  char preconditioner_to_string[BUFFER_SIZE];

  if (iterative_method_parameters->method == CGNIter) {
    strncpy (method_to_string, "CGN", BUFFER_SIZE);
  } else if (iterative_method_parameters->method == CGSIter) {
    strncpy (method_to_string, "CGS", BUFFER_SIZE);
  } else if (iterative_method_parameters->method == BiCGSTABIter) {
    strncpy (method_to_string, "BiCGStab", BUFFER_SIZE);
  } else if (iterative_method_parameters->method == QMRIter) {
    strncpy (method_to_string, "QMR", BUFFER_SIZE);
  } else if (iterative_method_parameters->method == GMRESIter) {
    strncpy (method_to_string, "GMRES", BUFFER_SIZE);
  } else {
    // nothing to do here
    strncpy (method_to_string, "nothing to do here", BUFFER_SIZE);
  }

  if (iterative_method_parameters->preconditioner_type == JacobiPrecond) {
    strncpy (preconditioner_to_string, "Jacobi", BUFFER_SIZE);
  } else if (iterative_method_parameters->preconditioner_type == SSORPrecond) {
    strncpy (preconditioner_to_string, "SSOR", BUFFER_SIZE);
  } else {
    // nothing to do here
    strncpy (preconditioner_to_string, "nothing to do here", BUFFER_SIZE);
  }

  printf ("Using %s iterative method with %s preconditioner\n", method_to_string, preconditioner_to_string);
*/

  return;
}
