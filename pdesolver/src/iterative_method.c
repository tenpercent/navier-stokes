#include "iterative_method.h"
#include <stdio.h>
#include <string.h>

double scalar_product (
    double const * v1,
    double const * v2,
    unsigned size) {

  register double result = 0;
  register unsigned i;
  for (i = 0; i < size; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

void Precond_Jacobi (
    Sparse_matrix const * matrix,
    double const        * c,
    double              * y,
    double                omega) {

  register unsigned i;
  for (i = 0; i < matrix->size; ++i) {
    y[i] = omega * c[i] / matrix->elements[i];
  }

}

void Precond_Null (
    Sparse_matrix const * matrix,
    double const        * c,
    double              * y,
    double                omega) {

  (void) omega;

  memcpy (y, c, matrix->size * sizeof(double));
}

/* Preconditioned BiCGSTAB method, described in: */
/* http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method */

void Iterative_method_BiCGSTAB (
    Sparse_matrix const * matrix,
    double              * x,
    double const        * b,
    unsigned              max_iter,
    Precond_type          precond,
    double                omega_precond, /* default = 1 */
    double                accuracy,
    double              * buffer /* we need 8 * size */
  ) {

  register unsigned size   = matrix->size;
  register double alpha    = 1;
  register double rho      = 1;
  register double omega    = 1;
  double *r       = buffer;
  double *rcap    = buffer + size;
  double *v       = buffer + size * 2;
  double *p       = buffer + size * 3;
  double *y       = buffer + size * 4;
  double *z       = buffer + size * 5;
  double *s       = buffer + size * 6;
  double *t       = buffer + size * 7;

  register unsigned i, iter;
  register double beta, rhoold;

  register const double bnorm2 = scalar_product(b, b, size);
  register const double accuracy2 = accuracy * accuracy;

  Sparse_matrix_Apply_to_vector (matrix, x, r);
  memset (v, 0, size * sizeof(double));
  memset (p, 0, size * sizeof(double));
  for (i = 0; i < size; ++i) {
    r[i] = b[i] - r[i];
    rcap[i] = r[i];
  }

  for (iter = 1; iter < max_iter; ++iter) {
    rhoold = rho;
    rho = scalar_product (rcap, r, size);
    beta = (rho * alpha) / (rhoold * omega);
    for (i = 0; i < size; ++i) {
      p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
    (*precond)(matrix, p, y, omega_precond);
    Sparse_matrix_Apply_to_vector (matrix, y, v);
    alpha = rho / scalar_product (rcap, v, size);
    for (i = 0; i < size; ++i) {
      s[i] = r[i] - alpha * v[i];
    }
    (*precond)(matrix, s, z, omega_precond);
    Sparse_matrix_Apply_to_vector (matrix, z, t);
    omega = scalar_product (t, s, size) /
            scalar_product (t, t, size);
    for (i = 0; i < size; ++i) {
      x[i] += (alpha * y[i] + omega * z[i]);
    }
    if (scalar_product(r, r, size) < accuracy2 * bnorm2) {
      return;
    }
    for (i = 0; i < size; ++i) {
      r[i] = s[i] - omega * t[i];
    }
  }
  printf ("\nWarning: maximum number of iterations reached: %u.\n", max_iter);
}
