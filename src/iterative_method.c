#include "iterative_method.h"
#include <string.h>

double scalar_product (
    double * v1,
    double * v2,
    unsigned size) {

  double result = 0;
  unsigned i;
  for (i = 0; i < size; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

/* Preconditioned BiCGSTAB method, described in: */
/* http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method */

void Precond_Jacobi (
    Sparse_matrix * matrix,
    double        * c,
    double        * y,
    double          omega) {

  unsigned i;
  for (i = 0; i < matrix->size; ++i) {
    y[i] = omega * c[i] / matrix->elements[i];
  }

}

void Iterative_method_BiCGSTAB (
    Sparse_matrix * matrix,
    double        * x,
    double        * b,
    unsigned        max_iter,
    Precond_type    precond,
    double          omega_precond, /* default = 1 */
    double          accuracy,
    double        * buffer /* we need 8 * size */
  ) {

  unsigned size   = matrix->size;
  double alpha    = 1;
  double rho      = 1;
  double omega    = 1;
  double *r       = buffer;
  double *rcap    = buffer + size;
  double *v       = buffer + size * 2;
  double *p       = buffer + size * 3;
  double *y       = buffer + size * 4;
  double *z       = buffer + size * 5;
  double *s       = buffer + size * 6;
  double *t       = buffer + size * 7;
/*double *K1rev_t = buffer + size * 8;
  double *K1rev_s = buffer + size * 9; */
  double *diff    = y;

  unsigned i, iter;
  double beta, rhoold, residual2;

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
    beta = (rho / rhoold) * (alpha / omega);
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
 /* Sparse_matrix_Apply_to_vector (K1_rev, t, K1rev_t);
    Sparse_matrix_Apply_to_vector (K1_rev, s, K1rev_s);
    omega = scalar_product (K1rev_t, K1rev_s, size) /
            scalar_product (K1rev_t, K1rev_t, size); */
    omega = scalar_product (t, s, size) /
            scalar_product (t, t, size);
    for (i = 0; i < size; ++i) {
      x[i] += (alpha * p[i] + alpha * y[i] + omega * z[i]);
    }
    Sparse_matrix_Apply_to_vector (matrix, x, diff);
    for (i = 0; i < size; ++i) {
      diff[i] -= b[i];
    }
    residual2 = scalar_product (diff, diff, size);
    if (residual2 < accuracy * accuracy) {
      break;
    }
    for (i = 0; i < size; ++i) {
      r[i] = s[i] - omega * t[i];
    }
  }
}
