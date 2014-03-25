#include "sparse_matrix.h"

/* Some typedefs forked from Laspack */

typedef void *(*Precond_Type) (
    Sparse_matrix *,   /* preconditioner matrix */
    double *,          /* y (FIXME: describe it) */
    double *,          /* c (FIXME: describe it) */
    double,            /* omega (relaxation constant) */
    double *           /* result */
);

typedef void *(*Iterative_method_type) (
    Sparse_matrix *,   /* matrix */
    double *,          /* unknown_vector */
    double *,          /* rightcol */
    unsigned,          /* maximal number of iterations */
 /* Precond_type,         preconditioner type */
    Sparse_matrix *,   /* K_rev (reverse matrix to preconditioner K = K1 * K2) */
    Sparse_matrix *,   /* K1_rev (reverse matrix to K1) */
    double,            /* omega (relaxation constant) */
    double,            /* accuracy */
    double *          /* buffer (10 * size) */
);
