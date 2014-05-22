#ifndef _ITERATIVE_METHOD_H
#define _ITERATIVE_METHOD_H

#ifndef NO_LASPACK
#include <itersolv.h>
#endif /* NO_LASPACK */

#include "sparse_matrix.h"

/* Some typedefs forked from Laspack */

typedef void (*Precond_type) (
    Sparse_matrix const *,  /* preconditioner matrix */
    double const *,         /* c (source rightcol) */
    double       *,         /* y (target rightcol) */
    double                  /* omega (relaxation constant) */
);

typedef void (*Iterative_Method_type) (
    Sparse_matrix const *,  /* matrix */
    double *,               /* unknown_vector */
    double const *,         /* rightcol */
    unsigned,               /* maximal number of iterations */
    Precond_type,           /* preconditioner type */
    double,                 /* omega (relaxation constant) */
    double,                 /* accuracy */
    double *                /* buffer (8 * size) */
);

typedef enum {
    Implementation_Native,
    Implementation_Laspack
} Implementation_type;

/* structure to store the information about used
 * preconditioner type and iterative method */
typedef struct {
    Implementation_type implementation;
    Precond_type preconditioner_type;
    Iterative_Method_type method;
#ifndef NO_LASPACK
    PrecondProcType preconditioner_type_laspack;
    IterProcType method_laspack;
#endif /* NO_LASPACK */
    double relaxation_constant;
    double accuracy;
} Iterative_Method_parameters;

void Precond_Jacobi (
    Sparse_matrix const * matrix,
    double const        * c,
    double              * y,
    double                omega);

void Precond_Null (
    Sparse_matrix const * matrix,
    double const        * c,
    double              * y,
    double                omega);

void Iterative_method_BiCGSTAB (
    Sparse_matrix const * matrix,
    double              * x,
    double        const * b,
    unsigned              max_iter,
    Precond_type          precond,
    double                omega,
    double                accuracy,
    double              * buffer);

#endif /* _ITERATIVE_METHOD_H */
